# These functions are
# Copyright (C) 1998-2011 T.W. Yee, University of Auckland.
# All rights reserved.





















getMaxMin = function(vov, objfun, y, x, w, extraargs = NULL, maximize = TRUE,
                     abs.arg = FALSE) {
    if (!is.vector(vov)) stop("'vov' must be a vector")
    objvals = vov
    for(ii in 1:length(vov))
        objvals[ii] = objfun(vov[ii], y = y, x = x, w = w, extraargs=extraargs)
    try.this = if (abs.arg) {
                   if (maximize) vov[abs(objvals) == max(abs(objvals))] else
                   vov[abs(objvals) == min(abs(objvals))]
               } else {
                   if (maximize) vov[objvals == max(objvals)] else
                   vov[objvals == min(objvals)]
               }
    if (!length(try.this)) stop("something has gone wrong!")
    if (length(try.this) == 1) try.this else sample(try.this, size=1)
}



 mccullagh89 = function(ltheta = "rhobit", lnu = "logoff",
               itheta = NULL, inu = NULL,
               etheta = list(),
               enu = if (lnu == "logoff") list(offset = 0.5) else list(),
               zero = NULL)
{
    if (mode(ltheta) != "character" && mode(ltheta) != "name")
        ltheta = as.character(substitute(ltheta))
    if (mode(lnu) != "character" && mode(lnu) != "name")
        lnu = as.character(substitute(lnu))
    if (length(zero) && !is.Numeric(zero, integer = TRUE, posit = TRUE))
        stop("bad input for argument 'zero'")
    if (!is.list(etheta)) etheta = list()
    if (!is.list(enu)) enu = list()

    new("vglmff",
    blurb = c("McCullagh (1989)'s distribution \n",
    "f(y) = (1-2*theta*y+theta^2)^(-nu) * [1 - y^2]^(nu-1/2) /\n",
            "       Beta[nu+1/2, 1/2], ",
            "  -1 < y < 1, -1 < theta < 1, nu > -1/2\n",
            "Links:     ",
            namesof("theta", ltheta, earg = etheta), ", ",
            namesof("nu", lnu, earg = enu),
            "\n",
            "\n",
            "Mean:     nu*theta/(1+nu)"),
    constraints = eval(substitute(expression({
        constraints = cm.zero.vgam(constraints, x, .zero, M)
    }), list( .zero = zero ))),
    initialize = eval(substitute(expression({
        y = as.numeric(y)
        if (any(y <= -1 | y >= 1))
            stop("all y values must be in (-1,1)")

        predictors.names =
          c(namesof("theta", .ltheta, earg = .etheta, tag = FALSE),
            namesof("nu",    .lnu,    earg = .enu,    tag = FALSE))
        if (!length(etastart)) {
            theta.init = if (length( .itheta)) rep( .itheta, length = n) else {
                mccullagh89.aux = function(thetaval, y, x, w, extraargs)
                mean((y-thetaval)*(thetaval^2-1)/(1-2*thetaval*y+thetaval^2))
                theta.grid = seq(-0.9, 0.9, by=0.05)
                try.this = getMaxMin(theta.grid, objfun=mccullagh89.aux,
                                     y = y,  x = x, w = w, maximize = FALSE,
                                     abs.arg = TRUE)
                try.this = rep(try.this, len = n)
                try.this
            }
            tmp = y / (theta.init-y)
            tmp[tmp < -0.4] = -0.4
            tmp[tmp > 10.0] = 10.0
            nu.init = rep(if (length( .inu)) .inu else tmp, length = n)
            nu.init[!is.finite(nu.init)] = 0.4
            etastart = cbind(theta2eta(theta.init, .ltheta, earg = .etheta ),
                             theta2eta(nu.init, .lnu, earg = .enu ))
        }
    }), list( .ltheta=ltheta, .lnu=lnu, .inu=inu, .itheta=itheta,
              .etheta = etheta, .enu=enu ))),
    inverse = eval(substitute(function(eta, extra = NULL) {
        Theta = eta2theta(eta[,1], .ltheta, earg = .etheta )
        nu = eta2theta(eta[,2], .lnu, earg = .enu )
        nu*Theta/(1+nu)
    }, list( .ltheta=ltheta, .lnu=lnu,
             .etheta = etheta, .enu=enu ))),
    last = eval(substitute(expression({
        misc$link =    c("theta" = .ltheta, "nu" = .lnu)
        misc$earg = list("theta" = .etheta, "nu" = .enu )
    }), list( .ltheta=ltheta, .lnu=lnu, .etheta = etheta, .enu=enu ))),
    loglikelihood = eval(substitute(
        function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
        Theta = eta2theta(eta[,1], .ltheta, earg = .etheta )
        nu = eta2theta(eta[,2], .lnu, earg = .enu )
        if (residuals) stop("loglikelihood residuals not ",
                            "implemented yet") else
            sum(w * ((nu-0.5)*log1p(-y^2) - nu * log1p(-2*Theta*y + Theta^2) -
                    lbeta(nu+0.5,0.5 )))
    }, list( .ltheta=ltheta, .lnu=lnu, .etheta = etheta, .enu=enu ))),
    vfamily = c("mccullagh89"),
    deriv = eval(substitute(expression({
        Theta = eta2theta(eta[,1], .ltheta, earg = .etheta )
        nu = eta2theta(eta[,2], .lnu, earg = .enu )
        dTheta.deta = dtheta.deta(Theta, .ltheta, earg = .etheta )
        dnu.deta = dtheta.deta(nu, .lnu, earg = .enu )
        dl.dTheta = 2 * nu * (y-Theta) / (1 -2*Theta*y + Theta^2)
        dl.dnu = log1p(-y^2) - log1p(-2*Theta*y + Theta^2) -
                 digamma(nu+0.5) + digamma(nu+1)
        c(w) * cbind(dl.dTheta * dTheta.deta,
                     dl.dnu * dnu.deta)
    }), list( .ltheta=ltheta, .lnu=lnu, .etheta = etheta, .enu=enu ))),
    weight = eval(substitute(expression({
        d2l.dTheta2 = (2 * nu^2 / (1+nu)) / (1-Theta^2)
        d2l.dnu2 = trigamma(nu+0.5) - trigamma(nu+1)
        wz = matrix(as.numeric(NA), n, M)  #diagonal matrix
        wz[,iam(1,1,M)] = d2l.dTheta2 * dTheta.deta^2
        wz[,iam(2,2,M)] = d2l.dnu2 * dnu.deta^2
        c(w) * wz
    }), list( .ltheta=ltheta, .lnu=lnu ))))
}




hzeta.control <- function(save.weight = TRUE, ...)
{
    list(save.weight = save.weight)
}



 hzeta = function(link = "loglog", earg = list(), ialpha = NULL, nsimEIM = 100)
{

    stopifnot(ialpha > 0)
    stopifnot(nsimEIM > 10, length(nsimEIM) == 1, nsimEIM == round(nsimEIM))

    if (mode(link) != "character" && mode(link) != "name")
        link = as.character(substitute(link))
    if (!is.list(earg)) earg = list()

    new("vglmff",
    blurb = c(
    "Haight's Zeta distribution f(y) = (2y-1)^(-alpha) - (2y+1)^(-alpha),\n",
            "    alpha>0, y = 1,2,....\n\n",
            "Link:    ",
            namesof("alpha", link, earg = earg), "\n\n",
            "Mean:     (1-2^(-alpha)) * zeta(alpha) if alpha>1",
            "\n",
            "Variance: (1-2^(1-alpha)) * zeta(alpha-1) - mean^2 if alpha>2"),
    initialize = eval(substitute(expression({
        y = as.numeric(y)
        if (ncol(cbind(y)) != 1)
            stop("response must be a vector or a one-column matrix")
        if (any(y < 1))
            stop("all y values must be in 1,2,3,....")
        predictors.names = namesof("alpha", .link, earg = .earg, tag = FALSE)
        if (!length(etastart)) {
            a.init = if (length( .ialpha)) .ialpha else {
                if ((meany <- weighted.mean(y,w)) < 1.5) 3.0 else
                if (meany < 2.5) 1.4 else 1.1 
            }
            a.init = rep(a.init, length = n) 
            etastart = theta2eta(a.init, .link, earg = .earg )
        }
    }), list( .link = link, .earg = earg, .ialpha=ialpha ))),
    inverse = eval(substitute(function(eta, extra = NULL) {
        alpha = eta2theta(eta, .link, earg = .earg )
        mu = (1-2^(-alpha)) * zeta(alpha)
        mu[alpha <= 1] = Inf
        mu
    }, list( .link = link, .earg = earg ))),
    last = eval(substitute(expression({
        misc$link =    c(alpha = .link)
        misc$earg = list(alpha = .earg)
        misc$nsimEIM = .nsimEIM
    }), list( .link = link, .earg = earg, .nsimEIM = nsimEIM ))),
    loglikelihood = eval(substitute(
        function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
        alpha = eta2theta(eta, .link, earg = .earg )
        if (residuals) stop("loglikelihood residuals not ",
                            "implemented yet") else {
            sum(w * dhzeta(x = y, alpha=alpha, log = TRUE))
        }
    }, list( .link = link, .earg = earg ))),
    vfamily = c("hzeta"),
    deriv = eval(substitute(expression({
        alpha = eta2theta(eta, .link, earg = .earg ) 
        dalpha.deta = dtheta.deta(alpha, .link, earg = .earg )
        d3 = deriv3(~ log((2*y-1)^(-alpha) - (2*y+1)^(-alpha)),
                    "alpha", hessian = FALSE)
        eval.d3 = eval(d3)
        dl.dalpha =  attr(eval.d3, "gradient")
        c(w) * dl.dalpha * dalpha.deta
    }), list( .link = link, .earg = earg ))),
    weight = eval(substitute(expression({
        sd3 = deriv3(~ log((2*ysim-1)^(-alpha) - (2*ysim+1)^(-alpha)),
                     "alpha", hessian = FALSE)
        run.var = 0
        for(ii in 1:( .nsimEIM )) {
            ysim = rhzeta(n, alpha=alpha)
            eval.sd3 = eval(sd3)
            dl.dalpha =  attr(eval.d3, "gradient")
            rm(ysim)
            temp3 = dl.dalpha
            run.var = ((ii-1) * run.var + temp3^2) / ii
        }
        wz = if (intercept.only)
            matrix(colMeans(cbind(run.var)),
                   n, dimm(M), byrow = TRUE) else cbind(run.var)

        wz = wz * dalpha.deta^2
        c(w) * wz
    }), list( .nsimEIM = nsimEIM ))))
}




dhzeta = function(x, alpha, log = FALSE)
{
    if (!is.logical(log.arg <- log))
        stop("bad input for argument 'log'")
    rm(log)

    if (!is.Numeric(alpha, posit = TRUE))
        stop("'alpha' must be numeric and have positive values")
    nn = max(length(x), length(alpha))
    x = rep(x, len = nn); alpha = rep(alpha, len = nn)
    ox = !is.finite(x)
    zero = ox | round(x) != x | x < 1
    ans = rep(0, len = nn)
    ans[!zero] = (2*x[!zero]-1)^(-alpha[!zero]) - (2*x[!zero]+1)^(-alpha[!zero])
    if (log.arg) log(ans) else ans
}


phzeta = function(q, alpha) 
{
    if (!is.Numeric(alpha, posit = TRUE))
        stop("'alpha' must be numeric and have positive values")
    nn = max(length(q), length(alpha))
    q = rep(q, len = nn)
    alpha = rep(alpha, len = nn)
    oq = !is.finite(q)
    zero = oq | q < 1
    q = floor(q)
    ans = 0 * q
    ans[!zero] = 1 - (2*q[!zero]+1)^(-alpha[!zero])
    ans
}


qhzeta = function(p, alpha) 
{
    if (!is.Numeric(alpha, posit = TRUE))
        stop("'alpha' must be numeric and have positive values")
    if (!is.Numeric(p, posit = TRUE) || any(p >= 1))
        stop("argument 'p' must have values inside the interval (0,1)")
    nn = max(length(p), length(alpha))
    p = rep(p, len = nn)
    alpha = rep(alpha, len = nn)
    ans = (((1 - p)^(-1/alpha) - 1) / 2) # p is in (0,1)
    floor(ans+1)
}

rhzeta = function(n, alpha) 
{
    if (!is.Numeric(alpha, posit = TRUE))
        stop("'alpha' must be numeric and have positive values")
    if (!is.Numeric(n, posit = TRUE, integ = TRUE, allow = 1))
        stop("argument 'n' must be a positive integer")
    ans = ((runif(n)^(-1/alpha) - 1) / 2)
    floor(ans+1)
}






 dirmultinomial <- function(lphi = "logit", ephi = list(),
                            iphi = 0.10, parallel = FALSE, zero = "M")
{

  if (mode(lphi) != "character" && mode(lphi) != "name")
    lphi <- as.character(substitute(lphi))
  if (length(zero) && 
     !(is.Numeric(zero, integer = TRUE, posit = TRUE) || is.character(zero )))
    stop("bad input for argument 'zero'")
  if (!is.Numeric(iphi, positive = TRUE) || max(iphi) >= 1.0)
    stop("bad input for argument 'iphi'")
  if (!is.list(ephi)) ephi <- list()

  new("vglmff",
  blurb = c("Dirichlet-multinomial distribution\n\n",
            "Links:    ",
            "log(prob[1]/prob[M]), ..., log(prob[M-1]/prob[M]), ",
            namesof("phi", lphi, earg = ephi), "\n", "\n",
            "Mean:     shape_j / sum_j(shape_j)"),
  constraints = eval(substitute(expression({
    .ZERO <- .zero
    if (is.character( .ZERO)) .ZERO <- eval(parse(text = .ZERO))
    .PARALLEL <- .parallel
    if (is.logical( .PARALLEL) && .PARALLEL) {
      mycmatrix <- if (length( .ZERO))
          stop("can only handle parallel = TRUE when zero = NULL") else
          cbind(rbind(matrix(1, M - 1, 1), 0), rbind(matrix(0, M - 1, 1), 1))
    } else
      mycmatrix <- if (M == 1) diag(1) else diag(M)
      constraints <- cm.vgam(mycmatrix, x, .PARALLEL,
                             constraints, int = TRUE)
      constraints <- cm.zero.vgam(constraints, x, .ZERO, M)
  }), list( .parallel = parallel, .zero = zero ))),
  initialize = eval(substitute(expression({
    delete.zero.colns <- TRUE
    eval(process.categorical.data.vgam)

    y <- as.matrix(y)
    ycount <- as.matrix(y * c(w))
    M <- ncol(y)
    if (max(abs(ycount - round(ycount ))) > 1.0e-6)
      warning("there appears to be non-integer responses")
    if (min(ycount) < 0)
      stop("all values of the response (matrix) must be non-negative")
    predictors.names <-
      c(paste("log(prob[,", 1:(M-1), "]/prob[,", M, "])", sep = ""),
        namesof("phi", .lphi, short = TRUE))
    extra$n2 <- w  # aka omega, must be integer # as.vector(apply(y, 1, sum))
    if (!length(etastart)) {
      prob.init <- colSums(ycount)
      prob.init <- prob.init / sum(prob.init)
      prob.init <- matrix(prob.init, n, M, byrow = TRUE)
      phi.init <- rep( .iphi, len = n)
      etastart <- cbind(log(prob.init[,-M]/prob.init[,M]),
                        theta2eta(phi.init, .lphi, earg = .ephi ))
    }
  }), list( .lphi = lphi, .ephi = ephi, .iphi=iphi ))),
  inverse = eval(substitute(function(eta, extra = NULL) {
    M <- if (is.matrix(eta)) ncol(eta) else 1
    temp <- cbind(exp(eta[,-M]), 1)
    temp / as.vector(temp %*% rep(1, M))
  }, list( .ephi = ephi, .lphi = lphi ))),
  last = eval(substitute(expression({
      misc$link <- c(rep("noLinkFunction", length = M-1), .lphi)
      names(misc$link) <- c(paste("prob", 1:(M-1), sep = ""), "phi")
      misc$earg <- vector("list", M)
      names(misc$earg) <- names(misc$link)
      for(ii in 1:(M-1)) misc$earg[[ii]] <- list()
      misc$earg[[M]] <- .ephi
      misc$expected <- TRUE
      if (intercept.only) {
        misc$shape<-probs[1,]*(1/phi[1]-1) # phi & probs computed in @deriv
      }
  }), list( .ephi = ephi, .lphi = lphi ))),
  loglikelihood = eval(substitute(
      function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
      M = if (is.matrix(eta)) ncol(eta) else 1
      probs <- cbind(exp(eta[,-M]), 1)
      probs <- probs / as.vector(probs %*% rep(1, M))
      phi <- eta2theta(eta[,M], .lphi, earg = .ephi )
      n <- length(phi)
      ycount <- as.matrix(y * c(w))
      if (residuals) stop("loglikelihood residuals not ",
                          "implemented yet") else {
        ans <- rep(0.0, len = n)
        omega <- extra$n2
        for(jay in 1:M) {
          maxyj <- max(ycount[,jay])
          loopOveri <- n < maxyj
          if (loopOveri) {
            for(iii in 1:n) {
                rrr <- 1:ycount[iii,jay] # a vector
                if (ycount[iii,jay] > 0)
                ans[iii] <- ans[iii] + sum(log((1-phi[iii]) *
                           probs[iii,jay] + (rrr-1)*phi[iii]))
            }
          } else {
            for(rrr in 1:maxyj) {
                index <- (rrr <= ycount[,jay]) & (ycount[,jay] > 0)
                if (any(index))
                    ans[index] <- ans[index] + log((1-phi[index]) *
                                 probs[index,jay] + (rrr-1)*phi[index])
            }
          }
        } # end of jay loop

        maxomega <- max(omega)
        loopOveri <- n < maxomega
        if (loopOveri) {
          for(iii in 1:n) {
            rrr <- 1:omega[iii]
            ans[iii]<- ans[iii] - sum(log1p(-phi[iii] + (rrr-1)*phi[iii]))
          }
        } else {
          for(rrr in 1:maxomega) {
            ind8 <- rrr <= omega
            ans[ind8] <- ans[ind8] - log1p(-phi[ind8] + (rrr-1)*phi[ind8])
          }
        }
        sum(ans)
    }
  }, list( .ephi = ephi, .lphi = lphi ))),
  vfamily = c("dirmultinomial "),
  deriv = eval(substitute(expression({
    probs <- cbind(exp(eta[,-M]), 1)
    probs <- probs / as.vector(probs %*% rep(1, M))
    phi <- eta2theta(eta[,M], .lphi, earg = .ephi )
    dl.dprobs <- matrix(0.0, n, M-1)
    dl.dphi <- rep(0.0, len = n)
    omega <- extra$n2
    ycount <- as.matrix(y * c(w))
    for(jay in 1:M) {
        maxyj <- max(ycount[,jay])
        loopOveri <- n < maxyj
        if (loopOveri) {
          for(iii in 1:n) {
            rrr <- 1:ycount[iii,jay]
            if (ycount[iii,jay] > 0) {
              PHI <- phi[iii]
              dl.dphi[iii] <- dl.dphi[iii] +
 sum((rrr-1-probs[iii,jay]) / ((1-PHI)*probs[iii,jay] + (rrr-1)*PHI))

              tmp9 <- (1-PHI) / ((1-PHI)*probs[iii,jay] + (rrr-1)*PHI)
              if (jay < M) {
                  dl.dprobs[iii,jay] <- dl.dprobs[iii,jay] + sum(tmp9)
              } else {
                  for(jay2 in 1:(M-1))
                     dl.dprobs[iii,jay2]<-dl.dprobs[iii,jay2]-sum(tmp9)
              }
            }
          }
        } else {
          for(rrr in 1:maxyj) {
            index <- (rrr <= ycount[,jay]) & (ycount[,jay] > 0)
            PHI <- phi[index]
            dl.dphi[index] <- dl.dphi[index] +
              (rrr-1-probs[index,jay]) / ((1-PHI)*probs[index,jay] +
              (rrr-1)*PHI)
            tmp9 <- (1-PHI) / ((1-PHI)*probs[index,jay] + (rrr-1)*PHI)
            if (jay < M) {
                dl.dprobs[index,jay] <- dl.dprobs[index,jay] + tmp9
            } else {
                for(jay2 in 1:(M-1))
                    dl.dprobs[index,jay2] <- dl.dprobs[index,jay2] - tmp9
            }
          }
        }
    } # end of jay loop
    maxomega <- max(omega)
    loopOveri <- n < maxomega
    if (loopOveri) {
      for(iii in 1:n) {
        rrr <- 1:omega[iii]
        dl.dphi[iii]<-dl.dphi[iii] - sum((rrr-2)/(1 + (rrr-2)*phi[iii]))
      }
    } else {
      for(rrr in 1:maxomega) {
        index <- rrr <= omega
        dl.dphi[index]<-dl.dphi[index] - (rrr-2)/(1 + (rrr-2)*phi[index])
      }
    }
    dprobs.deta <- probs[,-M] * (1 - probs[,-M])    # n x (M-1)
    dphi.deta <- dtheta.deta(phi, .lphi, earg = .ephi )
    ans <- cbind(dl.dprobs * dprobs.deta,
                 dl.dphi   * dphi.deta)
    ans
  }), list( .ephi = ephi, .lphi = lphi ))),
    weight = eval(substitute(expression({
      wz <- matrix(0, n, dimm(M))
      loopOveri <- n < maxomega
      if (loopOveri) {
          for(iii in 1:n) {
              rrr <- 1:omega[iii]  # A vector
              PHI <- phi[iii]
              pYiM.ge.rrr <- 1 - pbetabin.ab(q=rrr-1, size=omega[iii],
                  shape1<-probs[iii,M]*(1/PHI-1),
                  shape2<-(1-probs[iii,M])*(1/PHI-1))  # A vector
              denomM <- ((1-PHI)*probs[iii,M] + (rrr-1)*PHI)^2  # A vector
              wz[iii,iam(M,M,M)] <- wz[iii,iam(M,M,M)] +
                      sum(probs[iii,M]^2 * pYiM.ge.rrr / denomM) -
                      sum(1 / (1 + (rrr-2)*PHI)^2)
              for(jay in 1:(M-1)) {
                  denomj <- ((1-PHI)*probs[iii,jay] + (rrr-1)*PHI)^2
                  pYij.ge.rrr <- 1 - pbetabin.ab(q=rrr-1, size=omega[iii],
                      shape1<-probs[iii,jay]*(1/PHI-1),
                      shape2<-(1-probs[iii,jay])*(1/PHI-1))
                  wz[iii,iam(jay,jay,M)] <- wz[iii,iam(jay,jay,M)] + 
                      sum(pYij.ge.rrr / denomj) + 
                      sum(pYiM.ge.rrr / denomM)
                  for(kay in jay:(M-1)) if (kay > jay) {
                      wz[iii,iam(jay,kay,M)] <- wz[iii,iam(jay,kay,M)] + 
                          sum(pYiM.ge.rrr / denomM)
                  }
                  wz[iii,iam(jay,M,M)] <- wz[iii,iam(jay,M,M)] +
                          sum(probs[iii,jay] * pYij.ge.rrr / denomj) -
                          sum(probs[iii,M]   * pYiM.ge.rrr / denomM)
                  wz[iii,iam(M,M,M)] <- wz[iii,iam(M,M,M)] +
                          sum(probs[iii,jay]^2 * pYij.ge.rrr / denomj)
              } # end of jay loop
          } # end of iii loop
      } else {
          for(rrr in 1:maxomega) {
              ind5 <- rrr <= omega
              PHI <- phi[ind5]
              pYiM.ge.rrr <- 1 - pbetabin.ab(q=rrr-1, size=omega[ind5],
                  shape1<-probs[ind5,M]*(1/PHI-1),
                  shape2<-(1-probs[ind5,M])*(1/PHI-1))
              denomM <- ((1-PHI)*probs[ind5,M] + (rrr-1)*PHI)^2
              wz[ind5,iam(M,M,M)] <- wz[ind5,iam(M,M,M)] +
                      probs[ind5,M]^2 * pYiM.ge.rrr / denomM -
                      1 / (1 + (rrr-2)*PHI)^2
              for(jay in 1:(M-1)) {
                  denomj <- ((1-PHI)*probs[ind5,jay] + (rrr-1)*PHI)^2
                  pYij.ge.rrr <- 1 - pbetabin.ab(q=rrr-1, size=omega[ind5],
                      shape1<-probs[ind5,jay]*(1/PHI-1),
                      shape2<-(1-probs[ind5,jay])*(1/PHI-1))
                  wz[ind5,iam(jay,jay,M)] <- wz[ind5,iam(jay,jay,M)] + 
                      pYij.ge.rrr / denomj + pYiM.ge.rrr / denomM 
                  for(kay in jay:(M-1)) if (kay > jay) {
                      wz[ind5,iam(jay,kay,M)] <- wz[ind5,iam(jay,kay,M)] + 
                          pYiM.ge.rrr / denomM 
                  }
                  wz[ind5,iam(jay,M,M)] <- wz[ind5,iam(jay,M,M)] +
                          probs[ind5,jay] * pYij.ge.rrr / denomj -
                          probs[ind5,M]   * pYiM.ge.rrr / denomM
                  wz[ind5,iam(M,M,M)] <- wz[ind5,iam(M,M,M)] +
                          probs[ind5,jay]^2 * pYij.ge.rrr / denomj
              } # end of jay loop
          } # end of rrr loop
      }

      for(jay in 1:(M-1))
          for(kay in jay:(M-1))
              wz[,iam(jay,kay,M)] <- wz[,iam(jay,kay,M)] * (1-phi)^2
      for(jay in 1:(M-1))
          wz[,iam(jay,M,M)] <- wz[,iam(jay,M,M)] * (phi-1) / phi
      wz[,iam(M,M,M)] <- wz[,iam(M,M,M)] / phi^2

      d1Thetas.deta <- cbind(dprobs.deta, dphi.deta)
      index <- iam(NA, NA, M, both = TRUE, diag = TRUE)
      wz <- wz * d1Thetas.deta[,index$row] * d1Thetas.deta[,index$col]
      wz
  }), list( .ephi = ephi, .lphi = lphi ))))
}





dirmul.old = function(link = "loge", earg = list(), init.alpha = 0.01,
                      parallel = FALSE, zero = NULL)
{

    if (mode(link) != "character" && mode(link) != "name")
        link = as.character(substitute(link))
    if (length(zero) && !is.Numeric(zero, integer = TRUE, posit = TRUE))
        stop("bad input for argument 'zero'")
    if (!is.Numeric(init.alpha, posit = TRUE))
        stop("'init.alpha' must contain positive values only")
    if (!is.list(earg)) earg = list()

    new("vglmff",
    blurb = c("Dirichlet-Multinomial distribution\n\n",
              "Links:     ",
              namesof("shape1", link, earg = earg), ", ..., ",
              namesof("shapeM", link, earg = earg), "\n\n",
            "Posterior mean:    (n_j + shape_j)/(2*sum(n_j) + sum(shape_j))\n"),
    constraints = eval(substitute(expression({
        constraints = cm.vgam(matrix(1, M, 1), x, .parallel,
                              constraints, int = TRUE)
        constraints = cm.zero.vgam(constraints, x, .zero, M)
    }), list( .parallel = parallel, .zero = zero ))),
    initialize = eval(substitute(expression({
        y = as.matrix(y)
        M = ncol(y)
        if (any(y != round(y )))
            stop("all y values must be integer-valued")

        predictors.names = namesof(paste("shape", 1:M, sep = ""),
                                   .link, earg = .earg, short = TRUE)
        extra$n2 = rowSums(y)  # Nb. don't multiply by 2
        extra$y  = y
        if (!length(etastart)) {
            yy = if (is.numeric( .init.alpha)) 
                matrix( .init.alpha, n, M, byrow= TRUE) else
                matrix(runif(n*M), n, M)
            etastart = theta2eta(yy, .link, earg = .earg)
        }
    }), list( .link = link, .earg = earg, .init.alpha=init.alpha ))),
    inverse = eval(substitute(function(eta, extra = NULL) {
        shape = eta2theta(eta, .link, earg = .earg)
        M = if (is.matrix(eta)) ncol(eta) else 1
        sumshape = as.vector(shape %*% rep(1, len = M))
        (extra$y + shape) / (extra$n2 + sumshape)
    }, list( .link = link, .earg = earg ))),
    last = eval(substitute(expression({
        misc$link = rep( .link, length = M)
        names(misc$link) = paste("shape", 1:M, sep = "")
        misc$earg = vector("list", M)
        names(misc$earg) = names(misc$link)
        for(ii in 1:M) misc$earg[[ii]] = .earg
        misc$pooled.weight = pooled.weight
    }), list( .link = link, .earg = earg ))),
    loglikelihood = eval(substitute(
        function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
        shape = eta2theta(eta, .link, earg = .earg)
        M = if (is.matrix(eta)) ncol(eta) else 1
        sumshape = as.vector(shape %*% rep(1, len = M))
        if (residuals) stop("loglikelihood residuals not ",
                            "implemented yet") else
        sum(w*(lgamma(sumshape) - lgamma(extra$n2 + sumshape ))) +
            sum(w * (lgamma(y + shape) - lgamma(shape )))
    }, list( .link = link, .earg = earg ))),
    vfamily = c("dirmul.old"),
    deriv = eval(substitute(expression({
        shape = eta2theta(eta, .link, earg = .earg)
        sumshape = as.vector(shape %*% rep(1, len = M))
        dl.dsh = digamma(sumshape) - digamma(extra$n2 + sumshape) +
                 digamma(y + shape) - digamma(shape)
        dsh.deta = dtheta.deta(shape, .link, earg = .earg)
        c(w) * dl.dsh * dsh.deta
    }), list( .link = link, .earg = earg ))),
    weight = eval(substitute(expression({
        index = iam(NA, NA, M, both = TRUE, diag = TRUE)
        wz = matrix(trigamma(sumshape)-trigamma(extra$n2 + sumshape),
                    nrow=n, ncol=dimm(M))
        wz[,1:M] = wz[,1:M] + trigamma(y + shape) - trigamma(shape)
        wz = -wz * dsh.deta[, index$row] * dsh.deta[, index$col]


        if (TRUE && intercept.only) {
            sumw = sum(w)
            for(ii in 1:ncol(wz))
                wz[,ii] = sum(wz[,ii]) / sumw
            pooled.weight = TRUE
            wz = c(w) * wz   # Put back the weights
        } else
            pooled.weight = FALSE

        wz
    }), list( .link = link, .earg = earg ))))
}






rdiric = function(n, shape, dimension = NULL) {
    if (!is.numeric(dimension))
        dimension = length(shape)
    shape = rep(shape, len=dimension)

    ans = rgamma(n*dimension, rep(shape, rep(n, dimension)))
    dim(ans) = c(n, dimension) 


    ans = ans / rowSums(ans)
    ans
}




 dirichlet = function(link = "loge", earg = list(),
                      parallel = FALSE, zero = NULL)
{
    if (mode(link) != "character" && mode(link) != "name")
        link = as.character(substitute(link))
    if (length(zero) &&
        !is.Numeric(zero, integer = TRUE, posit = TRUE))
        stop("bad input for argument 'zero'")
    if (!is.list(earg)) earg = list()

    new("vglmff",
    blurb = c("Dirichlet distribution\n\n",
              "Links:     ",
              namesof("shapej", link, earg = earg), "\n\n",
              "Mean:     shape_j/(1 + sum(shape_j)), j = 1,..,ncol(y)"),
    constraints = eval(substitute(expression({
        constraints = cm.vgam(matrix(1, M, 1), x, .parallel, constraints, int= TRUE)
        constraints = cm.zero.vgam(constraints, x, .zero, M)
    }), list( .parallel = parallel, .zero = zero ))),
    initialize = eval(substitute(expression({
        y = as.matrix(y)
        M = ncol(y)
        if (any(y <= 0) || any(y>=1))
            stop("all y values must be > 0 and < 1")
        predictors.names = namesof(paste("shape", 1:M, sep = ""), .link,
                                   earg = .earg, short = TRUE)
        if (!length(etastart)) {
            yy = matrix(t(y) %*% rep(1/nrow(y), nrow(y)), nrow(y), M, byrow= TRUE)
            etastart = theta2eta(yy, .link, earg = .earg )
        }
    }), list( .link = link, .earg = earg ))),
    inverse = eval(substitute(function(eta, extra = NULL) {
        shape = eta2theta(eta, .link, earg = .earg )
        M = if (is.matrix(eta)) ncol(eta) else 1
        sumshape = rowSums(shape)
        shape / sumshape
    }, list( .link = link, .earg = earg ))),
    last = eval(substitute(expression({
        misc$link = c(shape= .link)
        temp.names = paste("shape", 1:M, sep = "")
        misc$link = rep( .link, len = M)
        names(misc$link) = temp.names
        misc$earg = vector("list", M)
        names(misc$earg) = names(misc$link)
        for(ii in 1:M) misc$earg[[ii]] = .earg
        misc$expected = TRUE
    }), list( .link = link, .earg = earg ))),
    loglikelihood = eval(substitute(
        function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
        shape = eta2theta(eta, .link, earg = .earg )
        M = if (is.matrix(eta)) ncol(eta) else 1
        sumshape = as.vector(shape %*% rep(1, len = M))
        if (residuals) stop("loglikelihood residuals not ",
                            "implemented yet") else {
          sum(c(w) * lgamma(sumshape)) -
          sum(c(w) * lgamma(shape)) +
          sum(c(w) * (shape-1) * log(y))
        }
    }, list( .link = link, .earg = earg ))),
    vfamily = c("dirichlet"),
    deriv = eval(substitute(expression({
        shape = eta2theta(eta, .link, earg = .earg )
        sumshape = as.vector(shape %*% rep(1, len = M))
        dl.dsh = digamma(sumshape) - digamma(shape) + log(y)
        dsh.deta = dtheta.deta(shape, .link, earg = .earg )
        c(w) * dl.dsh * dsh.deta
    }), list( .link = link, .earg = earg ))),
    weight = expression({
        index = iam(NA, NA, M, both = TRUE, diag = TRUE)
        wz = matrix(trigamma(sumshape), nrow=n, ncol=dimm(M))
        wz[,1:M] = wz[,1:M] - trigamma(shape)
        wz = -c(w) * wz * dsh.deta[, index$row] * dsh.deta[, index$col]
        wz
    }))
}




 zeta = function(x, deriv = 0) {



    deriv.arg = deriv
    rm(deriv)
    if (!is.Numeric(deriv.arg, allow = 1, integer = TRUE))
        stop("'deriv' must be a single non-negative integer")
    if (deriv.arg < 0 || deriv.arg > 2)
        stop("'deriv' must be 0, 1, or 2")


    if (deriv.arg > 0)
        return(Zeta.derivative(x, deriv.arg = deriv.arg))



    if (any(special <- Re(x) <= 1)) {
        ans <- x
        ans[special] <- Inf   # For Re(x) == 1

        special3 <- Re(x) < 1
        ans[special3] <- NA # For 0 < Re(x) < 1

        special4 <- (0 < Re(x)) & (Re(x) < 1) & (Im(x) == 0)
        ans[special4] <- Zeta.derivative(x[special4], deriv.arg = deriv.arg)


        special2 <- Re(x) < 0
        if (any(special2)) {
            x2 = x[special2]
            cx = 1-x2
            ans[special2] = 2^(x2) * pi^(x2-1) * sin(pi*x2/2) * gamma(cx) * Recall(cx)
        }

        if (any(!special)) {
            ans[!special] <- Recall(x[!special])
        }
        return(ans)
    }

    a = 12; k = 8
    B = c(1/6, -1/30,1/42,-1/30,5/66,-691/2730,7/6,-3617/510)
    ans = 0
    for(ii in 1:(a-1))
       ans = ans + 1.0 / ii^x
    ans = ans + 1.0 / ((x-1.0)* a^(x-1.0)) + 1.0 / (2.0 * a^x)

    term = (x/2) / a^(x+1)
    ans = ans + term * B[1]

    for(mm in 2:k) {
        term = term * (x+2*mm-2) * (x+2*mm-3) / (a * a * 2 * mm * (2*mm-1))
        ans = ans + term * B[mm]
    }
    ans
}



 Zeta.derivative = function(x, deriv.arg = 0)
{


    if (!is.Numeric(deriv.arg, allow = 1, integer = TRUE))
        stop("'deriv.arg' must be a single non-negative integer")
    if (deriv.arg < 0 || deriv.arg > 2)
        stop("'deriv.arg' must be 0, 1, or 2")

    if (any(Im(x) != 0))
        stop("Sorry, currently can only handle x real, not complex")
    if (any(x < 0))
        stop("Sorry, currently cannot handle x < 0")

    ok = is.finite(x) & x > 0 & x != 1   # Handles NAs
    ans = rep(as.numeric(NA), length(x))
    nn = sum(ok)  # Effective length (excludes x < 0 and x = 1 values)
    if (nn)
        ans[ok] = dotC(name = "vzetawr", as.double(x[ok]), ans = double(nn),
                  as.integer(deriv.arg), as.integer(nn))$ans



    if (deriv.arg == 0)
        ans[is.finite(x) & abs(x) < 1.0e-12] = -0.5

    ans
}



dzeta = function(x, p, log = FALSE)
{
    if (!is.logical(log.arg <- log))
        stop("bad input for argument 'log'")
    rm(log)

    if (!is.Numeric(p, posit = TRUE)) # || min(p) <= 1
        stop("'p' must be numeric and > 0")
    LLL = max(length(p), length(x))
    x = rep(x, len = LLL); p = rep(p, len = LLL)

    ox = !is.finite(x)
    zero = ox | round(x) != x | x < 1
    if (any(zero)) warning("non-integer x and/or x < 1 or NAs")
    ans = rep(if (log.arg) log(0) else 0, len = LLL)
    if (any(!zero)) {
        if (log.arg) {
            ans[!zero] = (-p[!zero]-1)*log(x[!zero]) - log(zeta(p[!zero]+1))
        } else {
            ans[!zero] = x[!zero]^(-p[!zero]-1) / zeta(p[!zero]+1)
        }
    }
    if (any(ox)) ans[ox] = NA
    ans
}

 zetaff = function(link = "loge", earg = list(), init.p = NULL)
{

    if (length(init.p) && !is.Numeric(init.p, positi = TRUE))
        stop("argument 'init.p' must be > 0")
    if (mode(link) != "character" && mode(link) != "name")
        link = as.character(substitute(link))
    if (!is.list(earg)) earg = list()

    new("vglmff",
    blurb = c("Zeta distribution ",
              "f(y) = 1/(y^(p+1) zeta(p+1)), p>0, y = 1,2,..\n\n",
              "Link:    ",
              namesof("p", link, earg = earg), "\n\n",
              "Mean:     zeta(p) / zeta(p+1), provided p>1\n",
              "Variance: zeta(p-1) / zeta(p+1) - mean^2, provided p>2"),
    initialize = eval(substitute(expression({
        y = as.numeric(y)
        if (any(y < 1))
            stop("all y values must be in 1,2,3,...")
        if (any(y != round(y )))
            warning("'y' should be integer-valued")
        if (ncol(cbind(y)) != 1)
            stop("response must be a vector or a one-column matrix")

        predictors.names = namesof("p", .link, earg = .earg, tag = FALSE)

        if (!length(etastart)) {
            zetaff.Loglikfun = function(pp, y, x, w, extraargs) {
                sum(w * dzeta(x = y, p=pp, log = TRUE))
            }
            p.grid = seq(0.1, 3.0, len=19)
            pp.init = if (length( .init.p )) .init.p else
                      getMaxMin(p.grid, objfun=zetaff.Loglikfun, y = y,  x = x, w = w)
            pp.init = rep(pp.init, length=length(y))
            if ( .link == "loglog") pp.init[pp.init <= 1] = 1.2
            etastart = theta2eta(pp.init, .link, earg = .earg)
        }
    }), list( .link = link, .earg = earg, .init.p=init.p ))),
    inverse = eval(substitute(function(eta, extra = NULL) {
        ans <- pp <- eta2theta(eta, .link, earg = .earg)
        ans[pp > 1] <- zeta(pp[pp > 1]) / zeta(pp[pp > 1] + 1)
        ans[pp <= 1] <- NA
        ans
    }, list( .link = link, .earg = earg ))),
    last = eval(substitute(expression({
        misc$link <-    c(pp = .link)
        misc$earg <- list(pp = .earg)
    }), list( .link = link, .earg = earg ))),
    loglikelihood = eval(substitute(
        function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
        pp = eta2theta(eta, .link, earg = .earg)
        if (residuals) stop("loglikelihood residuals not ",
                            "implemented yet") else {
            sum(w * dzeta(x = y, p=pp, log = TRUE))
        }
    }, list( .link = link, .earg = earg ))),
    vfamily = c("zetaff"),
    deriv = eval(substitute(expression({
        pp = eta2theta(eta, .link, earg = .earg)
        fred1 = zeta(pp+1)
        fred2 = zeta(pp+1, deriv=1)
        dl.dpp = -log(y) - fred2 / fred1
        dpp.deta = dtheta.deta(pp, .link, earg = .earg)
        c(w) * dl.dpp * dpp.deta
    }), list( .link = link, .earg = earg ))),
    weight = expression({
        ed2l.dpp2 = zeta(pp+1, deriv=2) / fred1 - (fred2/fred1)^2
        wz = c(w) * dpp.deta^2 * ed2l.dpp2
        wz
    }))
}



gharmonic = function(n, s = 1, lognexponent=0) {

    if (!is.Numeric(n, integ = TRUE, posit = TRUE))
        stop("bad input for argument 'n'")
    if (!is.Numeric(lognexponent, allow = 1))
        stop("bad input for argument 'lognexponent'")
    if (length(n) == 1 && length(s) == 1) {
        if (lognexponent != 0) sum(log(1:n)^lognexponent * (1:n)^(-s)) else
            sum((1:n)^(-s))
    } else {
        LEN = max(length(n), length(s))
        n = rep(n, len = LEN)
        ans = s = rep(s, len = LEN)
        if (lognexponent != 0) {
            for(ii in 1:LEN)
                ans[ii] = sum(log(1:n[ii])^lognexponent * (1:n[ii])^(-s[ii]))
        } else
            for(ii in 1:LEN)
                ans[ii] = sum((1:n[ii])^(-s[ii]))
        ans
    }
}

dzipf = function(x, N, s, log = FALSE)
{
    if (!is.logical(log.arg <- log))
        stop("bad input for argument 'log'")
    rm(log)

    if (!is.Numeric(x))
        stop("bad input for argument 'x'")
    if (!is.Numeric(N, integ = TRUE, posit = TRUE))
        stop("bad input for argument 'N'")
    if (!is.Numeric(s, posit = TRUE))
        stop("bad input for argument 's'")
    nn = max(length(x), length(N), length(s))
    x = rep(x, len = nn); N = rep(N, len = nn); s = rep(s, len = nn);
    ox = !is.finite(x)
    zero = ox | round(x) != x | x < 1 | x > N
    ans = (if (log.arg) log(0) else 0) * x
    if (any(!zero))
        if (log.arg) {
            ans[!zero] = (-s[!zero]) * log(x[!zero]) -
                         log(gharmonic(N[!zero], s[!zero]))
        } else {
            ans[!zero] = x[!zero]^(-s[!zero]) / gharmonic(N[!zero], s[!zero])
        }
    ans
}



pzipf = function(q, N, s) {
    if (!is.Numeric(q))
        stop("bad input for argument 'q'")
    if (!is.Numeric(N, integ = TRUE, posit = TRUE))
        stop("bad input for argument 'N'")
    if (!is.Numeric(s, posit = TRUE))
        stop("bad input for argument 's'")

    nn = max(length(q), length(N), length(s))
    q = rep(q, len = nn); N = rep(N, len = nn); s = rep(s, len = nn);
    oq = !is.finite(q)
    zeroOR1 = oq | q < 1 | q >= N
    floorq = floor(q)
    ans = 0 * floorq
    ans[oq | q >= N] = 1
    if (any(!zeroOR1))
        ans[!zeroOR1] = gharmonic(floorq[!zeroOR1], s[!zeroOR1]) /
                        gharmonic(N[!zeroOR1], s[!zeroOR1])
    ans
}


 zipf = function(N = NULL, link = "loge", earg = list(), init.s = NULL)
{
    if (length(N) &&
      (!is.Numeric(N, positi = TRUE, integ = TRUE, allow = 1) || N <= 1))
        stop("bad input for argument 'N'")
    enteredN = length(N)
    if (length(init.s) && !is.Numeric(init.s, positi = TRUE))
        stop("argument 'init.s' must be > 0")

    if (mode(link) != "character" && mode(link) != "name")
        link = as.character(substitute(link))
    if (!is.list(earg)) earg = list()

    new("vglmff",
    blurb = c("Zipf distribution f(y;s) = y^(-s) / sum((1:N)^(-s)),",
              " s>0, y = 1,2,...,N", ifelse(enteredN, paste(" = ",N,sep = ""), ""),
              "\n\n",
              "Link:    ",
              namesof("s", link, earg = earg),
              "\n\n",
              "Mean:    gharmonic(N,s-1) / gharmonic(N,s)"),
    initialize = eval(substitute(expression({
        if (ncol(cbind(y)) != 1)
            stop("response must be a vector or a one-column matrix")
        y = as.numeric(y)
        if (any(y != round(y )))
            stop("y must be integer-valued")
        predictors.names = namesof("s", .link, earg = .earg, tag = FALSE)
        NN = .N
        if (!is.Numeric(NN, allow = 1, posit = TRUE, integ = TRUE))
            NN = max(y)
        if (max(y) > NN)
            stop("maximum of the response is greater than argument 'N'")
        if (any(y < 1))
            stop("all response values must be in 1,2,3,...,N( = ", NN,")")
        extra$N = NN
        if (!length(etastart)) {
            llfun = function(ss, y, N, w) {
                sum(w * dzipf(x = y, N=extra$N, s=ss, log = TRUE))
            }
            ss.init = if (length( .init.s )) .init.s else
                getInitVals(gvals=seq(0.1, 3.0, len=19), llfun=llfun,
                            y = y, N=extra$N, w = w)
            ss.init = rep(ss.init, length=length(y))
            if ( .link == "loglog") ss.init[ss.init <= 1] = 1.2
            etastart = theta2eta(ss.init, .link, earg = .earg)
        }
    }), list( .link = link, .earg = earg, .init.s = init.s, .N = N ))),
    inverse = eval(substitute(function(eta, extra = NULL) {
        ss = eta2theta(eta, .link, earg = .earg)
        gharmonic(extra$N, s=ss - 1) / gharmonic(extra$N, s=ss)
    }, list( .link = link, .earg = earg ))),
    last = eval(substitute(expression({
        misc$expected = FALSE
        misc$link =    c(s = .link)
        misc$earg = list(s = .earg )
        misc$N = extra$N
    }), list( .link = link, .earg = earg ))),
    loglikelihood = eval(substitute(
        function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
        ss = eta2theta(eta, .link, earg = .earg)
        if (residuals) stop("loglikelihood residuals not ",
                            "implemented yet") else {
            sum(w * dzipf(x = y, N=extra$N, s=ss, log = TRUE))
        }
    }, list( .link = link, .earg = earg ))),
    vfamily = c("zipf"),
    deriv = eval(substitute(expression({
        ss = eta2theta(eta, .link, earg = .earg)
        fred1 = gharmonic(extra$N, ss)
        fred2 = gharmonic(extra$N, ss, lognexp=1)
        dl.dss = -log(y) + fred2 / fred1
        dss.deta = dtheta.deta(ss, .link, earg = .earg)
        d2ss.deta2 = d2theta.deta2(ss, .link, earg = .earg)
        c(w) * dl.dss * dss.deta
    }), list( .link = link, .earg = earg ))),
    weight = expression({
        d2l.dss = gharmonic(extra$N, ss, lognexp=2) / fred1 - (fred2/fred1)^2
        wz = c(w) * (dss.deta^2 * d2l.dss - d2ss.deta2 * dl.dss)
        wz
    }))
}



cauchy.control <- function(save.weight = TRUE, ...)
{
    list(save.weight = save.weight)
}

 cauchy = function(llocation = "identity", lscale = "loge",
                  elocation = list(), escale = list(),
                  ilocation = NULL, iscale = NULL,
                  iprobs = seq(0.2, 0.8, by=0.2),
                  imethod = 1, nsimEIM = NULL, zero = 2)
{
    if (mode(llocation) != "character" && mode(llocation) != "name")
        llocation = as.character(substitute(llocation))
    if (mode(lscale) != "character" && mode(lscale) != "name")
        lscale = as.character(substitute(lscale))
    if (!is.Numeric(imethod, allow = 1, integ = TRUE, posit = TRUE) ||
       imethod > 3)
        stop("argument 'imethod' must be 1 or 2 or 3")
    if (!is.list(elocation)) elocation = list()
    if (!is.list(escale)) escale = list()
    if (length(zero) && !is.Numeric(zero, integer = TRUE, posit = TRUE))
        stop("bad input for argument 'zero'")
    if (length(nsimEIM) &&
       (!is.Numeric(nsimEIM, allow = 1, integ = TRUE) || nsimEIM <= 50))
        stop("argument 'nsimEIM' should be an integer greater than 50")
    if (length(iscale) && !is.Numeric(iscale, posit = TRUE))
        stop("bad input for argument 'iscale'")
    if (!is.Numeric(iprobs, posit = TRUE) || max(iprobs) >= 1)
        stop("bad input for argument 'iprobs'")

    new("vglmff",
    blurb = c("Two parameter Cauchy distribution (location & scale unknown)\n\n",
              "Link:    ",
              namesof("location", llocation, earg = elocation), "\n",
              namesof("scale",    lscale,    earg = escale), "\n\n",
              "Mean:     NA\n",
              "Variance: NA"),
    constraints = eval(substitute(expression({
        constraints = cm.zero.vgam(constraints, x, .zero, M)
    }), list( .zero = zero ))),
    initialize = eval(substitute(expression({
        predictors.names = c(
          namesof("location", .llocation, earg = .elocation, tag = FALSE),
          namesof("scale",    .lscale,    earg = .escale,    tag = FALSE))
        if (ncol(cbind(y)) != 1)
            stop("response must be a vector or a one-column matrix")

        if (!length(etastart)) {
            loc.init = if (length( .ilocation)) .ilocation else {
                if ( .imethod == 2) median(rep(y,w)) else 
                if ( .imethod == 3) y else {
                    cauchy2.Loglikfun = function(loc, y, x, w, extraargs) {
                         iprobs = .iprobs
                         qy = quantile(rep(y,w), probs=iprobs)
                         ztry = tan(pi*(iprobs-0.5))
                         btry = (qy - loc) / ztry
                         scal = median(btry, na.rm = TRUE)
                         if (scal <= 0) scal = 0.1
                         sum(w * dcauchy(x = y, loc=loc, scale=scal, log = TRUE))
                     }
                     loc.grid = c(quantile(y, probs=seq(0.1, 0.9, by=0.05)))
                     try.this = getMaxMin(loc.grid, objfun=cauchy2.Loglikfun,
                                          y = y,  x = x, w = w)
                    try.this = rep(c(try.this), len = n)
                    try.this
                }
            }
            loc.init = rep(c(loc.init), len = n)


            sca.init = if (length( .iscale)) .iscale else {
                iprobs = .iprobs
                qy = quantile(rep(y,w), probs=iprobs)
                ztry = tan(pi*(iprobs-0.5))
                btry = (qy - loc.init[1]) / ztry
                sca.init = median(btry, na.rm = TRUE)
                if (sca.init <= 0) sca.init = 0.01
                sca.init
            }

            sca.init = rep(c(sca.init), len = n)
            if ( .llocation == "loge") loc.init = abs(loc.init)+0.01
            etastart = cbind(theta2eta(loc.init, .llocation, earg = .elocation),
                             theta2eta(sca.init, .lscale,    earg = .escale))
        }
    }), list( .ilocation = ilocation, .elocation = elocation, .llocation = llocation,
              .iscale = iscale, .escale = escale, .lscale = lscale,
              .iprobs=iprobs, .imethod = imethod ))),
    inverse = eval(substitute(function(eta, extra = NULL) {
        eta2theta(eta[,1], .llocation, earg = .elocation)
    }, list( .llocation = llocation,
             .elocation = elocation ))),
    last = eval(substitute(expression({
        misc$expected = TRUE
        misc$link =    c("location" = .llocation, "scale" =.lscale)
        misc$earg = list("location" = .elocation, "scale" = .escale)
        misc$imethod = .imethod
    }), list( .escale = escale, .elocation = elocation,
              .imethod = imethod,
              .llocation = llocation, .lscale = lscale ))),
    loglikelihood = eval(substitute(
        function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
        location = eta2theta(eta[,1], .llocation, earg = .elocation)
        myscale  = eta2theta(eta[,2], .lscale,    earg = .escale)
        if (residuals) stop("loglikelihood residuals not ",
                            "implemented yet") else {
            sum(w * dcauchy(x = y, loc=location, sc=myscale, log = TRUE))
        }
    }, list( .escale = escale, .lscale = lscale,
             .elocation = elocation, .llocation = llocation ))),
    vfamily = c("cauchy"),
    deriv = eval(substitute(expression({
        location = eta2theta(eta[,1], .llocation, earg = .elocation)
        myscale = eta2theta(eta[,2], .lscale, earg = .escale)
        dlocation.deta = dtheta.deta(location, .llocation, earg = .elocation)
        dscale.deta = dtheta.deta(myscale, .lscale, earg = .escale)
        Z = (y-location) / myscale
        dl.dlocation = 2 * Z / ((1 + Z^2) * myscale)
        dl.dscale = (Z^2 - 1) / ((1 + Z^2) * myscale)
        c(w) * cbind(dl.dlocation * dlocation.deta,
                     dl.dscale * dscale.deta)
    }), list( .escale = escale, .lscale = lscale,
              .elocation = elocation, .llocation = llocation ))),
    weight = eval(substitute(expression({
        run.varcov = 0
        ind1 = iam(NA, NA, M = M, both = TRUE, diag = TRUE)
        dthetas.detas = cbind(dlocation.deta, dscale.deta)
        if (length( .nsimEIM )) {
            for(ii in 1:( .nsimEIM )) {
                ysim = rcauchy(n, loc=location, scale=myscale)
                Z = (ysim-location) / myscale
                dl.dlocation = 2 * Z / ((1 + Z^2) * myscale)
                dl.dscale = (Z^2 - 1) / ((1 + Z^2) * myscale)
                rm(ysim)
                temp3 = matrix(c(dl.dlocation, dl.dscale), n, 2)
                run.varcov = ((ii-1) * run.varcov +
                           temp3[,ind1$row.index]*temp3[,ind1$col.index]) / ii
            }
            wz = if (intercept.only)
                matrix(colMeans(run.varcov),
                       n, ncol(run.varcov), byrow = TRUE) else run.varcov

            wz = wz * dthetas.detas[,ind1$row] * dthetas.detas[,ind1$col]
            wz = c(w) * matrix(wz, n, dimm(M))
        } else {
            wz = cbind(matrix(0.5 / myscale^2,n,2), matrix(0,n,1)) *
                 dthetas.detas[,ind1$row] * dthetas.detas[,ind1$col]
            wz = c(w) * wz[,1:M]  # diagonal wz
        }

        wz
    }), list( .escale = escale, .lscale = lscale, .nsimEIM = nsimEIM,
              .elocation = elocation, .llocation = llocation ))))
}







 cauchy1 = function(scale.arg = 1, llocation = "identity",
                    elocation = list(),
                    ilocation = NULL, imethod = 1)
{
    if (mode(llocation) != "character" && mode(llocation) != "name")
        llocation = as.character(substitute(llocation))
    if (!is.Numeric(scale.arg, posit = TRUE)) stop("bad input for 'scale.arg'")
    if (!is.Numeric(imethod, allow = 1, integ = TRUE, posit = TRUE) ||
       imethod > 3)
        stop("argument 'imethod' must be 1 or 2 or 3")
    if (!is.list(elocation)) elocation = list()

    new("vglmff",
    blurb = c("One-parameter Cauchy distribution ",
              "(location unknown, scale known)\n\n",
              "Link:    ",
              namesof("location", llocation, earg = elocation), "\n\n",
              "Mean:     NA\n",
              "Variance: NA"),
    initialize = eval(substitute(expression({
        predictors.names = namesof("location", .llocation,
                                   earg = .elocation, tag = FALSE)
        if (ncol(cbind(y)) != 1)
            stop("response must be a vector or a one-column matrix")

        if (!length(etastart)) {
            loc.init = if (length( .ilocation)) .ilocation else {
                if ( .imethod == 2) median(rep(y,w)) else 
                if ( .imethod == 3) y else {
                    cauchy1.Loglikfun = function(loc, y, x, w, extraargs) {
                         scal = extraargs
                         sum(w * dcauchy(x = y, loc=loc, scale=scal, log = TRUE))
                     }
                     loc.grid = quantile(y, probs=seq(0.1, 0.9, by=0.05))
                     try.this = getMaxMin(loc.grid, objfun=cauchy1.Loglikfun,
                                          y = y,  x = x, w = w, extraargs= .scale.arg)
                    try.this = rep(try.this, len = n)
                    try.this
                }
            }
            loc.init = rep(loc.init, len = n)
            if ( .llocation == "loge") loc.init = abs(loc.init)+0.01
            etastart = theta2eta(loc.init, .llocation, earg = .elocation)
        }
    }), list( .scale.arg=scale.arg, .ilocation = ilocation,
              .elocation = elocation, .llocation = llocation,
              .imethod = imethod ))),
    inverse = eval(substitute(function(eta, extra = NULL) {
        eta2theta(eta, .llocation, earg = .elocation)
    }, list( .llocation = llocation,
             .elocation = elocation ))),
    last = eval(substitute(expression({
        misc$expected = TRUE
        misc$link =    c("location" = .llocation)
        misc$earg = list("location" = .elocation )
        misc$scale.arg = .scale.arg 
    }), list( .scale.arg=scale.arg, .elocation = elocation,
             .llocation = llocation ))),
    loglikelihood = eval(substitute(
        function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
        location = eta2theta(eta, .llocation, earg = .elocation)
        if (residuals) stop("loglikelihood residuals not ",
                            "implemented yet") else {
            sum(w * dcauchy(x = y, loc=location, scale= .scale.arg, log = TRUE))
        }
    }, list( .scale.arg=scale.arg, .elocation = elocation,
             .llocation = llocation ))),
    vfamily = c("cauchy1"),
    deriv = eval(substitute(expression({
        location = eta2theta(eta, .llocation, earg = .elocation)
        temp = (y-location)/.scale.arg
        dl.dlocation = 2 * temp / ((1 + temp^2) * .scale.arg)
        dlocation.deta = dtheta.deta(location, .llocation, earg = .elocation)
        c(w) * dl.dlocation * dlocation.deta
    }), list( .scale.arg=scale.arg, .elocation = elocation,
              .llocation = llocation ))),
    weight = eval(substitute(expression({
        wz = c(w) * dlocation.deta^2 / ( .scale.arg^2 * 2)
        wz
    }), list( .scale.arg=scale.arg, .elocation = elocation,
              .llocation = llocation ))))
}






 logistic1 = function(llocation = "identity",
                     elocation = list(),
                     scale.arg = 1, imethod = 1)
{
    if (mode(llocation) != "character" && mode(llocation) != "name")
        llocation = as.character(substitute(llocation))
    if (!is.Numeric(scale.arg, allow = 1, posit = TRUE))
        stop("'scale.arg' must be a single positive number")
    if (!is.Numeric(imethod, allow = 1, integ = TRUE, posit = TRUE) ||
       imethod > 2)
        stop("argument 'imethod' must be 1 or 2")
    if (!is.list(elocation)) elocation = list()

    new("vglmff",
    blurb = c("One-parameter logistic distribution ",
            "(location unknown, scale known)\n\n",
            "Link:    ",
            namesof("location", llocation, earg = elocation), "\n\n",
            "Mean:     location", "\n",
            "Variance: (pi*scale)^2 / 3"),
    initialize = eval(substitute(expression({
        predictors.names = namesof("location", .llocation, 
                                   earg = .elocation, tag = FALSE)
        if (!length(etastart)) {
            location.init = if ( .imethod == 1) y else median(rep(y, w))
            location.init = rep(location.init, len = n)
            if ( .llocation == "loge") location.init = abs(location.init) + 0.001
            etastart = theta2eta(location.init, .llocation, earg = .elocation)
        }
    }), list( .imethod = imethod, .llocation = llocation,
              .elocation = elocation ))),
    inverse = eval(substitute(function(eta, extra = NULL) {
        eta2theta(eta, .llocation, earg = .elocation)
    }, list( .llocation = llocation,
             .elocation = elocation ))),
    last = eval(substitute(expression({
        misc$expected = TRUE
        misc$link =    c(location = .llocation)
        misc$earg = list(location = .elocation )
        misc$scale.arg = .scale.arg 
    }), list( .llocation = llocation, 
              .elocation = elocation, .scale.arg=scale.arg ))),
    loglikelihood = eval(substitute(
        function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
        location = eta2theta(eta, .llocation, earg = .elocation)
        zedd = (y-location)/.scale.arg
        if (residuals) stop("loglikelihood residuals not ",
                            "implemented yet") else {
            sum(w * dlogis(x = y, location = location,
                           scale = .scale.arg, log = TRUE))
        }
    }, list( .llocation = llocation,
             .elocation = elocation, .scale.arg=scale.arg ))),
    vfamily = c("logistic1"),
    deriv = eval(substitute(expression({
        location = eta2theta(eta, .llocation, earg = .elocation)
        ezedd = exp(-(y-location)/.scale.arg)
        dl.dlocation = (1 - ezedd) / ((1 + ezedd) * .scale.arg)
        dlocation.deta = dtheta.deta(location, .llocation, earg = .elocation)
        c(w) * dl.dlocation * dlocation.deta
    }), list( .llocation = llocation,
              .elocation = elocation, .scale.arg=scale.arg ))),
    weight = eval(substitute(expression({
        wz = c(w) * dlocation.deta^2 / ( .scale.arg^2 * 3) 
        wz
    }), list( .scale.arg=scale.arg ))))
}




 erlang = function(shape.arg, link = "loge", earg = list(), imethod = 1)
{

    if (!is.Numeric(shape.arg, allow = 1, integer = TRUE, positi = TRUE))
        stop("'shape' must be a positive integer")
    if (!is.Numeric(imethod, allow = 1, integer = TRUE, positi = TRUE) ||
       imethod > 2)
        stop("argument 'imethod' must be 1 or 2")

    if (mode(link) != "character" && mode(link) != "name")
        link = as.character(substitute(link))
    if (!is.list(earg)) earg = list()

    new("vglmff",
    blurb = c("Erlang distribution\n\n",
            "Link:    ", namesof("scale", link, earg = earg), "\n", "\n",
            "Mean:     shape * scale", "\n",
            "Variance: shape * scale^2"),
    initialize = eval(substitute(expression({
        if (ncol(y <- as.matrix(y)) > 1)
            stop("erlang cannot handle matrix responses yet")
        if (any(y < 0))
            stop("all y values must be >= 0")

        predictors.names =
          namesof("scale", .link, earg = .earg, tag = FALSE)

        if (!length(etastart)) {
            if ( .imethod == 1) 
                sc.init = y / .shape.arg
            if ( .imethod==2) {
                sc.init = median(y) / .shape.arg
                sc.init = rep(sc.init, length = n) 
            }
            etastart = theta2eta(sc.init, .link, earg = .earg)
        }
    }), list( .link = link, .earg = earg,
              .shape.arg=shape.arg, .imethod = imethod ))),
    inverse = eval(substitute(function(eta, extra = NULL) {
        sc = eta2theta(eta, .link, earg = .earg)
        .shape.arg * sc 
    }, list( .link = link, .earg = earg, .shape.arg=shape.arg ))),
    last = eval(substitute(expression({
        misc$expected = TRUE
        misc$link =    c(scale = .link)
        misc$earg = list(scale = .earg )
        misc$shape.arg = .shape.arg 
    }), list( .link = link, .earg = earg, .shape.arg=shape.arg ))),
    loglikelihood = eval(substitute(
        function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
        sc = eta2theta(eta, .link, earg = .earg)
        if (residuals) stop("loglikelihood residuals not ",
                            "implemented yet") else {
            sum(w * (( .shape.arg - 1) * log(y) - y / sc - .shape.arg * log(sc) -
                     lgamma( .shape.arg )))
        }
    }, list( .link = link, .earg = earg, .shape.arg=shape.arg ))),
    vfamily = c("erlang"),
    deriv = eval(substitute(expression({
        sc = eta2theta(eta, .link, earg = .earg)
        dl.dsc = (y / sc - .shape.arg) / sc
        dsc.deta = dtheta.deta(sc, .link, earg = .earg)
        c(w) * dl.dsc * dsc.deta
    }), list( .link = link, .earg = earg, .shape.arg=shape.arg ))),
    weight = eval(substitute(expression({
        ed2l.dsc2 = .shape.arg / sc^2
        wz = c(w) * dsc.deta^2 * ed2l.dsc2
        wz
    }), list( .earg = earg, .shape.arg=shape.arg ))))
}





dbort = function(x, Qsize = 1, a=0.5, log = FALSE) {
    if (!is.logical(log.arg <- log))
        stop("bad input for argument 'log'")
    rm(log)

    if (!is.Numeric(x)) stop("bad input for argument 'x'")
    if (!is.Numeric(Qsize, allow = 1, integ = TRUE, posit = TRUE))
        stop("bad input for argument 'Qsize'")
    if (!is.Numeric(a, posit = TRUE) || max(a) >= 1)
        stop("bad input for argument 'a'")
    N = max(length(x), length(Qsize), length(a))
    x = rep(x, len=N); Qsize = rep(Qsize, len=N); a = rep(a, len=N);

    xok = (x >= Qsize) & (x == round(x)) & (a > 0) & (a < 1)
    ans = rep(if (log.arg) log(0) else 0, len=N) # loglikelihood
    ans[xok] = lgamma(1 + Qsize[xok]) - lgamma(x[xok] + 1 - Qsize[xok]) +
               (x[xok] - 1 - Qsize[xok]) * log(x[xok]) +
               (x[xok] - Qsize[xok]) * log(a[xok]) - a[xok] * x[xok]
    if (!log.arg) {
        ans[xok] = exp(ans[xok])
    }
    ans
}


rbort = function(n, Qsize = 1, a=0.5) {
    if (!is.Numeric(n, integ = TRUE, posit = TRUE, allow = 1))
        stop("bad input for argument 'n'")
    if (!is.Numeric(Qsize, allow = 1, integ = TRUE, posit = TRUE))
        stop("bad input for argument 'Qsize'")
    if (!is.Numeric(a, posit = TRUE) || max(a) >= 1)
        stop("bad input for argument 'a'")
    N = n
    qsize = rep(Qsize, len=N); a = rep(a, len=N)
    totqsize = qsize
    fini = (qsize < 1)
    while(any(!fini)) {
        additions = rpois(sum(!fini), a[!fini])
        qsize[!fini] = qsize[!fini] + additions
        totqsize[!fini] = totqsize[!fini] + additions
        qsize = qsize - 1
        fini = fini | (qsize < 1)
    }
    totqsize
}


 borel.tanner = function(Qsize = 1, link = "logit", earg = list(), imethod = 1)
{
    if (!is.Numeric(Qsize, allow = 1, integ = TRUE, posit = TRUE))
        stop("bad input for argument 'Qsize'")
    if (mode(link) != "character" && mode(link) != "name")
        link = as.character(substitute(link))
    if (!is.list(earg)) earg = list()
    if (!is.Numeric(imethod, allow = 1, integ = TRUE, posit = TRUE) ||
       imethod > 4)
        stop("argument 'imethod' must be 1 or 2, 3 or 4")

    new("vglmff",
    blurb = c("Borel-Tanner distribution\n\n",
            "Link:    ",
            namesof("a", link, earg = earg), "\n\n",
            "Mean:     Qsize/(1-a)",
            "\n",
            "Variance: Qsize*a / (1-a)^3"),
    initialize = eval(substitute(expression({
        if (ncol(cbind(y)) != 1)
            stop("response must be a vector or a one-column matrix")
        if (any(y < .Qsize))
            stop("all y values must be >= ", .Qsize)
        if (any(y != round(y)))
            warning("response should be integer-valued")

        predictors.names = namesof("a", .link, earg = .earg, tag = FALSE)

        if (!length(etastart)) {
            a.init = switch(as.character( .imethod ),
                "1" = 1 - .Qsize / (y+1/8),
                "2" = rep(1 - .Qsize / weighted.mean(y,w), len = n),
                "3" = rep(1 - .Qsize / median(y), len = n),
                "4" = rep(0.5, len = n))
            etastart = theta2eta(a.init, .link, earg = .earg)
        }
    }), list( .link = link, .earg = earg, .Qsize=Qsize,
              .imethod = imethod ))),
    inverse = eval(substitute(function(eta, extra = NULL) {
        a = eta2theta(eta, .link, earg = .earg)
        .Qsize / (1 - a)
    }, list( .link = link, .earg = earg, .Qsize=Qsize ))),
    last = eval(substitute(expression({
        misc$expected = TRUE
        misc$link =    c(a = .link)
        misc$earg = list(a = .earg )
        misc$Qsize = .Qsize 
    }), list( .link = link, .earg = earg, .Qsize=Qsize ))),
    loglikelihood = eval(substitute(
        function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
        aa = eta2theta(eta, .link, earg = .earg)
        if (residuals) stop("loglikelihood residuals not ",
                            "implemented yet") else {
            sum(w * dbort(x = y, Qsize= .Qsize, a=aa, log = TRUE))
        }
    }, list( .link = link, .earg = earg, .Qsize=Qsize ))),
    vfamily = c("borel.tanner"),
    deriv = eval(substitute(expression({
        a = eta2theta(eta, .link, earg = .earg)
        dl.da = (y- .Qsize)/a - y 
        da.deta = dtheta.deta(a, .link, earg = .earg)
        c(w) * dl.da * da.deta
    }), list( .link = link, .earg = earg, .Qsize=Qsize ))),
    weight = eval(substitute(expression({
        ed2l.da2 = .Qsize / (a*(1-a))
        wz = c(w) * da.deta^2 * ed2l.da2
        wz
    }), list( .Qsize=Qsize ))))
}



dfelix = function(x, a=0.25, log = FALSE) {
    if (!is.logical(log.arg <- log))
        stop("bad input for argument 'log'")
    rm(log)

    if (!is.Numeric(x)) stop("bad input for argument 'x'")
    if (!is.Numeric(a, posit = TRUE)) stop("bad input for argument 'a'")
    N = max(length(x), length(a))
    x = rep(x, len=N); a = rep(a, len=N);

    xok = (x %% 2 == 1) & (x == round(x)) & (x >= 1) & (a > 0) & (a < 0.5)
    ans = rep(if (log.arg) log(0) else 0, len=N) # loglikelihood
    ans[xok] = ((x[xok]-3)/2) * log(x[xok]) + ((x[xok]-1)/2) * log(a[xok]) -
               lgamma(x[xok]/2 + 0.5) - a[xok] * x[xok]
    if (!log.arg) {
        ans[xok] = exp(ans[xok])
    }
    ans
}



 felix = function(link = "elogit",
            earg=if (link == "elogit") list(min=0, max=0.5) else list(),
            imethod = 1)
{
    if (mode(link) != "character" && mode(link) != "name")
        link = as.character(substitute(link))
    if (!is.list(earg)) earg = list()
    if (!is.Numeric(imethod, allow = 1, integ = TRUE, posit = TRUE) ||
       imethod > 4)
        stop("argument 'imethod' must be 1 or 2, 3 or 4")

    new("vglmff",
    blurb = c("Felix distribution\n\n",
            "Link:    ",
            namesof("a", link, earg = earg), "\n\n",
            "Mean:     1/(1-2*a)"),
    initialize = eval(substitute(expression({
        if (ncol(cbind(y)) != 1)
            stop("response must be a vector or a one-column matrix")
        if (any(y < 1) || any((y+1)/2 != round((y+1)/2)))
            warning("response should be positive, odd and integer-valued")

        predictors.names = namesof("a", .link, earg = .earg, tag = FALSE)

        if (!length(etastart)) {
            wymean = weighted.mean(y,w)
            a.init = switch(as.character( .imethod ),
                "1" = (y-1+1/8) / (2*(y+1/8)+1/8),
                "2" = rep((wymean-1+1/8) / (2*(wymean+1/8)+1/8), len = n),
                "3" = rep((median(y)-1+1/8) / (2*(median(y)+1/8)+1/8), len = n),
                "4" = rep(0.25, len = n))
            etastart = theta2eta(a.init, .link, earg = .earg)
        }
    }), list( .link = link, .earg = earg,
              .imethod = imethod ))),
    inverse = eval(substitute(function(eta, extra = NULL) {
        a = eta2theta(eta, .link, earg = .earg)
        1 / (1 - 2*a)
    }, list( .link = link, .earg = earg ))),
    last = eval(substitute(expression({
        misc$expected = TRUE
        misc$link =    c(a = .link)
        misc$earg = list(a = .earg )
    }), list( .link = link, .earg = earg ))),
    loglikelihood = eval(substitute(
        function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
        aa = eta2theta(eta, .link, earg = .earg)
        if (residuals) stop("loglikelihood residuals not ",
                            "implemented yet") else {
           sum(w * dfelix(x = y, a=aa, log = TRUE))
       }
    }, list( .link = link, .earg = earg ))),
    vfamily = c("felix"),
    deriv = eval(substitute(expression({
        a = eta2theta(eta, .link, earg = .earg)
        dl.da = (y- 1)/(2*a) - y 
        da.deta = dtheta.deta(a, .link, earg = .earg)
        c(w) * dl.da * da.deta
    }), list( .link = link, .earg = earg ))),
    weight = eval(substitute(expression({
        ed2l.da2 = 1 / (a*(1-2*a))
        wz = c(w) * da.deta^2 * ed2l.da2
        wz
    }), list( .link = link ))))
}





 betaff = function(A=0, B = 1,
          lmu = if (A == 0 & B == 1) "logit" else "elogit", lphi = "loge",
          emu = if (lmu == "elogit") list(min=A, max=B) else list(),
          ephi = list(),
          imu = NULL, iphi = NULL, imethod = 1, zero = NULL)
{
    if (!is.Numeric(A, allow = 1) || !is.Numeric(B, allow = 1) || A >= B)
        stop("A must be < B, and both must be of length one")
    stdbeta = (A == 0 && B == 1)

    if (mode(lmu) != "character" && mode(lmu) != "name")
        lmu = as.character(substitute(lmu))
    if (mode(lphi) != "character" && mode(lphi) != "name")
        lphi = as.character(substitute(lphi))
    if (length(zero) && !is.Numeric(zero, integer = TRUE, posit = TRUE))
        stop("bad input for argument 'zero'")
    if (length(imu) && (!is.Numeric(imu, posit = TRUE) ||
       any(imu <= A) || any(imu >= B)))
        stop("bad input for argument 'imu'")
    if (length(iphi) && !is.Numeric(iphi, posit = TRUE))
        stop("bad input for argument 'iphi'")
    if (!is.Numeric(imethod, allow = 1, integ = TRUE, posit = TRUE) ||
       imethod > 2)
        stop("argument 'imethod' must be 1 or 2")

    if (!is.list(emu)) emu = list()
    if (!is.list(ephi)) ephi = list()

    new("vglmff",
    blurb = c("Beta distribution parameterized by mu and a precision parameter\n",
            if (stdbeta) paste("f(y) = y^(mu*phi-1) * (1-y)^((1-mu)*phi-1)",
            "/ beta(mu*phi,(1-mu)*phi), 0<y<1, 0<mu<1, phi>0\n\n") else
            paste("f(y) = (y-",A,")^(mu1*phi-1) * (",B,
            "-y)^(((1-mu1)*phi)-1) / \n(beta(mu1*phi,(1-mu1)*phi) * (",
            B, "-", A, ")^(phi-1)),\n",
            A," < y < ",B, ", ", A," < mu < ",B,
            ", mu = ", A, " + ", (B-A), " * mu1",
            ", phi > 0\n\n", sep = ""),
            "Links:    ",
            namesof("mu",  lmu,  earg = emu),  ", ",
            namesof("phi", lphi, earg = ephi)),
    constraints = eval(substitute(expression({
        constraints = cm.zero.vgam(constraints, x, .zero, M)
    }), list( .zero = zero ))),
    initialize = eval(substitute(expression({
        if (min(y) <= .A || max(y) >= .B)
            stop("data not within (A, B)")
        if (ncol(cbind(y)) != 1)
            stop("response must be a vector or a one-column matrix")
        predictors.names = c(namesof("mu",  .lmu,  .emu,  short = TRUE),
                             namesof("phi", .lphi, .ephi, short = TRUE))
        if (!length(etastart)) {
          mu.init = if (is.Numeric( .imu)) .imu else
                    {if ( .imethod == 1) weighted.mean(y,w) else
                     median(rep(y,w))}
          mu1.init = (mu.init - .A) / ( .B - .A)  # In (0,1)
          phi.init = if (is.Numeric( .iphi)) .iphi else
                       max(0.01, -1 + ( .B-.A)^2 * mu1.init*(1-mu1.init)/var(y))
          etastart = matrix(0, n, 2)
          etastart[,1] = theta2eta(mu.init, .lmu, earg = .emu )
          etastart[,2] = theta2eta(phi.init, .lphi, earg = .ephi )
      }
    }), list( .lmu = lmu, .lphi = lphi, .imu=imu, .iphi=iphi,
              .A = A, .B = B, .emu = emu, .ephi = ephi, .imethod = imethod ))),
    inverse = eval(substitute(function(eta, extra = NULL) {
       mu = eta2theta(eta[,1], .lmu, .emu )
       mu
    }, list( .lmu = lmu, .emu = emu, .A = A, .B = B))),
    last = eval(substitute(expression({
        misc$link =    c(mu = .lmu, phi = .lphi)
        misc$earg = list(mu = .emu, phi = .ephi)
        misc$limits = c( .A, .B)
        misc$stdbeta = .stdbeta
    }), list( .lmu = lmu, .lphi = lphi, .A = A, .B = B, .emu = emu, .ephi = ephi,
              .stdbeta = stdbeta ))),
    loglikelihood = eval(substitute(
        function(mu, y, w, residuals = FALSE, eta, extra = NULL){
        mu = eta2theta(eta[,1], .lmu, .emu )
        m1u = if ( .stdbeta ) mu else (mu - .A) / ( .B - .A)
        phi = eta2theta(eta[,2], .lphi, .ephi )
        if (residuals) stop("loglikelihood residuals not ",
                            "implemented yet") else {
            shape1 = phi * m1u
            shape2 = (1 - m1u) * phi
            zedd = (y - .A) / ( .B - .A)
            sum(w * (dbeta(x=zedd, shape1=shape1, shape2=shape2, log = TRUE) -
                     log( abs( .B - .A ))))
        }
    }, list( .lmu = lmu, .lphi = lphi, .A = A, .B = B, .emu = emu, .ephi = ephi,
             .stdbeta = stdbeta ))),
    vfamily = "betaff",
    deriv = eval(substitute(expression({
        mu = eta2theta(eta[,1], .lmu, .emu )
        phi = eta2theta(eta[,2], .lphi, .ephi )
        m1u = if ( .stdbeta ) mu else (mu - .A) / ( .B - .A)
        dmu.deta = dtheta.deta(mu, .lmu, .emu )
        dmu1.dmu = 1 / ( .B - .A)
        dphi.deta = dtheta.deta(phi, .lphi, .ephi )
        temp1 = m1u*phi
        temp2 = (1-m1u)*phi
        if ( .stdbeta ) {
            dl.dmu1 = phi*(digamma(temp2) - digamma(temp1) + log(y) - log1p(-y))
            dl.dphi = digamma(phi) - mu*digamma(temp1) - (1-mu)*digamma(temp2) +
                mu*log(y) + (1-mu)*log1p(-y)
        } else {
            dl.dmu1 = phi*(digamma(temp2) - digamma(temp1) +
                           log(y-.A) - log( .B-y))
            dl.dphi = digamma(phi) - m1u*digamma(temp1) -
                      (1-m1u)*digamma(temp2) +
                      m1u*log(y-.A) + (1-m1u)*log( .B-y) - log( .B -.A)
        }
        c(w) * cbind(dl.dmu1 * dmu1.dmu * dmu.deta,
                     dl.dphi * dphi.deta)
    }), list( .lmu = lmu, .lphi = lphi,
              .emu = emu, .ephi = ephi,
              .A = A, .B = B,
              .stdbeta = stdbeta ))),
    weight = eval(substitute(expression({
        d2l.dmu12 = phi^2 * (trigamma(temp1) + trigamma(temp2))
        d2l.dphi2 = -trigamma(phi) + trigamma(temp1) * m1u^2 +
            trigamma(temp2) * (1-m1u)^2
        d2l.dmu1phi = temp1*trigamma(temp1) - temp2*trigamma(temp2)
        wz = matrix(as.numeric(NA), n, dimm(M))
        wz[,iam(1,1,M)] = d2l.dmu12 * dmu1.dmu^2 * dmu.deta^2
        wz[,iam(2,2,M)] = d2l.dphi2 * dphi.deta^2
        wz[,iam(1,2,M)] = d2l.dmu1phi * dmu1.dmu * dmu.deta * dphi.deta
        c(w) * wz
    }), list( .A = A, .B = B ))))
}





 beta.ab = function(lshape1 = "loge", lshape2 = "loge",
                    eshape1 = list(), eshape2 = list(),
                    i1 = NULL, i2 = NULL, trim = 0.05,
                    A = 0, B = 1, parallel = FALSE, zero = NULL)
{
    if (mode(lshape1) != "character" && mode(lshape1) != "name")
        lshape1 = as.character(substitute(lshape1))
    if (mode(lshape2) != "character" && mode(lshape2) != "name")
        lshape2 = as.character(substitute(lshape2))
    if (length(zero) && !is.Numeric(zero, integer = TRUE, posit = TRUE))
        stop("bad input for argument 'zero'")
    if (length( i1 ) && !is.Numeric( i1, posit = TRUE))
        stop("bad input for argument 'i1'")
    if (length( i2 ) && !is.Numeric( i2, posit = TRUE))
        stop("bad input for argument 'i2'")

    if (!is.Numeric(A, allow = 1) || !is.Numeric(B, allow = 1) || A >= B)
        stop("A must be < B, and both must be of length one")
    stdbeta = (A == 0 && B == 1)  # stdbeta == T iff standard beta distribution
    if (!is.list(eshape1)) eshape1 = list()
    if (!is.list(eshape2)) eshape2 = list()

    new("vglmff",
    blurb = c("Two-parameter Beta distribution ",
              "(shape parameters parameterization)\n",
              if (stdbeta)
              paste("y^(shape1-1) * (1-y)^(shape2-1) / B(shape1,shape2),",
              "0 <= y <= 1, shape1>0, shape2>0\n\n") else
              paste("(y-",A,")^(shape1-1) * (",B,
              "-y)^(shape2-1) / [B(shape1,shape2) * (",
              B, "-", A, ")^(shape1+shape2-1)], ",
               A," <= y <= ",B," shape1>0, shape2>0\n\n", sep = ""),
              "Links:    ",
              namesof("shape1", lshape1, earg = eshape1),  ", ",
              namesof("shape2", lshape2, earg = eshape2)),
    constraints = eval(substitute(expression({
        constraints = cm.vgam(matrix(1, M, 1), x, .parallel, constraints, int= TRUE)
        constraints = cm.zero.vgam(constraints, x, .zero, M)
    }), list( .parallel = parallel, .zero = zero ))),
    initialize = eval(substitute(expression({
        if (min(y) <= .A || max(y) >= .B)
            stop("data not within (A, B)")
        if (ncol(cbind(y)) != 1)
            stop("response must be a vector or a one-column matrix")
        predictors.names =
            c(namesof("shape1", .lshape1, earg = .eshape1, short = TRUE),
              namesof("shape2", .lshape2, earg = .eshape2, short = TRUE))

        if (!length(etastart)) {
            mu1d = mean(y, trim = .trim)
            uu = (mu1d - .A) / ( .B - .A) 
            DD = ( .B - .A)^2 
            pinit = max(0.01, uu^2 * (1 - uu) * DD / var(y) - uu)
            qinit = max(0.01, pinit * (1 - uu) / uu)
            etastart = matrix(0, n, 2)
            etastart[,1] = theta2eta( pinit, .lshape1, earg = .eshape1 )
            etastart[,2] = theta2eta( qinit, .lshape2, earg = .eshape2 )
        }
        if (is.Numeric( .i1 ))
            etastart[,1] = theta2eta( .i1, .lshape1, earg = .eshape1 )
        if (is.Numeric( .i2 ))
            etastart[,2] = theta2eta( .i2, .lshape2, earg = .eshape2 )
    }), list( .lshape1 = lshape1, .lshape2 = lshape2,
              .i1 = i1, .i2 = i2, .trim = trim, .A = A, .B = B,
              .eshape1 = eshape1, .eshape2 = eshape2 ))),
    inverse = eval(substitute(function(eta, extra = NULL) {
        shapes = cbind(eta2theta(eta[,1], .lshape1, earg = .eshape1 ),
                       eta2theta(eta[,2], .lshape2, earg = .eshape2 ))
        .A + ( .B-.A) * shapes[,1] / (shapes[,1] + shapes[,2])
    }, list( .lshape1 = lshape1, .lshape2 = lshape2, .A = A, .B = B, 
             .eshape1 = eshape1, .eshape2 = eshape2 ))),
    last = eval(substitute(expression({
        misc$link =    c(shape1 = .lshape1, shape2 = .lshape2)
        misc$earg = list(shape1 = .eshape1, shape2 = .eshape2)
        misc$limits = c( .A, .B)
    }), list( .lshape1 = lshape1, .lshape2 = lshape2, .A = A, .B = B, 
              .eshape1 = eshape1, .eshape2 = eshape2 ))),
    loglikelihood = eval(substitute(
        function(mu, y, w, residuals = FALSE, eta, extra = NULL){
        shapes = cbind(eta2theta(eta[,1], .lshape1, earg = .eshape1 ),
                       eta2theta(eta[,2], .lshape2, earg = .eshape2 ))
        if (residuals) stop("loglikelihood residuals not ",
                            "implemented yet") else {
            zedd = (y - .A) / ( .B - .A)
            sum(w * (dbeta(x=zedd, shape1=shapes[,1], shape2=shapes[,2],
                           log = TRUE) - log( abs( .B - .A ))))
        }
    }, list( .lshape1 = lshape1, .lshape2 = lshape2, .A = A, .B = B, 
             .eshape1 = eshape1, .eshape2 = eshape2 ))),
    vfamily = "beta.ab",
    deriv = eval(substitute(expression({
        shapes = cbind(eta2theta(eta[,1], .lshape1, earg = .eshape1 ),
                       eta2theta(eta[,2], .lshape2, earg = .eshape2 ))
        dshapes.deta = cbind(dtheta.deta(shapes[,1], .lshape1, earg = .eshape1),
                             dtheta.deta(shapes[,2], .lshape2, earg = .eshape2))
        dl.dshapes = cbind(log(y-.A), log( .B-y)) - digamma(shapes) +
                     digamma(shapes[,1] + shapes[,2]) - log( .B - .A)
        c(w) * dl.dshapes * dshapes.deta
    }), list( .lshape1 = lshape1, .lshape2 = lshape2, .A = A, .B = B, 
              .eshape1 = eshape1, .eshape2 = eshape2 ))),
    weight = expression({
        temp2 = trigamma(shapes[,1]+shapes[,2])
        d2l.dshape12 = temp2 - trigamma(shapes[,1])
        d2l.dshape22 = temp2 - trigamma(shapes[,2])
        d2l.dshape1shape2 = temp2

        wz = matrix(as.numeric(NA), n, dimm(M))   #3=dimm(M)
        wz[,iam(1,1,M)] = d2l.dshape12 * dshapes.deta[,1]^2
        wz[,iam(2,2,M)] = d2l.dshape22 * dshapes.deta[,2]^2
        wz[,iam(1,2,M)] = d2l.dshape1shape2 * dshapes.deta[,1] * dshapes.deta[,2]

        -c(w) * wz
    }))
}



 beta4 = function(link = "loge", earg = list(),
                  i1=2.3, i2=2.4, iA = NULL, iB = NULL)
{



    if (mode(link) != "character" && mode(link) != "name")
        link = as.character(substitute(link))
    if (!is.list(earg)) earg = list()

    new("vglmff",
    blurb = c("Four-parameter Beta distribution\n",
            "(y-A)^(shape1-1) * (B-y)^(shape2-1), A < y < B \n\n",
            "Links:    ",
            namesof("shape1", link, earg = earg),  ", ",
            namesof("shape2", link, earg = earg), ", ",
            " A, B"),
    initialize = eval(substitute(expression({
        if (!is.vector(y) || (is.matrix(y) && ncol(y) != 1))
            stop("y must be a vector or a one-column matrix")

        if (length( .iA) && any(y < .iA))
            stop("initial 'A' value out of range")
        if (length( .iB) && any(y > .iB))
            stop("initial 'B' value out of range")

        predictors.names = c(
          namesof("shape1", .link, earg = .earg, short = TRUE),
          namesof("shape2", .link, earg = .earg, short = TRUE), "A", "B")
        my.range = diff(range(y))
        if (!length(etastart)) {
            etastart = cbind(shape1= rep( .i1, len = length(y)),
                             shape2= .i2,
                             A = if (length( .iA)) .iA else min(y)-my.range/70,
                             B = if (length( .iB)) .iB else max(y)+my.range/70)
        }
    }), list( .i1=i1, .i2=i2, .iA=iA, .iB=iB, .link = link, .earg = earg ))), 
    inverse = eval(substitute(function(eta, extra = NULL) {
        shapes = eta2theta(eta[,1:2], .link, earg = .earg)
        .A = eta[,3]
        .B = eta[,4]
        .A + ( .B-.A) * shapes[,1] / (shapes[,1] + shapes[,2])
    }, list( .link = link, .earg = earg ))),
    last = eval(substitute(expression({
        misc$link =    c(shape1 = .link, shape2 = .link, 
                         A = "identity", B = "identity")
        misc$earg = list(shape1 = .earg, shape2 = .earg, 
                         A = list(),     B = list())
    }), list( .link = link, .earg = earg ))), 
    loglikelihood = eval(substitute(
        function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
        shapes = eta2theta(eta[,1:2], .link, earg = .earg)
        .A = eta[,3]
        .B = eta[,4]
        temp = lbeta(shapes[,1], shapes[,2])
        if (residuals) stop("loglikelihood residuals not ",
                            "implemented yet") else
        sum(w * ((shapes[,1]-1)*log(y-.A) + (shapes[,2]-1)*log( .B-y) - temp -
            (shapes[,1]+shapes[,2]-1)*log( .B-.A )))
    }, list( .link = link, .earg = earg ))), 
    vfamily = "beta4",
    deriv = eval(substitute(expression({
        shapes = eta2theta(eta[,1:2], .link, earg = .earg)
        .A = eta[,3]
        .B = eta[,4]
        dshapes.deta = dtheta.deta(shapes, .link, earg = .earg)
        rr1 = ( .B - .A)
        temp3 = (shapes[,1] + shapes[,2] - 1)
        temp1 = temp3 / rr1
        dl.dshapes = cbind(log(y-.A), log( .B-y)) - digamma(shapes) +
                     digamma(shapes[,1] + shapes[,2]) - log( .B - .A)
        dl.dA = -(shapes[,1]-1) / (y- .A)  + temp1
        dl.dB =  (shapes[,2]-1) / ( .B - y) - temp1
        c(w) * cbind(dl.dshapes * dshapes.deta, dl.dA, dl.dB)
    }), list( .link = link, .earg = earg ))), 
    weight = expression({

        temp2 = trigamma(shapes[,1]+shapes[,2])
        d2l.dshape12 = temp2 - trigamma(shapes[,1])
        d2l.dshape22 = temp2 - trigamma(shapes[,2])
        d2l.dshape1shape2 = temp2

        ed2l.dAA = -temp3 * shapes[,2] / ((shapes[,1]-2) * rr1^2)
        ed2l.dBB = -temp3 * shapes[,1] / ((shapes[,2]-2) * rr1^2)
        ed2l.dAB = -temp3 / (rr1^2)
        ed2l.dAshape1 = -shapes[,2] / ((shapes[,1]-1) * rr1)
        ed2l.dAshape2 = 1/rr1
        ed2l.dBshape1 = -1/rr1
        ed2l.dBshape2 = shapes[,1] / ((shapes[,2]-1) * rr1)

        wz = matrix(as.numeric(NA), n, dimm(M))   #10=dimm(M)
        wz[,iam(1,1,M)] = d2l.dshape12 * dshapes.deta[,1]^2
        wz[,iam(2,2,M)] = d2l.dshape22 * dshapes.deta[,2]^2
        wz[,iam(1,2,M)] = d2l.dshape1shape2 * dshapes.deta[,1] * dshapes.deta[,2]

        wz[,iam(3,3,M)] = ed2l.dAA
        wz[,iam(4,4,M)] = ed2l.dBB
        wz[,iam(4,3,M)] = ed2l.dAB

        wz[,iam(3,1,M)] = ed2l.dAshape1 * dshapes.deta[,1]
        wz[,iam(3,2,M)] = ed2l.dAshape2 * dshapes.deta[,2]
        wz[,iam(4,1,M)] = ed2l.dBshape1 * dshapes.deta[,1]
        wz[,iam(4,2,M)] = ed2l.dBshape2 * dshapes.deta[,2]


        -c(w) * wz
    }))
}



 simple.exponential = function()
{
  new("vglmff",
  blurb = c("Simple Exponential distribution\n",
          "Link:    log(rate)\n"),
  deviance= function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
      devy = -log(y) - 1
      devmu = -log(mu) - y/mu
      devi = 2 * (devy - devmu)
      if (residuals) sign(y - mu) * sqrt(abs(devi) * w) else sum(w * devi)
  },
  initialize=expression({
      predictors.names = "log(rate)"
      mustart = y + (y == 0) / 8
  }),
  inverse=function(eta, extra = NULL)
      exp(-eta),
  link=function(mu, extra = NULL)
      -log(mu),
  vfamily = "simple.exponential",
  deriv=expression({
      rate = 1 / mu
      dl.drate = mu - y
      drate.deta = dtheta.deta(rate, "loge")
      c(w) * dl.drate * drate.deta
  }),
  weight = expression({
      ed2l.drate2 = -1 / rate^2
      wz = -c(w) * drate.deta^2 * ed2l.drate2
      wz
  }))
}




 exponential <- function(link = "loge", earg = list(),
                         location = 0, expected = TRUE) {
  if (!is.Numeric(location, allow = 1))
      stop("bad input for argument 'location'")

  if (mode(link) != "character" && mode(link) != "name")
    link <- as.character(substitute(link))
  if (!is.list(earg)) earg = list()
  if (!is.logical(expected) || length(expected) != 1)
    stop("bad input for argument 'expected'")

  new("vglmff",
  blurb = c("Exponential distribution\n\n",
            "Link:     ", namesof("rate", link, tag = TRUE), "\n",
            "Mean:     ", "mu = ", 
             if (location == 0) "1/rate" else
             paste(location, "+ 1/rate"), "\n",
            "Variance: ",
            if (location == 0) "Exponential: mu^2" else
            paste("(mu - ", location, ")^2", sep = "")),
  deviance = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
    devy <- -log(y - .location) - 1
    devmu <- -log(mu - .location) - (y - .location) / (mu - .location)
    devi <- 2 * (devy - devmu)
    if (residuals) {
      sign(y - mu) * sqrt(abs(devi) * w)
    } else  {
      sum(w * devi)
    }
  }, list( .location = location, .earg = earg ))),
  initialize = eval(substitute(expression({
    if (ncol(cbind(y)) != 1)
      stop("response must be a vector or a one-column matrix")

    extra$loc <- .location # Passed into, e.g., @link, @deriv etc.

    if (any(y <= extra$loc))
      stop("all responses must be greater than ", extra$loc)

    predictors.names <- namesof("rate", .link, tag = FALSE)

    if (length(mustart) + length(etastart) == 0)
      mustart <- y + (y == extra$loc) / 8
    if (!length(etastart))
        etastart <- theta2eta(1 / (mustart - extra$loc),
                              .link, earg = .earg)
  }), list( .location = location, .link = link, .earg = earg ))),
  inverse = eval(substitute(function(eta, extra = NULL)
    extra$loc + 1 / eta2theta(eta, .link, earg = .earg),
  list( .link = link, .earg = earg ))),
  last = eval(substitute(expression({
    misc$location <- extra$loc
    misc$link <-    c(rate = .link)
    misc$earg <- list(rate = .earg)
    misc$expected <- .expected
  }), list( .link = link, .earg = earg, .expected = expected ))),
  link = eval(substitute(function(mu, extra = NULL) 
    theta2eta(1 / (mu - extra$loc), .link, earg = .earg),
  list( .link = link, .earg = earg ))),
  vfamily = c("exponential"),
  deriv = eval(substitute(expression({
    rate <- 1 / (mu - extra$loc)
    dl.drate <- mu - y
    drate.deta <- dtheta.deta(rate, .link, earg = .earg)
    c(w) * dl.drate * drate.deta
  }), list( .link = link, .earg = earg ))),
  weight = eval(substitute(expression({
    d2l.drate2 <- -((mu-extra$loc)^2)
    wz <- -(drate.deta^2) * d2l.drate2
    if (! .expected) {
      d2rate.deta2 <- d2theta.deta2(rate, .link, earg = .earg)
      wz <- wz - dl.drate * d2rate.deta2
    }
      c(w) * wz
  }), list( .link = link, .expected = expected, .earg = earg ))))
}




 gamma1 = function(link = "loge", earg = list())
{
    if (mode(link) != "character" && mode(link) != "name")
        link = as.character(substitute(link))
    if (!is.list(earg)) earg = list()

    new("vglmff",
    blurb = c("1-parameter Gamma distribution\n",
            "Link:     ",
            namesof("shape", link, earg = earg, tag = TRUE), "\n", 
            "Mean:       mu (=shape)\n",
            "Variance:   mu (=shape)"),
    initialize = eval(substitute(expression({
        if (any(y <= 0))
            stop("all responses must be positive")
        M = if (is.matrix(y)) ncol(y) else 1
        temp.names = if (M == 1) "shape" else paste("shape", 1:M, sep = "")
        predictors.names =
          namesof(temp.names, .link, earg = .earg, short = TRUE)
        if (!length(etastart))
            etastart = cbind(theta2eta(y + 1/8, .link, earg = .earg ))
    }), list( .link = link, .earg = earg ))), 
    inverse = eval(substitute(function(eta, extra = NULL)
        eta2theta(eta, .link, earg = .earg)),
    list( .link = link, .earg = earg )),
    last = eval(substitute(expression({
        temp.names = if (M == 1) "shape" else paste("shape", 1:M, sep = "")
        misc$link = rep( .link, length = M)
        names(misc$link) = temp.names
        misc$earg = vector("list", M)
        names(misc$earg) = names(misc$link)
        for(ii in 1:M) misc$earg[[ii]] = .earg
        misc$expected = TRUE
    }), list( .link = link, .earg = earg ))),
    link = eval(substitute(function(mu, extra = NULL)
        theta2eta(mu, .link, earg = .earg)),
    list( .link = link, .earg = earg )),
    loglikelihood= function(mu, y, w, residuals = FALSE, eta, extra = NULL)
        if (residuals) stop("loglikelihood residuals not ",
                            "implemented yet") else {
        sum(w * dgamma(x = y, shape=mu, scale = 1, log = TRUE))
    },
    vfamily = c("gamma1"),
    deriv = eval(substitute(expression({
        shape = mu
        dl.dshape = log(y) - digamma(shape)
        dshape.deta = dtheta.deta(shape, .link, earg = .earg)
        c(w) * dl.dshape * dshape.deta
    }), list( .link = link, .earg = earg ))),
    weight = expression({
        d2l.dshape = -trigamma(shape)
        wz = -(dshape.deta^2) * d2l.dshape
        c(w) * wz
    }))
}


 gamma2.ab = function(lrate = "loge", lshape = "loge",
                      erate = list(), eshape = list(),
                      irate = NULL, ishape = NULL, expected = TRUE, zero = 2)
{
    if (mode(lrate) != "character" && mode(lrate) != "name")
        lrate = as.character(substitute(lrate))
    if (mode(lshape) != "character" && mode(lshape) != "name")
        lshape = as.character(substitute(lshape))
    if (length( irate) && !is.Numeric(irate, posit = TRUE))
        stop("bad input for argument 'irate'")
    if (length( ishape) && !is.Numeric(ishape, posit = TRUE))
        stop("bad input for argument 'ishape'")
    if (length(zero) && !is.Numeric(zero, integer = TRUE, posit = TRUE))
        stop("bad input for argument 'zero'")
    if (!is.logical(expected) || length(expected) != 1)
        stop("bad input for argument 'expected'")
    if (!is.list(erate)) erate = list()
    if (!is.list(eshape)) eshape = list()

    new("vglmff",
    blurb = c("2-parameter Gamma distribution\n",
            "Links:    ",
            namesof("rate",  lrate,  earg = erate), ", ", 
            namesof("shape", lshape, earg = eshape), "\n",
            "Mean:     mu = shape/rate\n",
            "Variance: (mu^2)/shape = shape/rate^2"),
    constraints = eval(substitute(expression({
        constraints = cm.zero.vgam(constraints, x, .zero, M)
    }), list( .zero = zero ))),
    initialize = eval(substitute(expression({
        # Error check
        if (ncol(cbind(y)) != 1)
            stop("response must be a vector or a one-column matrix")
        if (any(y <= 0))
            stop("all responses must be positive")
        predictors.names =
            c(namesof("rate",  .lrate,  earg = .erate,  tag = FALSE),
              namesof("shape", .lshape, earg = .eshape, tag = FALSE))
        if (!length(etastart)) {
            mymu = y + 0.167 * (y == 0)
            junk = lsfit(x, y, wt = w, intercept = FALSE)
            var.y.est = sum(w * junk$resid^2) / (nrow(x) - length(junk$coef))
            init.shape =  if (length( .ishape)) .ishape else mymu^2 / var.y.est
            init.rate =  if (length( .irate)) .irate else init.shape / mymu
            init.rate = rep(init.rate, len = n)
            init.shape = rep(init.shape, len = n)
            if ( .lshape == "loglog")
                init.shape[init.shape <= 1] = 3.1 #Hopefully value is big enough
            etastart = cbind(theta2eta(init.rate, .lrate, earg = .erate),
                             theta2eta(init.shape, .lshape, earg = .eshape))
        }
    }), list( .lrate = lrate, .lshape = lshape, .irate=irate, .ishape = ishape,
              .erate = erate, .eshape = eshape ))),
    inverse = eval(substitute(function(eta, extra = NULL) {
        eta2theta(eta[,2], .lshape, earg = .eshape) / eta2theta(eta[,1], .lrate,
        earg = .erate)
    }, list( .lrate = lrate, .lshape = lshape,
             .erate = erate, .eshape = eshape ))),
    last = eval(substitute(expression({
        misc$link =    c(rate = .lrate, shape = .lshape)
        misc$earg = list(rate = .erate, shape = .eshape)
    }), list( .lrate = lrate, .lshape = lshape,
              .erate = erate, .eshape = eshape ))),
    loglikelihood = eval(substitute(
        function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
        rate = eta2theta(eta[,1], .lrate, earg = .erate)
        shape = eta2theta(eta[,2], .lshape, earg = .eshape)
        if (residuals) stop("loglikelihood residuals not ",
                            "implemented yet") else {
            sum(w * dgamma(x = y, shape = shape, rate=rate, log = TRUE))
        }
    }, list( .lrate = lrate, .lshape = lshape,
             .erate = erate, .eshape = eshape ))),
    vfamily = c("gamma2.ab"),
    deriv = eval(substitute(expression({
        rate = eta2theta(eta[,1], .lrate, earg = .erate)
        shape = eta2theta(eta[,2], .lshape, earg = .eshape)
        dl.drate = mu - y
        dl.dshape = log(y*rate) - digamma(shape)
        dratedeta = dtheta.deta(rate, .lrate, earg = .erate)
        dshape.deta = dtheta.deta(shape, .lshape, earg = .eshape)
        c(w) * cbind(dl.drate * dratedeta,
                     dl.dshape * dshape.deta)
    }), list( .lrate = lrate, .lshape = lshape,
              .erate = erate, .eshape = eshape ))),
    weight = eval(substitute(expression({
        d2l.dshape2 = -trigamma(shape)
        d2l.drate2 = -shape/(rate^2)
        d2l.drateshape = 1/rate
        wz = matrix(as.numeric(NA), n, dimm(M))  #3=dimm(M)
        wz[,iam(1,1,M)] = -d2l.drate2 * dratedeta^2
        wz[,iam(2,2,M)] = -d2l.dshape2 * dshape.deta^2
        wz[,iam(1,2,M)] = -d2l.drateshape * dratedeta * dshape.deta
        if (! .expected) {
            d2ratedeta2 = d2theta.deta2(rate, .lrate, earg = .erate)
            d2shapedeta2 = d2theta.deta2(shape, .lshape, earg = .eshape)
            wz[,iam(1,1,M)] = wz[,iam(1,1,M)] - dl.drate * d2ratedeta2
            wz[,iam(2,2,M)] = wz[,iam(2,2,M)] - dl.dshape * d2shapedeta2
        }
        c(w) * wz
    }), list( .lrate = lrate, .lshape = lshape,
              .erate = erate, .eshape = eshape, .expected = expected ))))
}



 gamma2 = function(lmu = "loge", lshape = "loge",
                   emu = list(), eshape = list(),
                   imethod = 1,
                   deviance.arg = FALSE, ishape = NULL, zero = -2)
{

    if (mode(lmu) != "character" && mode(lmu) != "name")
        lmu = as.character(substitute(lmu))
    if (mode(lshape) != "character" && mode(lshape) != "name")
        lshape = as.character(substitute(lshape))
    if (length(zero) && !is.Numeric(zero, integer = TRUE))
        stop("bad input for argument 'zero'")
    if (length( ishape) && !is.Numeric(ishape, posit = TRUE))
        stop("bad input for argument 'ishape'")
    if (!is.Numeric(imethod, allow = 1, integ = TRUE, posit = TRUE) ||
       imethod > 2)
        stop("argument 'imethod' must be 1 or 2")
    if (!is.list(emu)) emu = list()
    if (!is.list(eshape)) eshape = list()

    ans = 
    new("vglmff",
    blurb = c("2-parameter Gamma distribution",
            " (McCullagh and Nelder 1989 parameterization)\n",
            "Links:    ",
            namesof("mu",    lmu,    earg = emu), ", ", 
            namesof("shape", lshape, earg = eshape), "\n",
            "Mean:     mu\n",
            "Variance: (mu^2)/shape"),
    constraints = eval(substitute(expression({

        dotzero <- .zero
        Musual <- 2
        eval(negzero.expression)
        constraints = cm.zero.vgam(constraints, x, z_Index, M)
    }), list( .zero = zero ))),
    initialize = eval(substitute(expression({
        Musual <- 2

        assign("CQO.FastAlgorithm", ( .lmu == "loge" && .lshape == "loge"),
               envir = VGAM:::VGAMenv)
        if (any(function.name == c("cqo","cao")) &&
           is.Numeric( .zero, allow = 1) && .zero != -2)
            stop("argument zero=-2 is required")

        y = as.matrix(y)
        M = Musual * ncol(y)
        NOS = ncoly = ncol(y)  # Number of species
        temp1.names =
          if (NOS == 1) "mu"    else paste("mu",    1:NOS, sep = "")
        temp2.names =
          if (NOS == 1) "shape" else paste("shape", 1:NOS, sep = "")
        predictors.names =
            c(namesof(temp1.names, .lmu,    earg = .emu,    tag = FALSE),
              namesof(temp2.names, .lshape, earg = .eshape, tag = FALSE))
        predictors.names = predictors.names[interleave.VGAM(M, M = Musual)]


        # Error check
        if (any(y <= 0))
            stop("all responses must be positive") # see @loglikelihood
        if (!length(etastart)) {
            init.shape = matrix(1.0, n, NOS)
            mymu = y # + 0.167 * (y == 0)  # imethod == 1 (the default)
            if ( .imethod == 2) {
                for(ii in 1:ncol(y)) {
                    mymu[,ii] = weighted.mean(y[,ii], w = w)
                }
            }
            for(spp in 1:NOS) {
                junk = lsfit(x, y[,spp], wt = w, intercept = FALSE)
                var.y.est = sum(w * junk$resid^2) / (n - length(junk$coef))
                init.shape[,spp] = if (length( .ishape)) .ishape else
                    mymu[,spp]^2 / var.y.est
                if ( .lshape == "loglog") init.shape[init.shape[,spp] <=
                             1,spp] = 3.1 # Hopefully value is big enough
            }
            etastart = cbind(theta2eta(mymu, .lmu, earg = .emu ),
                             theta2eta(init.shape, .lshape, earg = .eshape ))
            etastart = etastart[,interleave.VGAM(M, M = Musual), drop = FALSE]
        }
    }), list( .lmu = lmu, .lshape = lshape, .ishape = ishape, .zero = zero,
              .emu = emu, .eshape = eshape,
              .imethod = imethod ))),
    inverse = eval(substitute(function(eta, extra = NULL) {
        Musual <- 2
        NOS = ncol(eta) / Musual
        eta2theta(eta[, 2*(1:NOS)-1, drop = FALSE], .lmu, earg = .emu )
    }, list( .lmu = lmu, .emu = emu ))),
    last = eval(substitute(expression({
        if (exists("CQO.FastAlgorithm", envir = VGAM:::VGAMenv))
            rm("CQO.FastAlgorithm", envir = VGAM:::VGAMenv)
        tmp34 = c(rep( .lmu,    length = NOS),
                  rep( .lshape, length = NOS))
        names(tmp34) =
           c(if (NOS == 1) "mu"    else paste("mu",    1:NOS, sep = ""), 
             if (NOS == 1) "shape" else paste("shape", 1:NOS, sep = ""))
        tmp34 = tmp34[interleave.VGAM(M, M = 2)]
        misc$link = tmp34 # Already named
        misc$earg = vector("list", M)
        misc$Musual <- Musual
        names(misc$earg) = names(misc$link)
        for(ii in 1:NOS) {
            misc$earg[[2*ii-1]] = .emu
            misc$earg[[2*ii  ]] = .eshape
        }
        misc$expected = TRUE
    }), list( .lmu = lmu, .lshape = lshape,
              .emu = emu, .eshape = eshape ))),
    link = eval(substitute(function(mu, extra = NULL) {
        temp = theta2eta(mu, .lmu, earg = .emu )
        temp = cbind(temp, NA * temp)
        temp[,interleave.VGAM(ncol(temp), M = 2), drop = FALSE]
    }, list( .lmu = lmu, .emu = emu ))),
    loglikelihood = eval(substitute(
        function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
        Musual <- 2
        NOS = ncol(eta) / Musual
        mymu = mu  # eta2theta(eta[,2*(1:NOS)-1], .lmu, earg = .emu )
        shapemat = eta2theta(eta[,2*(1:NOS), drop = FALSE], .lshape, earg = .eshape )
        if (residuals) stop("loglikelihood residuals not ",
                            "implemented yet") else {
            sum(w * dgamma(x = y, shape = c(shapemat), scale = c(mymu/shapemat),
                           log = TRUE))
        }
    }, list( .lmu = lmu, .lshape = lshape,
             .emu = emu, .eshape = eshape ))),
    vfamily = c("gamma2"),
    deriv = eval(substitute(expression({
        Musual <- 2
        NOS = ncol(eta) / Musual

        mymu  = eta2theta(eta[,2*(1:NOS)-1], .lmu,    earg = .emu )
        shape = eta2theta(eta[,2*(1:NOS)],   .lshape, earg = .eshape )

        dl.dmu = shape * (y / mymu - 1) / mymu
        dl.dshape = log(y) + log(shape) - log(mymu) + 1 - digamma(shape) -
                    y / mymu

        dmu.deta    = dtheta.deta(mymu,  .lmu,    earg = .emu )
        dshape.deta = dtheta.deta(shape, .lshape, earg = .eshape )

        myderiv = c(w) * cbind(dl.dmu * dmu.deta,
                               dl.dshape * dshape.deta)
        myderiv[, interleave.VGAM(M, M = Musual)]
    }), list( .lmu = lmu, .lshape = lshape,
              .emu = emu, .eshape = eshape ))),
    weight = eval(substitute(expression({
        ed2l.dmu2 = shape / (mymu^2)
        ed2l.dshape2 = trigamma(shape) - 1 / shape
        wz = matrix(as.numeric(NA), n, M)  # 2 = M; diagonal!

        wz[,2*(1:NOS)-1] = ed2l.dmu2 * dmu.deta^2
        wz[,2*(1:NOS)] = ed2l.dshape2 * dshape.deta^2
        c(w) * wz
    }), list( .lmu = lmu ))))

    if (deviance.arg) ans@deviance = eval(substitute(
        function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
        NOS = ncol(eta) / 2
        temp300 =  eta[,2*(1:NOS), drop = FALSE]
        if ( .lshape == "loge") {
            bigval = 28
            temp300[temp300 >  bigval] =  bigval
            temp300[temp300 < -bigval] = -bigval
        } else stop("can only handle the 'loge' link")
        shape =  eta2theta(temp300, .lshape, earg = .eshape )
        devi = -2 * (log(y/mu) - y/mu + 1)
        if (residuals) {
           warning("not 100% sure about these deviance residuals!")
           sign(y - mu) * sqrt(abs(devi) * w)
        } else
           sum(w * devi)
    }, list( .lshape = lshape )))
    ans
}



 geometric =function(link = "logit", earg = list(), expected = TRUE,
                     imethod= 1)
{
    if (!is.logical(expected) || length(expected) != 1)
        stop("bad input for argument 'expected'")
    if (mode(link) != "character" && mode(link) != "name")
        link = as.character(substitute(link))
    if (!is.list(earg)) earg = list()
    if (!is.Numeric(imethod, allow = 1, integ = TRUE, posit = TRUE) ||
       imethod > 3)
        stop("argument 'imethod' must be 1 or 2 or 3")

    new("vglmff",
    blurb = c("Geometric distribution (P[Y=y] = prob*(1-prob)^y, y=0,1,2,...)\n",
            "Link:     ",
            namesof("prob", link, earg = earg), "\n",
            "Mean:     mu = (1-prob)/prob\n",
            "Variance: mu*(1+mu) = (1-prob)/prob^2"),
    initialize = eval(substitute(expression({
        if (ncol(cbind(y)) != 1)
            stop("response must be a vector or a 1-column matrix")
        if (any(y < 0)) stop("all responses must be >= 0")
        if (any(y!=round(y ))) stop("response should be integer-valued")
        predictors.names = namesof("prob", .link, earg = .earg, tag = FALSE)
        if (!length(etastart)) {
            prob.init = if ( .imethod == 3)
                            1 / (1 + y + 1/16) else
                        if ( .imethod == 1)
                            1 / (1 + median(rep(y,w)) + 1/16) else
                        1 / (1 + weighted.mean(y,w) + 1/16)
            etastart = theta2eta(prob.init, .link, earg = .earg)
        }
    }), list( .link = link, .earg = earg, .imethod = imethod ))),
    inverse = eval(substitute(function(eta, extra = NULL) {
        prob = eta2theta(eta, .link, earg = .earg)
        (1-prob)/prob 
    }, list( .link = link, .earg = earg ))),
    last = eval(substitute(expression({
        misc$link =    c(prob = .link)
        misc$earg = list(prob = .earg )
        misc$expected = .expected
        misc$imethod = .imethod
    }), list( .link = link, .earg = earg, .expected = expected, .imethod = imethod ))),
    loglikelihood = eval(substitute(
        function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
        prob = eta2theta(eta, .link, earg = .earg)
        if (residuals) stop("loglikelihood residuals not ",
                            "implemented yet") else {
            sum(w * dgeom(x = y, prob=prob, log = TRUE))
        }
    }, list( .link = link, .earg = earg ))),
    vfamily = c("geometric"),
    deriv = eval(substitute(expression({
        prob = eta2theta(eta, .link, earg = .earg)
        dl.dprob = -y/(1-prob) + 1/prob 
        dprobdeta = dtheta.deta(prob, .link, earg = .earg)
        c(w) * cbind(dl.dprob * dprobdeta)
    }), list( .link = link, .earg = earg, .expected = expected ))),
    weight = eval(substitute(expression({
        ed2l.dprob2 = if ( .expected ) 1 / (prob^2 * (1-prob)) else
            y / (1-prob)^2 + 1 / prob^2
        wz = ed2l.dprob2 * dprobdeta^2
        if ( !( .expected )) wz = wz - dl.dprob * d2theta.deta2(prob, .link, earg = .earg)
        c(w) * wz
    }), list( .link = link, .earg = earg, .expected = expected ))))
}


dbetageom = function(x, shape1, shape2, log = FALSE) {
    if (!is.logical(log.arg <- log))
        stop("bad input for argument 'log'")
    rm(log)

    if (!is.Numeric(x)) stop("bad input for argument 'x'")
    if (!is.Numeric(shape1, pos = TRUE)) stop("bad input for argument 'shape1'")
    if (!is.Numeric(shape2, pos = TRUE)) stop("bad input for argument 'shape2'")
    N = max(length(x), length(shape1), length(shape2))
    x = rep(x, len=N); shape1 = rep(shape1, len=N); shape2 = rep(shape2, len=N)
    loglik = lbeta(1+shape1, shape2+abs(x)) - lbeta(shape1, shape2)
    xok = (x == round(x) & x >= 0)
    loglik[!xok] = log(0)
    if (log.arg) {
        loglik
    } else {
        exp(loglik)
    }
}


pbetageom = function(q, shape1, shape2, log.p = FALSE) {
    if (!is.Numeric(q)) stop("bad input for argument 'q'")
    if (!is.Numeric(shape1, pos = TRUE)) stop("bad input for argument 'shape1'")
    if (!is.Numeric(shape2, pos = TRUE)) stop("bad input for argument 'shape2'")
    N = max(length(q), length(shape1), length(shape2))
    q = rep(q, len=N); shape1 = rep(shape1, len=N); shape2 = rep(shape2, len=N)
    ans = q * 0  # Retains names(q)
    if (max(abs(shape1-shape1[1])) < 1.0e-08 &&
       max(abs(shape2-shape2[1])) < 1.0e-08) {
        qstar = floor(q)
        temp = if (max(qstar) >= 0) dbetageom(x=0:max(qstar), 
               shape1=shape1[1], shape2=shape2[1]) else 0*qstar
        unq = unique(qstar)
        for(i in unq) {
            index = qstar == i
            ans[index] = if (i >= 0) sum(temp[1:(1+i)]) else 0
        }
    } else
    for(ii in 1:N) {
        qstar = floor(q[ii])
        ans[ii] = if (qstar >= 0) sum(dbetageom(x=0:qstar, 
                 shape1=shape1[ii], shape2=shape2[ii])) else 0
    }
    if (log.p) log(ans) else ans
}

rbetageom = function(n, shape1, shape2) {
    if (!is.Numeric(n, integ = TRUE,allow = 1)) stop("bad input for argument 'n'")
    if (!is.Numeric(shape1, pos = TRUE)) stop("bad input for argument 'shape1'")
    if (!is.Numeric(shape2, pos = TRUE)) stop("bad input for argument 'shape2'")
    rgeom(n = n, prob = rbeta(n = n, shape1=shape1, shape2=shape2))
}





interleave.VGAM = function(L, M) c(matrix(1:L, nrow=M, byrow = TRUE))



negbinomial.control <- function(save.weight = TRUE, ...)
{
    list(save.weight = save.weight)
}



 negbinomial = function(lmu = "loge", lsize = "loge",
                        emu = list(), esize = list(),
                        imu = NULL,   isize = NULL,
                        quantile.probs = 0.75,
                        nsimEIM = 100, cutoff = 0.995, Maxiter = 5000,
                        deviance.arg = FALSE, imethod = 1,
                        parallel = FALSE,
                        shrinkage.init = 0.95, zero = -2)
{






  lmuuu = lmu
  emuuu = emu
  imuuu = imu


  if (length(imuuu) && !is.Numeric(imuuu, posit = TRUE))
    stop("bad input for argument 'imu'")
  if (length(isize) && !is.Numeric(isize, posit = TRUE))
    stop("bad input for argument 'isize'")

  if (!is.Numeric(cutoff, allow = 1) || cutoff < 0.8 || cutoff >= 1)
    stop("range error in the argument 'cutoff'")
  if (!is.Numeric(Maxiter, integ = TRUE, allow = 1) || Maxiter < 100)
    stop("bad input for argument 'Maxiter'")
  if (!is.Numeric(imethod, allow = 1, integ = TRUE, posit = TRUE) ||
     imethod > 3)
    stop("argument 'imethod' must be 1 or 2 or 3")
  if (!is.Numeric(shrinkage.init, allow = 1) || shrinkage.init < 0 ||
     shrinkage.init > 1)
    stop("bad input for argument 'shrinkage.init'")

  if (!is.null(nsimEIM)) {
    if (!is.Numeric(nsimEIM, allow = 1, integ = TRUE))
      stop("bad input for argument 'nsimEIM'")
    if (nsimEIM <= 10)
      warning("argument 'nsimEIM' should be an integer ",
               "greater than 10, say")
  }

  if (mode(lmuuu) != "character" && mode(lmuuu) != "name")
    lmuuu = as.character(substitute(lmuuu))
  if (mode(lsize) != "character" && mode(lsize) != "name")
    lsize = as.character(substitute(lsize))
  if (!is.list(emuuu)) emuuu = list()
  if (!is.list(esize)) esize = list()


    if (!is.logical( parallel ) || length( parallel ) != 1)
      stop("argument 'parallel' must be TRUE or FALSE")
    if ( parallel  && length(zero))
      stop("need to set 'zero = NULL' when parallel = TRUE")



  ans = 
  new("vglmff",
  blurb = c("Negative-binomial distribution\n\n",
            "Links:    ",
            namesof("mu",   lmuuu, earg = emuuu), ", ",
            namesof("size", lsize, earg = esize), "\n",
            "Mean:     mu\n",
            "Variance: mu * (1 + mu / size)"),
  constraints = eval(substitute(expression({

    dotzero <- .zero
    Musual <- 2
    eval(negzero.expression)

    if ( .parallel && ncol(cbind(y)) > 1)
      stop("univariate responses needed if parallel = TRUE")
    constraints = cm.vgam(matrix(1, M, 1), x, .parallel, constraints)
  }), list( .parallel = parallel, .zero = zero ))),

  infos = eval(substitute(function(...) {
    list(Musual = 2,
         zero = .zero)
  }, list( .zero = zero ))),

  initialize = eval(substitute(expression({
    Musual <- 2

    assign("CQO.FastAlgorithm",
          ( .lmuuu == "loge") && ( .lsize == "loge"),
           envir = VGAM:::VGAMenv)
    if (any(function.name == c("cqo","cao")) &&
        is.Numeric( .zero, allow = 1) && .zero != -2)
        stop("argument zero = -2 is required")

    if (any(y < 0))
      stop("negative values not allowed for the 'negbinomial' family")
    if (any(round(y) != y))
      stop("integer-values only allowed for the 'negbinomial' family")

    y = as.matrix(y) 
    M = 2 * ncol(y) 
    NOS = ncoly = ncol(y)  # Number of species
    predictors.names =
      c(namesof(if (NOS == 1) "mu" else paste("mu",   1:NOS, sep = ""),
                .lmuuu, earg = .emuuu, tag = FALSE),
       namesof(if (NOS == 1) "size" else paste("size", 1:NOS, sep = ""),
                .lsize, earg = .esize, tag = FALSE))
    predictors.names = predictors.names[interleave.VGAM(M, M = 2)]

    if (is.null( .nsimEIM )) {
       save.weight <- control$save.weight <- FALSE
    }

    if (is.numeric( .mu.init ))
      MU.INIT <- matrix( .mu.init, nrow(y), ncol(y), byrow = TRUE)

    if (!length(etastart)) {
      mu.init = y
      for(iii in 1:ncol(y)) {
        use.this = if ( .imethod == 1) {
            weighted.mean(y[, iii], w) + 1/16
        } else if ( .imethod == 3) {
            c(quantile(y[, iii], probs = .quantile.probs) + 1/16)
        } else {
            median(y[, iii]) + 1/16
        }

        if (is.numeric( .mu.init )) {
            mu.init[, iii] = MU.INIT[, iii]
        } else {
          medabsres = median(abs(y[, iii] - use.this)) + 1/32
          allowfun = function(z, maxtol=1) sign(z)*pmin(abs(z), maxtol)
          mu.init[, iii] = use.this + (1 - .sinit) * allowfun(y[, iii] -
                          use.this, maxtol=medabsres)

          mu.init[, iii] = abs(mu.init[, iii]) + 1 / 1024
        }
      } # of for(iii)

      if ( is.Numeric( .k.init )) {
        kay.init = matrix( .k.init, nrow = n, ncol = NOS, byrow = TRUE)
      } else {
        negbinomial.Loglikfun = function(kmat, y, x, w, extraargs) {
            mu = extraargs
            sum(w * dnbinom(x = y, mu = mu, size = kmat, log = TRUE))
        }
        k.grid = 2^((-7):7)
        k.grid = 2^(seq(-8, 8, length = 40))
        kay.init = matrix(0, nr=n, nc=NOS)
        for(spp. in 1:NOS) {
          kay.init[,spp.] = getMaxMin(k.grid,
                            objfun=negbinomial.Loglikfun,
                            y = y[,spp.], x = x, w = w,
                            extraargs= mu.init[,spp.])
        }
      }
      etastart = cbind(theta2eta(mu.init,  .lmuuu, earg = .emuuu),
                       theta2eta(kay.init, .lsize, earg = .esize))
      etastart = etastart[,interleave.VGAM(M, M = 2), drop = FALSE]
      }
  }), list( .lmuuu = lmuuu, .lsize = lsize,
            .emuuu = emuuu, .esize = esize,
            .mu.init = imu,
            .k.init = isize, .quantile.probs = quantile.probs,
            .sinit = shrinkage.init, .nsimEIM = nsimEIM,
            .zero = zero, .imethod = imethod ))),
  inverse = eval(substitute(function(eta, extra = NULL) {
      NOS = ncol(eta) / 2
      eta2theta(eta[, 2*(1:NOS)-1, drop = FALSE],
                .lmuuu, earg = .emuuu)
  }, list( .lmuuu = lmuuu, .emuuu = emuuu,
           .esize = esize ))),
  last = eval(substitute(expression({
      if (exists("CQO.FastAlgorithm", envir = VGAM:::VGAMenv))
          rm("CQO.FastAlgorithm", envir = VGAM:::VGAMenv)

      temp0303 = c(rep( .lmuuu, length = NOS),
                   rep( .lsize, length = NOS))
      names(temp0303) = c(if (NOS == 1) "mu"   else
                          paste("mu",   1:NOS, sep = ""),
                          if (NOS == 1) "size" else
                          paste("size", 1:NOS, sep = ""))
      temp0303 = temp0303[interleave.VGAM(M, M = 2)]
      misc$link = temp0303 # Already named
      misc$earg = vector("list", M)
      names(misc$earg) = names(misc$link)
      for(ii in 1:NOS) {
          misc$earg[[2*ii-1]] = .emuuu
          misc$earg[[2*ii  ]] = .esize
      }
      misc$cutoff = .cutoff 
      misc$imethod = .imethod 
      misc$nsimEIM = .nsimEIM
      misc$expected = TRUE
      misc$shrinkage.init = .sinit
  }), list( .lmuuu = lmuuu, .lsize = lsize,
            .emuuu = emuuu, .esize = esize,
            .cutoff = cutoff,
            .nsimEIM = nsimEIM,
            .sinit = shrinkage.init,
            .imethod = imethod ))),
  link = eval(substitute(function(mu, extra = NULL) {
    temp = theta2eta(mu, .lmuuu, earg = .emuuu)
    kayy = theta2eta(if (is.numeric( .isize)) .isize else 1.0,
                     .lsize, earg = .esize)
    kayy = 0 * temp + kayy  # Right dimension now.
    temp = cbind(temp, kayy)
    temp[, interleave.VGAM(ncol(temp), M = 2), drop = FALSE]
  }, list( .lmuuu = lmuuu, .emuuu = emuuu,
           .lsize = lsize, .esize = esize,
           .isize = isize ))),
  loglikelihood = eval(substitute(
      function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
      NOS = ncol(eta) / 2
      temp300 = eta[, 2*(1:NOS), drop = FALSE]
      if ( .lsize == "loge") {
          bigval = 28
          temp300 = ifelse(temp300 >  bigval,  bigval, temp300)
          temp300 = ifelse(temp300 < -bigval, -bigval, temp300)
      }
      kmat = eta2theta(temp300, .lsize, earg = .esize)
      if (residuals) stop("loglikelihood residuals not ",
                          "implemented yet") else
        sum(w * dnbinom(x = y, mu = mu, size = kmat, log = TRUE))
  }, list( .lsize = lsize, .emu = emu, .esize = esize ))),
  vfamily = c("negbinomial"),
  deriv = eval(substitute(expression({
    NOS = ncol(eta) / 2
    M = ncol(eta)
    temp3 = eta[, 2*(1:NOS), drop = FALSE]
    bigval = 28
    temp3 = ifelse(temp3 >  bigval,  bigval, temp3)
    temp3 = ifelse(temp3 < -bigval, -bigval, temp3)
    kmat = eta2theta(temp3, .lsize, earg = .esize)

    dl.dmu = y/mu - (y+kmat)/(kmat+mu)
    dl.dk = digamma(y+kmat) - digamma(kmat) - (y+kmat)/(mu+kmat) + 1 +
            log(kmat/(kmat+mu))

    dmu.deta = dtheta.deta(mu, .lmu, earg = .emu)
    dk.deta = dtheta.deta(kmat, .lsize, earg = .esize)

    dthetas.detas = cbind(dmu.deta, dk.deta)
    myderiv = c(w) * cbind(dl.dmu, dl.dk) * dthetas.detas
    myderiv[, interleave.VGAM(M, M = 2)]
  }), list( .lmu = lmu, .lsize = lsize, .emu = emu, .esize = esize ))),
  weight = eval(substitute(expression({
    wz = matrix(as.numeric(NA), n, M)  # wz is 'diagonal' 
    if (is.null( .nsimEIM)) {
      fred2 = dotFortran(name = "enbin9", ans = double(n*NOS),
                  as.double(kmat), as.double(mu), as.double( .cutoff ),
                  as.integer(n), ok = as.integer(1), as.integer(NOS),
                  sumpdf = double(1), as.double( .Machine$double.eps ),
                  as.integer( .Maxiter ))
      if (fred2$ok != 1)
        stop("error in Fortran subroutine exnbin9")
      dim(fred2$ans) = c(n, NOS)
      ed2l.dk2 = -fred2$ans - 1/kmat + 1/(kmat+mu)
      wz[,2*(1:NOS)] = dk.deta^2 * ed2l.dk2
    } else {
      run.varcov = matrix(0, n, NOS)
      ind1 = iam(NA, NA, M = M, both = TRUE, diag = TRUE)
      for(ii in 1:( .nsimEIM )) {
        ysim = rnbinom(n = n*NOS, mu = c(mu), size = c(kmat))
        if (NOS > 1) dim(ysim) = c(n, NOS)
          dl.dk = digamma(ysim+kmat) - digamma(kmat) -
                  (ysim+kmat)/(mu+kmat) + 1 + log(kmat/(kmat+mu))
              run.varcov = run.varcov + dl.dk^2
          }
          run.varcov = cbind(run.varcov / .nsimEIM)
          wz[, 2*(1:NOS)] = if (intercept.only)
              matrix(colMeans(run.varcov),
                     n, ncol(run.varcov), byrow = TRUE) else run.varcov

          wz[, 2*(1:NOS)] = wz[, 2*(1:NOS)] * dk.deta^2
      }
      ed2l.dmu2 = 1/mu - 1/(mu+kmat)
      wz[, 2*(1:NOS)-1] = dmu.deta^2 * ed2l.dmu2
      c(w) * wz
  }), list( .cutoff = cutoff, .Maxiter = Maxiter,
            .nsimEIM = nsimEIM ))))



  if (deviance.arg) ans@deviance = eval(substitute(
      function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
      NOS = ncol(eta) / 2
      temp300 =  eta[,2*(1:NOS), drop = FALSE]
      if ( .lsize == "loge") {
          bigval = 28
          temp300[temp300 >  bigval] =  bigval
          temp300[temp300 < -bigval] = -bigval
      } else stop("can only handle the 'loge' link")
      k =  eta2theta(temp300, .lsize, earg = .esize)
      devi = 2 * (y*log(ifelse(y < 1, 1, y)/mu) + (y+k)*log((mu+k)/(k+y)))
      if (residuals)
         sign(y - mu) * sqrt(abs(devi) * w) else
         sum(w * devi)
  }, list( .lsize = lsize, .emu = emu, .esize = esize )))

  ans
}







polya.control <- function(save.weight = TRUE, ...)
{
    list(save.weight = save.weight)
}



 polya <-
  function(lprob = "logit", lsize = "loge",
           eprob = list(),  esize = list(),
           iprob = NULL,    isize = NULL,
           quantile.probs = 0.75,
           nsimEIM = 100,
           deviance.arg = FALSE, imethod = 1,
           shrinkage.init = 0.95, zero = -2)
{



  if (length(iprob) && !is.Numeric(iprob, posit = TRUE))
    stop("bad input for argument 'iprob'")
  if (length(isize) && !is.Numeric(isize, posit = TRUE))
    stop("bad input for argument 'isize'")

  if (!is.Numeric(imethod, allow = 1, integ = TRUE, posit = TRUE) ||
     imethod > 3)
     stop("argument 'imethod' must be 1 or 2 or 3")
  if (!is.Numeric(shrinkage.init, allow = 1) || shrinkage.init < 0 ||
     shrinkage.init > 1)
     stop("bad input for argument 'shrinkage.init'")

  if (!is.Numeric(nsimEIM, allow = 1, integ = TRUE))
    stop("bad input for argument 'nsimEIM'")
  if (nsimEIM <= 10)
    warning("argument 'nsimEIM' should be an integer ",
            "greater than 10, say")

  if (mode(lprob) != "character" && mode(lprob) != "name")
    lprob = as.character(substitute(lprob))
  if (mode(lsize) != "character" && mode(lsize) != "name")
    lsize = as.character(substitute(lsize))
  if (!is.list(eprob)) eprob = list()
  if (!is.list(esize)) esize = list()


  ans = 
  new("vglmff",
  blurb = c("Polya (negative-binomial) distribution\n\n",
            "Links:    ",
            namesof("prob", lprob, earg = eprob), ", ",
            namesof("size", lsize, earg = esize), "\n",
            "Mean:     size * (1 - prob) / prob\n",
            "Variance: mean / prob"),
  constraints = eval(substitute(expression({

    dotzero <- .zero
    Musual <- 2
    eval(negzero.expression)

  }), list( .zero = zero ))),

  infos = eval(substitute(function(...) {
    list(Musual = 2,
         zero = .zero)
  }, list( .zero = zero ))),

  initialize = eval(substitute(expression({
    Musual = 2
    if (any(function.name == c("cqo", "cao")))
      stop("polya() does not work with cqo() or cao(). ",
           "Try negbinomial()")

    if (any(y < 0))
      stop("negative values not allowed for the 'polya' family")
    if (any(round(y) != y))
      stop("integer-values only allowed for the 'polya' family")

    y = as.matrix(y)
    M = 2 * ncol(y)
    NOS = ncoly = ncol(y)  # Number of species

    predictors.names =
      c(namesof(if (NOS == 1) "prob" else
                paste("prob", 1:NOS, sep = ""),
               .lprob, earg = .eprob, tag = FALSE),
        namesof(if (NOS == 1) "size" else
                paste("size", 1:NOS, sep = ""),
               .lsize,  earg = .esize,  tag = FALSE))
    predictors.names = predictors.names[interleave.VGAM(M, M = 2)]

    if (is.null( .nsimEIM )) {
       save.weight <- control$save.weight <- FALSE
    }

    
    PROB.INIT <- if (is.numeric( .pinit )) {
      matrix( .pinit, nrow(y), ncol(y), byrow = TRUE)
    } else {
      NULL
    }

    if (!length(etastart)) {
      mu.init = y
      for(iii in 1:ncol(y)) {
        use.this = if ( .imethod == 1) {
          weighted.mean(y[, iii], w) + 1/16
        } else if ( .imethod == 3) {
          c(quantile(y[, iii], probs = .quantile.probs) + 1/16)
        } else {
          median(y[, iii]) + 1/16
        }

        if (FALSE) {
          mu.init[, iii] = MU.INIT[, iii]
        } else {
          medabsres = median(abs(y[, iii] - use.this)) + 1/32
          allowfun = function(z, maxtol = 1) sign(z) * pmin(abs(z), maxtol)
          mu.init[, iii] = use.this + (1 - .sinit) * allowfun(y[, iii] -
                          use.this, maxtol = medabsres)

          mu.init[, iii] = abs(mu.init[, iii]) + 1 / 1024
        }
      }



      if ( is.Numeric( .kinit )) {
        kayy.init = matrix( .kinit, nrow = n, ncol = NOS, byrow = TRUE)
      } else {
        negbinomial.Loglikfun = function(kmat, y, x, w, extraargs) {
            mu = extraargs
            sum(w * dnbinom(x = y, mu = mu, size = kmat, log = TRUE))
        }
        k.grid = 2^((-7):7)
        k.grid = 2^(seq(-8, 8, length = 40))
        kayy.init = matrix(0, nrow = n, ncol = NOS)
        for(spp. in 1:NOS) {
          kayy.init[,spp.] = getMaxMin(k.grid,
                             objfun=negbinomial.Loglikfun,
                             y = y[,spp.], x = x, w = w,
                             extraargs= mu.init[,spp.])
        }
      }

      prob.init = if (length(PROB.INIT)) PROB.INIT else
                  kayy.init / (kayy.init + mu.init)


      etastart = cbind(theta2eta(prob.init, .lprob, earg = .eprob),
                       theta2eta(kayy.init, .lsize, earg = .esize))
      etastart = etastart[, interleave.VGAM(M, M = Musual), drop = FALSE]
      }
  }), list( .lprob = lprob, .lsize = lsize,
            .eprob = eprob, .esize = esize,
            .pinit = iprob, .kinit = isize,
            .quantile.probs = quantile.probs,
            .sinit = shrinkage.init, .nsimEIM = nsimEIM, .zero = zero,
            .imethod = imethod ))),
  inverse = eval(substitute(function(eta, extra = NULL) {
    Musual = 2
    NOS = ncol(eta) / Musual
    pmat = eta2theta(eta[, Musual*(1:NOS) - 1, drop = FALSE],
                     .lprob, earg = .eprob)
    kmat = eta2theta(eta[, Musual*(1:NOS)-  0, drop = FALSE],
                     .lsize, earg = .esize)
    kmat / (kmat + pmat)
  }, list( .lprob = lprob, .eprob = eprob,
           .lsize = lsize, .esize = esize ))),
  last = eval(substitute(expression({
    temp0303 = c(rep( .lprob, length = NOS),
                 rep( .lsize, length = NOS))
    names(temp0303) =
      c(if (NOS == 1) "prob" else paste("prob", 1:NOS, sep = ""),
        if (NOS == 1) "size" else paste("size", 1:NOS, sep = ""))
    temp0303 = temp0303[interleave.VGAM(M, M = 2)]
    misc$link = temp0303 # Already named

    misc$earg = vector("list", M)
    names(misc$earg) = names(misc$link)
    for(ii in 1:NOS) {
        misc$earg[[2*ii-1]] = .eprob
        misc$earg[[2*ii  ]] = .esize
    }

    misc$isize = .isize  
    misc$imethod = .imethod 
    misc$nsimEIM = .nsimEIM
    misc$expected = TRUE
    misc$shrinkage.init = .sinit
    misc$Musual = 2
  }), list( .lprob = lprob, .lsize = lsize,
            .eprob = eprob, .esize = esize,
            .isize = isize,
            .nsimEIM = nsimEIM,
            .sinit = shrinkage.init, .imethod = imethod ))),


  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
    Musual = 2
    NOS = ncol(eta) / Musual
    pmat  = eta2theta(eta[, Musual*(1:NOS) - 1, drop = FALSE],
                      .lprob, earg = .eprob)
    temp300 =         eta[, Musual*(1:NOS)    , drop = FALSE]
    if ( .lsize == "loge") {
      bigval = 28
      temp300 = ifelse(temp300 >  bigval,  bigval, temp300)
      temp300 = ifelse(temp300 < -bigval, -bigval, temp300)
    }
    kmat = eta2theta(temp300, .lsize, earg = .esize)
    if (residuals) stop("loglikelihood residuals not ",
                        "implemented yet") else
      sum(w * dnbinom(x = y, prob = pmat, size = kmat, log = TRUE))
  }, list( .lsize = lsize, .lprob = lprob,
           .esize = esize, .eprob = eprob ))),
  vfamily = c("polya"),
  deriv = eval(substitute(expression({
    Musual = 2
    NOS = ncol(eta) / Musual
    M = ncol(eta)

    pmat  = eta2theta(eta[, Musual*(1:NOS) - 1, drop = FALSE],
                      .lprob, earg = .eprob)
    temp3 =           eta[, Musual*(1:NOS)    , drop = FALSE]
    if ( .lsize == "loge") {
      bigval = 28
      temp3 = ifelse(temp3 >  bigval,  bigval, temp3)
      temp3 = ifelse(temp3 < -bigval, -bigval, temp3)
    }
    kmat = eta2theta(temp3, .lsize, earg = .esize)

    dl.dprob = kmat / pmat - y / (1.0 - pmat)
    dl.dkayy = digamma(y + kmat) - digamma(kmat) + log(pmat)

    dprob.deta = dtheta.deta(pmat, .lprob, earg = .eprob)
    dkayy.deta = dtheta.deta(kmat, .lsize, earg = .esize)
    dthetas.detas = cbind(dprob.deta, dkayy.deta)
    dThetas.detas = dthetas.detas[, interleave.VGAM(M, M = Musual)]
    myderiv = c(w) * cbind(dl.dprob, dl.dkayy) * dthetas.detas
    myderiv[, interleave.VGAM(M, M = Musual)]
  }), list( .lprob = lprob, .lsize = lsize,
            .eprob = eprob, .esize = esize ))),
  weight = eval(substitute(expression({
    wz = matrix(0.0, n, M + M - 1)  # wz is 'tridiagonal' 

    ind1 = iam(NA, NA, M = Musual, both = TRUE, diag = TRUE)
    mumat = as.matrix(mu)


    for(spp. in 1:NOS) {
      run.varcov = 0
      kvec = kmat[, spp.]
      pvec = pmat[, spp.]

      for(ii in 1:( .nsimEIM )) {
        ysim = rnbinom(n = n, prob = pvec, size = kvec)

        dl.dprob = kvec / pvec - ysim / (1.0 - pvec)
        dl.dkayy = digamma(ysim + kvec) - digamma(kvec) + log(pvec)
        temp3 = cbind(dl.dprob, dl.dkayy)
        run.varcov = run.varcov +
                     temp3[, ind1$row.index] *
                     temp3[, ind1$col.index]
      }
      run.varcov = cbind(run.varcov / .nsimEIM)

      wz1 = if (intercept.only)
          matrix(colMeans(run.varcov),
                 nrow = n, ncol = ncol(run.varcov), byrow = TRUE) else
          run.varcov

      wz1 = wz1 * dThetas.detas[, Musual * (spp. - 1) + ind1$row] *
                  dThetas.detas[, Musual * (spp. - 1) + ind1$col]


      for(jay in 1:Musual)
          for(kay in jay:Musual) {
              cptr = iam((spp. - 1) * Musual + jay,
                         (spp. - 1) * Musual + kay,
                         M = M)
              wz[, cptr] = wz1[, iam(jay, kay, M = Musual)]
          }
    } # End of for(spp.) loop


    c(w) * wz
  }), list( .nsimEIM = nsimEIM ))))




  if (deviance.arg) ans@deviance = eval(substitute(
      function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
    NOS = ncol(eta) / 2
    temp300 =  eta[, 2*(1:NOS), drop = FALSE]
    if ( .lsize == "loge") {
      bigval = 28
      temp300[temp300 >  bigval] =  bigval
      temp300[temp300 < -bigval] = -bigval
    } else {
      stop("can only handle the 'loge' link")
    }
    kayy =  eta2theta(temp300, .lsize, earg = .esize)
    devi = 2 * (y * log(ifelse(y < 1, 1, y) / mu) +
               (y + kayy) * log((mu + kayy) / (kayy + y)))
    if (residuals)
      sign(y - mu) * sqrt(abs(devi) * w) else
      sum(w * devi)
    }, list( .lsize = lsize, .eprob = eprob,
             .esize = esize )))

    ans
} # End of polya()




 simple.poisson = function()
{
    new("vglmff",
    blurb = c("Poisson distribution\n\n",
            "Link:     log(lambda)",
            "\n",
            "Variance: lambda"),
    deviance= function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
        nz = y > 0
        devi =  - (y - mu)
        devi[nz] = devi[nz] + y[nz] * log(y[nz]/mu[nz])
        if (residuals) sign(y - mu) * sqrt(2 * abs(devi) * w) else
            2 * sum(w * devi)
    },
    initialize=expression({
        if (ncol(cbind(y)) != 1)
            stop("response must be a vector or a one-column matrix")
        predictors.names = "log(lambda)"
        mu = y + 0.167 * (y == 0)
        if (!length(etastart))
            etastart = log(mu)
    }), 
    inverse=function(eta, extra = NULL)
        exp(eta),
    last = expression({
        misc$link = c(lambda = "loge")
        misc$earg = list(lambda = list())
    }),
    link=function(mu, extra = NULL)
        log(mu),
    vfamily = "simple.poisson",
    deriv=expression({
        lambda = mu
        dl.dlambda = -1 + y/lambda
        dlambda.deta = dtheta.deta(theta=lambda, link = "loge", earg = list())
        c(w) * dl.dlambda * dlambda.deta
    }),
    weight = expression({
        d2l.dlambda2 = 1 / lambda
        c(w) * d2l.dlambda2 * dlambda.deta^2
    }))
}













 studentt <-  function(ldf = "loglog", edf = list(), idf = NULL,
                       tol1 = 0.1,
                       imethod = 1)
{

  ldof <- ldf
  edof <- edf
  idof <- idf

  if (mode(ldof) != "character" && mode(ldof) != "name")
    ldof <- as.character(substitute(ldof))
  if (!is.list(edof)) edof <- list()

  if (length(idof))
    if (!is.Numeric(idof) || any(idof <= 1))
      stop("argument 'idf' should be > 1")

  if (!is.Numeric(tol1, posit  = TRUE))
    stop("argument 'tol1' should be positive")

  if (!is.Numeric(imethod, allow = 1, integ = TRUE, posit = TRUE) ||
     imethod > 3)
      stop("argument 'imethod' must be 1 or 2 or 3")


  new("vglmff",
  blurb = c("Student t-distribution\n\n",
            "Link:     ",
            namesof("df", ldof, earg = edof), "\n",
            "Variance: df/(df-2) if df > 2\n"),
  infos = eval(substitute(function(...) {
    list(Musual = 1,
         tol1 = .tol1 )
  }, list( .tol1 = tol1 ))),
  initialize = eval(substitute(expression({
    if (ncol(cbind(y)) != 1)
        stop("response must be a vector or a one-column matrix")

    predictors.names <- namesof("df", .ldof, earg = .edof, tag = FALSE)

    if (!length(etastart)) {

      init.df <- if (length( .idof )) .idof else {
        VarY = var(y)
        MadY = mad(y)
        if (VarY <= (1 + .tol1 )) VarY = 1.12
        if ( .imethod == 1) {
          2 * VarY / (VarY - 1)
        } else if ( .imethod == 2) {
          ifelse(MadY < 1.05, 30, ifelse(MadY > 1.2, 2, 5))
        } else
          10
      }


      etastart <- rep(theta2eta(init.df, .ldof , earg = .edof ),
                      len = length(y))
    }
  }), list( .ldof = ldof, .edof = edof, .idof = idof,
            .tol1 = tol1, .imethod = imethod ))), 
  inverse = eval(substitute(function(eta, extra = NULL) {
    Dof <- eta2theta(eta, .ldof, earg = .edof)
    ans <- 0 * eta
    ans[Dof <= 1] <- NA
    ans
  }, list( .ldof = ldof, .edof = edof ))),
  last = eval(substitute(expression({
    misc$link <-    c(df = .ldof )
    misc$earg <- list(df = .edof )
    misc$imethod <- .imethod
    misc$expected = TRUE
  }), list( .ldof = ldof,
            .edof = edof, .imethod = imethod ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
    Dof <-  eta2theta(eta, .ldof, earg = .edof)
    if (residuals) stop("loglikelihood residuals not ",
                        "implemented yet") else {
        sum(w * dt(x = y, df = Dof, log = TRUE))
    }
  }, list( .ldof = ldof, .edof = edof ))), 
  vfamily = c("studentt"),
  deriv = eval(substitute(expression({
    Dof <- eta2theta(eta, .ldof, earg = .edof)
    ddf.deta <-  dtheta.deta(theta = Dof, .ldof, earg = .edof)

    DDS  <- function(df)          digamma((df + 1) / 2) -  digamma(df / 2)
    DDSp <- function(df)  0.5 * (trigamma((df + 1) / 2) - trigamma(df / 2))

    temp0 <- 1 / Dof
    temp1 <-  temp0 * y^2
    dl.ddf <- 0.5 * (-temp0 - log1p(temp1) +
              (Dof + 1) * y^2 / (Dof^2 * (1 + temp1)) + DDS(Dof))
    c(w) * dl.ddf * ddf.deta
  }), list( .ldof = ldof, .edof = edof ))),
  weight = eval(substitute(expression({

    const2 = (Dof + 0) / (Dof + 3)
    const2[!is.finite(Dof)] <- 1  # Handles Inf

    tmp6 = DDS(Dof)
    edl2.dnu2 <- 0.5 * (tmp6 * (const2 * tmp6 - 2 / (Dof + 1)) - DDSp(Dof))
 
    wz <- c(w) * edl2.dnu2 * ddf.deta^2
    wz
  }), list( .ldof = ldof, .edof = edof ))))
}






    Kayfun.studentt <- function(df, bigno = .Machine$double.eps^(-0.46)) {
      ind1 <- is.finite(df)

      const4 = dnorm(0)
      ans <- df

      if (any(ind1))
        ans[ind1] <- exp(lgamma((df[ind1] + 1) / 2) -
                         lgamma( df[ind1]      / 2)) / sqrt(pi * df[ind1])
      ans[df <= 0] = NaN
      ind2 <- (df >= bigno)
      if (any(ind2)) {
        dff = df[ind2]
        ans[ind2] <- const4 # 1 / const3  # for handling df = Inf
      }
      ans[!ind1] <- const4 # 1 / const3  # for handling df = Inf

      ans
    }




 studentt3 <- function(llocation = "identity", elocation = list(),
                       lscale    = "loge",     escale   = list(),
                       ldf       = "loglog",   edf      = list(),
                       ilocation = NULL, iscale = NULL, idf = NULL,
                       imethod = 1,
                       zero = -(2:3))
{



  lloc <- llocation; lsca <- lscale; ldof <- ldf
  eloc <- elocation; esca <- escale; edof <- edf
  iloc <- ilocation; isca <- iscale; idof <- idf

  if (mode(lloc) != "character" && mode(lloc) != "name")
    lloc <- as.character(substitute(lloc))
  if (!is.list(eloc)) eloc <- list()

  if (mode(lsca) != "character" && mode(lsca) != "name")
    lsca <- as.character(substitute(lsca))
  if (!is.list(esca)) esca <- list()

  if (mode(ldof) != "character" && mode(ldof) != "name")
    ldof <- as.character(substitute(ldof))
  if (!is.list(edof)) edof <- list()

  if (!is.Numeric(imethod, allow = 1, integ = TRUE, posit = TRUE) ||
     imethod > 3)
      stop("argument 'imethod' must be 1 or 2 or 3")

  if (length(iloc))
    if (!is.Numeric(iloc))
      stop("bad input in argument 'ilocation'")
  if (length(isca))
    if (!is.Numeric(isca, posit = TRUE))
      stop("argument 'iscale' should be positive")
  if (length(idof))
    if (!is.Numeric(idof) || any(idof <= 1))
      stop("argument 'idf' should be > 1")

  new("vglmff",
  blurb = c("Student t-distribution\n\n",
            "Link:     ",
            namesof("location", lloc, earg = eloc), ", ",
            namesof("scale",    lsca, earg = esca), ", ",
            namesof("df",       ldof, earg = edof), "\n",
            "Variance: scale^2 * df / (df - 2) if df > 2\n"),
  constraints = eval(substitute(expression({

    dotzero <- .zero
    Musual <- 3
    eval(negzero.expression)
  }), list( .zero = zero ))),
    infos = eval(substitute(function(...) {
      list(Musual = 3,
           zero = .zero)
    }, list( .zero = zero ))),
  initialize = eval(substitute(expression({
    Musual <- 3
    if (ncol(cbind(w)) != 1)
      stop("prior weights must be a vector or a one-column matrix")

    y <- as.matrix(y)
    extra$NOS <- NOS <- ncoly <- ncol(y)  # Number of species
    extra$Musual <- Musual
    M <- Musual * ncoly #

    mynames1 <- paste("location", if (NOS > 1) 1:NOS else "", sep = "")
    mynames2 <- paste("scale",    if (NOS > 1) 1:NOS else "", sep = "")
    mynames3 <- paste("df",       if (NOS > 1) 1:NOS else "", sep = "")
    predictors.names <-
        c(namesof(mynames1, .lloc, earg = .eloc, tag = FALSE),
          namesof(mynames2, .lsca, earg = .esca, tag = FALSE),
          namesof(mynames3, .ldof, earg = .edof, tag = FALSE))
    predictors.names <-
      predictors.names[interleave.VGAM(Musual * NOS, M = Musual)]

    if (!length(etastart)) {

      init.loc <- if (length( .iloc )) .iloc else {
        if ( .imethod == 2) apply(y, 2, median) else
        if ( .imethod == 3) y else {
           colSums(w * y) / sum(w)
        }
      }

      sdvec <- apply(y, 2, sd)
      init.sca <- if (length( .isca )) .isca else
                  sdvec / 2.3

      sdvec    <- rep(sdvec,    len = max(length(sdvec), length(init.sca)))
      init.sca <- rep(init.sca, len = max(length(sdvec), length(init.sca)))
      ind9 <- (sdvec / init.sca <= (1 + 0.12))
      sdvec[ind9] <- sqrt(1.12) * init.sca[ind9]
      init.dof <- if (length( .idof )) .idof else
                  (2 * (sdvec / init.sca)^2) / ((sdvec / init.sca)^2  - 1)
      if (!is.Numeric(init.dof) || init.dof <= 1)
        init.dof <- rep(3, len = ncoly)

      mat1 <- matrix(theta2eta(init.loc, .lloc, earg = .eloc), n, NOS,
                     byrow = TRUE)
      mat2 <- matrix(theta2eta(init.sca, .lsca, earg = .esca), n, NOS,
                     byrow = TRUE)
      mat3 <- matrix(theta2eta(init.dof, .ldof, earg = .edof), n, NOS,
                     byrow = TRUE)
      etastart <- cbind(mat1, mat2, mat3)
      etastart <- etastart[, interleave.VGAM(ncol(etastart), M = Musual)]
    }
  }), list( .lloc = lloc, .eloc = eloc, .iloc = iloc,
            .lsca = lsca, .esca = esca, .isca = isca,
            .ldof = ldof, .edof = edof, .idof = idof,
            .imethod = imethod ))), 
  inverse = eval(substitute(function(eta, extra = NULL) {
    NOS    <- extra$NOS
    Musual <- extra$Musual
    Loc <-  eta2theta(eta[, Musual*(1:NOS)-2], .lloc, earg = .eloc)
    Dof <-  eta2theta(eta[, Musual*(1:NOS)-0], .ldof, earg = .edof)
    Loc[Dof <= 1] <- NA
    Loc
  }, list( .lloc = lloc, .eloc = eloc,
           .lsca = lsca, .esca = esca,
           .ldof = ldof, .edof = edof ))),
  last = eval(substitute(expression({
    Musual <- extra$Musual
    misc$link <- c(rep( .lloc, length = NOS),
                   rep( .lsca, length = NOS),
                   rep( .ldof, length = NOS))
    misc$link <- misc$link[interleave.VGAM(Musual * NOS, M = Musual)]
    temp.names <- c(mynames1, mynames2, mynames3)
    temp.names <- temp.names[interleave.VGAM(Musual * NOS, M = Musual)]
    names(misc$link) <- temp.names

    misc$earg <- vector("list", Musual * NOS)
    names(misc$earg) <- temp.names
    for(ii in 1:NOS) {
        misc$earg[[Musual*ii-2]] <- .eloc
        misc$earg[[Musual*ii-1]] <- .esca
        misc$earg[[Musual*ii  ]] <- .edof
    }
 
    misc$Musual <- Musual
    misc$imethod <- .imethod
    misc$expected = TRUE
  }), list( .lloc = lloc, .eloc = eloc,
            .lsca = lsca, .esca = esca,
            .ldof = ldof, .edof = edof,
            .imethod = imethod ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
    NOS <- extra$NOS
    Musual <- extra$Musual
    Loc <-  eta2theta(eta[, Musual*(1:NOS)-2], .lloc, earg = .eloc)
    Sca <-  eta2theta(eta[, Musual*(1:NOS)-1], .lsca, earg = .esca)
    Dof <-  eta2theta(eta[, Musual*(1:NOS)-0], .ldof, earg = .edof)
    zedd <- (y - Loc) / Sca
    if (residuals) stop("loglikelihood residuals not ",
                        "implemented yet") else {
        sum(w * (dt(x = zedd, df = Dof, log = TRUE) - log(Sca)))
    }
  }, list(  .lloc = lloc, .eloc = eloc,
            .lsca = lsca, .esca = esca,
            .ldof = ldof, .edof = edof ))), 
  vfamily = c("studentt3"),
  deriv = eval(substitute(expression({
    Musual <- extra$Musual
    NOS <- extra$NOS
    Loc <- eta2theta(eta[, Musual*(1:NOS)-2], .lloc, earg = .eloc)
    Sca <- eta2theta(eta[, Musual*(1:NOS)-1], .lsca, earg = .esca)
    Dof <- eta2theta(eta[, Musual*(1:NOS)-0], .ldof, earg = .edof)

    dloc.deta <- cbind(dtheta.deta(theta = Loc, .lloc, earg = .eloc))
    dsca.deta <- cbind(dtheta.deta(theta = Sca, .lsca, earg = .esca))
    ddof.deta <- cbind(dtheta.deta(theta = Dof, .ldof, earg = .edof))

    zedd  <- (y - Loc) / Sca
    temp0 <- 1 / Dof
    temp1 <- temp0 * zedd^2
    dl.dloc <- (Dof + 1) * zedd / (Sca * (Dof + zedd^2))
    dl.dsca <- zedd * dl.dloc - 1 / Sca
    dl.ddof <- 0.5 * (-temp0 - log1p(temp1) +
                     (Dof+1) * zedd^2 / (Dof^2 * (1 + temp1)) +
                     digamma((Dof+1)/2) - digamma(Dof/2))
 
    ans <- c(w) * cbind(dl.dloc * dloc.deta,
                        dl.dsca * dsca.deta,
                        dl.ddof * ddof.deta)
    ans <- ans[, interleave.VGAM(ncol(ans), M = Musual)]
    ans
  }), list( .lloc = lloc, .eloc = eloc,
            .lsca = lsca, .esca = esca,
            .ldof = ldof, .edof = edof ))),
  weight = eval(substitute(expression({

    const1 = (Dof + 1) / (Dof + 3)
    const2 = (Dof + 0) / (Dof + 3)
    const1[!is.finite(Dof)] <- 1  # Handles Inf
    const2[!is.finite(Dof)] <- 1  # Handles Inf

    const4 = dnorm(0)
    ed2l.dlocat2 =      const1 / (Sca * (Kayfun.studentt(Dof) / const4))^2
    ed2l.dscale2 = 2  * const2 /  Sca^2

    DDS  <- function(df)          digamma((df + 1) / 2) -  digamma(df / 2)
    DDSp <- function(df)  0.5 * (trigamma((df + 1) / 2) - trigamma(df / 2))


    tmp6 = DDS(Dof)
    edl2.dnu2 <- 0.5 * (tmp6 * (const2 * tmp6 - 2 / (Dof + 1)) - DDSp(Dof))
    ed2l.dshape2 <- cbind(edl2.dnu2)  # cosmetic name change

    ed2l.dshape.dlocat = cbind(0 * Sca)
    ed2l.dshape.dscale = cbind((-1 / (Dof + 1.0) + const2 * DDS(Dof)) / Sca)

    wz = matrix(0.0, n, dimm(M))
    wz[, Musual*(1:NOS) - 2] = ed2l.dlocat2 * dloc.deta^2
    wz[, Musual*(1:NOS) - 1] = ed2l.dscale2 * dsca.deta^2
    wz[, Musual*(1:NOS) - 0] = ed2l.dshape2 * ddof.deta^2

    for (ii in ((1:NOS) - 1)) {
      ind3 = 1 + ii
      wz[, iam(ii*Musual + 1, ii*Musual + 3, M = M)] <-
           ed2l.dshape.dlocat[, ind3] *
           dloc.deta[, ind3] * ddof.deta[, ind3]
      wz[, iam(ii*Musual + 2, ii*Musual + 3, M = M)] <-
           ed2l.dshape.dscale[, ind3] *
           dsca.deta[, ind3] * ddof.deta[, ind3]
    }

  while (all(wz[, ncol(wz)] == 0))
    wz <- wz[, -ncol(wz)]

    c(w) * wz
  }), list( .lloc = lloc, .eloc = eloc,
            .lsca = lsca, .esca = esca,
            .ldof = ldof, .edof = edof ))))
}





 studentt2 <- function(df = Inf,
                       llocation = "identity", elocation = list(),
                       lscale    = "loge",     escale   = list(),
                       ilocation = NULL, iscale = NULL,
                       imethod = 1,
                       zero = -2)
{



  lloc <- llocation; lsca <- lscale
  eloc <- elocation; esca <- escale
  iloc <- ilocation; isca <- iscale
  doff <- df

  if (mode(lloc) != "character" && mode(lloc) != "name")
    lloc <- as.character(substitute(lloc))
  if (!is.list(eloc)) eloc <- list()

  if (mode(lsca) != "character" && mode(lsca) != "name")
    lsca <- as.character(substitute(lsca))
  if (!is.list(esca)) esca <- list()

  if (is.finite(doff))
    if (!is.Numeric(doff, posit = TRUE))
    stop("argument 'df' must be positive")

  if (!is.Numeric(imethod, allow = 1, integ = TRUE, posit = TRUE) ||
     imethod > 3)
      stop("argument 'imethod' must be 1 or 2 or 3")

  if (length(iloc))
    if (!is.Numeric(iloc))
      stop("bad input in argument 'ilocation'")
  if (length(isca))
    if (!is.Numeric(isca, posit = TRUE))
      stop("argument 'iscale' should be positive")


  new("vglmff",
  blurb = c("Student t-distribution\n\n",
            "Link:     ",
            namesof("location", lloc, earg = eloc), ", ",
            namesof("scale",    lsca, earg = esca), "\n",
            "Variance: scale^2 * df / (df - 2) if df > 2\n"),
  constraints = eval(substitute(expression({

    dotzero <- .zero
    Musual <- 2
    eval(negzero.expression)
  }), list( .zero = zero ))),
    infos = eval(substitute(function(...) {
      list(Musual = 2,
           zero = .zero)
    }, list( .zero = zero ))),
  initialize = eval(substitute(expression({
    Musual <- 2
    if (ncol(cbind(w)) != 1)
      stop("prior weights must be a vector or a one-column matrix")

    y <- as.matrix(y)
    extra$NOS <- NOS <- ncoly <- ncol(y)  # Number of species
    extra$Musual <- Musual
    M <- Musual * ncoly #

    mynames1 <- paste("location", if (NOS > 1) 1:NOS else "", sep = "")
    mynames2 <- paste("scale",    if (NOS > 1) 1:NOS else "", sep = "")
    predictors.names <-
        c(namesof(mynames1, .lloc, earg = .eloc, tag = FALSE),
          namesof(mynames2, .lsca, earg = .esca, tag = FALSE))
    predictors.names <-
      predictors.names[interleave.VGAM(Musual * NOS, M = Musual)]

    if (!length(etastart)) {

      init.loc <- if (length( .iloc )) .iloc else {
        if ( .imethod == 2) apply(y, 2, median) else
        if ( .imethod == 3) y else {
           colSums(w * y) / sum(w)
        }
      }

      sdvec <- apply(y, 2, sd)
      init.sca <- if (length( .isca )) .isca else
                  sdvec / 2.3

      mat1 <- matrix(theta2eta(init.loc, .lloc, earg = .eloc), n, NOS,
                     byrow = TRUE)
      mat2 <- matrix(theta2eta(init.sca, .lsca, earg = .esca), n, NOS,
                     byrow = TRUE)
      etastart <- cbind(mat1, mat2)
      etastart <- etastart[, interleave.VGAM(ncol(etastart), M = Musual)]
    }
  }), list( .lloc = lloc, .eloc = eloc, .iloc = iloc,
            .lsca = lsca, .esca = esca, .isca = isca,
            .doff = doff,
            .imethod = imethod ))), 
  inverse = eval(substitute(function(eta, extra = NULL) {
    NOS <- extra$NOS
    Musual <- extra$Musual
    Loc <-  eta2theta(eta[, Musual*(1:NOS) - 1], .lloc, earg = .eloc)
    Dof <- matrix( .doff , nrow(cbind(Loc)), NOS, byrow = TRUE)
    Loc[Dof <= 1] <- NA
    Loc
  }, list( .lloc = lloc, .eloc = eloc,
           .lsca = lsca, .esca = esca,
           .doff = doff ))),
  last = eval(substitute(expression({
    Musual <- extra$Musual
    misc$link <- c(rep( .lloc, length = NOS),
                   rep( .lsca, length = NOS))
    temp.names <- c(mynames1, mynames2)
    temp.names <- temp.names[interleave.VGAM(Musual * NOS, M = Musual)]
    names(misc$link) <- temp.names
    misc$earg <- vector("list", Musual * NOS)
    names(misc$earg) <- temp.names
    for(ii in 1:NOS) {
        misc$earg[[Musual*ii-1]] <- .eloc
        misc$earg[[Musual*ii-0]] <- .esca
    }
 
    misc$Musual <- Musual
    misc$simEIM <- TRUE
    misc$df <- .doff
    misc$imethod <- .imethod
    misc$expected = TRUE
  }), list( .lloc = lloc, .eloc = eloc,
            .lsca = lsca, .esca = esca,
            .doff = doff,
            .imethod = imethod ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
    NOS <- extra$NOS
    Musual <- extra$Musual
    Loc <- eta2theta(eta[, Musual*(1:NOS)-1], .lloc, earg = .eloc)
    Sca <- eta2theta(eta[, Musual*(1:NOS)-0], .lsca, earg = .esca)
    Dof <- matrix( .doff , nrow(cbind(Loc)), NOS, byrow = TRUE)
    zedd <- (y - Loc) / Sca
    if (residuals) stop("loglikelihood residuals not ",
                        "implemented yet") else {
        sum(w * (dt(x = zedd, df = Dof, log = TRUE) - log(Sca)))
    }
  }, list(  .lloc = lloc, .eloc = eloc,
            .lsca = lsca, .esca = esca,
            .doff = doff ))), 
  vfamily = c("studentt2"),
  deriv = eval(substitute(expression({
    Musual <- extra$Musual
    NOS <- extra$NOS
    Loc <- eta2theta(eta[, Musual*(1:NOS)-1], .lloc, earg = .eloc)
    Sca <- eta2theta(eta[, Musual*(1:NOS)-0], .lsca, earg = .esca)
    Dof <- matrix( .doff , n, NOS, byrow = TRUE)

    dlocat.deta <- dtheta.deta(theta = Loc, .lloc, earg = .eloc)
    dscale.deta <- dtheta.deta(theta = Sca, .lsca, earg = .esca)

    zedd  <- (y - Loc) / Sca
    temp0 <- 1 / Dof
    temp1 <- temp0 * zedd^2
    dl.dlocat <- (Dof + 1) * zedd / (Sca * (Dof + zedd^2))
    dl.dlocat[!is.finite(Dof)] <- zedd / Sca  # Adjust for df=Inf
    dl.dscale <- zedd * dl.dlocat - 1 / Sca
 
    ans <- c(w) * cbind(dl.dlocat * dlocat.deta,
                        dl.dscale * dscale.deta)
    ans <- ans[, interleave.VGAM(ncol(ans), M = Musual)]
    ans
  }), list( .lloc = lloc, .eloc = eloc,
            .lsca = lsca, .esca = esca,
            .doff = doff ))),
  weight = eval(substitute(expression({

    const1 = (Dof + 1) / (Dof + 3)
    const2 = (Dof + 0) / (Dof + 3)
    const1[!is.finite( Dof )] <- 1  # Handles Inf
    const2[!is.finite( Dof )] <- 1  # Handles Inf

    const4 = dnorm(0)
    ed2l.dlocat2 =        const1 / (Sca * (Kayfun.studentt(Dof) / const4))^2

    ed2l.dscale2 = 2.0  * const2 /  Sca^2                 # 2.0 seems to work

    wz = matrix(as.numeric(NA), n, M)  #2=M; diagonal!
    wz[, Musual*(1:NOS) - 1] = ed2l.dlocat2 * dlocat.deta^2
    wz[, Musual*(1:NOS)    ] = ed2l.dscale2 * dscale.deta^2
    c(w) * wz
  }), list( .lloc = lloc, .eloc = eloc,
            .lsca = lsca, .esca = esca,
            .doff = doff  ))))
}





 
 chisq <- function(link = "loge", earg = list())
{
  if (mode(link) != "character" && mode(link) != "name")
      link <- as.character(substitute(link))
  if (!is.list(earg)) earg <- list()

  new("vglmff",
  blurb = c("Chi-squared distribution\n\n",
            "Link:     ",
            namesof("df", link, earg = earg, tag = FALSE)),
  initialize = eval(substitute(expression({
    if (ncol(cbind(w)) != 1)
      stop("argument 'weights' must be a vector or a one-column matrix")

    y <- as.matrix(y)

    extra$ncoly <- NOS <- ncol(y)  # Number of species
    mynames1 <- paste("df", if (NOS > 1) 1:NOS else "", sep = "")
    predictors.names <- namesof(mynames1, .link, earg = .earg, tag = FALSE)

    if (!length(mustart) && !length(etastart))
      mustart <- y + (1 / 8) * (y == 0)
  }), list( .link = link, .earg = earg ))),
  inverse = eval(substitute(function(eta, extra = NULL) {
    eta2theta(eta, .link, earg = .earg)
  }, list( .link = link, .earg = earg ))),
  last = eval(substitute(expression({
    misc$link <-    c(df = .link)
    misc$earg <- list(df = .earg )
  }), list( .link = link, .earg = earg ))),
  link = eval(substitute(function(mu, extra = NULL) {
    theta2eta(mu, .link, earg = .earg)
  }, list( .link = link, .earg = earg ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
    mydf <- eta2theta(eta, .link, earg = .earg)
    if (residuals) stop("loglikelihood residuals not ",
                        "implemented yet") else
        sum(w * dchisq(x = y, df = mydf, ncp = 0, log = TRUE))
  }, list( .link = link, .earg = earg ))),
  vfamily = "chisq",
  deriv = eval(substitute(expression({
    mydf <- eta2theta(eta, .link, earg = .earg)
    dl.dv <- (log(y / 2) - digamma(mydf / 2)) / 2
    dv.deta <- dtheta.deta(mydf, .link, earg = .earg)
    c(w) * dl.dv * dv.deta
  }), list( .link = link, .earg = earg ))),
  weight = eval(substitute(expression({
    ed2l.dv2 <- -trigamma(mydf / 2) / 4
    wz <- -ed2l.dv2 * dv.deta^2
    c(w) * wz
  }), list( .link = link, .earg = earg ))))
}







dsimplex = function(x, mu = 0.5, dispersion = 1, log = FALSE) {
  if (!is.logical(log.arg <- log))
      stop("bad input for argument 'log'")
  rm(log)
  sigma = dispersion 

  deeFun = function(y, mu)
      (((y - mu) / (mu * (1 - mu)))^2) / (y * (1 - y))
  logpdf = (-0.5 * log(2 * pi) - log(sigma) - 1.5 * log(x) -
            1.5 * log1p(-x) - 0.5 * deeFun(x, mu) / sigma^2)
  logpdf[x     <= 0.0] = -Inf # log(0.0)
  logpdf[x     >= 1.0] = -Inf # log(0.0)
  logpdf[mu    <= 0.0] = NaN
  logpdf[mu    >= 1.0] = NaN
  logpdf[sigma <= 0.0] = NaN
  if (log.arg) logpdf else exp(logpdf)
}


rsimplex = function(n, mu = 0.5, dispersion = 1) {
  use.n = if ((length.n <- length(n)) > 1) length.n else
          if (!is.Numeric(n, integ = TRUE, allow = 1, posit = TRUE))
              stop("bad input for argument 'n'") else n

  oneval <- (length(mu) == 1 && length(dispersion) == 1)
  answer = rep(0.0, len = use.n)
  mu = rep(mu, len = use.n); dispersion = rep(dispersion, len = use.n)
  Kay1 = 3 * (dispersion * mu * (1-mu))^2

  if (oneval) {
    Kay1 = Kay1[1] # Since oneval means there is only one unique value
    mymu =   mu[1]
    myroots = polyroot(c(-mymu^2, Kay1+2*mymu^2, -3*Kay1+1-2*mymu, 2*Kay1))
    myroots = myroots[abs(Im(myroots)) < 0.00001]
    myroots = Re(myroots)
    myroots = myroots[myroots >= 0.0]
    myroots = myroots[myroots <= 1.0]
    pdfmax = dsimplex(myroots, mymu, dispersion[1])
    pdfmax = rep(max(pdfmax), len = use.n) # For multiple peaks
  } else {
    pdfmax = numeric(use.n)
    for (ii in 1:use.n) {
      myroots = polyroot(c(-mu[ii]^2, Kay1[ii]+2*mu[ii]^2,
                           -3*Kay1[ii]+1-2*mu[ii], 2*Kay1[ii]))
      myroots = myroots[abs(Im(myroots)) < 0.00001]
      myroots = Re(myroots)
      myroots = myroots[myroots >= 0.0]
      myroots = myroots[myroots <= 1.0]
      pdfmax[ii] = max(dsimplex(myroots, mu[ii], dispersion[ii]))
    }
  }

  index = 1:use.n
  nleft = length(index)
  while (nleft > 0) {
    xx = runif(nleft) # , 0, 1
    yy = runif(nleft, max = pdfmax[index])
    newindex = (1:nleft)[yy < dsimplex(xx, mu[index], dispersion[index])]
    if (length(newindex)) {
      answer[index[newindex]] = xx[newindex]
      index = setdiff(index, index[newindex])
      nleft = nleft - length(newindex)
    }
  }
  answer
}






 simplex = function(lmu = "logit", lsigma = "loge",
                    emu = list(), esigma = list(),
                    imu = NULL, isigma = NULL,
                    imethod = 1, shrinkage.init = 0.95,
                    zero = 2) {


  if (mode(lmu) != "character" && mode(lmu) != "name")
      lmu = as.character(substitute(lmu))
  if (mode(lsigma) != "character" && mode(lsigma) != "name")
      lsigma = as.character(substitute(lsigma))

  if (!is.list(emu)) emu = list()
  if (!is.list(esigma)) esigma = list()

  if (!is.Numeric(imethod, allow = 1, integ = TRUE, posit = TRUE) ||
       imethod > 3)
      stop("argument 'imethod' must be 1 or 2 or 3")
  if (!is.Numeric(shrinkage.init, allow = 1) || shrinkage.init < 0 ||
       shrinkage.init > 1) stop("bad input for argument 'shrinkage.init'")

    if (length(zero) && !is.Numeric(zero, integer = TRUE, posit = TRUE))
        stop("bad input for argument 'zero'")

  new("vglmff",
  blurb = c("Univariate Simplex distribution\n\n",
          "f(y) = [2*pi*sigma^2*(y*(1-y))^3]^(-0.5) * \n",
          "       exp[-0.5*(y-mu)^2 / (sigma^2 * y*(1-y)*mu^2*(1-mu)^2)],\n",
          "   0 < y < 1, 0 < mu < 1, sigma > 0\n\n",
          "Links:     ",
          namesof("mu", lmu, earg = emu), ", ",
          namesof("sigma", lsigma, earg = esigma), "\n\n",
          "Mean:              mu\n",
          "Variance function: V(mu) = mu^3 * (1 - mu)^3"),
    constraints = eval(substitute(expression({
        constraints = cm.zero.vgam(constraints, x, .zero, M)
    }), list( .zero = zero ))),
  initialize = eval(substitute(expression({
      y = as.numeric(y)
      if (ncol(y <- cbind(y)) != 1)
        stop("response must be a vector or a one-column matrix")
      if (any(y <= 0.0 | y >= 1.0))
        stop("all 'y' values must be in (0,1)")

      predictors.names = c(
          namesof("mu",    .lmu,    earg = .emu,    tag = FALSE),
          namesof("sigma", .lsigma, earg = .esigma, tag = FALSE))

      deeFun = function(y, mu)
          (((y - mu) / (mu * (1 - mu)))^2) / (y * (1 - y))

      if (!length(etastart)) {
          use.this = if ( .imethod == 3) weighted.mean(y, w) else
                     if ( .imethod == 1) median(y) else
                                             mean(y, trim = 0.1)
          init.mu = (1 - .sinit) * y + .sinit * use.this
          mu.init = rep(if (length( .imu )) .imu else init.mu, length = n)
          sigma.init = if (length( .isigma )) rep( .isigma, leng = n) else {
          use.this = deeFun(y, mu=init.mu)
          rep(sqrt( if ( .imethod == 3) weighted.mean(use.this, w) else
                    if ( .imethod == 1) median(use.this) else
                                            mean(use.this, trim = 0.1)),
              length = n)
          }
          etastart = cbind(theta2eta(mu.init,    .lmu,    earg = .emu),
                           theta2eta(sigma.init, .lsigma, earg = .esigma))
      }
  }), list( .lmu = lmu, .lsigma = lsigma,
            .emu = emu, .esigma = esigma,
            .imu = imu, .isigma = isigma,
            .sinit = shrinkage.init, .imethod = imethod ))),
  inverse = eval(substitute(function(eta, extra = NULL) {
      eta2theta(eta[,1], .lmu, earg = .emu)
  }, list( .lmu = lmu, .emu = emu ))),
  last = eval(substitute(expression({
      misc$link = c(mu = .lmu, sigma = .lsigma)
      misc$earg = list(mu = .emu, sigma = .esigma)
      misc$imu    = .imu
      misc$isigma = .isigma
      misc$imethod = .imethod
      misc$shrinkage.init = .sinit
  }), list( .lmu = lmu, .lsigma = lsigma,
            .imu = imu, .isigma = isigma,
            .emu = emu, .esigma = esigma,
            .sinit = shrinkage.init, .imethod = imethod ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
      sigma = eta2theta(eta[,2], .lsigma, earg = .esigma)
      if (residuals) stop("loglikelihood residuals not ",
                          "implemented yet") else {
        sum(w * dsimplex(x = y, mu = mu, dispersion = sigma, log = TRUE))
      }
  }, list( .lsigma = lsigma, .emu = emu,
           .esigma = esigma ))),
  vfamily = c("simplex"),
  deriv = eval(substitute(expression({
      deeFun = function(y, mu)
          (((y - mu) / (mu * (1 - mu)))^2) / (y * (1 - y))
      sigma       = eta2theta(eta[,2], .lsigma, earg = .esigma)
      dmu.deta    = dtheta.deta(mu,    .lmu,    earg = .emu)
      dsigma.deta = dtheta.deta(sigma, .lsigma, earg = .esigma)

      dl.dmu = (y - mu) * (deeFun(y, mu) +
               1 / (mu * (1 - mu))^2) / (mu * (1 - mu) * sigma^2)

      dl.dsigma = (deeFun(y, mu) / sigma^2 - 1) / sigma
      cbind(dl.dmu * dmu.deta, dl.dsigma * dsigma.deta)
  }), list( .lmu = lmu, .lsigma = lsigma,
            .emu = emu, .esigma = esigma ))),
  weight = eval(substitute(expression({
      wz = matrix(0.0, n, M)  # Diagonal!!
      eim11 = 3 / (mu * (1 - mu)) + 1 / (sigma^2 * (mu * (1 - mu))^3)
      wz[, iam(1, 1, M)] = eim11 * dmu.deta^2
      wz[, iam(2, 2, M)] = (2 / sigma^2) * dsigma.deta^2
      c(w) * wz
  }), list( .lmu = lmu, .lsigma = lsigma,
            .emu = emu, .esigma = esigma ))))
}








 rig = function(lmu = "identity", llambda = "loge",
               emu = list(), elambda = list(), imu = NULL, ilambda=1)
{

    if (mode(lmu) != "character" && mode(lmu) != "name")
        lmu = as.character(substitute(lmu))
    if (mode(llambda) != "character" && mode(llambda) != "name")
        llambda = as.character(substitute(llambda))
    if (!is.Numeric(ilambda, posit = TRUE))
        stop("bad input for 'ilambda'")
    if (!is.list(emu)) emu = list()
    if (!is.list(elambda)) elambda = list()

    new("vglmff",
    blurb = c("Reciprocal inverse Gaussian distribution \n",
            "f(y) = [lambda/(2*pi*y)]^(0.5) * \n",
            "       exp[-0.5*(lambda/y) * (y-mu)^2], ",
            "  0 < y,\n",
            "Links:     ",
            namesof("mu", lmu, earg = emu), ", ",
            namesof("lambda", llambda, earg = elambda), "\n\n",
            "Mean:     mu"),
    initialize = eval(substitute(expression({
        if (ncol(cbind(y)) != 1)
            stop("response must be a vector or a one-column matrix")
        y = as.numeric(y)
        if (any(y <= 0))
            stop("all y values must be > 0")
        predictors.names = 
        c(namesof("mu", .lmu, earg = .emu, tag = FALSE),
          namesof("lambda", .llambda, earg = .elambda, tag = FALSE))
        if (!length(etastart)) {
            mu.init = rep(if (length( .imu)) .imu else
                           median(y), length = n)
            lambda.init = rep(if (length( .ilambda )) .ilambda else
                           sqrt(var(y)), length = n)
            etastart = cbind(theta2eta(mu.init, .lmu, earg = .emu),
                             theta2eta(lambda.init, .llambda, earg = .elambda))
        }
    }), list( .lmu = lmu, .llambda = llambda,
              .emu = emu, .elambda = elambda,
              .imu=imu, .ilambda=ilambda ))),
    inverse = eval(substitute(function(eta, extra = NULL) {
        eta2theta(eta[,1], .lmu, earg = .emu)
    }, list( .lmu = lmu,
             .emu = emu, .elambda = elambda ))),
    last = eval(substitute(expression({
        misc$d3 = d3    # because save.weights = FALSE
        misc$link = c(mu= .lmu, lambda= .llambda)
        misc$earg = list(mu= .emu, lambda= .elambda)
        misc$pooled.weight = pooled.weight
    }), list( .lmu = lmu, .llambda = llambda,
              .emu = emu, .elambda = elambda ))),
    loglikelihood = eval(substitute(
                  function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
        lambda = eta2theta(eta[,2], .llambda, earg = .elambda)
        if (residuals) stop("loglikelihood residuals not ",
                            "implemented yet") else
        sum(w * (-0.5*log(y) + 0.5*log(lambda) - (0.5*lambda/y) * (y-mu)^2))
    }, list( .llambda = llambda,
             .emu = emu, .elambda = elambda ))),
    vfamily = c("rig"),
    deriv = eval(substitute(expression({
        if (iter == 1) {
            d3 = deriv3(~ w * 
                 (-0.5*log(y) + 0.5*log(lambda) - (0.5*lambda/y) * (y-mu)^2),
                        c("mu", "lambda"), hessian= TRUE)
        }

        lambda = eta2theta(eta[,2], .llambda, earg = .elambda)

        eval.d3 = eval(d3)
        dl.dthetas =  attr(eval.d3, "gradient")

        dmu.deta = dtheta.deta(mu, .lmu, earg = .emu)
        dlambda.deta = dtheta.deta(lambda, .llambda, earg = .elambda)
        dtheta.detas = cbind(dmu.deta, dlambda.deta)

        dl.dthetas * dtheta.detas
    }), list( .lmu = lmu, .llambda = llambda,
              .emu = emu, .elambda = elambda ))),
    weight = eval(substitute(expression({
        d2l.dthetas2 =  attr(eval.d3, "hessian")

        wz = matrix(as.numeric(NA), n, dimm(M))  #3=dimm(M)
        wz[,iam(1,1,M)] = -d2l.dthetas2[,1,1] * dtheta.detas[,1]^2
        wz[,iam(2,2,M)] = -d2l.dthetas2[,2,2] * dtheta.detas[,2]^2
        wz[,iam(1,2,M)] = -d2l.dthetas2[,1,2] * dtheta.detas[,1] *
                                                 dtheta.detas[,2]
        if (!.expected) {
            d2mudeta2 = d2theta.deta2(mu, .lmu, earg = .emu)
            d2lambda = d2theta.deta2(lambda, .llambda, earg = .elambda)
            wz[,iam(1,1,M)] = wz[,iam(1,1,M)] - dl.dthetas[,1] * d2mudeta2
            wz[,iam(2,2,M)] = wz[,iam(2,2,M)] - dl.dthetas[,2] * d2lambda
        }

        if (intercept.only) {
            sumw = sum(w)
            for(ii in 1:ncol(wz))
                wz[,ii] = sum(wz[,ii]) / sumw
            pooled.weight = TRUE
            wz = c(w) * wz   # Put back the weights
        } else
            pooled.weight = FALSE

        wz
    }), list( .lmu = lmu, .llambda = llambda, .expected = FALSE,
              .emu = emu, .elambda = elambda ))))
}



 hypersecant = function(link.theta = "elogit",
    earg = if (link.theta == "elogit") list(min=-pi/2, max=pi/2) else list(),
    init.theta = NULL)
{

    if (mode(link.theta) != "character" && mode(link.theta) != "name")
        link.theta = as.character(substitute(link.theta))
    if (!is.list(earg)) earg = list()

    new("vglmff",
    blurb = c("Hyperbolic Secant distribution \n",
            "f(y) = exp(theta*y + log(cos(theta ))) / (2*cosh(pi*y/2))\n",
            "  for all y,\n",
            "Link:     ",
            namesof("theta", link.theta, earg = earg), "\n\n",
            "Mean:     tan(theta)"),
    initialize = eval(substitute(expression({
        if (ncol(cbind(y)) != 1)
            stop("response must be a vector or a one-column matrix")
        predictors.names = namesof("theta", .link.theta, earg = .earg, tag = FALSE)
        if (!length(etastart)) {
            theta.init = rep(if (length( .init.theta)) .init.theta else
                             median(y), length = n)
            etastart = theta2eta(theta.init, .link.theta, earg = .earg)
        }
    }), list( .link.theta = link.theta, .earg = earg,
              .init.theta=init.theta ))),
    inverse = eval(substitute(function(eta, extra = NULL) {
        theta = eta2theta(eta, .link.theta, earg = .earg)
        tan(theta)
    }, list( .link.theta = link.theta, .earg = earg ))),
    last = eval(substitute(expression({
        misc$link = c(theta= .link.theta )
        misc$earg = list(theta= .earg )
        misc$expected = TRUE
    }), list( .link.theta = link.theta, .earg = earg ))),
    loglikelihood = eval(substitute(function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
        theta = eta2theta(eta, .link.theta, earg = .earg)
        if (residuals) stop("loglikelihood residuals not ",
                            "implemented yet") else
        sum(w * (theta*y + log(cos(theta)) - log(cosh(pi*y/2 ))))
    }, list( .link.theta = link.theta, .earg = earg ))),
    vfamily = c("hypersecant"),
    deriv = eval(substitute(expression({
        theta = eta2theta(eta, .link.theta, earg = .earg)
        dl.dthetas =  y - tan(theta)
        dparam.deta = dtheta.deta(theta, .link.theta, earg = .earg)
        c(w) * dl.dthetas * dparam.deta
    }), list( .link.theta = link.theta, .earg = earg ))),
    weight = expression({
        d2l.dthetas2 =  1 / cos(theta)^2
        wz = c(w) * d2l.dthetas2 * dparam.deta^2
        wz
    }))
}



 hypersecant.1 = function(link.theta = "elogit",
    earg=if (link.theta == "elogit") list(min=-pi/2, max=pi/2) else list(),
    init.theta = NULL)
{

    if (mode(link.theta) != "character" && mode(link.theta) != "name")
        link.theta = as.character(substitute(link.theta))
    if (!is.list(earg)) earg = list()

    new("vglmff",
    blurb = c("Hyperbolic Secant distribution \n",
            "f(y) = (cos(theta)/pi) * y^(-0.5+theta/pi) * \n",
            "       (1-y)^(-0.5-theta/pi), ",
            "  0 < y < 1,\n",
            "Link:     ",
            namesof("theta", link.theta, earg = earg), "\n\n",
            "Mean:     0.5 + theta/pi", "\n",
            "Variance: (pi^2 - 4*theta^2) / (8*pi^2)"),
    initialize = eval(substitute(expression({
        if (ncol(cbind(y)) != 1)
            stop("response must be a vector or a one-column matrix")
        y = as.numeric(y)
        if (any(y <= 0 | y >= 1))
            stop("all y values must be in (0,1)")
        predictors.names = namesof("theta", .link.theta, earg = .earg, tag = FALSE)
        if (!length(etastart)) {
            theta.init = rep(if (length( .init.theta)) .init.theta else
                           median(y), length = n)

            etastart = theta2eta(theta.init, .link.theta, earg = .earg)
        }
    }), list( .link.theta = link.theta, .earg = earg,
              .init.theta=init.theta ))),
    inverse = eval(substitute(function(eta, extra = NULL) {
        theta = eta2theta(eta, .link.theta, earg = .earg)
        0.5 + theta/pi
    }, list( .link.theta = link.theta, .earg = earg ))),
    last = eval(substitute(expression({
        misc$link = c(theta= .link.theta)
        misc$earg = list(theta= .earg )
        misc$expected = TRUE
    }), list( .link.theta = link.theta, .earg = earg ))),
    loglikelihood = eval(substitute(function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
        theta = eta2theta(eta, .link.theta, earg = .earg)
        if (residuals) stop("loglikelihood residuals not ",
                            "implemented yet") else
        sum(w * (log(cos(theta)) + (-0.5+theta/pi)*log(y) +
                (-0.5-theta/pi)*log1p(-y )))
    }, list( .link.theta = link.theta, .earg = earg ))),
    vfamily = c("hypersecant.1"),
    deriv = eval(substitute(expression({
        theta = eta2theta(eta, .link.theta, earg = .earg)
        dl.dthetas =  -tan(theta) + log(y/(1-y)) / pi 
        dparam.deta = dtheta.deta(theta, .link.theta, earg = .earg)
        c(w) * dl.dthetas * dparam.deta
    }), list( .link.theta = link.theta, .earg = earg ))),
    weight = expression({
        d2l.dthetas2 =  1 / cos(theta)^2
        wz = c(w) * d2l.dthetas2 * dparam.deta^2
        wz
    }))
}



 leipnik = function(lmu = "logit", llambda = "loge",
                    emu = list(), elambda = list(), imu = NULL, ilambda = NULL)
{


    if (mode(lmu) != "character" && mode(lmu) != "name")
        lmu = as.character(substitute(lmu))
    if (mode(llambda) != "character" && mode(llambda) != "name")
        llambda = as.character(substitute(llambda))
    if (is.Numeric(ilambda) && any(ilambda <= -1))
        stop("ilambda must be > -1")
    if (!is.list(emu)) emu = list()
    if (!is.list(elambda)) elambda = list()

    new("vglmff",
    blurb = c("Leipnik's distribution \n",
    "f(y) = (y(1-y))^(-1/2) * [1 + (y-mu)^2 / (y*(1-y))]^(-lambda/2) /\n",
            "       Beta[(lambda+1)/2, 1/2], ",
            "  0 < y < 1,  lambda > -1\n",
            "Links:     ",
            namesof("mu", lmu, earg = emu), ", ",
            namesof("lambda", llambda, earg = elambda), "\n\n",
            "Mean:     mu\n",
            "Variance: mu*(1-mu)"),
    initialize = eval(substitute(expression({
        if (ncol(cbind(y)) != 1)
            stop("response must be a vector or a one-column matrix")
        y = as.numeric(y)
        if (any(y <= 0 | y >= 1))
            stop("all y values must be in (0,1)")
        predictors.names =
        c(namesof("mu", .lmu, earg = .emu, tag = FALSE),
          namesof("lambda", .llambda, earg = .elambda, tag = FALSE))
        if (!length(etastart)) {
            mu.init = rep(if (length( .imu)) .imu else
                          (y), length = n)
            lambda.init = rep(if (length( .ilambda)) .ilambda else
                           1/var(y), length = n)
            etastart = cbind(theta2eta(mu.init, .lmu, earg = .emu),
                             theta2eta(lambda.init, .llambda, earg = .elambda))
        }
    }), list( .lmu = lmu, .llambda = llambda, .imu=imu, .ilambda=ilambda,
              .emu = emu, .elambda = elambda ))),
    inverse = eval(substitute(function(eta, extra = NULL) {
        eta2theta(eta[,1], .lmu, earg = .emu)
    }, list( .lmu = lmu,
             .emu = emu, .elambda = elambda ))),
    last = eval(substitute(expression({
        misc$link = c(mu= .lmu, lambda= .llambda)
        misc$earg = list(mu= .emu, lambda= .elambda)
        misc$pooled.weight = pooled.weight
        misc$expected = FALSE
    }), list( .lmu = lmu, .llambda = llambda,
              .emu = emu, .elambda = elambda ))),
    loglikelihood = eval(substitute(
        function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
        lambda = eta2theta(eta[,2], .llambda, earg = .elambda)
        if (residuals) stop("loglikelihood residuals not ",
                            "implemented yet") else
        sum(w * (-0.5*log(y*(1-y)) - 0.5 * lambda *
                log1p((y-mu)^2 / (y*(1-y ))) - lgamma((lambda+1)/2) +
                lgamma(1+ lambda/2 )))
    }, list( .llambda = llambda,
             .emu = emu, .elambda = elambda ))),
    vfamily = c("leipnik"),
    deriv = eval(substitute(expression({
        lambda = eta2theta(eta[,2], .llambda, earg = .elambda)
        dl.dthetas =
          c(w) * cbind(dl.dmu = lambda*(y-mu) / (y*(1-y)+(y-mu)^2),
                       dl.dlambda= -0.5 * log1p((y-mu)^2 / (y*(1-y))) -
                         0.5*digamma((lambda+1)/2) +
                         0.5*digamma(1+lambda/2))
        dmu.deta = dtheta.deta(mu, .lmu, earg = .emu)
        dlambda.deta = dtheta.deta(lambda, .llambda, earg = .elambda)
        dtheta.detas = cbind(dmu.deta, dlambda.deta)
        dl.dthetas * dtheta.detas
    }), list( .lmu = lmu, .llambda = llambda,
              .emu = emu, .elambda = elambda ))),
    weight = eval(substitute(expression({
        if (is.R()) {
            denominator = y*(1-y) + (y-mu)^2
            d2l.dthetas2 =  array(NA, c(n,2,2))
            d2l.dthetas2[,1,1] = c(w) * lambda*(-y*(1-y)+(y-mu)^2)/denominator^2
            d2l.dthetas2[,1,2] = 
            d2l.dthetas2[,2,1] = c(w) * (y-mu) / denominator
            d2l.dthetas2[,2,2] = c(w) * (-0.25*trigamma((lambda+1)/2) +
                                       0.25*trigamma(1+lambda/2))
        } else {
            d2l.dthetas2 =  attr(eval.d3, "hessian")
        }

        wz = matrix(as.numeric(NA), n, dimm(M))  #3=dimm(M)
        wz[,iam(1,1,M)] = -d2l.dthetas2[,1,1] * dtheta.detas[,1]^2
        wz[,iam(2,2,M)] = -d2l.dthetas2[,2,2] * dtheta.detas[,2]^2
        wz[,iam(1,2,M)] = -d2l.dthetas2[,1,2] * dtheta.detas[,1] *
                                                dtheta.detas[,2]
        if (!.expected) {
            d2mudeta2 = d2theta.deta2(mu, .lmu, earg = .emu)
            d2lambda = d2theta.deta2(lambda, .llambda, earg = .elambda)
            wz[,iam(1,1,M)] = wz[,iam(1,1,M)] - dl.dthetas[,1] * d2mudeta2
            wz[,iam(2,2,M)] = wz[,iam(2,2,M)] - dl.dthetas[,2] * d2lambda
        }

        if (intercept.only) {
            sumw = sum(w)
            for(ii in 1:ncol(wz))
                wz[,ii] = sum(wz[,ii]) / sumw
            pooled.weight = TRUE
            wz = c(w) * wz   # Put back the weights
        } else
            pooled.weight = FALSE

        wz
    }), list( .lmu = lmu, .llambda = llambda, .expected = FALSE,
              .emu = emu, .elambda = elambda ))))
}





 invbinomial = function(lrho = "elogit", llambda = "loge",
          erho=if (lrho == "elogit") list(min = 0.5, max = 1) else list(),
          elambda = list(),
          irho = NULL,
          ilambda = NULL,
          zero = NULL)
{

    if (mode(lrho) != "character" && mode(lrho) != "name")
        lrho = as.character(substitute(lrho))
    if (mode(llambda) != "character" && mode(llambda) != "name")
        llambda = as.character(substitute(llambda))
    if (length(zero) && !is.Numeric(zero, integer = TRUE, posit = TRUE))
        stop("bad input for argument 'zero'")
    if (!is.list(erho)) erho = list()
    if (!is.list(elambda)) elambda = list()

    new("vglmff",
    blurb = c("Inverse binomial distribution\n\n",
            "Links:    ",
            namesof("rho", lrho, earg = erho), ", ", 
            namesof("lambda", llambda, earg = elambda), "\n", 
            "Mean:     lambda*(1-rho)/(2*rho-1)\n",
            "Variance: lambda*rho*(1-rho)/(2*rho-1)^3\n"),
    constraints = eval(substitute(expression({
        constraints = cm.zero.vgam(constraints, x, .zero, M)
    }), list( .zero = zero ))),
    initialize = eval(substitute(expression({
        if (ncol(cbind(y)) != 1)
            stop("response must be a vector or a one-column matrix")
        predictors.names =
        c(namesof("rho", .lrho, earg = .erho, tag = FALSE),
          namesof("lambda", .llambda, earg = .elambda, tag = FALSE))
        if (!length(etastart)) {
            covarn = sd(y)^2 / weighted.mean(y, w)
            temp1 = 0.5 + (1 + sqrt(1+8*covarn)) / (8*covarn)
            temp2 = 0.5 + (1 - sqrt(1+8*covarn)) / (8*covarn)
            init.rho = rep(if (length( .irho)) .irho else {
                ifelse(temp1 > 0.5 && temp1 < 1, temp1, temp2)
            }, length = n)
            init.lambda = rep(if (length( .ilambda)) .ilambda else {
                (2*init.rho-1) * weighted.mean(y, w) / (1-init.rho)
            }, length = n)
            etastart = cbind(theta2eta(init.rho, .lrho, earg = .erho),
                             theta2eta(init.lambda, .llambda, earg = .elambda))
        }
    }), list( .llambda = llambda, .lrho=lrho,
              .elambda = elambda, .erho=erho,
              .ilambda=ilambda, .irho=irho ))),
    inverse = eval(substitute(function(eta, extra = NULL) {
        rho = eta2theta(eta[,1], .lrho, earg = .erho)
        lambda = eta2theta(eta[,2], .llambda, earg = .elambda)
        ifelse(rho > 0.5, lambda*(1-rho)/(2*rho-1), NA)
    }, list( .llambda = llambda, .lrho=lrho,
             .elambda = elambda, .erho=erho ))),
    last = eval(substitute(expression({
        misc$link = c(rho= .lrho, lambda= .llambda)
        misc$earg = list(rho= .erho, lambda= .elambda)
        misc$pooled.weight = pooled.weight
    }), list( .llambda = llambda, .lrho=lrho,
              .elambda = elambda, .erho=erho ))),
    loglikelihood = eval(substitute(
        function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
        rho = eta2theta(eta[,1], .lrho, earg = .erho)
        lambda = eta2theta(eta[,2], .llambda, earg = .elambda)
        if (residuals) stop("loglikelihood residuals not ",
                            "implemented yet") else
        sum(w*(log(lambda) - lgamma(2*y+lambda) - lgamma(y+1) -
               lgamma(y+lambda+1) + y*log(rho) + y*log1p(-rho) +
               lambda*log(rho)))
    }, list( .llambda = llambda, .lrho=lrho,
             .elambda = elambda, .erho=erho ))),
    vfamily = c("invbinomial"),
    deriv = eval(substitute(expression({
        rho = eta2theta(eta[,1], .lrho, earg = .erho)
        lambda = eta2theta(eta[,2], .llambda, earg = .elambda)
        dl.drho = (y + lambda)/rho - y/(1-rho)
        dl.dlambda = 1/lambda - digamma(2*y+lambda) - digamma(y+lambda+1) +
                     log(rho)
        drho.deta = dtheta.deta(rho, .lrho, earg = .erho)
        dlambda.deta = dtheta.deta(lambda, .llambda, earg = .elambda)
        c(w) * cbind(dl.drho * drho.deta,
                     dl.dlambda * dlambda.deta )
    }), list( .llambda = llambda, .lrho=lrho,
              .elambda = elambda, .erho=erho ))),
    weight = eval(substitute(expression({
        ed2l.drho2 = (mu+lambda) / rho^2 + mu / (1-rho)^2
        d2l.dlambda2 = 1/(lambda^2) + trigamma(2*y+lambda)+trigamma(y+lambda+1)
        ed2l.dlambdarho = -1/rho
        wz = matrix(as.numeric(NA), n, dimm(M))  #3=dimm(M)
        wz[,iam(1,1,M)] = ed2l.drho2 * drho.deta^2
        wz[,iam(1,2,M)] = ed2l.dlambdarho * dlambda.deta * drho.deta
        wz[,iam(2,2,M)] =  d2l.dlambda2 * dlambda.deta^2

        d2rhodeta2 = d2theta.deta2(rho, .lrho, earg = .erho)
        d2lambda.deta2 = d2theta.deta2(lambda, .llambda, earg = .elambda)
        wz = c(w) * wz

        if (intercept.only) {
            pooled.weight = TRUE

            wz[,iam(2,2,M)] =  sum(wz[,iam(2,2,M)]) / sum(w)

        } else
            pooled.weight = FALSE

        wz
    }), list( .llambda = llambda, .lrho=lrho,
              .elambda = elambda, .erho=erho ))))
}



 genpoisson = function(llambda = "elogit", ltheta = "loge",
                  elambda=if (llambda == "elogit") list(min=-1,max=1) else list(),
                      etheta = list(),
                      ilambda = NULL, itheta = NULL,
                      use.approx = TRUE,
                      imethod = 1, zero=1)
{


    if (mode(llambda) != "character" && mode(llambda) != "name")
        llambda = as.character(substitute(llambda))
    if (mode(ltheta) != "character" && mode(ltheta) != "name")
        ltheta = as.character(substitute(ltheta))
    if (length(zero) && !is.Numeric(zero, integer = TRUE, posit = TRUE))
        stop("bad input for argument 'zero'")
    if (!is.list(elambda)) elambda = list()
    if (!is.list(etheta)) etheta = list()
    if (!is.Numeric(imethod, allow = 1, integ = TRUE, posit = TRUE) ||
       imethod > 2)
        stop("argument 'imethod' must be 1 or 2")
    if (!is.logical(use.approx) || length(use.approx) != 1)
        stop("'use.approx' must be logical value")

    new("vglmff",
    blurb = c("Generalized Poisson distribution\n\n",
            "Links:    ",
            namesof("lambda", llambda, earg = elambda), ", ", 
            namesof("theta", ltheta, earg = etheta), "\n", 
            "Mean:     theta / (1-lambda)\n",
            "Variance: theta / (1-lambda)^3"),
    constraints = eval(substitute(expression({
        constraints = cm.zero.vgam(constraints, x, .zero, M)
    }), list( .zero = zero ))),
    initialize = eval(substitute(expression({
        if (ncol(cbind(y)) != 1)
            stop("response must be a vector or a one-column matrix")
        predictors.names =
           c(namesof("lambda", .llambda, earg = .elambda, tag = FALSE),
             namesof("theta",  .ltheta,  earg = .etheta,  tag = FALSE))
        init.lambda = if ( .imethod == 1)
            1 - sqrt(weighted.mean(y,w) / var(y)) else 0.5
        init.theta  = if ( .imethod == 1)
            sqrt((0.01+weighted.mean(y,w)^3)/var(y)) else
            median(y)*(1-init.lambda)
        if (init.theta <= 0)
            init.theta = 0.1
        cutpt = if (init.lambda < 0) {
            mmm = max(trunc(-init.theta / init.lambda), 4)
            max(-1, -init.theta /mmm)
        } else -1
        if (init.lambda <= cutpt)
            init.lambda = cutpt + 0.1
        if (init.lambda >= 1)
            init.lambda = 0.9
        if (!length(etastart)) {
            lambda = rep(if (length( .ilambda)) .ilambda else
                       init.lambda, length = n)
            theta = rep(if (length( .itheta)) .itheta else init.theta, length = n)
            etastart = cbind(theta2eta(lambda, .llambda, earg = .elambda),
                             theta2eta(theta,  .ltheta,  earg = .etheta))
        }
    }), list( .ltheta=ltheta, .llambda = llambda,
              .etheta=etheta, .elambda = elambda,
              .imethod = imethod,
              .itheta=itheta, .ilambda=ilambda )) ),
    inverse = eval(substitute(function(eta, extra = NULL) {
        lambda = eta2theta(eta[,1], .llambda, earg = .elambda)
        theta = eta2theta(eta[,2], .ltheta, earg = .etheta)
        theta/(1-lambda)
    }, list( .ltheta=ltheta, .llambda = llambda,
             .etheta=etheta, .elambda = elambda ))),
    last = eval(substitute(expression({
        misc$link = c(lambda=.llambda, theta=.ltheta)
        misc$earg = list(lambda=.elambda, theta=.etheta)
        if (! .use.approx )
            misc$pooled.weight = pooled.weight
    }), list( .ltheta=ltheta, .llambda = llambda,
              .use.approx = use.approx,
              .etheta=etheta, .elambda = elambda ))),
    loglikelihood = eval(substitute(
        function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
        lambda = eta2theta(eta[,1], .llambda, earg = .elambda)
        theta = eta2theta(eta[,2], .ltheta, earg = .etheta)
        index = (y == 0)
        if (residuals) stop("loglikelihood residuals not ",
                            "implemented yet") else
        sum(w[index]*(-theta[index])) + 
        sum(w[!index] * (-y[!index]*lambda[!index]-theta[!index] +
            (y[!index]-1)*log(theta[!index]+y[!index]*lambda[!index]) +
            log(theta[!index]) - lgamma(y[!index]+1)) )
    }, list( .ltheta=ltheta, .llambda = llambda,
             .etheta=etheta, .elambda = elambda ))),
    vfamily = c("genpoisson"),
    deriv = eval(substitute(expression({
        lambda = eta2theta(eta[,1], .llambda, earg = .elambda)
        theta = eta2theta(eta[,2], .ltheta, earg = .etheta)
        dl.dlambda = -y + y*(y-1)/(theta+y*lambda)
        dl.dtheta = -1 + (y-1)/(theta+y*lambda) + 1/theta
        dTHETA.deta = dtheta.deta(theta, .ltheta, earg = .etheta)
        dlambda.deta = dtheta.deta(lambda, .llambda, earg = .elambda)
        c(w) * cbind(dl.dlambda * dlambda.deta,
                     dl.dtheta * dTHETA.deta )
    }), list( .ltheta=ltheta, .llambda = llambda,
              .etheta=etheta, .elambda = elambda ))),
    weight = eval(substitute(expression({
        wz = matrix(as.numeric(NA), n, dimm(M))  #3=dimm(M)
        if ( .use.approx ) {
            BBB = (theta+2)*(theta+2*lambda-theta*lambda)-(theta^2)*(1-lambda)
            d2l.dlambda2 = 2 * theta * (theta+2) / ((1-lambda) * BBB)
            d2l.dtheta2 = 2 * (1 + lambda * (2/theta - 1)) / BBB
            d2l.dthetalambda =  2 * theta / BBB
            wz[,iam(1,1,M)] = d2l.dlambda2 * dlambda.deta^2
            wz[,iam(2,2,M)] = d2l.dtheta2 * dTHETA.deta^2
            wz[,iam(1,2,M)] = d2l.dthetalambda * dTHETA.deta * dlambda.deta
            wz = c(w) * wz
        } else {
            d2l.dlambda2 = -y^2 * (y-1) / (theta+y*lambda)^2
            d2l.dtheta2 = -(y-1)/(theta+y*lambda)^2 - 1 / theta^2
            d2l.dthetalambda =  -y * (y-1) / (theta+y*lambda)^2 
            wz[,iam(1,1,M)] = -d2l.dlambda2 * dlambda.deta^2
            wz[,iam(2,2,M)] = -d2l.dtheta2 * dTHETA.deta^2
            wz[,iam(1,2,M)] = -d2l.dthetalambda * dTHETA.deta * dlambda.deta

            d2THETA.deta2 = d2theta.deta2(theta, .ltheta, earg = .etheta)
            d2lambdadeta2 = d2theta.deta2(lambda, .llambda, earg = .elambda)
            wz[,iam(1,1,M)] = wz[,iam(1,1,M)] - dl.dlambda * d2lambdadeta2
            wz[,iam(2,2,M)] = wz[,iam(2,2,M)] - dl.dtheta * d2THETA.deta2
            wz = c(w) * wz

            if (intercept.only) {
                sumw = sum(w)
                for(ii in 1:ncol(wz))
                    wz[,ii] = sum(wz[,ii]) / sumw
                pooled.weight = TRUE
                wz = c(w) * wz   # Put back the weights
            } else
                pooled.weight = FALSE
            }
        wz
    }), list( .ltheta=ltheta, .llambda = llambda,
              .use.approx = use.approx,
              .etheta=etheta, .elambda = elambda ))))
}






dlgamma = function(x, location = 0, scale = 1, k = 1, log = FALSE) {
  if (!is.logical(log.arg <- log))
    stop("bad input for argument 'log'")
  rm(log)

  if (!is.Numeric(scale, posit = TRUE))
    stop("bad input for argument 'scale'")
  if (!is.Numeric(k, posit = TRUE))
    stop("bad input for argument 'k'")
  z = (x-location) / scale
  if (log.arg) {
    k * z - exp(z) - log(scale) - lgamma(k)
  } else {
    exp(k * z - exp(z)) / (scale * gamma(k))
  }
}
plgamma = function(q, location = 0, scale = 1, k=1) {
  if (!is.Numeric(scale, posit = TRUE))
  stop("bad input for argument 'scale'")
  if (!is.Numeric(k, posit = TRUE))
  stop("bad input for argument 'k'")
  z = (q-location)/scale
  pgamma(exp(z), k)
}
qlgamma = function(p, location = 0, scale = 1, k=1) {
  if (!is.Numeric(scale, posit = TRUE))
    stop("bad input for argument 'scale'")
  if (!is.Numeric(k, posit = TRUE))
    stop("bad input for argument 'k'")
  q = qgamma(p, k)
  location + scale * log(q)
}
rlgamma = function(n, location = 0, scale = 1, k=1) {
  if (!is.Numeric(n, posit = TRUE, integ = TRUE, allow = 1)) 
    stop("bad input for argument 'n'")
  if (!is.Numeric(scale, posit = TRUE))
    stop("bad input for argument 'scale'")
  if (!is.Numeric(k, posit = TRUE))
    stop("bad input for argument 'k'")
  y = rgamma(n, k)
  location + scale * log(y)
}



 lgammaff = function(link = "loge", earg = list(), init.k = NULL)
{
    if (mode(link) != "character" && mode(link) != "name")
        link = as.character(substitute(link))
    if (!is.list(earg)) earg = list()

    new("vglmff",
    blurb = c("Log-gamma distribution f(y) = exp(ky - e^y)/gamma(k)), k>0\n\n",
            "Link:    ",
            namesof("k", link, earg = earg), "\n", "\n",
            "Mean:    digamma(k)", "\n"),
    initialize = eval(substitute(expression({
        if (ncol(cbind(y)) != 1)
            stop("response must be a vector or a one-column matrix")
        predictors.names = namesof("k", .link, earg = .earg, tag = FALSE) 
        if (!length(etastart)) {
            k.init = if (length( .init.k)) rep( .init.k, len = length(y)) else {
                medy = median(y)
                if (medy < 2) 5 else if (medy < 4) 20 else exp(0.7 * medy)
            }
            etastart = theta2eta(k.init, .link, earg = .earg)
        }
    }), list( .link = link, .earg = earg, .init.k=init.k ))),
    inverse = eval(substitute(function(eta, extra = NULL) {
        k = eta2theta(eta, .link, earg = .earg)
        digamma(k)
    }, list( .link = link, .earg = earg ))),
    last = eval(substitute(expression({
        misc$link = c(k= .link )
        misc$earg = list(k= .earg )
        misc$expected = TRUE
    }), list( .link = link, .earg = earg ))),
    loglikelihood = eval(substitute(
        function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
        kk = eta2theta(eta, .link, earg = .earg)
        if (residuals) stop("loglikelihood residuals not ",
                            "implemented yet") else {
            sum(w * dlgamma(x = y, location = 0, scale = 1, k=kk, log = TRUE))
        }
    }, list( .link = link, .earg = earg ))),
    vfamily = c("lgammaff"),
    deriv = eval(substitute(expression({
        k = eta2theta(eta, .link, earg = .earg) 
        dl.dk = y - digamma(k)
        dk.deta = dtheta.deta(k, .link, earg = .earg)
        c(w) * dl.dk * dk.deta
    }), list( .link = link, .earg = earg ))),
    weight = eval(substitute(expression({
        ed2l.dk2 = trigamma(k)
        wz = c(w) * dk.deta^2 * ed2l.dk2
        wz
    }), list( .link = link, .earg = earg ))))
}







 lgamma3ff = function(llocation = "identity", lscale = "loge", lshape = "loge",
                     elocation = list(), escale = list(), eshape = list(),
                     ilocation = NULL, iscale = NULL, ishape = 1, zero = NULL)
{
    if (mode(llocation) != "character" && mode(llocation) != "name")
        llocation = as.character(substitute(llocation))
    if (mode(lscale) != "character" && mode(lscale) != "name")
        lscale = as.character(substitute(lscale))
    if (mode(lshape) != "character" && mode(lshape) != "name")
        lshape = as.character(substitute(lshape))
    if (length(zero) && !is.Numeric(zero, integer = TRUE, posit = TRUE))
        stop("bad input for argument 'zero'")
    if (length(iscale) && !is.Numeric(iscale, posit = TRUE))
        stop("bad input for argument 'iscale'")
    if (!is.list(elocation)) elocation = list()
    if (!is.list(escale)) escale = list()
    if (!is.list(eshape)) eshape = list()

    new("vglmff",
    blurb = c("Log-gamma distribution",
            " f(y) = exp(k(y-a)/b - e^((y-a)/b))/(b*gamma(k)), ",
            "location=a, scale=b>0, shape=k>0\n\n",
            "Links:    ",
            namesof("location", llocation, earg = elocation), ", ",
            namesof("scale", lscale, earg = escale), ", ",
            namesof("shape", lshape, earg = eshape), "\n\n",
            "Mean:     a + b*digamma(k)", "\n"),
    constraints = eval(substitute(expression({
        constraints = cm.zero.vgam(constraints, x, .zero, M)
    }), list( .zero = zero ))),
    initialize = eval(substitute(expression({
        if (ncol(cbind(y)) != 1)
            stop("response must be a vector or a one-column matrix")
        predictors.names =
        c(namesof("location", .llocation, earg = .elocation, tag = FALSE),
          namesof("scale", .lscale, earg = .escale, tag = FALSE),
          namesof("shape", .lshape, earg = .eshape, tag = FALSE))
        if (!length(etastart)) {
            k.init = if (length( .ishape)) rep( .ishape, len = length(y)) else {
                rep(exp(median(y)), len = length(y))
            }
            scale.init = if (length( .iscale)) rep( .iscale, len = length(y)) else {
                rep(sqrt(var(y) / trigamma(k.init)), len = length(y))
            }
            loc.init = if (length( .iloc)) rep( .iloc, len = length(y)) else {
                rep(median(y) - scale.init * digamma(k.init), len = length(y))
            }
            etastart = cbind(theta2eta(loc.init, .llocation, earg = .elocation),
                             theta2eta(scale.init, .lscale, earg = .escale),
                             theta2eta(k.init, .lshape, earg = .eshape))
        }
    }), list( .llocation = llocation, .lscale = lscale, .lshape = lshape,
              .elocation = elocation, .escale = escale, .eshape = eshape,
              .iloc=ilocation, .iscale = iscale, .ishape = ishape ))),
    inverse = eval(substitute(function(eta, extra = NULL) {
        eta2theta(eta[,1], .llocation, earg = .elocation) +
        eta2theta(eta[,2], .lscale, earg = .escale) *
        digamma(eta2theta(eta[,3], .lshape, earg = .eshape))
    }, list( .llocation = llocation, .lscale = lscale, .lshape = lshape,
             .elocation = elocation, .escale = escale, .eshape = eshape ))),
    last = eval(substitute(expression({
        misc$link = c(location= .llocation, scale= .lscale, shape= .lshape)
        misc$earg = list(location= .elocation, scale= .escale, shape= .eshape)
    }), list( .llocation = llocation, .lscale = lscale, .lshape = lshape,
              .elocation = elocation, .escale = escale, .eshape = eshape ))),
    loglikelihood = eval(substitute(
        function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
        aa = eta2theta(eta[,1], .llocation, earg = .elocation)
        bb = eta2theta(eta[,2], .lscale, earg = .escale)
        kk = eta2theta(eta[,3], .lshape, earg = .eshape)
        if (residuals) stop("loglikelihood residuals not ",
                            "implemented yet") else {
            sum(w * dlgamma(x = y, location=aa, scale=bb, k=kk, log = TRUE))
        }
    }, list( .llocation = llocation, .lscale = lscale, .lshape = lshape,
             .elocation = elocation, .escale = escale, .eshape = eshape ))),
    vfamily = c("lgamma3ff"),
    deriv = eval(substitute(expression({
        a = eta2theta(eta[,1], .llocation, earg = .elocation)
        b = eta2theta(eta[,2], .lscale, earg = .escale)
        k = eta2theta(eta[,3], .lshape, earg = .eshape)
        zedd = (y-a)/b
        dl.da = (exp(zedd) - k) / b
        dl.db = (zedd * (exp(zedd) - k) - 1) / b
        dl.dk = zedd - digamma(k)
        da.deta = dtheta.deta(a, .llocation, earg = .elocation)
        db.deta = dtheta.deta(b, .lscale, earg = .escale)
        dk.deta = dtheta.deta(k, .lshape, earg = .eshape)
        c(w) * cbind(dl.da * da.deta,
                     dl.db * db.deta,
                     dl.dk * dk.deta)
    }), list( .llocation = llocation, .lscale = lscale, .lshape = lshape,
              .elocation = elocation, .escale = escale, .eshape = eshape ))),
    weight = eval(substitute(expression({
        ed2l.da2 = k / b^2
        ed2l.db2 = (1 + k*(trigamma(k+1) + (digamma(k+1))^2)) / b^2
        ed2l.dk2 = trigamma(k)
        ed2l.dadb = (1 + k*digamma(k)) / b^2
        ed2l.dadk = 1 / b
        ed2l.dbdk = digamma(k) / b
        wz = matrix(as.numeric(NA), n, dimm(M))
        wz[,iam(1,1,M)] = ed2l.da2 * da.deta^2
        wz[,iam(2,2,M)] = ed2l.db2 * db.deta^2
        wz[,iam(3,3,M)] = ed2l.dk2 * dk.deta^2
        wz[,iam(1,2,M)] = ed2l.dadb * da.deta * db.deta
        wz[,iam(1,3,M)] = ed2l.dadk * da.deta * dk.deta
        wz[,iam(2,3,M)] = ed2l.dbdk * db.deta * dk.deta
        wz = c(w) * wz
        wz
    }), list( .llocation = llocation, .lscale = lscale, .lshape = lshape,
              .elocation = elocation, .escale = escale, .eshape = eshape ))))
}

 prentice74 = function(llocation = "identity", lscale = "loge", lshape = "identity",
                      elocation = list(), escale = list(), eshape = list(),
                      ilocation = NULL, iscale = NULL, ishape = NULL, zero = 2:3)
{
    if (mode(llocation) != "character" && mode(llocation) != "name")
        llocation = as.character(substitute(llocation))
    if (mode(lscale) != "character" && mode(lscale) != "name")
        lscale = as.character(substitute(lscale))
    if (mode(lshape) != "character" && mode(lshape) != "name")
        lshape = as.character(substitute(lshape))
    if (length(zero) && !is.Numeric(zero, integer = TRUE, posit = TRUE))
        stop("bad input for argument 'zero'")
    if (length(iscale) && !is.Numeric(iscale, posit = TRUE))
        stop("bad input for argument 'iscale'")
    if (!is.list(elocation)) elocation = list()
    if (!is.list(escale)) escale = list()
    if (!is.list(eshape)) eshape = list()

    new("vglmff",
    blurb = c("Log-gamma distribution (Prentice, 1974)",
            " f(y) = |q| * exp(w/q^2 - e^w) / (b*gamma(1/q^2)) ,\n",
            "w=(y-a)*q/b + digamma(1/q^2), location=a, scale=b>0, shape=q\n\n",
            "Links:    ",
            namesof("location", llocation, earg = elocation), ", ",
            namesof("scale", lscale, earg = escale), ", ",
            namesof("shape", lshape, earg = eshape), "\n", "\n",
            "Mean:     a", "\n"),
    constraints = eval(substitute(expression({
        constraints = cm.zero.vgam(constraints, x, .zero, M)
    }), list( .zero = zero ))),
    initialize = eval(substitute(expression({
        if (ncol(cbind(y)) != 1)
            stop("response must be a vector or a one-column matrix")
        predictors.names =
        c(namesof("location", .llocation, earg = .elocation, tag = FALSE),
          namesof("scale", .lscale, earg = .escale, tag = FALSE),
          namesof("shape", .lshape, earg = .eshape, tag = FALSE))
        if (!length(etastart)) {
            sdy = sqrt(var(y))
            k.init = if (length( .ishape)) rep( .ishape, len = length(y)) else {
                skewness = mean((y-mean(y))^3) / sdy^3 # <0 Left Skewed
                rep(-skewness, len = length(y))
            }
            scale.init = if (length( .iscale)) rep( .iscale, len = length(y)) else {
                rep(sdy, len = length(y))
            }
            loc.init = if (length( .iloc)) rep( .iloc, len = length(y)) else {
                rep(median(y), len = length(y))
            }
            etastart = cbind(theta2eta(loc.init, .llocation, earg = .elocation),
                             theta2eta(scale.init, .lscale, earg = .escale),
                             theta2eta(k.init, .lshape, earg = .eshape))
        }
    }), list( .llocation = llocation, .lscale = lscale, .lshape = lshape,
              .elocation = elocation, .escale = escale, .eshape = eshape,
              .iloc=ilocation, .iscale = iscale, .ishape = ishape ))),
    inverse = eval(substitute(function(eta, extra = NULL) {
        eta2theta(eta[,1], .llocation, earg = .elocation)
    }, list( .llocation = llocation, .lscale = lscale, .lshape = lshape,
             .elocation = elocation, .escale = escale, .eshape = eshape ))),
    last = eval(substitute(expression({
        misc$link = c(location= .llocation, scale= .lscale, shape= .lshape)
        misc$earg = list(location= .elocation, scale= .escale, shape= .eshape)
    }), list( .llocation = llocation, .lscale = lscale, .lshape = lshape,
              .elocation = elocation, .escale = escale, .eshape = eshape ))),
    loglikelihood = eval(substitute(
        function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
        a = eta2theta(eta[,1], .llocation, earg = .elocation)
        b = eta2theta(eta[,2], .lscale, earg = .escale)
        k = eta2theta(eta[,3], .lshape, earg = .eshape)
        tmp55 = k^(-2)
        doubw = (y-a)*k/b + digamma(tmp55)
        if (residuals) stop("loglikelihood residuals not ",
                            "implemented yet") else
        sum(w*(log(abs(k)) -log(b) -lgamma(tmp55) + doubw*tmp55 -exp(doubw )))
    }, list( .llocation = llocation, .lscale = lscale, .lshape = lshape,
             .elocation = elocation, .escale = escale, .eshape = eshape ))),
    vfamily = c("prentice74"),
    deriv = eval(substitute(expression({
        a = eta2theta(eta[,1], .llocation, earg = .elocation)
        b = eta2theta(eta[,2], .lscale, earg = .escale)
        k = eta2theta(eta[,3], .lshape, earg = .eshape)
        tmp55 = k^(-2)
        mustar = digamma(tmp55)
        doubw = (y-a)*k/b + mustar
        sigmastar2 = trigamma(tmp55)
        dl.da = k*(exp(doubw) - tmp55) / b
        dl.db = ((doubw - mustar) * (exp(doubw) - tmp55) - 1) / b
        dl.dk = 1/k - 2 * (doubw - mustar) / k^3 - (exp(doubw) - tmp55) *
                ((doubw - mustar) / k - 2 * sigmastar2 / k^3)
        da.deta = dtheta.deta(a, .llocation, earg = .elocation)
        db.deta = dtheta.deta(b, .lscale, earg = .escale)
        dk.deta = dtheta.deta(k, .lshape, earg = .eshape)
        c(w) * cbind(dl.da * da.deta,
                     dl.db * db.deta,
                     dl.dk * dk.deta)
    }), list( .llocation = llocation, .lscale = lscale, .lshape = lshape,
              .elocation = elocation, .escale = escale, .eshape = eshape ))),
    weight = eval(substitute(expression({
        ed2l.da2 = 1 / b^2
        ed2l.db2 = (1 + sigmastar2*tmp55) / b^2
        ed2l.dk2 = tmp55 - 3*sigmastar2*tmp55^2 + 4*sigmastar2*tmp55^4 *
                   (sigmastar2 - k^2)
        ed2l.dadb = k / b^2
        ed2l.dadk = (2*(sigmastar2*tmp55^2 - tmp55) - 1) / b
        ed2l.dbdk = (sigmastar2*tmp55 - 1) / (b*k)
        wz = matrix(as.numeric(NA), n, dimm(M))
        wz[,iam(1,1,M)] = ed2l.da2 * da.deta^2
        wz[,iam(2,2,M)] = ed2l.db2 * db.deta^2
        wz[,iam(3,3,M)] = ed2l.dk2 * dk.deta^2
        wz[,iam(1,2,M)] = ed2l.dadb * da.deta * db.deta
        wz[,iam(1,3,M)] = ed2l.dadk * da.deta * dk.deta
        wz[,iam(2,3,M)] = ed2l.dbdk * db.deta * dk.deta
        wz = c(w) * wz
        wz
    }), list( .llocation = llocation, .lscale = lscale, .lshape = lshape,
              .elocation = elocation, .escale = escale, .eshape = eshape ))))
}



dgengamma = function(x, scale = 1, d = 1, k = 1, log = FALSE) {
    if (!is.logical(log.arg <- log))
        stop("bad input for argument 'log'")
    rm(log)

    if (!is.Numeric(scale, posit = TRUE))
      stop("bad input for argument 'scale'")
    if (!is.Numeric(d, posit = TRUE))
      stop("bad input for argument 'd'")
    if (!is.Numeric(k, posit = TRUE))
      stop("bad input for argument 'k'")
    N = max(length(x), length(scale), length(d), length(k))
    x = rep(x, len=N); scale = rep(scale, len=N);
    d = rep(d, len=N); k = rep(k, len=N); 

    Loglik = rep(log(0), len=N)
    xok = x > 0
    if (any(xok)) {
        zedd = (x[xok]/scale[xok])^d[xok]
        Loglik[xok] = log(d[xok]) + (-d[xok]*k[xok]) * log(scale[xok]) +
                   (d[xok]*k[xok]-1) * log(x[xok]) - zedd - lgamma(k[xok])
    }
    if (log.arg) {
        Loglik
    } else {
        exp(Loglik)
    }
}




pgengamma = function(q, scale = 1, d = 1, k=1) {
    if (!is.Numeric(scale, posit = TRUE))
      stop("bad input for argument 'scale'")
    if (!is.Numeric(d, posit = TRUE))
      stop("bad input for argument 'd'")
    if (!is.Numeric(k, posit = TRUE))
      stop("bad input for argument 'k'")
    z = (q/scale)^d
    pgamma(z, k)
}


qgengamma = function(p, scale = 1, d = 1, k=1) {
    if (!is.Numeric(scale, posit = TRUE))
      stop("bad input for argument 'scale'")
    if (!is.Numeric(d, posit = TRUE))
      stop("bad input for argument 'd'")
    if (!is.Numeric(k, posit = TRUE))
      stop("bad input for argument 'k'")
    q = qgamma(p, k)
    scale * q^(1/d)
}


rgengamma = function(n, scale = 1, d = 1, k=1) {
    if (!is.Numeric(n, posit = TRUE, integ = TRUE, allow = 1)) 
        stop("bad input for 'n'")
    if (!is.Numeric(scale, posit = TRUE))
      stop("bad input for 'scale'")
    if (!is.Numeric(d, posit = TRUE))
      stop("bad input for 'd'")
    if (!is.Numeric(k, posit = TRUE))
      stop("bad input for 'k'")
    y = rgamma(n, k)
    scale * y^(1/d)
}



 gengamma = function(lscale = "loge", ld = "loge", lk = "loge",
                  escale = list(), ed = list(), ek = list(),
                  iscale = NULL, id = NULL, ik = NULL, zero = NULL)
{
    if (mode(lscale) != "character" && mode(lscale) != "name")
        lscale = as.character(substitute(lscale))
    if (mode(ld) != "character" && mode(ld) != "name")
        ld = as.character(substitute(ld))
    if (mode(lk) != "character" && mode(lk) != "name")
        lk = as.character(substitute(lk))
    if (length(zero) && !is.Numeric(zero, integer = TRUE, posit = TRUE))
        stop("bad input for argument 'zero'")
    if (length(iscale) && !is.Numeric(iscale, posit = TRUE))
        stop("bad input for argument 'iscale'")
    if (!is.list(escale)) escale = list()
    if (!is.list(ed)) ed = list()
    if (!is.list(ek)) ek = list()

    new("vglmff",
    blurb = c("Generalized gamma distribution",
         " f(y) = d * b^(-d*k) * y^(d*k-1) * exp(-(y/b)^d) /  gamma(k),\n",
         "scale=b>0, d>0, k>0, y>0\n\n",
         "Links:    ",
         namesof("scale", lscale, earg = escale), ", ",
         namesof("d", ld, earg = ed), ", ",
         namesof("k", lk, earg = ek), "\n", "\n",
         "Mean:     b*k", "\n"),
    constraints = eval(substitute(expression({
        constraints = cm.zero.vgam(constraints, x, .zero, M)
    }), list( .zero = zero ))),
    initialize = eval(substitute(expression({
        if (ncol(cbind(y)) != 1)
            stop("response must be a vector or a one-column matrix")
        if (any(y <= 0)) stop("response must be have positive values only")
        predictors.names = 
            c(namesof("scale", .lscale, earg = .escale, tag = FALSE),
              namesof("d", .ld, earg = .ed, tag = FALSE),
              namesof("k", .lk, earg = .ek, tag = FALSE))
        if (!length(etastart)) {
            b.init = if (length( .iscale)) rep( .iscale, len = length(y)) else {
                rep(mean(y^2) / mean(y), len = length(y))
            }
            k.init = if (length( .ik)) rep( .ik, len = length(y)) else {
                rep(mean(y) / b.init, len = length(y))
            }
            d.init = if (length( .id)) rep( .id, len = length(y)) else {
                rep(digamma(k.init) / mean(log(y/b.init)), len = length(y))
            }
            etastart = cbind(theta2eta(b.init, .lscale, earg = .escale),
                             theta2eta(d.init, .ld, earg = .ed),
                             theta2eta(k.init, .lk, earg = .ek))
        }
    }), list( .lscale = lscale, .ld = ld, .lk = lk,
              .escale = escale, .ed = ed, .ek = ek,
              .iscale = iscale, .id=id, .ik=ik ))),
    inverse = eval(substitute(function(eta, extra = NULL) {
        b = eta2theta(eta[,1], .lscale, earg = .escale)
        k = eta2theta(eta[,3], .lk, earg = .ek)
        b * k
    }, list( .ld = ld, .lscale = lscale, .lk = lk,
             .escale = escale, .ed = ed, .ek = ek ))),
    last = eval(substitute(expression({
        misc$link = c(scale= .lscale, d= .ld, k= .lk)
        misc$earg = list(scale= .escale, d= .ed, k= .ek)
    }), list( .lscale = lscale, .ld = ld, .lk = lk,
              .escale = escale, .ed = ed, .ek = ek ))),
    loglikelihood = eval(substitute(
        function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
        b = eta2theta(eta[,1], .lscale, earg = .escale)
        d = eta2theta(eta[,2], .ld, earg = .ed)
        k = eta2theta(eta[,3], .lk, earg = .ek)
        if (residuals) stop("loglikelihood residuals not ",
                            "implemented yet") else {
            sum(w * dgengamma(x = y, scale=b, d=d, k=k, log = TRUE))
        }
    }, list( .lscale = lscale, .ld = ld, .lk = lk,
             .escale = escale, .ed = ed, .ek = ek ))),
    vfamily = c("gengamma"),
    deriv = eval(substitute(expression({
        b = eta2theta(eta[,1], .lscale, earg = .escale)
        d = eta2theta(eta[,2], .ld, earg = .ed)
        k = eta2theta(eta[,3], .lk, earg = .ek)
        tmp22 = (y/b)^d
        tmp33 = log(y/b)
        dl.db = d * (tmp22 - k) / b
        dl.dd = 1/d + tmp33 * (k - tmp22)
        dl.dk = d * tmp33 - digamma(k)
        db.deta = dtheta.deta(b, .lscale, earg = .escale)
        dd.deta = dtheta.deta(d, .ld, earg = .ed)
        dk.deta = dtheta.deta(k, .lk, earg = .ek)
        c(w) * cbind(dl.db * db.deta,
                     dl.dd * dd.deta,
                     dl.dk * dk.deta)
    }), list( .lscale = lscale, .ld = ld, .lk = lk,
              .escale = escale, .ed = ed, .ek = ek ))),
    weight = eval(substitute(expression({
        ed2l.db2 = k * (d/b)^2
        ed2l.dd2 = (1 + k * (trigamma(k+1) + (digamma(k+1))^2)) / d^2 
        ed2l.dk2 = trigamma(k)
        ed2l.dbdd = -(1 + k*digamma(k)) / b
        ed2l.dbdk = d / b
        ed2l.dddk = -digamma(k) / d
        wz = matrix(as.numeric(NA), n, dimm(M))
        wz[,iam(1,1,M)] = ed2l.db2 * db.deta^2
        wz[,iam(2,2,M)] = ed2l.dd2 * dd.deta^2
        wz[,iam(3,3,M)] = ed2l.dk2 * dk.deta^2
        wz[,iam(1,2,M)] = ed2l.dbdd * db.deta * dd.deta
        wz[,iam(1,3,M)] = ed2l.dbdk * db.deta * dk.deta
        wz[,iam(2,3,M)] = ed2l.dddk * dd.deta * dk.deta
        wz = c(w) * wz
        wz
    }), list( .lscale = lscale, .ld = ld, .lk = lk,
              .escale = escale, .ed = ed, .ek = ek ))))
}


dlog = function(x, prob, log = FALSE) {
    if (!is.logical(log.arg <- log))
        stop("bad input for argument 'log'")
    rm(log)

    if (!is.Numeric(prob, posit = TRUE) || max(prob) >= 1)
        stop("bad input for argument 'prob'")
    N = max(length(x), length(prob))
    if (length(x) != N)
        x = rep(x, len=N)
    if (length(prob) != N)
        prob = rep(prob, len=N)
    ox = !is.finite(x)
    zero = ox | round(x) != x | x < 1
    ans = rep(0.0, len = length(x))
        if (log.arg) {
            ans[ zero] = log(0.0)
            ans[!zero] = x[!zero] * log(prob[!zero]) - log(x[!zero]) -
                         log(-log1p(-prob[!zero]))
        } else {
            ans[!zero] = -(prob[!zero]^(x[!zero])) / (x[!zero] *
                         log1p(-prob[!zero]))
        }
    if (any(ox))
        ans[ox] = NA
    ans
}



plog  = function(q, prob, log.p = FALSE) {
    if (!is.Numeric(q)) stop("bad input for argument 'q'")
    if (!is.Numeric(prob, posit = TRUE) || max(prob) >= 1)
        stop("bad input for argument 'prob'")
    N = max(length(q), length(prob))
    q = rep(q, len=N); prob = rep(prob, len=N);

    bigno = 10
    owen1965 = (q * (1 - prob) > bigno)
    if (specialCase <- any(owen1965)) {
        qqq = q[owen1965]
        ppp = prob[owen1965]
        pqp = qqq * (1 - ppp)
        bigans = (ppp^(1+qqq) / (1-ppp)) * (1/qqq -
                 1 / (            pqp * (qqq-1)) +
                 2 / ((1-ppp)   * pqp * (qqq-1) * (qqq-2)) -
                 6 / ((1-ppp)^2 * pqp * (qqq-1) * (qqq-2) * (qqq-3)) +
                24 / ((1-ppp)^3 * pqp * (qqq-1) * (qqq-2) * (qqq-3) * (qqq-4)))
        bigans = 1 + bigans / log1p(-ppp)
    }

    floorq = pmax(1, floor(q)) # Ensures at least one element per q value
    floorq[owen1965] = 1
    seqq = sequence(floorq)
    seqp = rep(prob, floorq)
    onevector = (seqp^seqq / seqq) / (-log1p(-seqp))
    rlist =  dotC(name = "tyee_C_cum8sum",
                  as.double(onevector), answer = double(N),
                  as.integer(N), as.double(seqq),
                  as.integer(length(onevector)), notok=integer(1))
    if (rlist$notok != 0) stop("error in 'cum8sum'")
    ans = if (log.p) log(rlist$answer) else rlist$answer
    if (specialCase)
        ans[owen1965] = if (log.p) log(bigans) else bigans
    ans[q < 1] = if (log.p) log(0.0) else 0.0
    ans
}





 if (FALSE)
plog = function(q, prob, log.p = FALSE) {
    if (!is.Numeric(q)) stop("bad input for argument 'q'")
    if (!is.Numeric(prob, posit = TRUE) || max(prob) >= 1)
        stop("bad input for argument 'prob'")
    N = max(length(q), length(prob))
    q = rep(q, len=N); prob = rep(prob, len=N);
    ans = q * 0  # Retains names(q)
    if (max(abs(prob-prob[1])) < 1.0e-08) {
        qstar = floor(q)
        temp = if (max(qstar) >= 1) dlog(x=1:max(qstar), 
               prob=prob[1]) else 0*qstar
        unq = unique(qstar)
        for(ii in unq) {
            index = qstar == ii
            ans[index] = if (ii >= 1) sum(temp[1:ii]) else 0
        }
    } else
    for(ii in 1:N) {
        qstar = floor(q[ii])
        ans[ii] = if (qstar >= 1) sum(dlog(x=1:qstar, prob=prob[ii])) else 0
    }
    if (log.p) log(ans) else ans
}







rlog = function(n, prob, Smallno=1.0e-6) {
    if (!is.Numeric(n, posit = TRUE, integ = TRUE))
        stop("bad input for argument 'n'")
    if (!is.Numeric(prob, allow = 1, posit = TRUE) || max(prob) >= 1)
        stop("bad input for argument 'prob'")
    if (!is.Numeric(Smallno, posit = TRUE, allow = 1) || Smallno > 0.01 ||
       Smallno < 2 * .Machine$double.eps)
        stop("bad input for argument 'Smallno'")
    ans = rep(0.0, len = n)

    ptr1 = 1; ptr2 = 0
    a = -1 / log1p(-prob)
    mean = a*prob/(1-prob)    # E(Y)
    sigma = sqrt(a*prob*(1-a*prob)) / (1-prob)   # sd(Y)
    ymax = dlog(x = 1, prob)
    while(ptr2 < n) {
        Lower = 0.5 # A continuity correction is used = 1 - 0.5.
        Upper = mean + 5 * sigma
        while(plog(q=Upper, prob) < 1-Smallno)
            Upper = Upper + sigma
        Upper = Upper + 0.5
        x = round(runif(2*n, min=Lower, max=Upper))
        index = runif(2*n, max=ymax) < dlog(x,prob)
        sindex = sum(index)
        if (sindex) {
            ptr2 = min(n, ptr1 + sindex - 1)
            ans[ptr1:ptr2] = (x[index])[1:(1+ptr2-ptr1)]
            ptr1 = ptr2 + 1
        }
    }
    ans
}








 logff = function(link = "logit", earg = list(), init.c = NULL)
{
    if (length(init.c) &&
       (!is.Numeric(init.c, posit = TRUE) || max(init.c) >= 1))
        stop("init.c must be in (0,1)")
    if (mode(link) != "character" && mode(link) != "name")
        link = as.character(substitute(link))
    if (!is.list(earg)) earg = list()

    new("vglmff",
    blurb = c("Logarithmic distribution f(y) = a * c^y / y, y = 1,2,3,...,\n",
            "            0 < c < 1, a = -1 / log(1-c)  \n\n",
            "Link:    ", namesof("c", link, earg = earg), "\n", "\n",
            "Mean:    a * c / (1 - c)", "\n"),
    initialize = eval(substitute(expression({
        if (ncol(cbind(y)) != 1)
            stop("response must be a vector or a one-column matrix")
        predictors.names = namesof("c", .link, earg = .earg, tag = FALSE) 
        if (!length(etastart)) {
            llfun = function(cc, y, w) {
                a = -1 / log1p(-cc)
                sum(w * (log(a) + y * log(cc) - log(y)))
            }
            c.init = if (length( .init.c )) .init.c else
                getInitVals(gvals=seq(0.05, 0.95, len=9), llfun=llfun, y = y, w = w)
            c.init = rep(c.init, length=length(y))
            etastart = theta2eta(c.init, .link, earg = .earg)
        }
    }), list( .link = link, .earg = earg, .init.c=init.c ))),
    inverse = eval(substitute(function(eta, extra = NULL) {
        cc = eta2theta(eta, .link, earg = .earg)
        a = -1 / log1p(-cc)
        a * cc / (1-cc)
    }, list( .link = link, .earg = earg ))),
    last = eval(substitute(expression({
        misc$link = c(c= .link)
        misc$earg = list(c= .earg)
        misc$expected = TRUE
    }), list( .link = link, .earg = earg ))),
    loglikelihood = eval(substitute(
        function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
        cc = eta2theta(eta, .link, earg = .earg)
        a = -1 / log1p(-cc)
        if (residuals) stop("loglikelihood residuals ",
                            "not implemented yet") else {
            sum(w * dlog(x = y, prob=-expm1(-1/a), log = TRUE))
        }
    }, list( .link = link, .earg = earg ))),
    vfamily = c("logff"),
    deriv = eval(substitute(expression({
        cc = eta2theta(eta, .link, earg = .earg)
        a = -1 / log1p(-cc)
        dl.dc = 1 / ((1-cc) * log1p(-cc)) + y / cc
        dc.deta = dtheta.deta(cc, .link, earg = .earg)
        c(w) * dl.dc * dc.deta
    }), list( .link = link, .earg = earg ))),
    weight = eval(substitute(expression({
        ed2l.dc2 = a * (1 - a * cc) / (cc * (1-cc)^2)
        wz = c(w) * dc.deta^2 * ed2l.dc2
        wz
    }), list( .link = link, .earg = earg ))))
}





 levy = function(delta = NULL, link.gamma = "loge",
                earg = list(), idelta = NULL, igamma = NULL)
{



    delta.known = is.Numeric(delta, allow = 1)
    if (mode(link.gamma) != "character" && mode(link.gamma) != "name")
        link.gamma = as.character(substitute(link.gamma))
    if (!is.list(earg)) earg = list()

    new("vglmff",
    blurb = c("Levy distribution f(y) = sqrt(gamma/(2*pi)) * ",
            "(y-delta)^(-3/2) * \n",
            "          exp(-gamma / (2*(y-delta ))),\n",
            "          delta < y, gamma > 0",
            if (delta.known) paste(", delta = ", delta, ",", sep = ""),
            "\n\n",
            if (delta.known) "Link:    " else "Links:   ",
            namesof("gamma", link.gamma, earg = earg),
            if (! delta.known) 
                c(", ", namesof("delta", "identity", earg = list())),
            "\n\n",
            "Mean:    NA", 
            "\n"),
    initialize = eval(substitute(expression({
        if (ncol(cbind(y)) != 1)
            stop("response must be a vector or a one-column matrix")
        predictors.names = 
            c(namesof("gamma", .link.gamma, earg = .earg, tag = FALSE),
              if ( .delta.known) NULL else 
              namesof("delta", "identity", earg = list(), tag = FALSE))

        if (!length(etastart)) {
            delta.init = if ( .delta.known) {
                           if (min(y,na.rm= TRUE) <= .delta)
                               stop("delta must be < min(y)")
                           .delta 
                         } else {
                           if (length( .idelta)) .idelta else
                               min(y,na.rm= TRUE) - 1.0e-4 *
                               diff(range(y,na.rm= TRUE))
                         }
            gamma.init = if (length( .igamma)) .igamma else
                         median(y - delta.init) # = 1/median(1/(y-delta.init))
            gamma.init = rep(gamma.init, length=length(y))
            etastart = cbind(theta2eta(gamma.init, .link.gamma, earg = .earg),
                             if ( .delta.known) NULL else delta.init)
                             
        }
    }), list( .link.gamma = link.gamma, .earg = earg,
             .delta.known=delta.known,
             .delta=delta,
             .idelta=idelta,
             .igamma=igamma ))),
    inverse = eval(substitute(function(eta, extra = NULL) {
        eta = as.matrix(eta)
        mygamma = eta2theta(eta[,1], .link.gamma, earg = .earg)
        delta = if ( .delta.known) .delta else eta[,2]


        NA * mygamma
    }, list( .link.gamma = link.gamma, .earg = earg,
             .delta.known=delta.known,
             .delta=delta ))),
    last = eval(substitute(expression({
        misc$link = if ( .delta.known) NULL else c(delta = "identity")
        misc$link = c(gamma = .link.gamma, misc$link)
        misc$earg = if ( .delta.known) list(gamma = .earg) else
                    list(gamma = .earg, delta = list())
        if ( .delta.known)
            misc$delta = .delta
    }), list( .link.gamma = link.gamma, .earg = earg,
             .delta.known=delta.known,
             .delta=delta ))),
    loglikelihood = eval(substitute(
        function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
        eta = as.matrix(eta)
        mygamma = eta2theta(eta[,1], .link.gamma, earg = .earg)
        delta = if ( .delta.known) .delta else eta[,2]
        if (residuals) stop("loglikelihood residuals ",
                            "not implemented yet") else
        sum(w * 0.5 * (log(mygamma) -3*log(y-delta) - mygamma / (y-delta )))
    }, list( .link.gamma = link.gamma, .earg = earg,
             .delta.known=delta.known,
             .delta=delta ))),
    vfamily = c("levy"),
    deriv = eval(substitute(expression({
        eta = as.matrix(eta)
        mygamma = eta2theta(eta[,1], .link.gamma, earg = .earg)
        delta = if ( .delta.known) .delta else eta[,2]
        if (! .delta.known)
            dl.ddelta  = (3 - mygamma / (y-delta)) / (2 * (y-delta))
        dl.dgamma = 0.5 * (1 / mygamma - 1 / (y-delta))
        dgamma.deta = dtheta.deta(mygamma, .link.gamma, earg = .earg)
        c(w) * cbind(dl.dgamma * dgamma.deta, 
                     if ( .delta.known) NULL else dl.ddelta)
    }), list( .link.gamma = link.gamma, .earg = earg,
             .delta.known=delta.known,
             .delta=delta ))),
    weight = eval(substitute(expression({
        wz = matrix(as.numeric(NA), n, dimm(M))   # M = if (delta is known) 1 else 2
        wz[,iam(1,1,M)] = 1 * dgamma.deta^2 
        if (! .delta.known) {
            wz[,iam(1,2,M)] =  3 * dgamma.deta
            wz[,iam(2,2,M)] =  21
        }
        wz = c(w) * wz / (2 * mygamma^2) 
        wz
    }), list( .link.gamma = link.gamma, .earg = earg,
             .delta.known=delta.known,
             .delta=delta ))))
}


        

 if (FALSE) 
 stoppa = function(y0,
                  link.alpha = "loge",
                  link.theta = "loge", ealpha = list(), etheta = list(),
                  ialpha = NULL,
                  itheta=1.0,
                  zero = NULL)
{
    if (!is.Numeric(y0, allo=1) || y0 <= 0)
        stop("y0 must be a positive value")

    if (mode(link.alpha) != "character" && mode(link.alpha) != "name")
        link.alpha = as.character(substitute(link.alpha))
    if (mode(link.theta) != "character" && mode(link.theta) != "name")
        link.theta = as.character(substitute(link.theta))
    if (length(zero) && !is.Numeric(zero, integer = TRUE, posit = TRUE))
        stop("bad input for argument 'zero'")
    if (!is.list(ealpha)) ealpha = list()
    if (!is.list(etheta)) etheta = list()

    new("vglmff",
    blurb = c("Stoppa distribution\n\n",
            "Links:    ",
            namesof("alpha", link.alpha, earg = ealpha), ", ", 
            namesof("theta", link.theta, earg = etheta), "\n", 
            "Mean:     theta*y0*beta(1-1/alpha, theta)"),
    constraints = eval(substitute(expression({
        constraints = cm.zero.vgam(constraints, x, .zero, M)
    }), list( .zero = zero ))),
    initialize = eval(substitute(expression({
        predictors.names = 
        c(namesof("alpha", .link.alpha, earg = .ealpha, tag = FALSE),
          namesof("theta", .link.theta, earg = .etheta, tag = FALSE))

        y0 = .y0 
        if (min(y) < y0) stop("y0 must lie in the interval (0, min(y))")
        if (!length( .ialpha) || !length( .itheta)) {
            qvec = c( .25, .5, .75)   # Arbitrary; could be made an argument
            init.theta = if (length( .itheta)) .itheta else 1
            xvec = log1p(-qvec^(1/init.theta))
            fit0 = lsfit(x = xvec, y=log(quantile(y, qvec))-log(y0), intercept = FALSE)
        }

        extra$y0 = y0
        if (!length(etastart)) {
            alpha = rep(if (length( .ialpha)) .ialpha else
                        -1/fit0$coef[1], length = n)
            theta = rep(if (length( .itheta)) .itheta else 1.0, length = n)
            etastart = cbind(theta2eta(alpha, .link.alpha, earg = .ealpha),
                             theta2eta(theta, .link.theta, earg = .etheta))
        }
    }), list( .link.theta = link.theta, .link.alpha = link.alpha,
            .y0=y0,
            .itheta=itheta, .ialpha=ialpha ))),
    inverse = eval(substitute(function(eta, extra = NULL) {
        alpha = eta2theta(eta[,1], .link.alpha, earg = .ealpha)
        theta = eta2theta(eta[,2], .link.theta, earg = .etheta)
        theta * extra$y0 * beta(1-1/alpha, theta)
    }, list( .link.theta = link.theta, .link.alpha = link.alpha ))),
    last = eval(substitute(expression({
        misc$link = c(alpha= .link.alpha, theta= .link.theta)
    }), list( .link.theta = link.theta, .link.alpha = link.alpha ))),
    loglikelihood = eval(substitute(
            function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
        alpha = eta2theta(eta[,1], .link.alpha, earg = .ealpha)
        theta = eta2theta(eta[,2], .link.theta, earg = .etheta)
        if (residuals) stop("loglikelihood residuals ",
                            "not implemented yet") else
        sum(w*(log(theta*alpha) + alpha*log(extra$y0) -(alpha+1)*log(y)+
               (theta-1) * log1p(-(y/extra$y0)^(-alpha))))
    }, list( .link.theta = link.theta, .link.alpha = link.alpha ))),
    vfamily = c("stoppa"),
    deriv = eval(substitute(expression({
        alpha = eta2theta(eta[,1], .link.alpha, earg = .ealpha)
        theta = eta2theta(eta[,2], .link.theta, earg = .etheta)
        temp8  = (y / extra$y0)^(-alpha)
        temp8a = log(temp8)
        temp8b = log1p(-temp8)
        dl.dalpha = 1/alpha - log(y/extra$y0) + (theta-1) * temp8 *
                    log(y / extra$y0) / (1-temp8)
        dl.dtheta = 1/theta + temp8b
        dalpha.deta = dtheta.deta(alpha, .link.alpha, earg = .ealpha)
        dTHETA.deta = dtheta.deta(theta, .link.theta, earg = .etheta)
        c(w) * cbind( dl.dalpha * dalpha.deta,
                      dl.dtheta * dTHETA.deta )
    }), list( .link.theta = link.theta, .link.alpha = link.alpha ))),
    weight = eval(substitute(expression({
        ed2l.dalpha = 1/alpha^2 + theta * (2 * log(extra$y0) * (digamma(2)-
                      digamma(theta+4)) - (trigamma(1) +
                      trigamma(theta+3)) / alpha^3) / (alpha *
                      (theta+1) * (theta+2) / n)
        ed2l.dtheta = 1 / theta^2
        ed2l.dalphatheta = (digamma(2)-digamma(theta+2)) / (alpha*(theta+1))
        wz = matrix(as.numeric(NA), n, dimm(M))  #3=dimm(M)
        wz[,iam(1,1,M)] = ed2l.dalpha * dalpha.deta^2
        wz[,iam(2,2,M)] = ed2l.dtheta * dTHETA.deta^2
        wz[,iam(1,2,M)] = ed2l.dalpha * dTHETA.deta * dalpha.deta
        wz = c(w) * wz
        wz
    }), list( .link.theta = link.theta, .link.alpha = link.alpha ))) )
}




dlino = function(x, shape1, shape2, lambda = 1, log = FALSE) {
    if (!is.logical(log.arg <- log))
        stop("bad input for argument 'log'")
    rm(log)

    loglik =  dbeta(x = x, shape1=shape1, shape2=shape2, log = TRUE) +
              shape1 * log(lambda) -
              (shape1+shape2) * log1p(-(1-lambda)*x)
    if (log.arg) loglik else exp(loglik)
}

plino = function(q, shape1, shape2, lambda=1) {
    if (!is.Numeric(q)) stop("bad input for 'q'")
    if (!is.Numeric(shape1, posit = TRUE)) 
        stop("bad input for argument 'shape1'")
    if (!is.Numeric(shape2, posit = TRUE)) 
        stop("bad input for argument 'shape2'")
    if (!is.Numeric(lambda, posit = TRUE)) 
        stop("bad input for argument 'lambda'")
    pbeta(q=lambda*q/(1 - (1-lambda)*q), shape1=shape1, shape2=shape2)
}

qlino = function(p, shape1, shape2, lambda=1) {
    if (!is.Numeric(p, posit = TRUE) || any(p >= 1)) 
        stop("bad input for argument 'p'")
    if (!is.Numeric(shape1, posit = TRUE)) 
        stop("bad input for argument 'shape1'")
    if (!is.Numeric(lambda, posit = TRUE)) 
        stop("bad input for argument 'lambda'")
    Y = qbeta(p=p, shape1=shape1, shape2=shape2)
    Y / (lambda + (1-lambda)*Y)
}


rlino = function(n, shape1, shape2, lambda=1) {
    if (!is.Numeric(n, posit = TRUE, integ = TRUE, allow = 1)) 
        stop("bad input for argument 'n'")
    if (!is.Numeric(shape1, posit = TRUE)) 
        stop("bad input for argument 'shape1'")
    if (!is.Numeric(shape2, posit = TRUE)) 
        stop("bad input for argument 'shape2'")
    if (!is.Numeric(lambda, posit = TRUE)) 
        stop("bad input for argument 'lambda'")
    Y = rbeta(n = n, shape1=shape1, shape2=shape2)
    Y / (lambda + (1-lambda)*Y)
}



 lino = function(lshape1 = "loge",
                 lshape2 = "loge",
                 llambda = "loge",
                 eshape1 = list(), eshape2 = list(), elambda = list(),
                 ishape1 = NULL, ishape2 = NULL, ilambda = 1, zero = NULL)
{
    if (mode(lshape1) != "character" && mode(lshape1) != "name")
        lshape1 = as.character(substitute(lshape1))
    if (mode(lshape2) != "character" && mode(lshape2) != "name")
        lshape2 = as.character(substitute(lshape2))
    if (mode(llambda) != "character" && mode(llambda) != "name")
        llambda = as.character(substitute(llambda))
    if (length(zero) && !is.Numeric(zero, integer = TRUE, posit = TRUE))
        stop("bad input for argument 'zero'")
    if (!is.Numeric(ilambda, positive = TRUE))
        stop("bad input for argument 'ilambda'")
    if (!is.list(eshape1)) eshape1 = list()
    if (!is.list(eshape2)) eshape2 = list()
    if (!is.list(elambda)) elambda = list()

    new("vglmff",
    blurb = c("Generalized Beta distribution (Libby and Novick, 1982)\n\n",
            "Links:    ",
            namesof("shape1", lshape1, earg = eshape1), ", ", 
            namesof("shape2", lshape2, earg = eshape2), ", ", 
            namesof("lambda", llambda, earg = elambda), "\n", 
            "Mean:     something complicated"),
    constraints = eval(substitute(expression({
        constraints = cm.zero.vgam(constraints, x, .zero, M)
    }), list( .zero = zero ))),
    initialize = eval(substitute(expression({
        predictors.names = 
        c(namesof("shape1", .lshape1, earg = .eshape1, tag = FALSE),
          namesof("shape2", .lshape2, earg = .eshape2, tag = FALSE),
          namesof("lambda", .llambda, earg = .elambda, tag = FALSE))
        if (min(y) <= 0 || max(y) >= 1)
            stop("values of the response must be between 0 and 1 (0,1)")
        if (ncol(cbind(y)) != 1)
            stop("response must be a vector or a one-column matrix")
        if (!length(etastart)) {
            lambda.init = rep(if (length( .ilambda )) .ilambda else 1,
                              length = n)
            sh1.init = if (length( .ishape1 ))
                         rep( .ishape1, length = n) else NULL
            sh2.init = if (length( .ishape2 ))
                         rep( .ishape2, length = n) else NULL
            txY.init = lambda.init * y / (1+lambda.init*y - y)
            mean1 = mean(txY.init)
            mean2 = mean(1/txY.init)
            if (!is.Numeric(sh1.init))
                sh1.init = rep((mean2 - 1) / (mean2 - 1/mean1), length = n)
            if (!is.Numeric(sh2.init))
                sh2.init = rep(sh1.init * (1-mean1) / mean1, length = n)
            etastart = cbind(theta2eta(sh1.init, .lshape1, earg = .eshape1),
                             theta2eta(sh2.init, .lshape2, earg = .eshape2),
                             theta2eta(lambda.init, .llambda, earg = .elambda))
        }
    }), list( .lshape1 = lshape1, .lshape2 = lshape2, .llambda = llambda,
              .eshape1 = eshape1, .eshape2 = eshape2, .elambda = elambda,
              .ishape1=ishape1, .ishape2=ishape2, .ilambda=ilambda ))),
    inverse = eval(substitute(function(eta, extra = NULL) {
        sh1 = eta2theta(eta[,1], .lshape1, earg = .eshape1)
        sh2 = eta2theta(eta[,2], .lshape2, earg = .eshape2)
        lambda = eta2theta(eta[,3], .llambda, earg = .elambda)
        rep(as.numeric(NA), length = nrow(eta))
    }, list( .lshape1 = lshape1, .lshape2 = lshape2, .llambda = llambda,
             .eshape1 = eshape1, .eshape2 = eshape2, .elambda = elambda ))),
    last = eval(substitute(expression({
        misc$link = c(shape1 = .lshape1, shape2 = .lshape2, lambda = .llambda)
        misc$earg = list(shape1 = .eshape1, shape2 = .eshape2, lambda = .elambda)
    }), list( .lshape1 = lshape1, .lshape2 = lshape2, .llambda = llambda,
              .eshape1 = eshape1, .eshape2 = eshape2, .elambda = elambda ))),
    loglikelihood = eval(substitute(
            function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
        sh1 = eta2theta(eta[,1], .lshape1, earg = .eshape1)
        sh2 = eta2theta(eta[,2], .lshape2, earg = .eshape2)
        lambda = eta2theta(eta[,3], .llambda, earg = .elambda)
        if (residuals) stop("loglikelihood residuals ",
                            "not implemented yet") else {
            sum(w * dlino(y, shape1=sh1, shape2=sh2, lambda=lambda, log = TRUE))
        }
    }, list( .lshape1 = lshape1, .lshape2 = lshape2, .llambda = llambda,
             .eshape1 = eshape1, .eshape2 = eshape2, .elambda = elambda ))),
    vfamily = c("lino"),
    deriv = eval(substitute(expression({
        sh1 = eta2theta(eta[,1], .lshape1, earg = .eshape1)
        sh2 = eta2theta(eta[,2], .lshape2, earg = .eshape2)
        lambda = eta2theta(eta[,3], .llambda, earg = .elambda)
        temp1 = log1p(-(1-lambda) * y)
        temp2 = digamma(sh1+sh2)
        dl.dsh1 = log(lambda) + log(y) - digamma(sh1) + temp2 - temp1
        dl.dsh2 = log1p(-y) - digamma(sh2) + temp2 - temp1
        dl.dlambda = sh1/lambda - (sh1+sh2) * y / (1 - (1-lambda) * y)
        dsh1.deta = dtheta.deta(sh1, .lshape1, earg = .eshape1)
        dsh2.deta = dtheta.deta(sh2, .lshape2, earg = .eshape2)
        dlambda.deta = dtheta.deta(lambda, .llambda, earg = .elambda)
        c(w) * cbind( dl.dsh1 * dsh1.deta,
                      dl.dsh2    * dsh2.deta,
                      dl.dlambda * dlambda.deta)
    }), list( .lshape1 = lshape1, .lshape2 = lshape2, .llambda = llambda,
              .eshape1 = eshape1, .eshape2 = eshape2, .elambda = elambda ))),
    weight = eval(substitute(expression({
        temp3 = trigamma(sh1+sh2)
        ed2l.dsh1 = trigamma(sh1) - temp3
        ed2l.dsh2 = trigamma(sh2) - temp3
        ed2l.dlambda2 = sh1 * sh2 / (lambda^2 * (sh1+sh2+1))
        ed2l.dsh1sh2 = -temp3
        ed2l.dsh1lambda = -sh2 / ((sh1+sh2)*lambda)
        ed2l.dsh2lambda =  sh1 / ((sh1+sh2)*lambda)
        wz = matrix(as.numeric(NA), n, dimm(M))  #M==3 means 6=dimm(M)
        wz[,iam(1,1,M)] = ed2l.dsh1 * dsh1.deta^2
        wz[,iam(2,2,M)] = ed2l.dsh2 * dsh2.deta^2
        wz[,iam(3,3,M)] = ed2l.dlambda2 * dlambda.deta^2
        wz[,iam(1,2,M)] = ed2l.dsh1sh2 * dsh1.deta * dsh2.deta
        wz[,iam(1,3,M)] = ed2l.dsh1lambda * dsh1.deta * dlambda.deta
        wz[,iam(2,3,M)] = ed2l.dsh2lambda * dsh2.deta * dlambda.deta
        wz = c(w) * wz
        wz
    }), list( .lshape1 = lshape1, .lshape2 = lshape2, .llambda = llambda,
              .eshape1 = eshape1, .eshape2 = eshape2, .elambda = elambda ))))
}


 genbetaII= function(link.a = "loge",
                     link.scale = "loge",
                     link.p = "loge",
                     link.q = "loge",
                     earg.a = list(), earg.scale = list(),
                     earg.p = list(), earg.q = list(),
                     init.a = NULL,
                     init.scale = NULL,
                     init.p=1.0,
                     init.q=1.0,
                     zero = NULL)
{

    if (mode(link.a) != "character" && mode(link.a) != "name")
        link.a = as.character(substitute(link.a))
    if (mode(link.scale) != "character" && mode(link.scale) != "name")
        link.scale = as.character(substitute(link.scale))
    if (mode(link.p) != "character" && mode(link.p) != "name")
        link.p = as.character(substitute(link.p))
    if (mode(link.q) != "character" && mode(link.q) != "name")
        link.q = as.character(substitute(link.q))
    if (length(zero) && !is.Numeric(zero, integer = TRUE, posit = TRUE))
        stop("bad input for argument 'zero'")
    if (!is.list(earg.a)) earg.a = list()
    if (!is.list(earg.scale)) earg.scale = list()
    if (!is.list(earg.p)) earg.p = list()
    if (!is.list(earg.q)) earg.q = list()

    new("vglmff",
    blurb = c("Generalized Beta II distribution\n\n",
            "Links:    ",
            namesof("a", link.a, earg = earg.a), ", ", 
            namesof("scale", link.scale, earg = earg.scale), ", ", 
            namesof("p", link.p, earg = earg.p), ", ", 
            namesof("q", link.q, earg = earg.q), "\n", 
            "Mean:     scale*gamma(p + 1/a)*gamma(q - 1/a)/(gamma(p)*gamma(q))"),
    constraints = eval(substitute(expression({
        constraints = cm.zero.vgam(constraints, x, .zero, M)
    }), list( .zero = zero ))),
    initialize = eval(substitute(expression({
        predictors.names = 
        c(namesof("a", .link.a, earg = .earg.a, tag = FALSE),
          namesof("scale", .link.scale, earg = .earg.scale, tag = FALSE),
          namesof("p", .link.p, earg = .earg.p, tag = FALSE),
          namesof("q", .link.q, earg = .earg.q, tag = FALSE))

        if (!length( .init.a) || !length( .init.scale )) {
            qvec = c( .25, .5, .75)   # Arbitrary; could be made an argument
            init.q = if (length( .init.q)) .init.q else 1
            xvec = log( (1-qvec)^(-1/ init.q ) - 1 )
            fit0 = lsfit(x = xvec, y=log(quantile(y, qvec )))
        }

        if (!length(etastart)) {
            aa = rep(if (length( .init.a)) .init.a else 1/fit0$coef[2],
                     length = n)
            scale = rep(if (length( .init.scale )) .init.scale else
                        exp(fit0$coef[1]), length = n)
            qq = rep(if (length( .init.q)) .init.q else 1.0, length = n)
            parg = rep(if (length( .init.p)) .init.p else 1.0, length = n)
            etastart = cbind(theta2eta(aa, .link.a, earg = .earg.a),
                             theta2eta(scale, .link.scale, earg = .earg.scale),
                             theta2eta(parg, .link.p, earg = .earg.p),
                             theta2eta(qq, .link.q, earg = .earg.q))
        }
    }), list( .link.a = link.a, .link.scale = link.scale,
              .link.p = link.p, .link.q = link.q,
              .earg.a = earg.a, .earg.scale = earg.scale, 
              .earg.p=earg.p, .earg.q=earg.q,
              .init.a = init.a, .init.scale = init.scale, 
              .init.p=init.p, .init.q=init.q ))),
    inverse = eval(substitute(function(eta, extra = NULL) {
        aa = eta2theta(eta[,1], .link.a, earg = .earg.a)
        scale = eta2theta(eta[,2], .link.scale, earg = .earg.scale)
        parg = eta2theta(eta[,3], .link.p, earg = .earg.p)
        qq = eta2theta(eta[,4], .link.q, earg = .earg.q)
        scale*gamma(parg + 1/aa)*gamma(qq-1/aa)/(gamma(parg)*gamma(qq))
    }, list( .link.a = link.a, .link.scale = link.scale,
              .earg.a = earg.a, .earg.scale = earg.scale, 
              .earg.p=earg.p, .earg.q=earg.q,
             .link.p = link.p, .link.q = link.q ))),
    last = eval(substitute(expression({
        misc$link = c(a= .link.a, scale= .link.scale,
                      p= .link.p, q= .link.q)
        misc$earg = list(a= .earg.a, scale= .earg.scale,
                      p= .earg.p, q= .earg.q)
    }), list( .link.a = link.a, .link.scale = link.scale,
              .earg.a = earg.a, .earg.scale = earg.scale, 
              .earg.p=earg.p, .earg.q=earg.q,
              .link.p = link.p, .link.q = link.q ))),
    loglikelihood = eval(substitute(
            function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
        aa = eta2theta(eta[,1], .link.a, earg = .earg.a)
        scale = eta2theta(eta[,2], .link.scale, earg = .earg.scale)
        parg = eta2theta(eta[,3], .link.p, earg = .earg.p)
        qq = eta2theta(eta[,4], .link.q, earg = .earg.q)
        if (residuals) stop("loglikelihood residuals ",
                            "not implemented yet") else {
            sum(w*(log(aa) + (aa*parg-1)*log(y) - aa*parg*log(scale) +
                   -lbeta(parg, qq) - (parg+qq)*log1p((y/scale)^aa)))
        }
    }, list( .link.a = link.a, .link.scale = link.scale,
              .earg.a = earg.a, .earg.scale = earg.scale, 
              .earg.p=earg.p, .earg.q=earg.q,
            .link.p = link.p, .link.q = link.q ))),
    vfamily = c("genbetaII"),
    deriv = eval(substitute(expression({
        aa = eta2theta(eta[,1], .link.a, earg = .earg.a)
        scale = eta2theta(eta[,2], .link.scale, earg = .earg.scale)
        parg = eta2theta(eta[,3], .link.p, earg = .earg.p)
        qq = eta2theta(eta[,4], .link.q, earg = .earg.q)

        temp1 = log(y/scale)
        temp2 = (y/scale)^aa
        temp3 = digamma(parg + qq)
        temp3a = digamma(parg)
        temp3b = digamma(qq)
        temp4 = log1p(temp2)

        dl.da = 1/aa + parg * temp1 - (parg+qq) * temp1 / (1+1/temp2)
        dl.dscale = (aa/scale) * (-parg + (parg+qq) / (1+1/temp2))
        dl.dp = aa * temp1 + temp3 - temp3a - temp4
        dl.dq = temp3 - temp3b - temp4
        da.deta = dtheta.deta(aa, .link.a, earg = .earg.a)
        dscale.deta = dtheta.deta(scale, .link.scale, earg = .earg.scale)
        dp.deta = dtheta.deta(parg, .link.p, earg = .earg.p)
        dq.deta = dtheta.deta(qq, .link.q, earg = .earg.q)
        c(w) * cbind( dl.da * da.deta,
                      dl.dscale * dscale.deta,
                      dl.dp * dp.deta,
                      dl.dq * dq.deta )
    }), list( .link.a = link.a, .link.scale = link.scale,
              .earg.a = earg.a, .earg.scale = earg.scale, 
              .earg.p=earg.p, .earg.q=earg.q,
              .link.p = link.p, .link.q = link.q ))),
    weight = eval(substitute(expression({
        temp5  = trigamma(parg + qq)
        temp5a = trigamma(parg)
        temp5b = trigamma(qq)
        ed2l.da = (1 + parg+qq + parg * qq * (temp5a + temp5b +
                  (temp3b - temp3a + (parg-qq)/(parg*qq))^2 - 
                  (parg^2 + qq^2) / (parg*qq)^2)) / (aa^2 * (1+parg+qq))
        ed2l.dscale = aa^2 * parg * qq / (scale^2 * (1+parg+qq))
        ed2l.dp = temp5a - temp5
        ed2l.dq = temp5b - temp5
        ed2l.dascale = (parg - qq - parg*qq*(temp3a -temp3b)) /
                       (scale*(1 + parg+qq))
        ed2l.dap= -(qq   * (temp3a -temp3b) -1) / (aa*(parg+qq))
        ed2l.daq= -(parg * (temp3b -temp3a) -1) / (aa*(parg+qq))
        ed2l.dscalep =  aa * qq   / (scale*(parg+qq))
        ed2l.dscaleq = -aa * parg / (scale*(parg+qq))
        ed2l.dpq = -temp5
        wz = matrix(as.numeric(NA), n, dimm(M))  #M==4 means 10=dimm(M)
        wz[,iam(1,1,M)] = ed2l.da * da.deta^2
        wz[,iam(2,2,M)] = ed2l.dscale * dscale.deta^2
        wz[,iam(3,3,M)] = ed2l.dp * dp.deta^2
        wz[,iam(4,4,M)] = ed2l.dq * dq.deta^2
        wz[,iam(1,2,M)] = ed2l.dascale * da.deta * dscale.deta
        wz[,iam(1,3,M)] = ed2l.dap * da.deta * dp.deta
        wz[,iam(1,4,M)] = ed2l.daq * da.deta * dq.deta
        wz[,iam(2,3,M)] = ed2l.dscalep * dscale.deta * dp.deta
        wz[,iam(2,4,M)] = ed2l.dscaleq * dscale.deta * dq.deta
        wz[,iam(3,4,M)] = ed2l.dpq * dp.deta * dq.deta
        wz = c(w) * wz
        wz
    }), list( .link.a = link.a, .link.scale = link.scale,
              .earg.a = earg.a, .earg.scale = earg.scale, 
              .earg.p=earg.p, .earg.q=earg.q,
              .link.p = link.p, .link.q = link.q ))))
}


rsinmad <- function(n, a, scale = 1, q.arg)
    qsinmad(runif(n), a, scale, q.arg)

rlomax <- function(n, scale = 1, q.arg)
    rsinmad(n, a = 1, scale, q.arg)

rfisk <- function(n, a, scale = 1)
    rsinmad(n, a, scale, q.arg=1)

rparalogistic <- function(n, a, scale = 1)
    rsinmad(n, a, scale, a)

rdagum <- function(n, a, scale = 1, p.arg)
    qdagum(runif(n), a, scale = 1, p.arg)

rinvlomax <- function(n, scale = 1, p.arg)
    rdagum(n, a = 1, scale, p.arg)

rinvparalogistic <- function(n, a, scale = 1)
    rdagum(n, a, scale, a)




qsinmad <- function(p, a, scale = 1, q.arg) {
    bad = (p < 0) | (p > 1)
    ans = NA * p
    a = rep(a, len = length(p))[!bad]
    scale = rep(scale, len = length(p))[!bad]
    q = rep(q.arg, len = length(p))[!bad]
    xx = p[!bad]
    ans[!bad] = scale* ((1 - xx)^(-1/q) - 1)^(1/a)
    ans
}

qlomax <- function(p, scale = 1, q.arg)
    qsinmad(p, a = 1, scale, q.arg)

qfisk <- function(p, a, scale = 1)
    qsinmad(p, a, scale, q.arg=1)

qparalogistic <- function(p, a, scale = 1)
    qsinmad(p, a, scale, a)

qdagum <- function(p, a, scale = 1, p.arg) {
    bad = (p < 0) | (p > 1)
    ans = NA * p
    a = rep(a, len = length(p))[!bad]
    scale = rep(scale, len = length(p))[!bad]
    p.arg = rep(p.arg, len = length(p))[!bad]
    xx = p[!bad]
    ans[!bad] = scale* (xx^(-1/p.arg) - 1)^(-1/a)
    ans
}

qinvlomax <- function(p, scale = 1, p.arg)
    qdagum(p, a = 1, scale, p.arg)

qinvparalogistic <- function(p, a, scale = 1)
    qdagum(p, a, scale, a)






psinmad <- function(q, a, scale = 1, q.arg) {
    zero = q <= 0
    a = rep(a, len = length(q))[!zero]
    scale = rep(scale, len = length(q))[!zero]
    q.arg = rep(q.arg, len = length(q))[!zero]
    ans = 0 * q
    xx = q[!zero]
    ans[!zero] = 1 - (1 + (xx/scale)^a)^(-q.arg)
    ans
}

plomax = function(q, scale = 1, q.arg)
    psinmad(q, a = 1, scale, q.arg)

pfisk = function(q, a, scale = 1)
    psinmad(q, a, scale, q.arg=1)

pparalogistic = function(q, a, scale = 1)
    psinmad(q, a, scale, a)



pdagum <- function(q, a, scale = 1, p.arg) {
    zero <- q <= 0
    a <- rep(a, len = length(q))[!zero]
    scale <- rep(scale, len = length(q))[!zero]
    p <- rep(p.arg, len = length(q))[!zero]
    ans <- 0 * q
    xx <- q[!zero]
    ans[!zero] <- (1 + (xx/scale)^(-a))^(-p)
    ans
}

pinvlomax <- function(q, scale = 1, p.arg)
    pdagum(q, a = 1, scale, p.arg)

pinvparalogistic <- function(q, a, scale = 1)
    pdagum(q, a, scale, a)



dsinmad <- function(x, a, scale = 1, q.arg, log = FALSE) {
    if (!is.logical(log.arg <- log))
        stop("bad input for argument 'log'")
    rm(log)
    LLL <- max(length(x), length(a), length(scale), length(q.arg))
    x     <- rep(x,     len = LLL);
    a     <- rep(a,     len = LLL)
    scale <- rep(scale, len = LLL);
    q.arg <- rep(q.arg, len = LLL)

    Loglik <- rep(log(0), len = LLL)
    xok <- (x > 0)  # Avoids evaluating log(x) if x is negative.
    Loglik[xok] <- log(a[xok]) + log(q.arg[xok]) + (a[xok]-1)*log(x[xok]) -
                  a[xok]*log(scale[xok]) -
             (1+q.arg[xok]) * log1p((x[xok]/scale[xok])^a[xok])
    if (log.arg) Loglik else exp(Loglik)
}

dlomax <- function(x, scale = 1, q.arg, log = FALSE)
    dsinmad(x, a = 1, scale, q.arg, log = log)

dfisk <- function(x, a, scale = 1, log = FALSE)
    dsinmad(x, a, scale, q.arg = 1, log = log)

dparalogistic <- function(x, a, scale = 1, log = FALSE)
    dsinmad(x, a, scale, a, log = log)



ddagum <- function(x, a, scale = 1, p.arg, log = FALSE) {
    if (!is.logical(log.arg <- log))
        stop("bad input for argument 'log'")
    rm(log)

    LLL = max(length(x), length(a), length(scale), length(p.arg))
    x = rep(x, len = LLL); a = rep(a, len = LLL)
    scale = rep(scale, len = LLL); p.arg = rep(p.arg, len = LLL)

    Loglik = rep(log(0), len = LLL)
    xok = (x > 0)  # Avoids evaluating log(x) if x is negative.
    Loglik[xok] = log(a[xok]) + log(p.arg[xok]) +
                  (a[xok]*p.arg[xok]-1)*log(x[xok]) -
                  a[xok]*p.arg[xok]*log(scale[xok]) -
             (1+p.arg[xok]) * log1p((x[xok]/scale[xok])^a[xok])
    Loglik[p.arg <= 0] = NaN
    if (log.arg) Loglik else exp(Loglik)

}

dinvlomax <- function(x, scale = 1, p.arg, log = FALSE)
    ddagum(x, a = 1, scale, p.arg, log = log)

dinvparalogistic <- function(x, a, scale = 1, log = FALSE)
    ddagum(x, a, scale, a, log = log)



 sinmad = function(link.a = "loge",
                  link.scale = "loge",
                  link.q = "loge",
                  earg.a = list(), earg.scale = list(), earg.q = list(),
                  init.a = NULL, 
                  init.scale = NULL,
                  init.q=1.0, 
                  zero = NULL)
{

    if (mode(link.a) != "character" && mode(link.a) != "name")
        link.a = as.character(substitute(link.a))
    if (mode(link.scale) != "character" && mode(link.scale) != "name")
        link.scale = as.character(substitute(link.scale))
    if (mode(link.q) != "character" && mode(link.q) != "name")
        link.q = as.character(substitute(link.q))
    if (length(zero) && !is.Numeric(zero, integer = TRUE, posit = TRUE))
        stop("bad input for argument 'zero'")
    if (!is.list(earg.a)) earg.a = list()
    if (!is.list(earg.scale)) earg.scale = list()
    if (!is.list(earg.q)) earg.q = list()

    new("vglmff",
    blurb = c("Singh-Maddala distribution\n\n",
            "Links:    ",
            namesof("a", link.a, earg = earg.a), ", ", 
            namesof("scale", link.scale, earg = earg.scale), ", ", 
            namesof("q", link.q, earg = earg.q), "\n", 
            "Mean:     scale*gamma(1 + 1/a)*gamma(q - 1/a)/gamma(q)"),
    constraints = eval(substitute(expression({
        constraints = cm.zero.vgam(constraints, x, .zero, M)
    }), list( .zero = zero ))),
    initialize = eval(substitute(expression({
        if (ncol(cbind(y)) != 1)
            stop("response must be a vector or a one-column matrix")
        predictors.names = 
            c(namesof("a", .link.a, earg = .earg.a, tag = FALSE),
              namesof("scale", .link.scale, earg = .earg.scale, tag = FALSE),
              namesof("q", .link.q, earg = .earg.q, tag = FALSE))
        parg = 1

        if (!length( .init.a) || !length( .init.scale )) {
            qvec = c( .25, .5, .75)   # Arbitrary; could be made an argument
            init.q = if (length( .init.q)) .init.q else 1
            xvec = log( (1-qvec)^(-1/ init.q ) - 1 )
            fit0 = lsfit(x = xvec, y=log(quantile(y, qvec )))
        }

        if (!length(etastart)) {
            aa = rep(if (length( .init.a)) .init.a else 1/fit0$coef[2],
                     length = n)
            scale = rep(if (length( .init.scale )) .init.scale else
                        exp(fit0$coef[1]), length = n)
            qq = rep(if (length( .init.q)) .init.q else 1.0, length = n)
            etastart = cbind(theta2eta(aa, .link.a, earg = .earg.a),
                             theta2eta(scale, .link.scale, earg = .earg.scale),
                             theta2eta(qq, .link.q, earg = .earg.q))
        }
    }), list( .link.a = link.a, .link.scale = link.scale,
              .link.q = link.q,
              .earg.a = earg.a, .earg.scale = earg.scale, 
              .earg.q=earg.q,
              .init.a = init.a, .init.scale = init.scale, 
              .init.q=init.q ))),
    inverse = eval(substitute(function(eta, extra = NULL) {
        aa = eta2theta(eta[,1], .link.a, earg = .earg.a)
        scale = eta2theta(eta[,2], .link.scale, earg = .earg.scale)
        qq = eta2theta(eta[,3], .link.q, earg = .earg.q)
        scale*gamma(1 + 1/aa)*gamma(qq-1/aa)/(gamma(qq))
    }, list( .link.a = link.a, .link.scale = link.scale,
              .earg.a = earg.a, .earg.scale = earg.scale, 
              .earg.q=earg.q,
             .link.q = link.q ))),
    last = eval(substitute(expression({
        misc$link = c(a= .link.a, scale= .link.scale, q= .link.q)
        misc$earg = list(a= .earg.a, scale= .earg.scale, q= .earg.q)
    }), list( .link.a = link.a, .link.scale = link.scale,
              .earg.a = earg.a, .earg.scale = earg.scale, 
              .earg.q=earg.q,
              .link.q = link.q ))),
    loglikelihood = eval(substitute(
            function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
        aa = eta2theta(eta[,1], .link.a, earg = .earg.a)
        scale = eta2theta(eta[,2], .link.scale, earg = .earg.scale)
        parg = 1
        qq = eta2theta(eta[,3], .link.q, earg = .earg)
        if (residuals) stop("loglikelihood residuals ",
                            "not implemented yet") else {
            sum(w * dsinmad(x = y, a=aa, scale = scale, q.arg=qq, log = TRUE))
        }
    }, list( .link.a = link.a, .link.scale = link.scale, .link.q = link.q,
             .earg.a = earg.a, .earg.scale = earg.scale, .earg.q=earg.q ))),
    vfamily = c("sinmad"),
    deriv = eval(substitute(expression({
        aa = eta2theta(eta[,1], .link.a, earg = .earg.a)
        scale = eta2theta(eta[,2], .link.scale, earg = .earg.scale)
        parg = 1
        qq = eta2theta(eta[,3], .link.q, earg = .earg.q)

        temp1 = log(y/scale)
        temp2 = (y/scale)^aa
        temp3a = digamma(parg)
        temp3b = digamma(qq)

        dl.da = 1/aa + parg * temp1 - (parg+qq) * temp1 / (1+1/temp2)
        dl.dscale = (aa/scale) * (-parg + (parg+qq) / (1+1/temp2))
        dl.dq = digamma(parg + qq) - temp3b - log1p(temp2)
        da.deta = dtheta.deta(aa, .link.a, earg = .earg.a)
        dscale.deta = dtheta.deta(scale, .link.scale, earg = .earg.scale)
        dq.deta = dtheta.deta(qq, .link.q, earg = .earg.q)
        c(w) * cbind( dl.da * da.deta,
                     dl.dscale * dscale.deta,
                     dl.dq * dq.deta )
    }), list( .link.a = link.a, .link.scale = link.scale,
              .earg.a = earg.a, .earg.scale = earg.scale, 
              .earg.q=earg.q,
              .link.q = link.q ))),
    weight = eval(substitute(expression({
        ed2l.da = (1 + parg+qq + parg * qq * (trigamma(parg) + trigamma(qq) +
                  (temp3b - temp3a + (parg-qq)/(parg*qq))^2 - 
                  (parg^2 + qq^2) / (parg*qq)^2)) / (aa^2 * (1+parg+qq))
        ed2l.dscale = aa^2 * parg * qq / (scale^2 * (1+parg+qq))
        ed2l.dq = 1/qq^2
        ed2l.dascale = (parg - qq - parg*qq*(temp3a -temp3b)) /
                       (scale*(1 + parg+qq))
        ed2l.daq= -(parg * (temp3b -temp3a) -1) / (aa*(parg+qq))
        ed2l.dscaleq = -aa * parg / (scale*(parg+qq))
        wz = matrix(as.numeric(NA), n, dimm(M))  #M==3 means 6=dimm(M)
        wz[,iam(1,1,M)] = ed2l.da * da.deta^2
        wz[,iam(2,2,M)] = ed2l.dscale * dscale.deta^2
        wz[,iam(3,3,M)] = ed2l.dq * dq.deta^2
        wz[,iam(1,2,M)] = ed2l.dascale * da.deta * dscale.deta
        wz[,iam(1,3,M)] = ed2l.daq * da.deta * dq.deta
        wz[,iam(2,3,M)] = ed2l.dscaleq * dscale.deta * dq.deta
        wz = c(w) * wz
        wz
    }), list( .link.a = link.a, .link.scale = link.scale,
              .earg.a = earg.a, .earg.scale = earg.scale, 
              .earg.q=earg.q,
              .link.q = link.q ))))
}


 dagum = function(link.a = "loge",
                  link.scale = "loge",
                  link.p = "loge",
                  earg.a = list(), earg.scale = list(), earg.p = list(),
                  init.a = NULL, 
                  init.scale = NULL,
                  init.p=1.0, 
                  zero = NULL)
{

    if (mode(link.a) != "character" && mode(link.a) != "name")
        link.a = as.character(substitute(link.a))
    if (mode(link.scale) != "character" && mode(link.scale) != "name")
        link.scale = as.character(substitute(link.scale))
    if (mode(link.p) != "character" && mode(link.p) != "name")
        link.p = as.character(substitute(link.p))
    if (!is.list(earg.a)) earg.a = list()
    if (!is.list(earg.scale)) earg.scale = list()
    if (!is.list(earg.p)) earg.p = list()

    if (length(zero) && !is.Numeric(zero, integer = TRUE, posit = TRUE))
        stop("bad input for argument 'zero'")

    new("vglmff",
    blurb = c("Dagum distribution\n\n",
            "Links:    ",
            namesof("a", link.a, earg = earg.a), ", ", 
            namesof("scale", link.scale, earg = earg.scale), ", ", 
            namesof("p", link.p, earg = earg.p), "\n", 
            "Mean:     scale*gamma(p + 1/a)*gamma(1 - 1/a)/gamma(p)"),
    constraints = eval(substitute(expression({
        constraints = cm.zero.vgam(constraints, x, .zero, M)
    }), list( .zero = zero ))),
    initialize = eval(substitute(expression({
        if (ncol(cbind(y)) != 1)
            stop("response must be a vector or a one-column matrix")

        predictors.names <-
          c(namesof("a",     .link.a,     earg = .earg.a,     tag = FALSE),
            namesof("scale", .link.scale, earg = .earg.scale, tag = FALSE),
            namesof("p",     .link.p,     earg = .earg.p,     tag = FALSE))

        if (!length( .init.a) || !length( .init.scale )) {
            qvec = c( .25, .5, .75)   # Arbitrary; could be made an argument
            init.p = if (length( .init.p)) .init.p else 1
            xvec = log( qvec^(-1/ init.p ) - 1 )
            fit0 = lsfit(x = xvec, y=log(quantile(y, qvec )))
        }

        if (!length(etastart)) {
            parg = rep(if (length( .init.p)) .init.p else 1.0, length = n)
            aa = rep(if (length( .init.a)) .init.a else -1/fit0$coef[2],
                     length = n)
            scale = rep(if (length( .init.scale )) .init.scale else
                        exp(fit0$coef[1]), length = n)
            etastart = cbind(theta2eta(aa, .link.a, earg = .earg.a),
                             theta2eta(scale, .link.scale, earg = .earg.scale),
                             theta2eta(parg, .link.p, earg = .earg.p))
        }
    }), list( .link.a = link.a, .link.scale = link.scale,
              .link.p = link.p,
              .earg.a = earg.a, .earg.scale = earg.scale, 
              .earg.p=earg.p,
              .init.a = init.a, .init.scale = init.scale, 
              .init.p=init.p ))),
    inverse = eval(substitute(function(eta, extra = NULL) {
        aa = eta2theta(eta[,1], .link.a, earg = .earg.a)
        scale = eta2theta(eta[,2], .link.scale, earg = .earg.scale)
        parg = eta2theta(eta[,3], .link.p, earg = .earg.p)
        qq = 1
        scale*gamma(parg + 1/aa)*gamma(qq-1/aa)/(gamma(parg)*gamma(qq))
    }, list( .link.a = link.a, .link.scale = link.scale,
              .earg.a = earg.a, .earg.scale = earg.scale, 
              .earg.p=earg.p,
             .link.p = link.p ))),
    last = eval(substitute(expression({
        misc$link = c(a= .link.a, scale= .link.scale, p= .link.p )
        misc$earg = list(a= .earg.a, scale= .earg.scale, p= .earg.p)
    }), list( .link.a = link.a, .link.scale = link.scale,
              .earg.a = earg.a, .earg.scale = earg.scale, 
              .earg.p=earg.p,
              .link.p = link.p ))),
    loglikelihood = eval(substitute(
            function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
        aa = eta2theta(eta[,1], .link.a, earg = .earg.a)
        scale = eta2theta(eta[,2], .link.scale, earg = .earg.scale)
        parg = eta2theta(eta[,3], .link.p, earg = .earg.p)
        qq = 1
        if (residuals) stop("loglikelihood residuals ",
                            "not implemented yet") else {
            sum(w * ddagum(x = y, a=aa, scale = scale, p.arg=parg, log = TRUE))
        }
    }, list( .link.a = link.a, .link.scale = link.scale, .link.p = link.p, 
             .earg.a = earg.a, .earg.scale = earg.scale, .earg.p=earg.p ))),
    vfamily = c("dagum"),
    deriv = eval(substitute(expression({
        aa = eta2theta(eta[,1], .link.a, earg = .earg.a)
        scale = eta2theta(eta[,2], .link.scale, earg = .earg.scale)
        parg = eta2theta(eta[,3], .link.p, earg = .earg.p)
        qq = 1

        temp1 = log(y/scale)
        temp2 = (y/scale)^aa
        temp3a = digamma(parg)
        temp3b = digamma(qq)

        dl.da = 1/aa + parg * temp1 - (parg+qq) * temp1 / (1+1/temp2)
        dl.dscale = (aa/scale) * (-parg + (parg+qq) / (1+1/temp2))
        dl.dp = aa * temp1 + digamma(parg + qq) - temp3a - log1p(temp2)
        da.deta = dtheta.deta(aa, .link.a, earg = .earg.a)
        dscale.deta = dtheta.deta(scale, .link.scale, earg = .earg.scale)
        dp.deta = dtheta.deta(parg, .link.p, earg = .earg.p)
        c(w) * cbind( dl.da * da.deta,
                     dl.dscale * dscale.deta,
                     dl.dp * dp.deta )
    }), list( .link.a = link.a, .link.scale = link.scale,
              .earg.a = earg.a, .earg.scale = earg.scale, 
              .earg.p=earg.p,
              .link.p = link.p ))),
    weight = eval(substitute(expression({
        ed2l.da = (1 + parg+qq + parg * qq * (trigamma(parg) + trigamma(qq) + 
                  (temp3b - temp3a + (parg-qq)/(parg*qq))^2 - 
                  (parg^2 + qq^2) / (parg*qq)^2)) / (aa^2 * (1+parg+qq))
        ed2l.dscale = aa^2 * parg * qq / (scale^2 * (1+parg+qq))
        ed2l.dp = 1/parg^2 
        ed2l.dascale = (parg - qq - parg*qq*(temp3a -temp3b)) /
                       (scale*(1 + parg+qq))
        ed2l.dap= -(qq   * (temp3a -temp3b) -1) / (aa*(parg+qq))
        ed2l.dscalep =  aa * qq   / (scale*(parg+qq))
        wz = matrix(as.numeric(NA), n, dimm(M))  #M==3 means 6=dimm(M)
        wz[,iam(1,1,M)] = ed2l.da * da.deta^2
        wz[,iam(2,2,M)] = ed2l.dscale * dscale.deta^2
        wz[,iam(3,3,M)] = ed2l.dp * dp.deta^2
        wz[,iam(1,2,M)] = ed2l.dascale * da.deta * dscale.deta
        wz[,iam(1,3,M)] = ed2l.dap * da.deta * dp.deta
        wz[,iam(2,3,M)] = ed2l.dscalep * dscale.deta * dp.deta
        wz = c(w) * wz
        wz
    }), list( .link.a = link.a, .link.scale = link.scale,
              .earg.a = earg.a, .earg.scale = earg.scale, 
              .earg.p=earg.p,
              .link.p = link.p ))))
}



 betaII = function(link.scale = "loge", link.p = "loge", link.q = "loge",
                   earg.scale = list(), earg.p = list(), earg.q = list(),
                   init.scale = NULL, init.p=1.0, init.q=1.0, zero = NULL)
{

    if (mode(link.scale) != "character" && mode(link.scale) != "name")
        link.scale = as.character(substitute(link.scale))
    if (mode(link.p) != "character" && mode(link.p) != "name")
        link.p = as.character(substitute(link.p))
    if (mode(link.q) != "character" && mode(link.q) != "name")
        link.q = as.character(substitute(link.q))
    if (length(zero) && !is.Numeric(zero, integer = TRUE, posit = TRUE))
        stop("bad input for argument 'zero'")
    if (!is.list(earg.scale)) earg.scale = list()
    if (!is.list(earg.p)) earg.p = list()
    if (!is.list(earg.q)) earg.q = list()

    new("vglmff",
    blurb = c("Beta II distribution\n\n",
            "Links:    ",
            namesof("scale", link.scale, earg = earg.scale), ", ", 
            namesof("p", link.p, earg = earg.p), ", ", 
            namesof("q", link.q, earg = earg.q), "\n", 
            "Mean:     scale*gamma(p + 1)*gamma(q - 1)/(gamma(p)*gamma(q))"),
    constraints = eval(substitute(expression({
        constraints = cm.zero.vgam(constraints, x, .zero, M)
    }), list( .zero = zero ))),
    initialize = eval(substitute(expression({
        if (ncol(cbind(y)) != 1)
            stop("response must be a vector or a one-column matrix")
        predictors.names = 
        c(namesof("scale", .link.scale, earg = .earg.scale, tag = FALSE),
          namesof("p", .link.p, earg = .earg.p, tag = FALSE),
          namesof("q", .link.q, earg = .earg.q, tag = FALSE))

        if (!length( .init.scale )) {
            qvec = c( .25, .5, .75)   # Arbitrary; could be made an argument
            init.q = if (length( .init.q)) .init.q else 1
            xvec = log( (1-qvec)^(-1/ init.q ) - 1 )
            fit0 = lsfit(x = xvec, y=log(quantile(y, qvec )))
        }

        if (!length(etastart)) {
          scale = rep(if (length( .init.scale )) .init.scale else
                      exp(fit0$coef[1]), length = n)
          qq = rep(if (length( .init.q)) .init.q else 1.0, length = n)
          parg = rep(if (length( .init.p)) .init.p else 1.0, length = n)
          etastart = cbind(theta2eta(scale, .link.scale, earg = .earg.scale),
                           theta2eta(parg, .link.p, earg = .earg.p),
                           theta2eta(qq, .link.q, earg = .earg.q))
        }
    }), list( .link.scale = link.scale,
              .link.p = link.p, .link.q = link.q,
              .earg.scale = earg.scale, 
              .earg.p=earg.p, .earg.q=earg.q,
              .init.scale = init.scale, 
              .init.p=init.p, .init.q=init.q ))),
    inverse = eval(substitute(function(eta, extra = NULL) {
        aa = 1
        scale = eta2theta(eta[,1], .link.scale, earg = .earg.scale)
        parg = eta2theta(eta[,2], .link.p, earg = .earg.p)
        qq = eta2theta(eta[,3], .link.q, earg = .earg.q)
        scale*gamma(parg + 1/aa)*gamma(qq-1/aa)/(gamma(parg)*gamma(qq))
    }, list( .link.scale = link.scale,
              .earg.scale = earg.scale, 
              .earg.p=earg.p, .earg.q=earg.q,
             .link.p = link.p, .link.q = link.q ))),
    last = eval(substitute(expression({
        misc$link = c(scale= .link.scale, p= .link.p, q= .link.q)
        misc$earg = list(scale= .earg.scale, p= .earg.p, q= .earg.q)
    }), list( .link.scale = link.scale,
              .earg.scale = earg.scale, 
              .earg.p=earg.p, .earg.q=earg.q,
              .link.p = link.p, .link.q = link.q ))),
    loglikelihood = eval(substitute(
            function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
        aa = 1
        scale = eta2theta(eta[,1], .link.scale, earg = .earg.scale)
        parg = eta2theta(eta[,2], .link.p, earg = .earg.p)
        qq = eta2theta(eta[,3], .link.q, earg = .earg.q)
        if (residuals) stop("loglikelihood residuals ",
                            "not implemented yet") else
            sum(w*(log(aa) + (aa*parg-1)*log(y) - aa*parg*log(scale) +
                  (-lbeta(parg, qq)) - (parg+qq)*log1p((y/scale)^aa)))
    }, list( .link.scale = link.scale,
              .earg.scale = earg.scale, 
              .earg.p=earg.p, .earg.q=earg.q,
             .link.p = link.p, .link.q = link.q ))),
    vfamily = c("betaII"),
    deriv = eval(substitute(expression({
        aa = 1
        scale = eta2theta(eta[,1], .link.scale, earg = .earg.scale)
        parg = eta2theta(eta[,2], .link.p, earg = .earg.p)
        qq = eta2theta(eta[,3], .link.q, earg = .earg.q)

        temp1 = log(y/scale)
        temp2 = (y/scale)^aa
        temp3 = digamma(parg + qq)
        temp3a = digamma(parg)
        temp3b = digamma(qq)
        temp4 = log1p(temp2)

        dl.dscale = (aa/scale) * (-parg + (parg+qq) / (1+1/temp2))
        dl.dp = aa * temp1 + temp3 - temp3a - temp4
        dl.dq = temp3 - temp3b - temp4
        dscale.deta = dtheta.deta(scale, .link.scale, earg = .earg.scale)
        dp.deta = dtheta.deta(parg, .link.p, earg = .earg.p)
        dq.deta = dtheta.deta(qq, .link.q, earg = .earg.q)
        c(w) * cbind( dl.dscale * dscale.deta,
                      dl.dp * dp.deta,
                      dl.dq * dq.deta )
    }), list( .link.scale = link.scale,
              .earg.scale = earg.scale, 
              .earg.p=earg.p, .earg.q=earg.q,
              .link.p = link.p, .link.q = link.q ))),
    weight = eval(substitute(expression({
        temp5  = trigamma(parg + qq)
        ed2l.dscale = aa^2 * parg * qq / (scale^2 * (1+parg+qq))
        ed2l.dp = trigamma(parg) - temp5
        ed2l.dq = trigamma(qq) - temp5
        ed2l.dscalep =  aa * qq   / (scale*(parg+qq))
        ed2l.dscaleq = -aa * parg / (scale*(parg+qq))
        ed2l.dpq = -temp5
        wz = matrix(as.numeric(NA), n, dimm(M))  #M==3 means 6=dimm(M)
        wz[,iam(1,1,M)] = ed2l.dscale * dscale.deta^2
        wz[,iam(2,2,M)] = ed2l.dp * dp.deta^2
        wz[,iam(3,3,M)] = ed2l.dq * dq.deta^2
        wz[,iam(1,2,M)] = ed2l.dscalep * dscale.deta * dp.deta
        wz[,iam(1,3,M)] = ed2l.dscaleq * dscale.deta * dq.deta
        wz[,iam(2,3,M)] = ed2l.dpq * dp.deta * dq.deta
        wz = c(w) * wz
        wz
    }), list( .link.scale = link.scale,
              .earg.scale = earg.scale, 
              .earg.p=earg.p, .earg.q=earg.q,
              .link.p = link.p, .link.q = link.q ))))
}



 lomax = function(link.scale = "loge",
                 link.q = "loge",
                 earg.scale = list(), earg.q = list(),
                 init.scale = NULL,
                 init.q=1.0, 
                 zero = NULL)
{

    if (mode(link.scale) != "character" && mode(link.scale) != "name")
        link.scale = as.character(substitute(link.scale))
    if (mode(link.q) != "character" && mode(link.q) != "name")
        link.q = as.character(substitute(link.q))
    if (length(zero) && !is.Numeric(zero, integer = TRUE, posit = TRUE))
        stop("bad input for argument 'zero'")
    if (!is.list(earg.scale)) earg.scale = list()
    if (!is.list(earg.q)) earg.q = list()

    new("vglmff",
    blurb = c("Lomax distribution\n\n",
            "Links:    ",
            namesof("scale", link.scale, earg = earg.scale), ", ", 
            namesof("q", link.q, earg = earg.q), "\n", 
            "Mean:     scale/(q-1)"),
    constraints = eval(substitute(expression({
        constraints = cm.zero.vgam(constraints, x, .zero, M)
    }), list( .zero = zero ))),
    initialize = eval(substitute(expression({
        if (ncol(cbind(y)) != 1)
            stop("response must be a vector or a one-column matrix")
        predictors.names =
        c(namesof("scale", .link.scale, earg = .earg.scale, tag = FALSE),
          namesof("q", .link.q, earg = .earg.q, tag = FALSE))
        aa = parg = 1

        if (!length( .init.scale )) {
            qvec = c( .25, .5, .75)   # Arbitrary; could be made an argument
            init.q = if (length( .init.q)) .init.q else 1
            xvec = log( (1-qvec)^(-1/ init.q ) - 1 )
            fit0 = lsfit(x = xvec, y=log(quantile(y, qvec )))
        }

        if (!length(etastart)) {
          qq = rep(if (length( .init.q)) .init.q else 1.0, length = n)
          scale = rep(if (length( .init.scale )) .init.scale else
                      exp(fit0$coef[1]), length = n)
          etastart = cbind(theta2eta(scale, .link.scale, earg = .earg.scale),
                           theta2eta(qq, .link.q, earg = .earg.q))
        }
    }), list( .link.scale = link.scale, .link.q = link.q,
              .earg.scale = earg.scale, .earg.q=earg.q,
              .init.scale = init.scale, .init.q=init.q ))),
    inverse = eval(substitute(function(eta, extra = NULL) {
        scale = eta2theta(eta[,1], .link.scale, earg = .earg.scale)
        qq = eta2theta(eta[,2], .link.q, earg = .earg.q)
        scale/(qq-1)
    }, list( .link.scale = link.scale, .link.q = link.q,
             .earg.scale = earg.scale, .earg.q=earg.q ))),
    last = eval(substitute(expression({
        misc$link = c(scale= .link.scale, q= .link.q)
        misc$earg = list(scale= .earg.scale, q= .earg.q)
    }), list( .link.scale = link.scale, .link.q = link.q,
              .earg.scale = earg.scale, .earg.q=earg.q ))),
    loglikelihood = eval(substitute(
            function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
        aa = 1
        scale = eta2theta(eta[,1], .link.scale, earg = .earg.scale)
        parg = 1
        qq = eta2theta(eta[,2], .link.q, earg = .earg.q)
        if (residuals) stop("loglikelihood residuals ",
                            "not implemented yet") else {
            sum(w * dlomax(x = y, scale = scale, q.arg=qq, log = TRUE))
        }
    }, list( .link.scale = link.scale, .link.q = link.q,
             .earg.scale = earg.scale, .earg.q=earg.q ))),
    vfamily = c("lomax"),
    deriv = eval(substitute(expression({
        aa = 1
        scale = eta2theta(eta[,1], .link.scale, earg = .earg.scale)
        parg = 1
        qq = eta2theta(eta[,2], .link.q, earg = .earg.q)
        temp2 = (y/scale)^aa

        dl.dscale = (aa/scale) * (-parg + (parg+qq) / (1+1/temp2))
        dl.dq = digamma(parg + qq) - digamma(qq) - log1p(temp2)
        dscale.deta = dtheta.deta(scale, .link.scale, earg = .earg.scale)
        dq.deta = dtheta.deta(qq, .link.q, earg = .earg.q)
        c(w) * cbind( dl.dscale * dscale.deta,
                      dl.dq * dq.deta )
    }), list( .link.scale = link.scale, .link.q = link.q,
              .earg.scale = earg.scale, .earg.q=earg.q ))),
    weight = eval(substitute(expression({
        ed2l.dscale = aa^2 * parg * qq / (scale^2 * (1+parg+qq))
        ed2l.dq = 1/qq^2 
        ed2l.dscaleq = -aa * parg / (scale*(parg+qq))
        wz = matrix(as.numeric(NA), n, dimm(M))  #M==2 means 3=dimm(M)
        wz[,iam(1,1,M)] = ed2l.dscale * dscale.deta^2
        wz[,iam(2,2,M)] = ed2l.dq * dq.deta^2
        wz[,iam(1,2,M)] = ed2l.dscaleq * dscale.deta * dq.deta
        wz = c(w) * wz
        wz
    }), list( .link.scale = link.scale, .link.q = link.q,
              .earg.scale = earg.scale, .earg.q=earg.q ))))
}


 fisk = function(link.a = "loge",
                 link.scale = "loge",
                 earg.a = list(), earg.scale = list(),
                 init.a = NULL, 
                 init.scale = NULL,
                 zero = NULL)
{

    if (mode(link.a) != "character" && mode(link.a) != "name")
        link.a = as.character(substitute(link.a))
    if (mode(link.scale) != "character" && mode(link.scale) != "name")
        link.scale = as.character(substitute(link.scale))
    if (length(zero) && !is.Numeric(zero, integer = TRUE, posit = TRUE))
        stop("bad input for argument 'zero'")
    if (!is.list(earg.a)) earg.a = list()
    if (!is.list(earg.scale)) earg.scale = list()

    new("vglmff",
    blurb = c("Fisk distribution\n\n",
            "Links:    ",
            namesof("a", link.a, earg = earg.a), ", ", 
            namesof("scale", link.scale, earg = earg.scale), "\n", 
                "Mean:     scale * gamma(1 + 1/a) * gamma(1 - 1/a)"),
    constraints = eval(substitute(expression({
        constraints = cm.zero.vgam(constraints, x, .zero, M)
    }), list( .zero = zero ))),
    initialize = eval(substitute(expression({
        predictors.names =
        c(namesof("a", .link.a, earg = .earg.a, tag = FALSE),
          namesof("scale", .link.scale, earg = .earg.scale, tag = FALSE))
        qq = parg = 1

        if (!length( .init.scale )) {
            qvec = c( .25, .5, .75)   # Arbitrary; could be made an argument
            xvec = log( 1/qvec - 1 )
            fit0 = lsfit(x = xvec, y=log(quantile(y, qvec )))
        }

        if (!length(etastart)) {
          aa = rep(if (length( .init.a)) .init.a else -1/fit0$coef[2],
                   length = n)
          scale = rep(if (length( .init.scale )) .init.scale else
                      exp(fit0$coef[1]), length = n)
          etastart = cbind(theta2eta(aa, .link.a, earg = .earg.a),
                           theta2eta(scale, .link.scale, earg = .earg.scale))
        }
    }), list( .link.a = link.a, .link.scale = link.scale,
              .earg.a = earg.a, .earg.scale = earg.scale, 
              .init.a = init.a, .init.scale = init.scale ))),
    inverse = eval(substitute(function(eta, extra = NULL) {
        aa = eta2theta(eta[,1], .link.a, earg = .earg.a)
        scale = eta2theta(eta[,2], .link.scale, earg = .earg.scale)
        qq = 1
        scale*gamma(1 + 1/aa)*gamma(1-1/aa)
    }, list( .link.a = link.a, .link.scale = link.scale,
              .earg.a = earg.a, .earg.scale = earg.scale ))),
    last = eval(substitute(expression({
        misc$link = c(a= .link.a, scale= .link.scale)
        misc$earg = list(a= .earg.a, scale= .earg.scale)
    }), list( .link.a = link.a, .link.scale = link.scale,
              .earg.a = earg.a, .earg.scale = earg.scale
            ))),
    loglikelihood = eval(substitute(
            function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
        aa = eta2theta(eta[,1], .link.a, earg = .earg)
        scale = eta2theta(eta[,2], .link.scale, earg = .earg.scale)
        parg = qq = 1
        if (residuals) stop("loglikelihood residuals ",
                            "not implemented yet") else {
            sum(w * dfisk(x = y, a=aa, scale = scale, log = TRUE))
        }
    }, list( .link.a = link.a, .link.scale = link.scale,
             .earg.a = earg.a, .earg.scale = earg.scale ))),
    vfamily = c("fisk"),
    deriv = eval(substitute(expression({
        aa = eta2theta(eta[,1], .link.a, earg = .earg.a)
        scale = eta2theta(eta[,2], .link.scale, earg = .earg.scale)
        parg = qq = 1

        temp1 = log(y/scale)
        temp2 = (y/scale)^aa
        temp3a = digamma(parg)
        temp3b = digamma(qq)

        dl.da = 1/aa + parg * temp1 - (parg+qq) * temp1 / (1+1/temp2)
        dl.dscale = (aa/scale) * (-parg + (parg+qq) / (1+1/temp2))
        da.deta = dtheta.deta(aa, .link.a, earg = .earg.a)
        dscale.deta = dtheta.deta(scale, .link.scale, earg = .earg.scale)
        c(w) * cbind( dl.da * da.deta,
                      dl.dscale * dscale.deta )
    }), list( .link.a = link.a, .link.scale = link.scale,
              .earg.a = earg.a, .earg.scale = earg.scale ))),
    weight = eval(substitute(expression({
        ed2l.da = (1 + parg+qq + parg * qq * (trigamma(parg) + trigamma(qq) + 
                  (temp3b - temp3a + (parg-qq)/(parg*qq))^2 - 
                  (parg^2 + qq^2) / (parg*qq)^2)) / (aa^2 * (1+parg+qq))
        ed2l.dscale = aa^2 * parg * qq / (scale^2 * (1+parg+qq))
        ed2l.dascale = (parg - qq - parg*qq*(temp3a -temp3b)) /
                       (scale*(1 + parg+qq))
        wz = matrix(as.numeric(NA), n, dimm(M))  #M==2 means 3=dimm(M)
        wz[,iam(1,1,M)] = ed2l.da * da.deta^2
        wz[,iam(2,2,M)] = ed2l.dscale * dscale.deta^2
        wz[,iam(1,2,M)] = ed2l.dascale * da.deta * dscale.deta
        wz = c(w) * wz
        wz
    }), list( .link.a = link.a, .link.scale = link.scale,
              .earg.a = earg.a, .earg.scale = earg.scale ))))
}


 invlomax = function(link.scale = "loge",
                     link.p = "loge",
                     earg.scale = list(), earg.p = list(),
                     init.scale = NULL,
                     init.p=1.0, 
                     zero = NULL)
{

    if (mode(link.scale) != "character" && mode(link.scale) != "name")
        link.scale = as.character(substitute(link.scale))
    if (mode(link.p) != "character" && mode(link.p) != "name")
        link.p = as.character(substitute(link.p))
    if (length(zero) && !is.Numeric(zero, integer = TRUE, posit = TRUE))
        stop("bad input for argument 'zero'")
    if (!is.list(earg.scale)) earg.scale = list()
    if (!is.list(earg.p)) earg.p = list()

    new("vglmff",
    blurb = c("Inverse Lomax distribution\n\n",
            "Links:    ",
            namesof("scale", link.scale, earg = earg.scale), ", ", 
            namesof("p", link.p, earg = earg.p), "\n", 
            "Mean:     does not exist"),
    constraints = eval(substitute(expression({
        constraints = cm.zero.vgam(constraints, x, .zero, M)
    }), list( .zero = zero ))),
    initialize = eval(substitute(expression({
        if (ncol(cbind(y)) != 1)
            stop("response must be a vector or a one-column matrix")
        predictors.names =
        c(namesof("scale", .link.scale, earg = .earg.scale, tag = FALSE),
          namesof("p", .link.p, earg = .earg.p, tag = FALSE))
        qq = aa = 1

        if (!length( .init.scale )) {
            qvec = c( .25, .5, .75)   # Arbitrary; could be made an argument
            init.p = if (length( .init.p)) .init.p else 1
            xvec = log( qvec^(-1/ init.p ) - 1 )
            fit0 = lsfit(x = xvec, y=log(quantile(y, qvec )))
        }
        if (!length(etastart)) {
          scale = rep(if (length( .init.scale )) .init.scale else
                      exp(fit0$coef[1]), length = n)
          parg = rep(if (length( .init.p)) .init.p else 1.0, length = n)
          etastart = cbind(theta2eta(scale, .link.scale, earg = .earg.scale),
                           theta2eta(parg, .link.p, earg = .earg.p))
        }
    }), list( .link.scale = link.scale,
              .link.p = link.p,
              .earg.scale = earg.scale, 
              .earg.p=earg.p,
              .init.scale = init.scale, 
              .init.p=init.p ))),
    inverse = eval(substitute(function(eta, extra = NULL) {
        rep(as.numeric(NA), len = nrow(eta))
    }, list( .link.scale = link.scale,
              .earg.scale = earg.scale, 
              .earg.p=earg.p,
             .link.p = link.p ))),
    last = eval(substitute(expression({
        misc$link = c(scale= .link.scale, p= .link.p )
        misc$earg = list(scale= .earg.scale, p= .earg.p )
    }), list( .link.scale = link.scale,
              .earg.scale = earg.scale, 
              .earg.p=earg.p,
              .link.p = link.p ))),
    loglikelihood = eval(substitute(
            function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
        aa = qq = 1
        scale = eta2theta(eta[,1], .link.scale, earg = .earg.scale)
        parg = eta2theta(eta[,2], .link.p, earg = .earg.p)
        if (residuals) stop("loglikelihood residuals ",
                            "not implemented yet") else {
           sum(w * dinvlomax(x = y, scale = scale, p.arg=parg, log = TRUE))
        }
    }, list( .link.scale = link.scale, .link.p = link.p,
             .earg.scale = earg.scale, .earg.p=earg.p ))),
    vfamily = c("invlomax"),
    deriv = eval(substitute(expression({
        aa = qq = 1 
        scale = eta2theta(eta[,1], .link.scale, earg = .earg.scale)
        parg = eta2theta(eta[,2], .link.p, earg = .earg.p)

        temp1 = log(y/scale)
        temp2 = (y/scale)^aa

        dl.dscale = (aa/scale) * (-parg + (parg+qq) / (1+1/temp2))
        dl.dp = aa * temp1 + digamma(parg + qq) - digamma(parg) - log1p(temp2)
        dscale.deta = dtheta.deta(scale, .link.scale, earg = .earg.scale)
        dp.deta = dtheta.deta(parg, .link.p, earg = .earg.p)
        c(w) * cbind( dl.dscale * dscale.deta,
                      dl.dp * dp.deta )
    }), list( .link.scale = link.scale, .link.p = link.p,
              .earg.scale = earg.scale, .earg.p=earg.p ))),
    weight = eval(substitute(expression({
        ed2l.dscale = aa^2 * parg * qq / (scale^2 * (1+parg+qq))
        ed2l.dp = 1/parg^2 
        ed2l.dscalep =  aa * qq   / (scale*(parg+qq))
        wz = matrix(as.numeric(NA), n, dimm(M))  #M==2 means 3=dimm(M)
        wz[,iam(1,1,M)] = ed2l.dscale * dscale.deta^2
        wz[,iam(2,2,M)] = ed2l.dp * dp.deta^2
        wz[,iam(1,2,M)] = ed2l.dscalep * dscale.deta * dp.deta
        wz = c(w) * wz
        wz
    }), list( .link.scale = link.scale, .link.p = link.p,
              .earg.scale = earg.scale, .earg.p=earg.p ))))
}


 paralogistic = function(link.a = "loge",
                         link.scale = "loge",
                         earg.a = list(), earg.scale = list(), 
                         init.a=1.0,
                         init.scale = NULL,
                         zero = NULL)
{

    if (mode(link.a) != "character" && mode(link.a) != "name")
        link.a = as.character(substitute(link.a))
    if (mode(link.scale) != "character" && mode(link.scale) != "name")
        link.scale = as.character(substitute(link.scale))
    if (length(zero) && !is.Numeric(zero, integer = TRUE, posit = TRUE))
        stop("bad input for argument 'zero'")
    if (!is.list(earg.a)) earg.a = list()
    if (!is.list(earg.scale)) earg.scale = list()

    new("vglmff",
    blurb = c("Paralogistic distribution\n\n",
            "Links:    ",
            namesof("a", link.a, earg = earg.a), ", ", 
            namesof("scale", link.scale, earg = earg.scale), "\n", 
            "Mean:     scale*gamma(1 + 1/a)*gamma(a - 1/a)/gamma(a)"),
    constraints = eval(substitute(expression({
        constraints = cm.zero.vgam(constraints, x, .zero, M)
    }), list( .zero = zero ))),
    initialize = eval(substitute(expression({
        if (ncol(cbind(y)) != 1)
            stop("response must be a vector or a one-column matrix")
        predictors.names =
        c(namesof("a", .link.a, earg = .earg.a, tag = FALSE),
          namesof("scale", .link.scale, earg = .earg.scale, tag = FALSE))
        parg = 1

        if (!length( .init.a) || !length( .init.scale )) {
            qvec = c( .25, .5, .75)   # Arbitrary; could be made an argument
            init.a = if (length( .init.a)) .init.a else 1
            xvec = log( (1-qvec)^(-1/ init.a ) - 1 )
            fit0 = lsfit(x = xvec, y=log(quantile(y, qvec )))
        }

        if (!length(etastart)) {
          aa = rep(if (length( .init.a)) .init.a else 1/fit0$coef[2],
                   length = n)
          scale = rep(if (length( .init.scale )) .init.scale else
                  exp(fit0$coef[1]), length = n)
          etastart = cbind(theta2eta(aa, .link.a, earg = .earg.a),
                           theta2eta(scale, .link.scale, earg = .earg.scale))
        }
    }), list( .link.a = link.a, .link.scale = link.scale,
              .earg.a = earg.a, .earg.scale = earg.scale, 
              .init.a = init.a, .init.scale = init.scale
              ))),
    inverse = eval(substitute(function(eta, extra = NULL) {
        aa = eta2theta(eta[,1], .link.a, earg = .earg.a)
        scale = eta2theta(eta[,2], .link.scale, earg = .earg.scale)
        qq = aa
        scale*gamma(1 + 1/aa)*gamma(qq-1/aa)/(gamma(qq))
    }, list( .link.a = link.a, .link.scale = link.scale,
             .earg.a = earg.a, .earg.scale = earg.scale ))),
    last = eval(substitute(expression({
        misc$link = c(a= .link.a, scale= .link.scale)
        misc$earg = list(a= .earg.a, scale= .earg.scale )
    }), list( .link.a = link.a, .link.scale = link.scale,
              .earg.a = earg.a, .earg.scale = earg.scale ))),
    loglikelihood = eval(substitute(
            function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
        aa = eta2theta(eta[,1], .link.a, earg = .earg.a)
        scale = eta2theta(eta[,2], .link.scale, earg = .earg.scale)
        parg = 1
        qq = aa
        if (residuals) stop("loglikelihood residuals ",
                            "not implemented yet") else {
            sum(w * dparalogistic(x = y, a=aa, scale = scale, log = TRUE))
        }
    }, list( .link.a = link.a, .link.scale = link.scale,
              .earg.a = earg.a, .earg.scale = earg.scale ))),
    vfamily = c("paralogistic"),
    deriv = eval(substitute(expression({
        aa = eta2theta(eta[,1], .link.a, earg = .earg.a)
        scale = eta2theta(eta[,2], .link.scale, earg = .earg.scale)
        parg = 1
        qq = aa

        temp1 = log(y/scale)
        temp2 = (y/scale)^aa
        temp3a = digamma(parg)
        temp3b = digamma(qq)

        dl.da = 1/aa + parg * temp1 - (parg+qq) * temp1 / (1+1/temp2)
        dl.dscale = (aa/scale) * (-parg + (parg+qq) / (1+1/temp2))
        da.deta = dtheta.deta(aa, .link.a, earg = .earg.a)
        dscale.deta = dtheta.deta(scale, .link.scale, earg = .earg.scale)
        c(w) * cbind( dl.da * da.deta,
                   dl.dscale * dscale.deta)
    }), list( .link.a = link.a, .link.scale = link.scale,
              .earg.a = earg.a, .earg.scale = earg.scale ))),
    weight = eval(substitute(expression({
        ed2l.da = (1 + parg+qq + parg * qq * (trigamma(parg) + trigamma(qq) +
                  (temp3b - temp3a + (parg-qq)/(parg*qq))^2 - 
                  (parg^2 + qq^2) / (parg*qq)^2)) / (aa^2 * (1+parg+qq))
        ed2l.dscale = aa^2 * parg * qq / (scale^2 * (1+parg+qq))
        ed2l.dascale = (parg - qq - parg*qq*(temp3a -temp3b)) /
                       (scale*(1 + parg+qq))
        wz = matrix(as.numeric(NA), n, dimm(M))  #M==2 means 3=dimm(M)
        wz[,iam(1,1,M)] = ed2l.da * da.deta^2
        wz[,iam(2,2,M)] = ed2l.dscale * dscale.deta^2
        wz[,iam(1,2,M)] = ed2l.dascale * da.deta * dscale.deta
        wz = c(w) * wz
        wz
    }), list( .link.a = link.a, .link.scale = link.scale,
              .earg.a = earg.a, .earg.scale = earg.scale ))))
}


 invparalogistic = function(link.a = "loge",
                            link.scale = "loge",
                    earg.a = list(), earg.scale = list(), 
                            init.a=1.0, 
                            init.scale = NULL,
                            zero = NULL)
{

    if (mode(link.a) != "character" && mode(link.a) != "name")
        link.a = as.character(substitute(link.a))
    if (mode(link.scale) != "character" && mode(link.scale) != "name")
        link.scale = as.character(substitute(link.scale))
    if (length(zero) && !is.Numeric(zero, integer = TRUE, posit = TRUE))
        stop("bad input for argument 'zero'")
    if (!is.list(earg.a)) earg.a = list()
    if (!is.list(earg.scale)) earg.scale = list()

    new("vglmff",
    blurb = c("Inverse paralogistic distribution\n\n",
            "Links:    ",
            namesof("a", link.a, earg = earg.a), ", ", 
            namesof("scale", link.scale, earg = earg.scale), "\n", 
               "Mean:     scale*gamma(a + 1/a)*gamma(1 - 1/a)/gamma(a)"),
    constraints = eval(substitute(expression({
        constraints = cm.zero.vgam(constraints, x, .zero, M)
    }), list( .zero = zero ))),
    initialize = eval(substitute(expression({
        if (ncol(cbind(y)) != 1)
            stop("response must be a vector or a one-column matrix")
        predictors.names =
        c(namesof("a", .link.a, earg = .earg.a, tag = FALSE),
          namesof("scale", .link.scale, earg = .earg.scale, tag = FALSE))

        if (!length( .init.a) || !length( .init.scale )) {
            qvec = c( .25, .5, .75)   # Arbitrary; could be made an argument
            init.p = if (length( .init.a)) .init.a else 1
            xvec = log( qvec^(-1/ init.p ) - 1 )
            fit0 = lsfit(x = xvec, y=log(quantile(y, qvec )))
        }

        qq = 1
        if (!length(etastart)) {
          aa = rep(if (length( .init.a)) .init.a else -1/fit0$coef[2],
                   length = n)
          scale = rep(if (length( .init.scale )) .init.scale else
                      exp(fit0$coef[1]), length = n)
          etastart = cbind(theta2eta(aa, .link.a, earg = .earg.a),
                         theta2eta(scale, .link.scale, earg = .earg.scale))
        }
    }), list( .link.a = link.a, .link.scale = link.scale,
              .earg.a = earg.a, .earg.scale = earg.scale,
              .init.a = init.a, .init.scale = init.scale ))),
    inverse = eval(substitute(function(eta, extra = NULL) {
        aa = eta2theta(eta[,1], .link.a, earg = .earg.a)
        scale = eta2theta(eta[,2], .link.scale, earg = .earg.scale)
        parg = aa
        qq = 1
        scale * gamma(parg + 1/aa) *
        gamma(qq - 1/aa) / (gamma(parg) * gamma(qq))
    }, list( .link.a = link.a, .link.scale = link.scale,
             .earg.a = earg.a, .earg.scale = earg.scale ))),
    last = eval(substitute(expression({
        misc$link = c(a= .link.a, scale= .link.scale )
        misc$earg = list(a= .earg.a, scale= .earg.scale )
    }), list( .link.a = link.a, .link.scale = link.scale,
              .earg.a = earg.a, .earg.scale = earg.scale ))),
    loglikelihood = eval(substitute(
            function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
        aa = eta2theta(eta[,1], .link.a, earg = .earg.a)
        scale = eta2theta(eta[,2], .link.scale, earg = .earg.scale)
        parg = aa
        qq = 1
        if (residuals) stop("loglikelihood residuals ",
                            "not implemented yet") else {
            sum(w * dinvparalogistic(x = y, a=aa, scale = scale, log = TRUE))
        }
    }, list( .link.a = link.a, .link.scale = link.scale,
             .earg.a = earg.a, .earg.scale = earg.scale ))),
    vfamily = c("invparalogistic"),
    deriv = eval(substitute(expression({
        aa = eta2theta(eta[,1], .link.a, earg = .earg.a)
        scale = eta2theta(eta[,2], .link.scale, earg = .earg.scale)
        parg = aa 
        qq = 1

        temp1 = log(y/scale)
        temp2 = (y/scale)^aa
        temp3a = digamma(parg)
        temp3b = digamma(qq)

        dl.da = 1/aa + parg * temp1 - (parg+qq) * temp1 / (1+1/temp2)
        dl.dscale = (aa/scale) * (-parg + (parg+qq) / (1+1/temp2))
        da.deta = dtheta.deta(aa, .link.a, earg = .earg.a)
        dscale.deta = dtheta.deta(scale, .link.scale, earg = .earg.scale)
        c(w) * cbind( dl.da * da.deta,
                      dl.dscale * dscale.deta )
    }), list( .link.a = link.a, .link.scale = link.scale,
              .earg.a = earg.a, .earg.scale = earg.scale ))),
    weight = eval(substitute(expression({
        ed2l.da = (1 + parg+qq + parg * qq * (trigamma(parg) + trigamma(qq) +
                  (temp3b - temp3a + (parg-qq)/(parg*qq))^2 - 
                  (parg^2 + qq^2) / (parg*qq)^2)) / (aa^2 * (1+parg+qq))
        ed2l.dscale = aa^2 * parg * qq / (scale^2 * (1+parg+qq))
        ed2l.dascale = (parg - qq - parg*qq*(temp3a -temp3b)) /
                       (scale*(1 + parg+qq))
        wz = matrix(as.numeric(NA), n, dimm(M))  #M==3 means 6=dimm(M)
        wz[,iam(1,1,M)] = ed2l.da * da.deta^2
        wz[,iam(2,2,M)] = ed2l.dscale * dscale.deta^2
        wz[,iam(1,2,M)] = ed2l.dascale * da.deta * dscale.deta
        wz = c(w) * wz
        wz
    }), list( .link.a = link.a, .link.scale = link.scale,
              .earg.a = earg.a, .earg.scale = earg.scale ))))
}



 if (FALSE)
 genlognormal = function(link.sigma = "loge", link.r = "loge",
                        esigma = list(), er = list(),
                        init.sigma = 1, init.r = 1, zero = NULL)
{
warning("2/4/04; doesn't work, possibly because first derivs are ",
        "not continuous (sign() is used). Certainly, the derivs wrt ",
        "mymu are problematic (run with maxit=4:9 and look at weight ",
        "matrices). Possibly fundamentally cannot be estimated by IRLS. ",
        "Pooling doesn't seem to help")

    if (mode(link.sigma) != "character" && mode(link.sigma) != "name")
        link.sigma = as.character(substitute(link.sigma))
    if (mode(link.r) != "character" && mode(link.r) != "name")
        link.r = as.character(substitute(link.r))
    if (length(zero) && !is.Numeric(zero, integer = TRUE, posit = TRUE))
        stop("bad input for argument 'zero'")
    if (!is.list(esigma)) esigma = list()
    if (!is.list(er)) er = list()

    new("vglmff",
    blurb = c("Three-parameter generalized lognormal distribution\n\n",
            "Links:    ",
            "loc; ",
            namesof("sigma", link.sigma, earg = esigma, tag = TRUE), ", ",
            namesof("r",     link.r,     earg = er,     tag = TRUE)),
    constraints = eval(substitute(expression({
        constraints = cm.zero.vgam(constraints, x, .zero, M)
    }), list( .zero = zero ))),
    initialize = eval(substitute(expression({
        if (ncol(cbind(y)) != 1)
            stop("response must be a vector or a one-column matrix")
        predictors.names = 
        c(namesof("loc", "identity", earg = list(), tag = FALSE),
          namesof("sigma", .link.sigma, earg = .esigma, tag = FALSE),
          namesof("r", .link.r, earg = .er, tag = FALSE))

        if (!length( .init.sigma) || !length( .init.r)) {
            init.r = if (length( .init.r)) .init.r else 1
            sigma.init = (0.5 * sum(abs(log(y) - mean(log(y )))^init.r))^(1/init.r)
        }
        if (any(y <= 0)) stop("y must be positive")

        if (!length(etastart)) {
            sigma.init = rep(if (length( .init.sigma)) .init.sigma else
                             sigma.init, len = n)
            r.init = if (length( .init.r)) .init.r else init.r
            etastart = cbind(mu=rep(log(median(y)), len = n),
                             sigma=sigma.init,
                             r = r.init)
        }
    }), list( .link.sigma = link.sigma, .link.r = link.r,
             .init.sigma=init.sigma, .init.r=init.r ))),
    inverse = eval(substitute(function(eta, extra = NULL) {
        mymu = eta2theta(eta[,1], "identity", earg = list())
        sigma = eta2theta(eta[,2], .link.sigma, earg = .esigma)
        r = eta2theta(eta[,3], .link.r, earg = .er)
        r
    }, list( .link.sigma = link.sigma, .link.r = link.r ))),
    last = eval(substitute(expression({
        misc$link = c(loc = "identity", "sigma" = .link.sigma, r = .link.r )
        misc$expected = TRUE
    }), list( .link.sigma = link.sigma, .link.r = link.r ))),
    loglikelihood = eval(substitute(
        function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
        mymu = eta2theta(eta[,1], "identity", earg = list())
        sigma = eta2theta(eta[,2], .link.sigma, earg = .esigma)
        r = eta2theta(eta[,3], .link.r, earg = .er)
        temp89 = (abs(log(y)-mymu)/sigma)^r
        if (residuals) stop("loglikelihood residuals ",
                            "not implemented yet") else
        sum(w * (-log(r^(1/r) * sigma) - lgamma(1+1/r) - temp89/r))
    }, list( .link.sigma = link.sigma, .link.r = link.r ))),
    vfamily = c("genlognormal3"),
    deriv = eval(substitute(expression({
        mymu = eta2theta(eta[,1], "identity", earg = list())
        sigma = eta2theta(eta[,2], .link.sigma, earg = .esigma)
        r = eta2theta(eta[,3], .link.r, earg = .er)
        ss = 1 + 1/r
        temp33 = (abs(log(y)-mymu)/sigma)
        temp33r1 = temp33^(r-1)
        dl.dmymu = temp33r1 * sign(log(y)-mymu) / sigma
        dl.dsigma = (temp33*temp33r1 - 1) / sigma
        dl.dr = (log(r) - 1 + digamma(ss) + temp33*temp33r1)/r^2 -
                temp33r1 * log(temp33r1) / r

        dmymu.deta = dtheta.deta(mymu, "identity", earg = list())
        dsigma.deta = dtheta.deta(sigma, .link.sigma, earg = .esigma)
        dr.deta = dtheta.deta(r, .link.r, earg = .er)
        c(w) * cbind(dl.dmymu * dmymu.deta, 
                     dl.dsigma * dsigma.deta, 
                     dl.dr * dr.deta)
    }), list( .link.sigma = link.sigma, .link.r = link.r ))),
    weight = expression({
        wz = matrix(0, n, 6)  # 5 will have small savings of 1 column
        B = log(r) + digamma(ss)
        ed2l.dmymu2 = (r-1) * gamma(1-1/r) / (sigma^2 * r^(2/r) * gamma(ss))
        ed2l.dsigma2 = r / sigma^2
        ed2l.dr2 = (ss * trigamma(ss) + B^2 - 1) / r^3 
        ed2l.dsigmar = -B / (r * sigma)
        wz[,iam(1,1,M)] = ed2l.dmymu2 * dmymu.deta^2
        wz[,iam(2,2,M)] = ed2l.dsigma2 * dsigma.deta^2
        wz[,iam(3,3,M)] = ed2l.dr2 * dr.deta^2
        wz[,iam(2,3,M)] = ed2l.dsigmar * dsigma.deta * dr.deta
        wz = c(w) * wz
        wz
    }))
}


 betaprime = function(link = "loge", earg = list(), i1=2, i2 = NULL, zero = NULL)
{
    if (mode(link) != "character" && mode(link) != "name")
        link = as.character(substitute(link))
    if (!is.list(earg)) earg = list()

    new("vglmff",
    blurb = c("Beta-prime distribution\n",
            "y^(shape1-1) * (1+y)^(-shape1-shape2) / Beta(shape1,shape2),",
            " y>0, shape1>0, shape2>0\n\n",
            "Links:    ",
            namesof("shape1", link, earg = earg),  ", ",
            namesof("shape2", link, earg = earg), "\n",
            "Mean:     shape1/(shape2-1) provided shape2>1"),
    constraints = eval(substitute(expression({
        constraints = cm.zero.vgam(constraints, x, .zero, M)
    }), list( .zero = zero ))),
    initialize = eval(substitute(expression({
        if (ncol(y <- as.matrix(y)) > 1)
            stop("betaprime cannot handle matrix responses yet")
        if (min(y) <= 0)
            stop("response must be positive")
        predictors.names = c(namesof("shape1", .link, earg = .earg, short = TRUE),
                             namesof("shape2", .link, earg = .earg, short = TRUE))
        if (is.numeric( .i1) && is.numeric( .i2)) {
            vec = c( .i1, .i2)
            vec = c(theta2eta(vec[1], .link, earg = .earg),
                    theta2eta(vec[2], .link, earg = .earg))
            etastart = matrix(vec, n, 2, byrow= TRUE)
        }
        if (!length(etastart)) {
            init1 = if (length( .i1)) rep( .i1, len = n) else rep(1, len = n)
            init2 = if (length( .i2)) rep( .i2, len = n) else 1 + init1 / (y + 0.1)
            etastart = matrix(theta2eta(c(init1, init2), .link, earg = .earg),
                              n,2,byrow = TRUE)
        }
    }), list( .link = link, .earg = earg, .i1=i1, .i2=i2 ))), 
    inverse = eval(substitute(function(eta, extra = NULL) {
        shapes = eta2theta(eta, .link, earg = .earg)
        ifelse(shapes[,2] > 1, shapes[,1]/(shapes[,2]-1), NA)
    }, list( .link = link, .earg = earg ))),
    last = eval(substitute(expression({
        misc$link = c(shape1 = .link, shape2 = .link)
        misc$earg = list(shape1 = .earg, shape2 = .earg)
    }), list( .link = link, .earg = earg ))),
    loglikelihood = eval(substitute(
        function(mu, y, w, residuals = FALSE, eta, extra = NULL){
        shapes = eta2theta(eta, .link, earg = .earg)
        if (residuals) stop("loglikelihood residuals ",
                            "not implemented yet") else {
            sum(w *((shapes[,1]-1) * log(y) - lbeta(shapes[,1], shapes[,2]) -
                   (shapes[,2]+shapes[,1]) * log1p(y)))
        }
    }, list( .link = link, .earg = earg ))),
    vfamily = "betaprime",
    deriv = eval(substitute(expression({
        shapes = eta2theta(eta, .link, earg = .earg)
        dshapes.deta = dtheta.deta(shapes, .link, earg = .earg)
        dl.dshapes = cbind(log(y) - log1p(y) - digamma(shapes[,1]) + 
                           digamma(shapes[,1]+shapes[,2]),
                           - log1p(y) - digamma(shapes[,2]) + 
                           digamma(shapes[,1]+shapes[,2]))
        c(w) * dl.dshapes * dshapes.deta
    }), list( .link = link, .earg = earg ))),
    weight = expression({
        temp2 = trigamma(shapes[,1]+shapes[,2])
        d2l.dshape12 = temp2 - trigamma(shapes[,1])
        d2l.dshape22 = temp2 - trigamma(shapes[,2])
        d2l.dshape1shape2 = temp2

        wz = matrix(as.numeric(NA), n, dimm(M))   #3=dimm(M)
        wz[,iam(1,1,M)] = d2l.dshape12 * dshapes.deta[,1]^2
        wz[,iam(2,2,M)] = d2l.dshape22 * dshapes.deta[,2]^2
        wz[,iam(1,2,M)] = d2l.dshape1shape2 * dshapes.deta[,1] * dshapes.deta[,2]

        -c(w) * wz
    }))
}






dmaxwell = function(x, a, log = FALSE) {
    if (!is.logical(log.arg <- log))
        stop("bad input for argument 'log'")
    rm(log)

    L = max(length(x), length(a))
    x = rep(x, len = L); a = rep(a, len = L);
    logdensity = rep(log(0), len = L)
    xok = (x > 0)
    logdensity[xok] = 0.5 * log(2/pi) + 1.5 * log(a[xok]) +
                      2 * log(x[xok]) - 0.5 * a[xok] * x[xok]^2
    if (log.arg) logdensity else exp(logdensity)
}

pmaxwell = function(q, a) {
    if (any(a <= 0)) stop("argument 'a' must be positive")
    L = max(length(q), length(a)) 
    q = rep(q, len = L); a = rep(a, len = L); 
    ifelse(q > 0, erf(q*sqrt(a/2)) - q*exp(-0.5*a*q^2) * sqrt(2*a/pi), 0)
}

rmaxwell = function(n, a) {
    if (!is.Numeric(n, posit = TRUE, allow = 1)) 
        stop("bad input for argument 'n'")
    if (any(a <= 0)) stop("argument 'a' must be positive")
    sqrt(2 * rgamma(n = n, 1.5) / a)
}

qmaxwell = function(p, a) {
    if (!is.Numeric(p, posit = TRUE) || any(p>=1)) 
        stop("bad input for argument 'p'")
    if (any(a <= 0)) stop("argument 'a' must be positive")
    N = max(length(p), length(a)); p = rep(p, len=N); a = rep(a, len=N)
    sqrt(2 * qgamma(p=p, 1.5) / a)
}


 maxwell = function(link = "loge", earg = list()) {
    if (mode(link) != "character" && mode(link) != "name")
        link = as.character(substitute(link))
    if (!is.list(earg)) earg = list()

    new("vglmff",
    blurb = c("Maxwell distribution f(y) = sqrt(2/pi) * a^(3/2) * y^2 *",
            " exp(-0.5*a*y^2), y>0, a>0\n",
            "Link:    ", namesof("a", link, earg = earg), "\n", "\n",
            "Mean:    sqrt(8 / (a * pi))"),
    initialize = eval(substitute(expression({
        if (ncol(cbind(y)) != 1)
            stop("response must be a vector or a one-column matrix")
        predictors.names = namesof("a", .link, earg = .earg, tag = FALSE) 
        if (!length(etastart)) {
            a.init = rep(8 / (pi*(y+0.1)^2), length=length(y))
            etastart = theta2eta(a.init, .link, earg = .earg)
        }
    }), list( .link = link, .earg = earg ))),
    inverse = eval(substitute(function(eta, extra = NULL) {
        a = eta2theta(eta, .link, earg = .earg)
        sqrt(8 / (a * pi))
    }, list( .link = link, .earg = earg ))),
    last = eval(substitute(expression({
        misc$link = c(a= .link)
        misc$earg = list(a = .earg)
    }), list( .link = link, .earg = earg ))),
    loglikelihood = eval(substitute(
        function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
        aa = eta2theta(eta, .link, earg = .earg)
        if (residuals) stop("loglikelihood residuals ",
                            "not implemented yet") else
            sum(w * dmaxwell(x = y, a=aa, log = TRUE))
    }, list( .link = link, .earg = earg ))),
    vfamily = c("maxwell"),
    deriv = eval(substitute(expression({
        a = eta2theta(eta, .link, earg = .earg)
        dl.da = 1.5 / a - 0.5 * y^2
        da.deta = dtheta.deta(a, .link, earg = .earg)
        c(w) * dl.da * da.deta
    }), list( .link = link, .earg = earg ))),
    weight = eval(substitute(expression({
        ed2l.da2 = 1.5 / a^2
        wz = c(w) * da.deta^2 * ed2l.da2
        wz
    }), list( .link = link, .earg = earg ))))
}




dnaka = function(x, shape, scale = 1, log = FALSE) {
    if (!is.logical(log.arg <- log))
        stop("bad input for argument 'log'")
    rm(log)
    L = max(length(x), length(shape), length(scale))
    x = rep(x, len = L); shape = rep(shape, len = L); scale = rep(scale, len = L);

    logdensity = rep(log(0), len = L)
    xok = (x > 0)
    logdensity[xok] = dgamma(x = x[xok]^2, shape = shape[xok],
                             scale = scale[xok]/shape[xok], log = TRUE) +
                      log(2) + log(x[xok])
    if (log.arg) logdensity else exp(logdensity)
}


pnaka = function(q, shape, scale = 1) {
    if (!is.Numeric(q))
        stop("bad input for argument 'q'")
    if (!is.Numeric(shape, posit = TRUE))
        stop("bad input for argument 'shape'")
    if (!is.Numeric(scale, posit = TRUE))
        stop("bad input for argument 'scale'")
    L = max(length(q), length(shape), length(scale))
    q = rep(q, len = L); shape = rep(shape, len = L); scale = rep(scale, len = L);
    ifelse(q <= 0, 0, pgamma(shape * q^2 / scale, shape))
}


qnaka = function(p, shape, scale = 1, ...) {
    if (!is.Numeric(p, posit = TRUE) || max(p) >= 1)
        stop("bad input for argument 'p'")
    if (!is.Numeric(shape, posit = TRUE))
        stop("bad input for argument 'shape'")
    if (!is.Numeric(scale, posit = TRUE))
        stop("bad input for argument 'scale'")
    L = max(length(p), length(shape), length(scale))
    p = rep(p, len = L); shape = rep(shape, len = L);
    scale = rep(scale, len = L);
    ans = rep(0.0, len = L)
    myfun = function(x, shape, scale = 1, p)
        pnaka(q=x, shape = shape, scale = scale) - p
    for(ii in 1:L) {
        EY = sqrt(scale[ii]/shape[ii]) * gamma(shape[ii]+0.5) / gamma(shape[ii])
        Upper = 5 * EY
        while(pnaka(q=Upper, shape = shape[ii], scale = scale[ii]) < p[ii])
            Upper = Upper + scale[ii]
        ans[ii] = uniroot(f = myfun, lower = 0, upper = Upper,
                          shape = shape[ii], scale = scale[ii],
                          p = p[ii], ...)$root
    }
    ans
}


rnaka = function(n, shape, scale = 1, Smallno=1.0e-6) {
    if (!is.Numeric(n, posit = TRUE, integ = TRUE))
        stop("bad input for argument 'n'")
    if (!is.Numeric(scale, posit = TRUE, allow = 1))
        stop("bad input for argument 'scale'")
    if (!is.Numeric(shape, posit = TRUE, allow = 1))
        stop("bad input for argument 'shape'")
    if (!is.Numeric(Smallno, posit = TRUE, allow = 1) || Smallno > 0.01 ||
       Smallno < 2 * .Machine$double.eps)
        stop("bad input for argument 'Smallno'")
    ans = rep(0.0, len = n)

    ptr1 = 1; ptr2 = 0
    ymax = dnaka(x = sqrt(scale * (1 - 0.5 / shape)),
                 shape = shape, scale = scale)
    while(ptr2 < n) {
        EY = sqrt(scale / shape) * gamma(shape + 0.5) / gamma(shape)
        Upper = EY + 5 * scale
        while(pnaka(q=Upper, shape = shape, scale = scale) < 1-Smallno)
            Upper = Upper + scale
        x = runif(2*n, min=0, max=Upper)
        index = runif(2*n, max=ymax) < dnaka(x, shape = shape,
                                             scale = scale)
        sindex = sum(index)
        if (sindex) {
            ptr2 = min(n, ptr1 + sindex - 1)
            ans[ptr1:ptr2] = (x[index])[1:(1+ptr2-ptr1)]
            ptr1 = ptr2 + 1
        }
    }
    ans
}






 nakagami = function(lshape = "loge", lscale = "loge",
                     eshape = list(), escale = list(), ishape = NULL, iscale = 1) {
    if (mode(lshape) != "character" && mode(lshape) != "name")
        lshape = as.character(substitute(lshape))
    if (mode(lscale) != "character" && mode(lscale) != "name")
        lscale = as.character(substitute(lscale))
    if (!is.null(iscale) && !is.Numeric(iscale, positi = TRUE))
        stop("argument 'iscale' must be a positive number or NULL")
    if (!is.list(eshape)) eshape = list()
    if (!is.list(escale)) escale = list()

    new("vglmff",
    blurb = c("Nakagami distribution f(y) = 2 * (shape/scale)^shape *\n",
            "                             ",
            "y^(2*shape-1) * exp(-shape*y^2/scale) / gamma(shape),\n",
            "                             ",
            "y>0, shape>0, scale>0\n",
            "Links:    ",
            namesof("shape", lshape, earg = eshape), ", ",
            namesof("scale", lscale, earg = escale),
            "\n",
            "\n",
            "Mean:    sqrt(scale/shape) * gamma(shape+0.5) / gamma(shape)"),
    initialize = eval(substitute(expression({
        if (ncol(cbind(y)) != 1)
            stop("response must be a vector or a one-column matrix")
        predictors.names = c(namesof("shape", .lshape, earg = .eshape, tag = FALSE),
                             namesof("scale", .lscale, earg = .escale, tag = FALSE))
        if (!length(etastart)) {
            init2 = if (is.Numeric( .iscale, posit = TRUE))
                        rep( .iscale, len = n) else rep(1, len = n)
            init1 = if (is.Numeric( .ishape, posit = TRUE))
                        rep( .ishape, len = n) else
                    rep(init2 / (y+1/8)^2, len = n)
            etastart = cbind(theta2eta(init1, .lshape, earg = .eshape),
                             theta2eta(init2, .lscale, earg = .escale))
        }
    }), list( .lscale = lscale, .lshape = lshape,
              .escale = escale, .eshape = eshape,
              .ishape = ishape, .iscale = iscale ))),
    inverse = eval(substitute(function(eta, extra = NULL) {
        shape = eta2theta(eta[,1], .lshape, earg = .eshape)
        scale = eta2theta(eta[,2], .lscale, earg = .escale)
        sqrt(scale/shape) * gamma(shape+0.5) / gamma(shape)
    }, list( .lscale = lscale, .lshape = lshape,
             .escale = escale, .eshape = eshape ))),
    last = eval(substitute(expression({
        misc$link = c(shape= .lshape, scale= .lscale)
        misc$earg = list(shape = .eshape, scale = .escale)
        misc$expected = TRUE
    }), list( .lscale = lscale, .lshape = lshape,
              .escale = escale, .eshape = eshape ))),
    loglikelihood = eval(substitute(
        function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
        shape = eta2theta(eta[,1], .lshape, earg = .eshape)
        scale = eta2theta(eta[,2], .lscale, earg = .escale)
        if (residuals) stop("loglikelihood residuals ",
                            "not implemented yet") else
            sum(w * dnaka(x = y, shape = shape, scale = scale, log = TRUE))
    }, list( .lscale = lscale, .lshape = lshape,
             .escale = escale, .eshape = eshape ))),
    vfamily = c("nakagami"),
    deriv = eval(substitute(expression({
        shape = eta2theta(eta[,1], .lshape, earg = .eshape)
        Scale = eta2theta(eta[,2], .lscale, earg = .escale)
        dl.dshape = 1 + log(shape/Scale) - digamma(shape) +
                    2 * log(y) - y^2 / Scale
        dl.dscale = -shape/Scale + shape * (y/Scale)^2
        dshape.deta = dtheta.deta(shape, .lshape, earg = .eshape)
        dscale.deta = dtheta.deta(Scale, .lscale, earg = .escale)
        c(w) * cbind(dl.dshape * dshape.deta,
                     dl.dscale * dscale.deta)
    }), list( .lscale = lscale, .lshape = lshape,
              .escale = escale, .eshape = eshape ))),
    weight = eval(substitute(expression({
        d2l.dshape2 = trigamma(shape) - 1/shape
        d2l.dscale2 = shape / Scale^2
        wz = matrix(as.numeric(NA), n, M)  # diagonal
        wz[,iam(1,1,M)] = d2l.dshape2 * dshape.deta^2
        wz[,iam(2,2,M)] = d2l.dscale2 * dscale.deta^2
        c(w) * wz
    }), list( .lscale = lscale, .lshape = lshape,
              .escale = escale, .eshape = eshape ))))
}



drayleigh = function(x, scale = 1, log = FALSE) {
  if (!is.logical(log.arg <- log))
        stop("bad input for argument 'log'")
  rm(log)

  L = max(length(x), length(scale))
  x = rep(x, len = L); scale = rep(scale, len = L);
  logdensity = rep(log(0), len = L)
  xok = (x > 0)
  logdensity[xok] = log(x[xok]) - 0.5 * (x[xok]/scale[xok])^2 -
                    2*log(scale[xok])
  if (log.arg) logdensity else exp(logdensity)
}


prayleigh = function(q, scale = 1) {
  if (any(scale <= 0))
    stop("argument 'scale' must be positive")
  L = max(length(q), length(scale)) 
  q = rep(q, len = L); scale = rep(scale, len = L);
  ifelse(q > 0,  -expm1(-0.5*(q/scale)^2), 0)
}


qrayleigh = function(p, scale = 1) {
  if (any(scale <= 0))
    stop("argument 'scale' must be positive")
  if (any(p <= 0) || any(p >= 1))
    stop("argument 'p' must be between 0 and 1")
  scale * sqrt(-2 * log1p(-p))
}


rrayleigh = function(n, scale = 1) {
  if (any(scale <= 0))
    stop("argument 'scale' must be positive")
  scale * sqrt(-2 * log(runif(n)))
}



 rayleigh = function(lscale = "loge", escale = list(), nrfs = 1 / 3 + 0.01) {
    if (mode(lscale) != "character" && mode(lscale) != "name")
        lscale = as.character(substitute(lscale))
    if (!is.list(escale)) escale = list()
    if (!is.Numeric(nrfs, allow = 1) || nrfs<0 || nrfs > 1)
        stop("bad input for 'nrfs'")

    new("vglmff",
    blurb = c("Rayleigh distribution\n\n",
            "f(y) = y*exp(-0.5*(y/scale)^2)/scale^2, y>0, scale>0\n\n",
            "Link:    ",
            namesof("scale", lscale, earg = escale), "\n\n",
            "Mean:    scale * sqrt(pi / 2)"),
    initialize = eval(substitute(expression({
        if (ncol(cbind(y)) != 1)
            stop("response must be a vector or a one-column matrix")

        predictors.names =
          namesof("scale", .lscale, earg = .escale, tag = FALSE) 

        if (!length(etastart)) {
            b.init = (y + 1/8) / sqrt(pi/2)
            etastart = theta2eta(b.init, .lscale, earg = .escale)
        }
    }), list( .lscale = lscale, .escale = escale ))),
    inverse = eval(substitute(function(eta, extra = NULL) {
        Scale = eta2theta(eta, .lscale, earg = .escale)
        Scale * sqrt(pi/2)
    }, list( .lscale = lscale, .escale = escale ))),
    last = eval(substitute(expression({
        misc$link =    c(scale = .lscale)
        misc$earg = list(scale = .escale)
        misc$nrfs = .nrfs
    }), list( .lscale = lscale, .escale = escale, .nrfs = nrfs  ))),
    loglikelihood = eval(substitute(
        function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
        Scale = eta2theta(eta, .lscale, earg = .escale)
        if (residuals) stop("loglikelihood residuals ",
                            "not implemented yet") else {
            sum(w * drayleigh(x = y, scale = Scale, log = TRUE))
        }
    }, list( .lscale = lscale, .escale = escale ))),
    vfamily = c("rayleigh"),
    deriv = eval(substitute(expression({
        Scale = eta2theta(eta, .lscale, earg = .escale)
        dl.dScale = ((y/Scale)^2 - 2) / Scale
        dScale.deta = dtheta.deta(Scale, .lscale, earg = .escale)
        c(w) * dl.dScale * dScale.deta
    }), list( .lscale = lscale, .escale = escale ))),
    weight = eval(substitute(expression({
        d2l.dScale2 = (3 * (y/Scale)^2 - 2) / Scale^2
        ed2l.dScale2 = 4 / Scale^2
        wz = c(w) * dScale.deta^2 *
             ((1 - .nrfs) * d2l.dScale2 + .nrfs * ed2l.dScale2)
        wz
    }), list( .lscale = lscale, .escale = escale, .nrfs = nrfs ))))
}





dparetoIV = function(x, location = 0, scale = 1, inequality = 1, shape = 1, log = FALSE) {
    if (!is.logical(log.arg <- log))
        stop("bad input for argument 'log'")
    rm(log)

    N = max(length(x), length(location), length(scale), length(inequality),
            length(shape))
    x = rep(x, len=N); location = rep(location, len=N)
    scale = rep(scale, len=N); inequality = rep(inequality, len=N)
    shape = rep(shape, len=N)

    logdensity = rep(log(0), len=N)
    xok = (x > location)
    zedd = (x - location) / scale
    logdensity[xok] = log(shape[xok]) - log(scale[xok]) -  log(inequality[xok])+
                      (1/inequality[xok]-1) * log(zedd[xok]) - 
                      (shape[xok]+1) * log1p(zedd[xok]^(1/inequality[xok]))
    if (log.arg) logdensity else exp(logdensity)
}

pparetoIV = function(q, location = 0, scale = 1, inequality = 1, shape=1) {
    if (!is.Numeric(q)) stop("bad input for argument 'q'")
    if (!is.Numeric(scale, posit = TRUE)) 
        stop("bad input for argument 'scale'")
    if (!is.Numeric(inequality, posi = TRUE)) 
        stop("bad input for argument 'inequality'")
    if (!is.Numeric(shape, posit = TRUE)) 
        stop("bad input for argument 'shape'")
    N = max(length(q), length(location), length(scale), length(inequality),
            length(shape))
    q = rep(q, len=N); location = rep(location, len=N)
    scale = rep(scale, len=N); inequality = rep(inequality, len=N)
    shape = rep(shape, len=N)
    answer = q * 0
    ii = q > location
    zedd = (q[ii] - location[ii]) / scale[ii]
    answer[ii] = 1 - (1 + zedd^(1/inequality[ii]))^(-shape[ii])
    answer
}

qparetoIV = function(p, location = 0, scale = 1, inequality = 1, shape=1) {
    if (!is.Numeric(p, posit = TRUE) || any(p >= 1)) 
        stop("bad input for argument 'p'")
    if (!is.Numeric(scale, posit = TRUE)) 
        stop("bad input for argument 'scale'")
    if (!is.Numeric(inequality, posi = TRUE)) 
        stop("bad input for argument 'inequality'")
    if (!is.Numeric(shape, posit = TRUE)) 
        stop("bad input for argument 'shape'")
    location + scale * (-1 + (1-p)^(-1/shape))^inequality
}

rparetoIV = function(n, location = 0, scale = 1, inequality = 1, shape=1) {
    if (!is.Numeric(n, posit = TRUE, integ = TRUE, allow = 1)) 
        stop("bad input for argument n")
    if (!is.Numeric(scale, posit = TRUE)) stop("bad input for argument 'scale'")
    if (!is.Numeric(inequality, posi = TRUE)) 
        stop("bad input for argument 'inequality'")
    if (!is.Numeric(shape, posit = TRUE)) stop("bad input for argument 'shape'")
    location + scale * (-1 + runif(n)^(-1/shape))^inequality
}


dparetoIII = function(x, location = 0, scale = 1, inequality = 1, log = FALSE)
    dparetoIV(x = x, location=location, scale = scale, inequality=inequality,
              shape = 1, log = log)

pparetoIII = function(q, location = 0, scale = 1, inequality=1)
    pparetoIV(q=q, location=location, scale = scale, inequality=inequality,
              shape=1)

qparetoIII = function(p, location = 0, scale = 1, inequality=1)
    qparetoIV(p=p, location=location, scale = scale, inequality=inequality,
              shape=1)

rparetoIII = function(n, location = 0, scale = 1, inequality=1)
    rparetoIV(n = n, location=location, scale = scale, inequality=inequality,
              shape=1)



dparetoII = function(x, location = 0, scale = 1, shape = 1, log = FALSE)
    dparetoIV(x = x, location=location, scale = scale,
              inequality = 1, shape = shape,
              log = log)

pparetoII = function(q, location = 0, scale = 1, shape=1)
    pparetoIV(q=q, location=location, scale = scale,
              inequality = 1, shape = shape)

qparetoII = function(p, location = 0, scale = 1, shape=1)
    qparetoIV(p=p, location=location, scale = scale,
              inequality = 1, shape = shape)

rparetoII = function(n, location = 0, scale = 1, shape=1)
    rparetoIV(n = n, location=location, scale = scale,
              inequality = 1, shape = shape)


dparetoI = function(x, scale = 1, shape=1)
    dparetoIV(x = x, location=scale, scale = scale, inequality = 1, shape = shape)

pparetoI = function(q, scale = 1, shape=1)
    pparetoIV(q=q, location=scale, scale = scale, inequality = 1, shape = shape)

qparetoI = function(p, scale = 1, shape=1)
    qparetoIV(p=p, location=scale, scale = scale, inequality = 1, shape = shape)

rparetoI = function(n, scale = 1, shape=1)
    rparetoIV(n = n, location=scale, scale = scale, inequality = 1, shape = shape)



 paretoIV = function(location = 0,
                    lscale = "loge",
                    linequality = "loge",
                    lshape = "loge",
                    escale = list(), einequality = list(), eshape = list(),
                    iscale = 1, iinequality = 1, ishape = NULL,
                    imethod = 1) {
    if (mode(lscale) != "character" && mode(lscale) != "name")
        lscale = as.character(substitute(lscale))
    if (mode(linequality) != "character" && mode(linequality) != "name")
        linequality = as.character(substitute(linequality))
    if (mode(lshape) != "character" && mode(lshape) != "name")
        lshape = as.character(substitute(lshape))
    if (!is.Numeric(location))
        stop("argument 'location' must be numeric")
    if (is.Numeric(iscale) && any(iscale <= 0))
        stop("argument 'iscale' must be positive")
    if (is.Numeric(iinequality) && any(iinequality <= 0))
        stop("argument 'iinequality' must be positive")
    if (is.Numeric(ishape) && any(ishape <= 0))
        stop("argument 'ishape' must be positive")
    if (!is.Numeric(imethod, allow = 1, integ = TRUE) || imethod>2)
        stop("bad input for argument 'imethod'")
    if (linequality == "nloge" && location != 0)
        warning("The Burr distribution has 'location = 0' and 'linequality=nloge'")
    if (!is.list(escale)) escale = list()
    if (!is.list(einequality)) einequality = list()
    if (!is.list(eshape)) eshape = list()

    new("vglmff",
    blurb = c("Pareto(IV) distribution F(y)=1-[1+((y - ", location,
            ")/scale)^(1/inequality)]^(-shape),",
            "\n", "         y > ",
            location, ", scale > 0, inequality > 0, shape > 0,\n",
            "Links:    ", namesof("scale", lscale, earg = escale ), ", ",
                          namesof("inequality", linequality, earg = einequality ), ", ",
                          namesof("shape", lshape, earg = eshape ), "\n",
            "Mean:    location + scale * NA"),
    initialize = eval(substitute(expression({
        if (ncol(cbind(y)) != 1)
            stop("response must be a vector or a one-column matrix")
        predictors.names = 
        c(namesof("scale", .lscale, earg = .escale, tag = FALSE),
          namesof("inequality", .linequality, earg = .einequality, tag = FALSE),
          namesof("shape", .lshape, earg = .eshape, tag = FALSE))
        extra$location = location = .location
        if (any(y <= location))
        stop("the response must have values > than the 'location' argument")
        if (!length(etastart)) {
            inequality.init = if (length( .iinequality)) .iinequality else  1
            scale.init = if (length( .iscale)) .iscale else 1
            shape.init = if (length( .ishape)) .ishape else NULL
            if (!length(shape.init)) {
                zedd = (y - location) / scale.init
                if ( .imethod == 1) {
                    A1 = weighted.mean(1/(1 + zedd^(1/inequality.init)), w = w)
                    A2 = weighted.mean(1/(1 + zedd^(1/inequality.init))^2, w = w)
                } else {
                    A1 = median(1/(1 + zedd^(1/inequality.init )))
                    A2 = median(1/(1 + zedd^(1/inequality.init))^2)
                }
                shape.init = max(0.01, (2*A2-A1)/(A1-A2))
            }
            etastart=cbind(
              theta2eta(rep(scale.init, len = n), .lscale, earg = .escale),
              theta2eta(rep(inequality.init, len = n), .linequality, earg = .einequality),
              theta2eta(rep(shape.init, len = n), .lshape, earg = .eshape))
        }
    }), list( .location = location, .lscale = lscale,
             .linequality = linequality, .lshape = lshape, .imethod = imethod,
             .escale = escale, .einequality = einequality, .eshape = eshape,
             .iscale = iscale, .iinequality=iinequality, .ishape = ishape ))),
    inverse = eval(substitute(function(eta, extra = NULL) {
        location = extra$location
        Scale = eta2theta(eta[,1], .lscale, earg = .escale)
        inequality = eta2theta(eta[,2], .linequality, earg = .einequality)
        shape = eta2theta(eta[,3], .lshape, earg = .eshape)
        location + Scale * NA
    }, list( .lscale = lscale, .linequality = linequality, .lshape = lshape,
             .escale = escale, .einequality = einequality, .eshape = eshape ))),
    last = eval(substitute(expression({
        misc$link = c("scale" = .lscale, "inequality" = .linequality,
                    "shape" = .lshape)
        misc$earg = list(scale = .escale, inequality= .einequality,
                         shape = .eshape)
        misc$location = extra$location # Use this for prediction
    }), list( .lscale = lscale,  .linequality = linequality, .lshape = lshape,
              .escale = escale, .einequality = einequality, .eshape = eshape ))),
    loglikelihood = eval(substitute(
        function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
        location = extra$location
        Scale = eta2theta(eta[,1], .lscale, earg = .escale)
        inequality = eta2theta(eta[,2], .linequality, earg = .einequality)
        shape = eta2theta(eta[,3], .lshape, earg = .eshape)
        zedd = (y - location) / Scale
        if (residuals) stop("loglikelihood residuals ",
                            "not implemented yet") else {
            sum(w * dparetoIV(x = y, location=location, scale=Scale,
                              inequality=inequality, shape = shape, log = TRUE))
        }
    }, list( .lscale = lscale,  .linequality = linequality, .lshape = lshape,
             .escale = escale, .einequality = einequality, .eshape = eshape ))),
    vfamily = c("paretoIV"),
    deriv = eval(substitute(expression({
        location = extra$location
        Scale = eta2theta(eta[,1], .lscale, earg = .escale)
        inequality = eta2theta(eta[,2], .linequality, earg = .einequality)
        shape = eta2theta(eta[,3], .lshape, earg = .eshape)
        zedd = (y - location) / Scale
        temp100 = 1 + zedd^(1/inequality)
        dl.dscale = (shape  - (1+shape) / temp100) / (inequality * Scale)
        dl.dinequality = ((log(zedd) * (shape - (1+shape)/temp100)) /
                         inequality - 1) / inequality
        dl.dshape = -log(temp100) + 1/shape
        dscale.deta = dtheta.deta(Scale, .lscale, earg = .escale)
        dinequality.deta = dtheta.deta(inequality, .linequality, earg = .einequality)
        dshape.deta = dtheta.deta(shape, .lshape, earg = .eshape)
        c(w) * cbind(dl.dscale * dscale.deta,
                     dl.dinequality * dinequality.deta, 
                     dl.dshape * dshape.deta)
    }), list( .lscale = lscale,  .linequality = linequality, .lshape = lshape,
              .escale = escale, .einequality = einequality, .eshape = eshape ))),
    weight = eval(substitute(expression({
        temp200 = digamma(shape) - digamma(1) - 1
        d2scale.deta2 = shape / ((inequality*Scale)^2 * (shape+2))
        d2inequality.deta2 = (shape * (temp200^2 + trigamma(shape) + trigamma(1)
                             ) + 2*(temp200+1)) / (inequality^2 * (shape+2))
        d2shape.deta2 = 1 / shape^2
        d2si.deta2 = (shape*(-temp200) -1) / (inequality^2 * Scale * (shape+2))
        d2ss.deta2 = -1 / ((inequality*Scale) * (shape+1))
        d2is.deta2 = temp200 / (inequality*(shape+1))
        wz = matrix(0, n, dimm(M))
        wz[,iam(1,1,M)] = dscale.deta^2 * d2scale.deta2
        wz[,iam(2,2,M)] = dinequality.deta^2 * d2inequality.deta2
        wz[,iam(3,3,M)] = dshape.deta^2 * d2shape.deta2
        wz[,iam(1,2,M)] = dscale.deta * dinequality.deta * d2si.deta2
        wz[,iam(1,3,M)] = dscale.deta * dshape.deta * d2ss.deta2
        wz[,iam(2,3,M)] = dinequality.deta * dshape.deta * d2is.deta2
        c(w) * wz
    }), list( .lscale = lscale,  .linequality = linequality, .lshape = lshape,
              .escale = escale, .einequality = einequality, .eshape = eshape ))))
}




 paretoIII = function(location = 0,
                      lscale = "loge",
                      linequality = "loge",
                      escale = list(), einequality = list(),
                      iscale = NULL, iinequality = NULL) {
    if (mode(lscale) != "character" && mode(lscale) != "name")
        lscale = as.character(substitute(lscale))
    if (mode(linequality) != "character" && mode(linequality) != "name")
        linequality = as.character(substitute(linequality))
    if (!is.Numeric(location))
        stop("argument 'location' must be numeric")
    if (is.Numeric(iscale) && any(iscale <= 0))
        stop("argument 'iscale' must be positive")
    if (is.Numeric(iinequality) && any(iinequality <= 0))
        stop("argument 'iinequality' must be positive")
    if (!is.list(escale)) escale = list()
    if (!is.list(einequality)) einequality = list()

    new("vglmff",
    blurb = c("Pareto(III) distribution F(y)=1-[1+((y - ", location,
            ")/scale)^(1/inequality)]^(-1),",
            "\n", "         y > ",
            location, ", scale > 0, inequality > 0, \n",
            "Links:    ",
            namesof("scale", lscale, earg = escale ), ", ",
            namesof("inequality", linequality, earg = einequality ), "\n",
            "Mean:    location + scale * NA"),
    initialize = eval(substitute(expression({
        if (ncol(cbind(y)) != 1)
            stop("the response must be a vector or a one-column matrix")
        predictors.names = 
        c(namesof("scale", .lscale, earg = .escale, tag = FALSE),
          namesof("inequality", .linequality, earg = .einequality, tag = FALSE))
        extra$location = location = .location
        if (any(y <= location))
        stop("the response must have values > than the 'location' argument")
        if (!length(etastart)) {
            inequality.init = if (length( .iinequality)) .iinequality else  NULL
            scale.init = if (length( .iscale)) .iscale else NULL
            if (!length(inequality.init) || !length(scale.init)) {
                probs = (1:4)/5
                ytemp = quantile(x=log(y-location), probs=probs)
                fittemp = lsfit(x=logit(probs), y = ytemp, int = TRUE)
                if (!length(inequality.init))
                    inequality.init = max(fittemp$coef["X"], 0.01)
                if (!length(scale.init))
                    scale.init = exp(fittemp$coef["Intercept"])
            }
            etastart=cbind(
            theta2eta(rep(scale.init, len = n), .lscale, earg = .escale),
            theta2eta(rep(inequality.init, len = n), .linequality,
                      earg = .einequality))
        }
    }), list( .location = location, .lscale = lscale, .linequality = linequality,
              .escale = escale, .einequality = einequality,
              .iscale = iscale, .iinequality=iinequality ))),
    inverse = eval(substitute(function(eta, extra = NULL) {
        location = extra$location
        Scale = eta2theta(eta[,1], .lscale, earg = .escale)
        inequality = eta2theta(eta[,2], .linequality, earg = .einequality)
        location + Scale * NA
    }, list( .lscale = lscale, .linequality = linequality,
             .escale = escale, .einequality = einequality ))),
    last = eval(substitute(expression({
        misc$link = c("scale" = .lscale, "inequality" = .linequality)
        misc$earg = list(scale = .escale, inequality= .einequality)
        misc$location = extra$location # Use this for prediction
    }), list( .lscale = lscale, .linequality = linequality,
              .escale = escale, .einequality = einequality ))),
    loglikelihood = eval(substitute(
        function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
        location = extra$location
        Scale = eta2theta(eta[,1], .lscale, earg = .escale)
        inequality = eta2theta(eta[,2], .linequality, earg = .einequality)
        zedd = (y - location) / Scale
        if (residuals) stop("loglikelihood residuals ",
                            "not implemented yet") else {
            sum(w * dparetoIII(x = y, location=location, scale=Scale,
                               inequality=inequality, log = TRUE))
        }
    }, list( .lscale = lscale, .linequality = linequality,
             .escale = escale, .einequality = einequality ))),
    vfamily = c("paretoIII"),
    deriv = eval(substitute(expression({
        location = extra$location
        Scale = eta2theta(eta[,1], .lscale, earg = .escale)
        inequality = eta2theta(eta[,2], .linequality, earg = .einequality)
        shape = 1
        zedd = (y - location) / Scale
        temp100 = 1 + zedd^(1/inequality)
        dl.dscale = (shape  - (1+shape) / temp100) / (inequality * Scale)
        dl.dinequality = ((log(zedd) * (shape - (1+shape)/temp100)) /
                         inequality - 1) / inequality
        dscale.deta = dtheta.deta(Scale, .lscale, earg = .escale)
        dinequality.deta = dtheta.deta(inequality, .linequality, earg = .einequality)
        c(w) * cbind(dl.dscale * dscale.deta,
                     dl.dinequality * dinequality.deta)
    }), list( .lscale = lscale, .linequality = linequality,
              .escale = escale, .einequality = einequality ))),
    weight = eval(substitute(expression({
        d2scale.deta2 = 1 / ((inequality*Scale)^2 * 3)
        d2inequality.deta2 = (1 + 2* trigamma(1)) / (inequality^2 * 3)
        wz = matrix(0, n, M) # It is diagonal
        wz[,iam(1,1,M)] = dscale.deta^2 * d2scale.deta2
        wz[,iam(2,2,M)] = dinequality.deta^2 * d2inequality.deta2
        c(w) * wz
    }), list( .lscale = lscale,  .linequality = linequality,
              .escale = escale, .einequality = einequality ))))
}





 paretoII = function(location = 0,
                     lscale = "loge",
                     lshape = "loge",
                     escale = list(), eshape = list(),
                     iscale = NULL, ishape = NULL) {
    if (mode(lscale) != "character" && mode(lscale) != "name")
        lscale = as.character(substitute(lscale))
    if (mode(lshape) != "character" && mode(lshape) != "name")
        lshape = as.character(substitute(lshape))
    if (!is.Numeric(location))
        stop("argument 'location' must be numeric")
    if (is.Numeric(iscale) && any(iscale <= 0))
        stop("argument 'iscale' must be positive")
    if (is.Numeric(ishape) && any(ishape <= 0))
        stop("argument 'ishape' must be positive")
    if (!is.list(escale)) escale = list()
    if (!is.list(eshape)) eshape = list()

    new("vglmff",
    blurb = c("Pareto(II) distribution F(y)=1-[1+(y - ", location,
            ")/scale]^(-shape),",
            "\n", "         y > ",
            location, ", scale > 0,  shape > 0,\n",
            "Links:    ", namesof("scale", lscale, earg = escale ), ", ",
                          namesof("shape", lshape, earg = eshape ), "\n",
            "Mean:    location + scale * NA"),
    initialize = eval(substitute(expression({
        if (ncol(cbind(y)) != 1)
            stop("the response must be a vector or a one-column matrix")
        predictors.names = 
        c(namesof("scale", .lscale, earg = .escale, tag = FALSE),
          namesof("shape", .lshape, earg = .eshape, tag = FALSE))
        extra$location = location = .location
        if (any(y <= location))
        stop("the response must have values > than the 'location' argument")
        if (!length(etastart)) {
            scale.init = if (length( .iscale)) .iscale else NULL
            shape.init = if (length( .ishape)) .ishape else  NULL
            if (!length(shape.init) || !length(scale.init)) {
                probs = (1:4)/5
                scale.init.0 = 1
                ytemp = quantile(x=log(y-location+scale.init.0), probs=probs)
                fittemp = lsfit(x=log1p(-probs), y = ytemp, int = TRUE)
                if (!length(shape.init))
                    shape.init = max(-1/fittemp$coef["X"], 0.01)
                if (!length(scale.init))
                    scale.init = exp(fittemp$coef["Intercept"])
            }
            etastart=cbind(
            theta2eta(rep(scale.init, len = n), .lscale, earg = .escale),
            theta2eta(rep(shape.init, len = n), .lshape, earg = .eshape))
        }
    }), list( .location = location, .lscale = lscale,
              .escale = escale, .eshape = eshape, 
              .lshape = lshape, .iscale = iscale, .ishape = ishape ))),
    inverse = eval(substitute(function(eta, extra = NULL) {
        location = extra$location
        Scale = eta2theta(eta[,1], .lscale, earg = .escale)
        shape = eta2theta(eta[,2], .lshape, earg = .eshape)
        location + Scale * NA
    }, list( .lscale = lscale, .lshape = lshape,
             .escale = escale, .eshape = eshape ))),
    last = eval(substitute(expression({
        misc$link =    c("scale" = .lscale, "shape" = .lshape)
        misc$earg = list("scale" = .escale, "shape" = .eshape)
        misc$location = extra$location # Use this for prediction
    }), list( .lscale = lscale, .lshape = lshape,
              .escale = escale, .eshape = eshape ))),
    loglikelihood = eval(substitute(
        function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
        location = extra$location
        Scale = eta2theta(eta[,1], .lscale, earg = .escale)
        shape = eta2theta(eta[,2], .lshape, earg = .eshape)
        zedd = (y - location) / Scale
        if (residuals) stop("loglikelihood residuals ",
                            "not implemented yet") else {
            sum(w * dparetoII(x = y, location=location, scale=Scale,
                              shape = shape, log = TRUE))
        }
    }, list( .lscale = lscale, .lshape = lshape,
             .escale = escale, .eshape = eshape ))),
    vfamily = c("paretoII"),
    deriv = eval(substitute(expression({
        location = extra$location
        Scale = eta2theta(eta[,1], .lscale, earg = .escale)
        shape = eta2theta(eta[,2], .lshape, earg = .eshape)
        zedd = (y - location) / Scale
        temp100 = 1 + zedd
        dl.dscale = (shape  - (1+shape) / temp100) / (1 * Scale)
        dl.dshape = -log(temp100) + 1/shape
        dscale.deta = dtheta.deta(Scale, .lscale, earg = .escale)
        dshape.deta = dtheta.deta(shape, .lshape, earg = .eshape)
        c(w) * cbind(dl.dscale * dscale.deta,
                     dl.dshape * dshape.deta)
    }), list( .lscale = lscale, .lshape = lshape,
              .escale = escale, .eshape = eshape ))),
    weight = eval(substitute(expression({
        d2scale.deta2 = shape / (Scale^2 * (shape+2))
        d2shape.deta2 = 1 / shape^2
        d2ss.deta2 = -1 / (Scale * (shape+1))
        wz = matrix(0, n, dimm(M))
        wz[,iam(1,1,M)] = dscale.deta^2 * d2scale.deta2
        wz[,iam(2,2,M)] = dshape.deta^2 * d2shape.deta2
        wz[,iam(1,2,M)] = dscale.deta * dshape.deta * d2ss.deta2
        c(w) * wz
    }), list( .lscale = lscale,  .lshape = lshape,
              .escale = escale, .eshape = eshape ))))
}





dpareto = function(x, location, shape, log = FALSE) {
    if (!is.logical(log.arg <- log))
        stop("bad input for argument 'log'")
    rm(log)

    L = max(length(x), length(location), length(shape)) 
    x = rep(x, len = L); location = rep(location, len = L); shape= rep(shape, len = L)

    logdensity = rep(log(0), len = L)
    xok = (x > location)
    logdensity[xok] = log(shape[xok]) + shape[xok] * log(location[xok]) -
                      (shape[xok]+1) * log(x[xok])
    if (log.arg) logdensity else exp(logdensity)
}

ppareto = function(q, location, shape) {
    if (any(location <= 0)) stop("argument 'location' must be positive")
    if (any(shape <= 0)) stop("argument 'shape' must be positive")
    L = max(length(q), length(location), length(shape))
    q = rep(q, len = L); location = rep(location, len = L); shape= rep(shape, len = L)
    ifelse(q > location, 1 - (location/q)^shape, 0)
}

qpareto = function(p, location, shape) {
    if (any(location <= 0)) stop("argument 'location' must be positive")
    if (any(shape <= 0)) stop("argument 'shape' must be positive")
   if (any(p <= 0) || any(p >= 1)) stop("argument 'p' must be between 0 and 1")
    location / (1 - p)^(1/shape)
}

rpareto = function(n, location, shape) {
    if (!is.Numeric(n, posit = TRUE, integ = TRUE, allow = 1)) 
        stop("bad input for argument 'n'")
    if (any(location <= 0)) stop("argument 'location' must be positive")
    if (any(shape <= 0)) stop("argument 'shape' must be positive")
    location / runif(n)^(1/shape)
}



 pareto1 = function(lshape = "loge", earg = list(), location = NULL) {
    if (mode(lshape) != "character" && mode(lshape) != "name")
        lshape = as.character(substitute(lshape))
    if (is.Numeric(location) && location <= 0)
        stop("argument 'location' must be positive")
    if (!is.list(earg)) earg = list()

    new("vglmff",
    blurb = c("Pareto distribution f(y) = shape * location^shape / y^(shape+1),",
            " 0<location<y, shape>0\n",
            "Link:    ", namesof("shape", lshape, earg = earg), "\n", "\n",
            "Mean:    location*shape/(shape-1) for shape>1"),
    initialize = eval(substitute(expression({
        if (ncol(cbind(y)) != 1)
            stop("response must be a vector or a one-column matrix")
        predictors.names = namesof("shape", .lshape, earg = .earg, tag = FALSE) 
        locationhat = if (!length( .location)) {
            locationEstimated = TRUE
            min(y) # - .smallno
        } else {
            locationEstimated = FALSE
            .location
        }
        if (any(y < locationhat))
            stop("the value of location is too high (requires 0 < location < min(y))")
        extra$location = locationhat
        extra$locationEstimated = locationEstimated
        if (!length(etastart)) {
            k.init = (y + 1/8) / (y - locationhat + 1/8)
            etastart = theta2eta(k.init, .lshape, earg = .earg)
        }
    }), list( .lshape = lshape, .earg = earg,
              .location = location ))),
    inverse = eval(substitute(function(eta, extra = NULL) {
        k = eta2theta(eta, .lshape, earg = .earg)
        location = extra$location
        ifelse(k > 1, k * location / (k-1), NA)
    }, list( .lshape = lshape, .earg = earg ))),
    last = eval(substitute(expression({
        misc$link =    c(k = .lshape)
        misc$earg = list(k = .earg)
        misc$location = extra$location # Use this for prediction
    }), list( .lshape = lshape, .earg = earg ))),
    loglikelihood = eval(substitute(
        function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
        k = eta2theta(eta, .lshape, earg = .earg)
        location = extra$location
        if (residuals) stop("loglikelihood residuals ",
                            "not implemented yet") else {


            sum(w * (log(k) + k * log(location) - (k+1) * log(y )))
        }
    }, list( .lshape = lshape, .earg = earg ))),
    vfamily = c("pareto1"),
    deriv = eval(substitute(expression({
        location = extra$location
        k = eta2theta(eta, .lshape, earg = .earg)
        dl.dk = 1/k + log(location/y)
        dk.deta = dtheta.deta(k, .lshape, earg = .earg)
        c(w) * dl.dk * dk.deta
    }), list( .lshape = lshape, .earg = earg ))),
    weight = eval(substitute(expression({
        ed2l.dk2 = 1 / k^2
        wz = c(w) * dk.deta^2 * ed2l.dk2
        wz
    }), list( .lshape = lshape, .earg = earg ))))
}





dtpareto = function(x, lower, upper, shape, log = FALSE) {

  if (!is.logical(log.arg <- log))
    stop("bad input for argument 'log'")
  rm(log)

  if (!is.Numeric(x))
    stop("bad input for argument 'x'")
  if (!is.Numeric(lower, pos = TRUE))
    stop("argument 'lower' must be positive")
  if (!is.Numeric(upper, pos = TRUE))
    stop("argument 'upper' must be positive")
  if (!is.Numeric(shape, pos = TRUE))
    stop("argument 'shape' must be positive")

  L = max(length(x), length(lower), length(upper), length(shape))
  x = rep(x, len = L); shape = rep(shape, len = L)
  lower = rep(lower, len = L); upper = rep(upper, len = L);


  logdensity <- rep(log(0), len = L)
  xok <- (0 < lower) & (lower < x) & (x < upper) & (shape > 0)

  logdensity[xok] <- log(shape[xok]) + shape[xok] * log(lower[xok]) -
                     (shape[xok] + 1) * log(x[xok]) -
                     log1p(-(lower[xok] / upper[xok])^(shape[xok]))

  logdensity[shape <= 0] <- NaN
  logdensity[upper < lower] <- NaN
  logdensity[0 > lower] <- NaN
  if (log.arg) logdensity else exp(logdensity)
}


ptpareto = function(q, lower, upper, shape) {
    if (!is.Numeric(q))
      stop("bad input for argument 'q'")
    if (!is.Numeric(lower, pos = TRUE))
      stop("argument 'lower' must be positive")
    if (!is.Numeric(upper, pos = TRUE))
      stop("argument 'upper' must be positive")
    if (!is.Numeric(shape, pos = TRUE))
      stop("argument 'shape' must be positive")

    L = max(length(q), length(lower), length(upper), length(shape)) 
    q = rep(q, len = L); lower = rep(lower, len = L);
    upper = rep(upper, len = L); shape= rep(shape, len = L)

    ans = q * 0
    xok <- (0 < lower) & (lower < q) & (q < upper) & (shape > 0)
    ans[xok] = (1 - (lower[xok]/q[xok])^shape[xok]) / (1 -
                    (lower[xok]/upper[xok])^shape[xok])
    ans[q >= upper] = 1
    ans[upper < lower] <- NaN
    ans
}


qtpareto = function(p, lower, upper, shape) {
    if (!is.Numeric(p, posit = TRUE))
      stop("bad input for argument 'p'")
    if (!is.Numeric(lower, pos = TRUE))
      stop("argument 'lower' must be positive")
    if (!is.Numeric(upper, pos = TRUE))
      stop("argument 'upper' must be positive")
    if (!is.Numeric(shape, pos = TRUE))
      stop("argument 'shape' must be positive")
    if (max(p) >= 1)
      stop("argument 'p' must be in (0, 1)")
    if (min(upper - lower, na.rm = TRUE) < 0)
      stop("argument 'upper' must be greater than 'lower' values")

    lower / (1 - p*(1-(lower/upper)^shape))^(1/shape)
}


rtpareto = function(n, lower, upper, shape) {
    if (!is.Numeric(lower, pos = TRUE))
      stop("argument 'lower' must be positive")
    if (!is.Numeric(upper, pos = TRUE))
      stop("argument 'upper' must be positive")
    if (!is.Numeric(shape, pos = TRUE))
      stop("argument 'shape' must be positive")

    qtpareto(p = runif(n), lower = lower, upper = upper, shape = shape)
}




 tpareto1 = function(lower, upper, lshape = "loge", earg = list(),
                     ishape = NULL, imethod = 1) {
    if (mode(lshape) != "character" && mode(lshape) != "name")
        lshape = as.character(substitute(lshape))

    if (!is.Numeric(lower, posit = TRUE, allow = 1))
        stop("bad input for argument 'lower'")
    if (!is.Numeric(upper, posit = TRUE, allow = 1))
        stop("bad input for argument 'upper'")
    if (lower >= upper)
        stop("lower < upper is required")

    if (length(ishape) && !is.Numeric(ishape, posit = TRUE))
        stop("bad input for argument 'ishape'")
    if (!is.list(earg)) earg = list()
    if (!is.Numeric(imethod, allow = 1, integ = TRUE, posit = TRUE) ||
       imethod > 2)
        stop("argument 'imethod' must be 1 or 2")

    new("vglmff",
    blurb = c("Truncated Pareto distribution f(y) = shape * lower^shape /",
            "(y^(shape+1) * (1-(lower/upper)^shape)),",
            " 0 < lower < y < upper < Inf, shape>0\n",
            "Link:    ", namesof("shape", lshape, earg = earg), "\n", "\n",
            "Mean:    shape*lower^shape*(upper^(1-shape)-lower^(1-shape)) /",
                      " ((1-shape) * (1-(lower/upper)^shape))"),
    initialize = eval(substitute(expression({
        if (ncol(cbind(y)) != 1)
            stop("response must be a vector or a one-column matrix")

        predictors.names = namesof("shape", .lshape, earg = .earg,
                                   tag = FALSE)
        if (any(y <= .lower))
            stop("the value of argument 'lower' is too high ",
                 "(requires '0 < lower < min(y)')")

        extra$lower = .lower
        if (any(y >= .upper))
            stop("the value of argument 'upper' is too low ",
                 "(requires 'max(y) < upper')")
        extra$upper = .upper

        if (!length(etastart)) {
            shape.init = if (is.Numeric( .ishape)) 0 * y + .ishape else
            if ( .imethod == 2) {
                0 * y + median(rep((y + 1/8) / (y - .lower + 1/8), times=w))
            } else {
                tpareto1.Loglikfun = function(shape, y, x, w, extraargs) {
                     myratio = .lower / .upper
                     sum(w * (log(shape) + shape * log( .lower) -
                              (shape+1) * log(y) - log1p(-myratio^shape)))
                 }
                 shape.grid = 2^((-4):4)
                 try.this = getMaxMin(shape.grid, objfun = tpareto1.Loglikfun,
                                      y = y,  x = x, w = w)
                try.this = rep(try.this, len = n)
                try.this
            }
            etastart = theta2eta(shape.init, .lshape, earg = .earg)
        }
    }), list( .lshape = lshape, .earg = earg,
              .ishape = ishape,
              .imethod = imethod,
              .lower = lower, .upper = upper ))),
    inverse = eval(substitute(function(eta, extra = NULL) {
        shape = eta2theta(eta, .lshape, earg = .earg)
        myratio = .lower / .upper
        constprop = shape * .lower^shape / (1 - myratio^shape)
        constprop * ( .upper^(1-shape) - .lower^(1-shape)) / (1-shape)
    }, list( .lshape = lshape, .earg = earg,
             .lower = lower, .upper = upper ))),
    last = eval(substitute(expression({
        misc$link =    c(shape = .lshape)
        misc$earg = list(shape = .earg)
        misc$lower = extra$lower
        misc$upper = extra$upper
        misc$expected = TRUE
    }), list( .lshape = lshape, .earg = earg,
              .lower = lower, .upper = upper ))),
    loglikelihood = eval(substitute(
        function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
        shape = eta2theta(eta, .lshape, earg = .earg)
        if (residuals) stop("loglikelihood residuals ",
                            "not implemented yet") else {
          ans = sum(w * dtpareto(x = y, lower = .lower , upper = .upper ,
                                 shape = shape, log = TRUE))
          ans
        }
    }, list( .lshape = lshape, .earg = earg,
             .lower = lower, .upper = upper ))),
    vfamily = c("tpareto1"),
    deriv = eval(substitute(expression({
        shape = eta2theta(eta, .lshape, earg = .earg)
        myratio = .lower / .upper
        myratio2 =  myratio^shape
        tmp330 = myratio2 * log(myratio) / (1 - myratio2)
        dl.dshape = 1 / shape + log( .lower) - log(y) + tmp330 
        dshape.deta = dtheta.deta(shape, .lshape, earg = .earg)
        c(w) * dl.dshape * dshape.deta
    }), list( .lshape = lshape, .earg = earg,
              .lower = lower, .upper = upper ))),
    weight = eval(substitute(expression({
        ed2l.dshape2 = 1 / shape^2 - tmp330^2 / myratio2
        wz = c(w) * dshape.deta^2 * ed2l.dshape2
        wz
    }), list( .lshape = lshape, .earg = earg,
              .lower = lower, .upper = upper ))))
}




erf = function(x)
    2 * pnorm(x * sqrt(2)) - 1

erfc = function(x)
    2 * pnorm(x * sqrt(2), lower = FALSE)



 wald <- function(link.lambda = "loge", earg = list(), init.lambda = NULL)
{
    if (mode(link.lambda) != "character" && mode(link.lambda) != "name")
        link.lambda = as.character(substitute(link.lambda))
    if (!is.list(earg)) earg = list()

    new("vglmff",
    blurb = c("Standard Wald distribution\n\n",
           "f(y) = sqrt(lambda/(2*pi*y^3)) * exp(-lambda*(y-1)^2/(2*y)), y&lambda>0",
           "\n", 
           "Link:     ", 
                         namesof("lambda", link.lambda, earg = earg), "\n",
           "Mean:     ", "1\n",
           "Variance: 1 / lambda"),
    initialize = eval(substitute(expression({
        if (ncol(cbind(y)) != 1)
            stop("response must be a vector or a one-column matrix")
        if (any(y <= 0)) stop("Require the response to have positive values")
        predictors.names = 
        namesof("lambda", .link.lambda, earg = .earg, short = TRUE)
        if (!length(etastart)) {
            initlambda = if (length( .init.lambda)) .init.lambda else
                         1 / (0.01 + (y-1)^2)
            initlambda = rep(initlambda, len = n)
            etastart = cbind(theta2eta(initlambda, link=.link.lambda, earg = .earg))
        }
    }), list( .link.lambda = link.lambda, .earg = earg,
             .init.lambda=init.lambda ))),
    inverse=function(eta, extra = NULL) {
        0*eta + 1
    },
    last = eval(substitute(expression({
        misc$link = c(lambda = .link.lambda )
        misc$earg = list(lambda = .earg )
    }), list( .link.lambda = link.lambda, .earg = earg ))),
    loglikelihood = eval(substitute(
             function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
        lambda = eta2theta(eta, link=.link.lambda, earg = .earg)
        if (residuals) stop("loglikelihood residuals ",
                            "not implemented yet") else
        sum(w * (0.5 * log(lambda/(2*pi*y^3)) - lambda * (y-1)^2 / (2*y)))
    }, list( .link.lambda = link.lambda, .earg = earg ))),
    vfamily = "wald",
    deriv = eval(substitute(expression({
        lambda = eta2theta(eta, link=.link.lambda, earg = .earg)
        dl.dlambda = 0.5 / lambda + 1 - 0.5 * (y + 1/y)
        dlambda.deta = dtheta.deta(theta=lambda, link=.link.lambda, earg = .earg)
        c(w) * cbind(dl.dlambda * dlambda.deta)
    }), list( .link.lambda = link.lambda, .earg = earg ))),
    weight = eval(substitute(expression({
        d2l.dlambda2 = 0.5 / (lambda^2)
        c(w) * cbind(dlambda.deta^2 * d2l.dlambda2)
    }), list( .link.lambda = link.lambda, .earg = earg ))))
}


 expexp = function(lshape = "loge", lscale = "loge",
                  eshape = list(), escale = list(),
                  ishape=1.1, iscale = NULL,  # ishape cannot be 1
                  tolerance = 1.0e-6,
                  zero = NULL) {

    if (mode(lshape) != "character" && mode(lshape) != "name")
        lshape = as.character(substitute(lshape))
    if (mode(lscale) != "character" && mode(lscale) != "name")
        lscale = as.character(substitute(lscale))
    if (length(zero) && !is.Numeric(zero, integer = TRUE, posit = TRUE))
        stop("bad input for argument 'zero'")
    if (!is.Numeric(tolerance, posit = TRUE, allow = 1) || tolerance>1.0e-2)
        stop("bad input for argument 'tolerance'")
    if (!is.Numeric(ishape, posit = TRUE))
        stop("bad input for argument 'ishape'")
    if (length(iscale) && !is.Numeric(iscale, posit = TRUE))
        stop("bad input for argument 'iscale'")
    ishape[ishape == 1] = 1.1   # Fails in @deriv
    if (!is.list(escale)) escale = list()
    if (!is.list(eshape)) eshape = list()

    new("vglmff",
    blurb = c("Exponentiated Exponential Distribution\n",
           "Links:    ",
           namesof("shape", lshape, earg = eshape), ", ",
           namesof("scale", lscale, earg = escale),"\n",
           "Mean:     (digamma(shape+1)-digamma(1))/scale"),
    constraints = eval(substitute(expression({
        constraints = cm.zero.vgam(constraints, x, .zero, M)
    }), list( .zero = zero ))),
    initialize = eval(substitute(expression({
        if (ncol(cbind(y)) != 1)
            stop("response must be a vector or a one-column matrix")
        predictors.names = 
        c(namesof("shape", .lshape, earg = .eshape, short = TRUE), 
          namesof("scale", .lscale, earg = .escale, short = TRUE))
        if (!length(etastart)) {
            shape.init = if (!is.Numeric( .ishape, posit = TRUE))
                   stop("argument 'ishape' must be positive") else
                   rep( .ishape, len = n)
            scale.init = if (length( .iscale)) rep( .iscale, len = n) else
                        (digamma(shape.init+1) - digamma(1)) / (y+1/8)
            scale.init = rep(weighted.mean(scale.init, w = w), len = n)
            etastart = cbind(theta2eta(shape.init, .lshape, earg = .eshape),
                             theta2eta(scale.init, .lscale, earg = .escale))
        }
    }), list( .lshape = lshape, .lscale = lscale, .iscale = iscale, .ishape = ishape,
              .eshape = eshape, .escale = escale ))),
    inverse = eval(substitute(function(eta, extra = NULL) {
        shape = eta2theta(eta[,1], .lshape, earg = .eshape)
        scale = eta2theta(eta[,2], .lscale, earg = .escale)
        (digamma(shape+1)-digamma(1)) / scale
    }, list( .lshape = lshape, .lscale = lscale,
             .eshape = eshape, .escale = escale ))),
    last = eval(substitute(expression({
        misc$link =    c("shape" = .lshape, "scale" = .lscale)
        misc$earg = list("shape" = .eshape, "scale" = .escale)
        misc$expected = TRUE
    }), list( .lshape = lshape, .lscale = lscale,
              .eshape = eshape, .escale = escale ))),
    loglikelihood= eval(substitute(
        function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
        shape = eta2theta(eta[,1], .lshape, earg = .eshape)
        scale = eta2theta(eta[,2], .lscale, earg = .escale)
        if (residuals) stop("loglikelihood residuals ",
                            "not implemented yet") else
        sum(w * (log(shape) + log(scale) + 
                 (shape-1)*log1p(-exp(-scale*y)) - scale*y))
    }, list( .lscale = lscale, .lshape = lshape,
             .eshape = eshape, .escale = escale ))),
    vfamily = c("expexp"),
    deriv = eval(substitute(expression({
        shape = eta2theta(eta[,1], .lshape, earg = .eshape)
        scale = eta2theta(eta[,2], .lscale, earg = .escale)
        dl.dscale = 1/scale + (shape-1)*y*exp(-scale*y) / (-expm1(-scale*y)) - y
        dl.dshape = 1/shape + log1p(-exp(-scale*y))
        dscale.deta = dtheta.deta(scale, .lscale, earg = .escale)
        dshape.deta = dtheta.deta(shape, .lshape, earg = .eshape)
        c(w) * cbind(dl.dshape * dshape.deta,
                     dl.dscale * dscale.deta)
    }), list( .lshape = lshape, .lscale = lscale,
              .eshape = eshape, .escale = escale ))),
    weight = eval(substitute(expression({
        d11 = 1 / shape^2  # True for all shape
        d22 = d12 = rep(as.numeric(NA), len = n)
        index2 = abs(shape - 2) > .tolerance  # index2 = shape != 1
        largeno = 10000
        if (any(index2)) {
            Shape = shape[index2]
            Shape[abs(Shape-1) < .tolerance] = 1.001 # digamma(0) is undefined
            Scale = scale[index2]
            tmp200 = trigamma(1)-trigamma(Shape-1) +
                  (digamma(Shape-1)-digamma(1))^2    # Fails when Shape == 1
            tmp300 = trigamma(1)-digamma(Shape)+(digamma(Shape)-digamma(1))^2
            d22[index2] = (1 + Shape*(Shape-1)*tmp200/(Shape-2)) / Scale^2 +
                          Shape*tmp300 / Scale^2
        }
        if (any(!index2)) {
            Scale = scale[!index2]
            d22[!index2] = (1 + 4 * sum(1/(2 + (0:largeno))^3)) / Scale^2
        }

        index1 = abs(shape - 1) > .tolerance  # index1 = shape != 1
        if (any(index1)) {
            Shape = shape[index1]
            Scale = scale[index1]
            d12[index1] = -(Shape*(digamma(Shape)-digamma(1))/(Shape-1) -
                          digamma(Shape+1) + digamma(1)) / Scale
        }
        if (any(!index1)) {
            Scale = scale[!index1]
            d12[!index1] = -sum(1/(2 + (0:largeno))^2) / Scale
        }
        wz = matrix(0, n, dimm(M))
        wz[,iam(1,1,M)] = dshape.deta^2 * d11
        wz[,iam(2,2,M)] = dscale.deta^2 * d22
        wz[,iam(1,2,M)] = dscale.deta * dshape.deta * d12
        c(w) * wz
    }), list( .tolerance=tolerance ))))
}





 expexp1 = function(lscale = "loge",
                   escale = list(),
                   iscale = NULL,
                   ishape=1) {
    if (mode(lscale) != "character" && mode(lscale) != "name")
        lscale = as.character(substitute(lscale))
    if (!is.list(escale)) escale = list()
    if (length(iscale) && !is.Numeric(iscale, posit = TRUE))
        stop("bad input for argument 'iscale'")

    new("vglmff",
    blurb = c("Exponentiated Exponential Distribution",
            " (profile likelihood estimation)\n",
           "Links:    ",
           namesof("scale", lscale, earg = escale), "\n",
           "Mean:     (digamma(shape+1)-digamma(1))/scale"),
    initialize = eval(substitute(expression({
        if (ncol(cbind(y)) != 1)
            stop("response must be a vector or a one-column matrix")
        predictors.names = namesof("scale", .lscale, earg = .escale, short = TRUE)
        if (length(w) != n || !is.Numeric(w, integer = TRUE, posit = TRUE))
            stop("weights must be a vector of positive integers")
        if (!intercept.only)
  stop("this family function only works for an intercept-only, i.e., y ~ 1")
        extra$yvector = y
        extra$sumw = sum(w)
        extra$w = w
        if (!length(etastart)) {
            shape.init = if (!is.Numeric( .ishape, posit = TRUE))
                   stop("argument 'ishape' must be positive") else
                   rep( .ishape, len = n)
            scaleinit = if (length( .iscale)) rep( .iscale, len = n) else
                        (digamma(shape.init+1) - digamma(1)) / (y+1/8)  
            etastart = cbind(theta2eta(scaleinit, .lscale, earg = .escale))
        }
    }), list( .lscale = lscale, .iscale = iscale, .ishape = ishape,
              .escale = escale ))),
    inverse = eval(substitute(function(eta, extra = NULL) {
        scale = eta2theta(eta, .lscale, earg = .escale)
        temp7 =  -expm1(-scale*extra$yvector)
        shape = -extra$sumw / sum(extra$w*log(temp7)) # \gamma(\theta)
        (digamma(shape+1)-digamma(1)) / scale
    }, list( .lscale = lscale,
             .escale = escale ))),
    last = eval(substitute(expression({
        misc$link =    c("scale" = .lscale)
        misc$earg = list("scale" = .escale)
        temp7 =  -expm1(-scale*y)
        shape = -extra$sumw / sum(w*log(temp7)) # \gamma(\theta)
        misc$shape = shape   # Store the ML estimate here
        misc$pooled.weight = pooled.weight
    }), list( .lscale = lscale, .escale = escale ))),
    loglikelihood= eval(substitute(
        function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
        scale = eta2theta(eta, .lscale, earg = .escale)
        temp7 =  -expm1(-scale*y)
        shape = -extra$sumw / sum(w*log(temp7)) # \gamma(\theta)
        if (residuals) stop("loglikelihood residuals ",
                            "not implemented yet") else
        sum(w * (log(shape) + log(scale) + 
                 (shape-1)*log1p(-exp(-scale*y)) - scale*y))
    }, list( .lscale = lscale, .escale = escale ))),
    vfamily = c("expexp1"),
    deriv = eval(substitute(expression({
        scale = eta2theta(eta, .lscale, earg = .escale)
        temp6 = exp(-scale*y)
        temp7 = 1-temp6
        shape = -extra$sumw / sum(w*log(temp7)) # \gamma(\theta)
        d1 = 1/scale + (shape-1)*y*temp6/temp7 - y
        c(w) * cbind(d1 * dtheta.deta(scale, .lscale, earg = .escale))
    }), list( .lscale = lscale, .escale = escale ))),
    weight = eval(substitute(expression({
        d11 = 1/scale^2  + y*(temp6/temp7^2) * ((shape-1) *
              (y*temp7+temp6) - y*temp6 / (log(temp7))^2)
        wz = matrix(0, n, dimm(M))
        wz[,iam(1,1,M)] = dtheta.deta(scale, .lscale, earg = .escale)^2 * d11 -
                          d2theta.deta2(scale, .lscale, earg = .escale) * d1

        if (FALSE && intercept.only) {
            sumw = sum(w)
            for(ii in 1:ncol(wz))
                wz[,ii] = sum(wz[,ii]) / sumw
            pooled.weight = TRUE
            wz = c(w) * wz   # Put back the weights
        } else
            pooled.weight = FALSE
        c(w) * wz
    }), list( .lscale = lscale, .escale = escale ))))
}



betaffqn.control <- function(save.weight = TRUE, ...)
{
    list(save.weight = save.weight)
}



 if (FALSE)
 betaffqn = function(link = "loge", earg = list(),
                    i1 = NULL, i2 = NULL, trim=0.05, A=0, B=1)
{
    if (mode(link) != "character" && mode(link) != "name")
        link = as.character(substitute(link))

    if (!is.Numeric(A, allow = 1) || !is.Numeric(B, allow = 1) || A >= B)
        stop("A must be < B, and both must be of length one")
    stdbeta = (A == 0 && B == 1)  # stdbeta==T iff standard beta distribution
    if (!is.list(earg)) earg = list()

    new("vglmff",
    blurb = c("Two-parameter Beta distribution\n",
            if (stdbeta)
            "y^(shape1-1) * (1-y)^(shape2-1), 0<=y <= 1, shape1>0, shape2>0\n\n"
            else
            paste("(y-",A,")^(shape1-1) * (",B,
            "-y)^(shape2-1), ",A,"<=y <= ",B," shape1>0, shape2>0\n\n", sep = ""),
            "Links:    ",
            namesof("shape1", link, earg = earg),  ", ",
            namesof("shape2", link, earg = earg)),
    initialize = eval(substitute(expression({
        if (ncol(cbind(y)) != 1)
            stop("response must be a vector or a one-column matrix")
        if (min(y) <= .A || max(y) >= .B)
            stop("data not within (A, B)")
        predictors.names = c(namesof("shape1", .link, earg = .earg, short = TRUE),
                             namesof("shape2", .link, earg = .earg, short = TRUE))
        if (is.numeric( .i1) && is.numeric( .i2)) {
            vec = c( .i1, .i2)
            vec = c(theta2eta(vec[1], .link, earg = .earg),
                    theta2eta(vec[2], .link, earg = .earg))
            etastart = matrix(vec, n, 2, byrow= TRUE)
        }

        # For QN update below
        if (length(w) != n || !is.Numeric(w, posit = TRUE))
            stop("weights must be a vector of positive weights")

        if (!length(etastart)) {
            mu1d = mean(y, trim=.trim)
            uu = (mu1d-.A) / ( .B - .A) 
            DD = ( .B - .A)^2 
            pinit = uu^2 * (1-uu)*DD/var(y) - uu   # But var(y) is not robust
            qinit = pinit * (1-uu) / uu
            etastart = matrix(theta2eta(c(pinit,qinit), .link, earg = .earg),
                              n,2,byrow = TRUE)
        }
    }), list( .link = link, .earg = earg, .i1=i1, .i2=i2, .trim=trim, .A = A, .B = B ))), 
    inverse = eval(substitute(function(eta, extra = NULL) {
        shapes = eta2theta(eta, .link, earg = .earg)
        .A + ( .B-.A) * shapes[,1] / (shapes[,1] + shapes[,2])
    }, list( .link = link, .earg = earg, .A = A, .B = B ))),
    last = eval(substitute(expression({
        misc$link = c(shape1 = .link, shape2 = .link)
        misc$earg = list(shape1 = .earg, shape2 = .earg)
        misc$limits = c( .A, .B)
        misc$expected = FALSE
        misc$BFGS = TRUE
    }), list( .link = link, .earg = earg, .A = A, .B = B ))),
    loglikelihood = eval(substitute(
         function(mu, y, w, residuals = FALSE, eta, extra = NULL){
        shapes = eta2theta(eta, .link, earg = .earg)
        temp = lbeta(shapes[,1], shapes[,2])
        if (residuals) stop("loglikelihood residuals ",
                            "not implemented yet") else {



        sum(w * ((shapes[,1]-1)*log(y-.A) + (shapes[,2]-1)*log( .B-y) - temp -
            (shapes[,1]+shapes[,2]-1)*log( .B-.A )))
        }
    }, list( .link = link, .earg = earg, .A = A, .B = B ))),
    vfamily = "betaffqn",
    deriv = eval(substitute(expression({
        shapes = eta2theta(eta, .link, earg = .earg)
        dshapes.deta = dtheta.deta(shapes, .link, earg = .earg)
        dl.dshapes = cbind(log(y-.A), log( .B-y)) - digamma(shapes) +
                     digamma(shapes[,1] + shapes[,2]) - log( .B - .A)
        if (iter == 1) {
            etanew = eta
        } else {
            derivold = derivnew
            etaold = etanew
            etanew = eta
        }
        derivnew = c(w) * dl.dshapes * dshapes.deta
        derivnew
    }), list( .link = link, .earg = earg, .A = A, .B = B ))),
    weight = expression({
        if (iter == 1) {
            wznew = cbind(matrix(w, n, M), matrix(0, n, dimm(M)-M))
        } else {
            wzold = wznew
            wznew = qnupdate(w = w, wzold=wzold, dderiv=(derivold - derivnew),
                             deta=etanew-etaold, M = M,
                             trace=trace)  # weights incorporated in args
        }
        wznew
    }))
}





 logistic2 = function(llocation = "identity",
                     lscale = "loge",
                     elocation = list(),
                     escale = list(),
                     ilocation = NULL, iscale = NULL,
                     imethod = 1, zero = NULL) {
    if (mode(llocation) != "character" && mode(llocation) != "name")
        llocation = as.character(substitute(llocation))
    if (mode(lscale) != "character" && mode(lscale) != "name")
        lscale = as.character(substitute(lscale))
    if (!is.Numeric(imethod, allow = 1, integ = TRUE, posit = TRUE) ||
       imethod > 2) stop("argument 'imethod' must be 1 or 2")
    if (length(zero) && !is.Numeric(zero, integer = TRUE, posit = TRUE))
        stop("bad input for argument 'zero'")
    if (length(iscale) && !is.Numeric(iscale, posit = TRUE))
        stop("bad input for argument 'iscale'")
    if (!is.list(elocation)) elocation = list()
    if (!is.list(escale)) escale = list()

    new("vglmff",
    blurb = c("Two-parameter logistic distribution\n\n",
            "Links:    ",
            namesof("location", llocation, earg = elocation), ", ",
            namesof("scale", lscale, earg = escale),
            "\n", "\n",
            "Mean:     location", "\n",
            "Variance: (pi*scale)^2 / 3"),
    constraints = eval(substitute(expression({
        constraints = cm.zero.vgam(constraints, x, .zero, M)
    }), list( .zero = zero ))),
    initialize = eval(substitute(expression({
        if (ncol(cbind(y)) != 1)
            stop("response must be a vector or a one-column matrix")
        predictors.names = 
        c(namesof("location", .llocation, earg = .elocation, tag = FALSE),
          namesof("scale", .lscale, earg = .escale, tag = FALSE))
        if (!length(etastart)) {
            if ( .imethod == 1) {
                location.init = y
                scale.init = sqrt(3) * sd(y) / pi
            } else {
                location.init = median(rep(y, w))
                scale.init = sqrt(3) * sum(w*(y-location.init)^2) / (sum(w)*pi)
            }
            location.init = if (length( .ilocation)) rep( .ilocation, len = n) else
                             rep(location.init, len = n)
            if ( .llocation == "loge") location.init = abs(location.init) + 0.001
            scale.init = if (length( .iscale)) rep( .iscale, len = n) else
                             rep(1, len = n)
            etastart = cbind(
            theta2eta(location.init, .llocation, earg = .elocation),
            theta2eta(scale.init, .lscale, earg = .escale))
        }
    }), list( .imethod = imethod,
              .elocation = elocation, .escale = escale,
              .llocation = llocation, .lscale = lscale,
              .ilocation = ilocation, .iscale = iscale ))),
    inverse = eval(substitute(function(eta, extra = NULL) {
        eta2theta(eta[,1], .llocation, earg = .elocation)
    }, list( .llocation = llocation,
             .elocation = elocation, .escale = escale ))),
    last = eval(substitute(expression({
        misc$link =    c(location = .llocation, scale = .lscale)
        misc$earg = list(location = .elocation, scale = .escale)
    }), list( .llocation = llocation, .lscale = lscale,
              .elocation = elocation, .escale = escale ))),
    loglikelihood = eval(substitute(
        function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
        location = eta2theta(eta[,1], .llocation, earg = .elocation)
        Scale = eta2theta(eta[,2], .lscale, earg = .escale)
        if (residuals) stop("loglikelihood residuals ",
                            "not implemented yet") else {
            sum(w * dlogis(x = y, location = location,
                           scale = Scale, log = TRUE))
        }
    }, list( .llocation = llocation, .lscale = lscale,
             .elocation = elocation, .escale = escale ))),
    vfamily = c("logistic2"),
    deriv = eval(substitute(expression({
        location = eta2theta(eta[,1], .llocation, earg = .elocation)
        Scale = eta2theta(eta[,2], .lscale, earg = .escale)
        zedd = (y-location) / Scale
        ezedd = exp(-zedd)
        dl.dlocation = (1-ezedd) / ((1 + ezedd) * Scale)
        dlocation.deta = dtheta.deta(location, .llocation, earg = .elocation)
        dl.dscale =  zedd * (1-ezedd) / ((1 + ezedd) * Scale) - 1/Scale
        dscale.deta = dtheta.deta(Scale, .lscale, earg = .escale)
        c(w) * cbind(dl.dlocation * dlocation.deta,
                     dl.dscale * dscale.deta)
    }), list( .llocation = llocation, .lscale = lscale,
              .elocation = elocation, .escale = escale ))),
    weight = eval(substitute(expression({
        d2l.location2 = 1 / (3*Scale^2)
        d2l.dscale2 = (3 + pi^2) / (9*Scale^2)
        wz = matrix(as.numeric(NA), nrow=n, ncol=M) # diagonal
        wz[,iam(1,1,M)] = d2l.location2 * dlocation.deta^2
        wz[,iam(2,2,M)] = d2l.dscale2 * dscale.deta^2
        c(w) * wz
    }), list( .llocation = llocation, .lscale = lscale,
              .elocation = elocation, .escale = escale ))))
}











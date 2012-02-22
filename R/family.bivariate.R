# These functions are
# Copyright (C) 1998-2012 T.W. Yee, University of Auckland.
# All rights reserved.












bilogistic4.control <- function(save.weight = TRUE, ...)
{
    list(save.weight=save.weight)
}


 bilogistic4 = function(llocation = "identity",
                        lscale = "loge",
                        iloc1 = NULL, iscale1 = NULL,
                        iloc2 = NULL, iscale2 = NULL,
                        imethod = 1, zero = NULL) {
    if (mode(llocation) != "character" && mode(llocation) != "name")
        llocation = as.character(substitute(llocation))
    if (mode(lscale) != "character" && mode(lscale) != "name")
        lscale = as.character(substitute(lscale))
    if (!is.Numeric(imethod, allowable.length = 1,
                    integer.valued = TRUE, positive = TRUE) ||
       imethod > 2) stop("imethod must be 1 or 2")

    new("vglmff",
    blurb = c("Bivariate logistic distribution\n\n",
            "Link:    ",
            namesof("location1", llocation), ", ",
            namesof("scale1",    lscale), ", ",
            namesof("location2", llocation), ", ",
            namesof("scale2",    lscale),
            "\n", "\n",
            "Means:     location1, location2"),
    constraints = eval(substitute(expression({
        constraints = cm.zero.vgam(constraints, x, .zero, M)
    }), list( .zero = zero))),
    initialize = eval(substitute(expression({
        if (!is.matrix(y) || ncol(y) != 2)
            stop("the response must be a 2-column matrix") 

        predictors.names = c(namesof("location1", .llocation, tag= FALSE),
                             namesof("scale1", .lscale, tag= FALSE),
                             namesof("location2", .llocation, tag= FALSE),
                             namesof("scale2", .lscale, tag= FALSE))

        if (!length(etastart)) {
            if ( .imethod == 1) {
                location.init1 = y[, 1]
                scale.init1 = sqrt(3) * sd(y[, 1]) / pi
                location.init2 = y[, 2]
                scale.init2 = sqrt(3) * sd(y[, 2]) / pi
            } else {
                location.init1 = median(rep(y[, 1], w))
                location.init2 = median(rep(y[, 2], w))
                scale.init1=sqrt(3)*sum(w*(y[, 1]-location.init1)^2)/(sum(w)*pi)
                scale.init2=sqrt(3)*sum(w*(y[, 2]-location.init2)^2)/(sum(w)*pi)
            }
            loc1.init = if (length(.iloc1)) rep(.iloc1, length.out = n) else
                             rep(location.init1, length.out = n)
            loc2.init = if (length(.iloc2)) rep(.iloc2, length.out = n) else
                             rep(location.init2, length.out = n)
            scale1.init = if (length(.iscale1)) rep(.iscale1, length.out = n) else
                             rep(1, length.out = n)
            scale2.init = if (length(.iscale2)) rep(.iscale2, length.out = n) else
                             rep(1, length.out = n)
            if (.llocation == "loge") location.init1 = abs(location.init1) + 0.001
            if (.llocation == "loge") location.init2 = abs(location.init2) + 0.001
            etastart = cbind(theta2eta(location.init1, .llocation),
                             theta2eta(scale1.init, .lscale),
                             theta2eta(location.init2, .llocation),
                             theta2eta(scale2.init, .lscale))
        }
    }), list(.imethod = imethod, .iloc1=iloc1, .iloc2=iloc2,
             .llocation=llocation,
             .iscale1=iscale1, .iscale2=iscale2, .lscale=lscale))),
    linkinv = function(eta, extra = NULL) {
        cbind(eta[, 1], eta[, 2])
    },
    last = eval(substitute(expression({
        misc$link = c(location1= .llocation, scale1= .lscale,
                      location2= .llocation, scale2= .lscale)
        misc$expected = FALSE
        misc$BFGS = TRUE
    }), list(.lscale=lscale, .llocation=llocation))),
    loglikelihood = eval(substitute(
        function(mu,y,w,residuals= FALSE,eta, extra = NULL) {
        loc1 = eta2theta(eta[, 1], .llocation)
        Scale1 = eta2theta(eta[, 2], .lscale)
        loc2 = eta2theta(eta[, 3], .llocation)
        Scale2 = eta2theta(eta[, 4], .lscale)
        zedd1 = (y[, 1]-loc1) / Scale1
        zedd2 = (y[, 2]-loc2) / Scale2
        if (residuals) stop("loglikelihood residuals not ",
                            "implemented yet") else
        sum(w * (-zedd1 - zedd2 - 3 * log1p(exp(-zedd1)+exp(-zedd2)) -
                 log(Scale1) - log(Scale2)))
    }, list(.lscale=lscale, .llocation=llocation))),
    vfamily = c("bilogistic4"),
    deriv = eval(substitute(expression({
        loc1 = eta2theta(eta[, 1], .llocation)
        Scale1 = eta2theta(eta[, 2], .lscale)
        loc2 = eta2theta(eta[, 3], .llocation)
        Scale2 = eta2theta(eta[, 4], .lscale)
        zedd1 = (y[, 1]-loc1) / Scale1
        zedd2 = (y[, 2]-loc2) / Scale2
        ezedd1 = exp(-zedd1)
        ezedd2 = exp(-zedd2)
        denom = 1 + ezedd1 + ezedd2
        dl.dloc1 = (1 - 3 * ezedd1 / denom) / Scale1
        dl.dloc2 = (1 - 3 * ezedd2 / denom) / Scale2
        dl.dscale1 = (zedd1 - 1 - 3 * ezedd1 * zedd1 / denom) / Scale1
        dl.dscale2 = (zedd2 - 1 - 3 * ezedd2 * zedd2 / denom) / Scale2
        dloc1.deta = dtheta.deta(loc1, .llocation) 
        dloc2.deta = dtheta.deta(loc2, .llocation) 
        dscale1.deta = dtheta.deta(Scale1, .lscale) 
        dscale2.deta = dtheta.deta(Scale2, .lscale) 
        if (iter == 1) {
            etanew = eta
        } else {
            derivold = derivnew
            etaold = etanew
            etanew = eta
        }
        derivnew = c(w) * cbind(dl.dloc1 * dloc1.deta,
                                dl.dscale1 * dscale1.deta,
                                dl.dloc2 * dloc2.deta,
                                dl.dscale2 * dscale2.deta)
        derivnew
    }), list(.lscale=lscale, .llocation=llocation))),
    weight = eval(substitute(expression({
        if (iter == 1) {
            wznew = cbind(matrix(w, n, M), matrix(0, n, dimm(M)-M))
        } else {
            wzold = wznew
            wznew = qnupdate(w=w, wzold=wzold, dderiv=(derivold - derivnew),
                             deta=etanew-etaold, M=M,
                             trace=trace)  # weights incorporated in args
        }
        wznew
    }), list(.lscale=lscale, .llocation=llocation))))
}






dbilogis4 = function(x1, x2, loc1 = 0, scale1 = 1,
                     loc2 = 0, scale2 = 1, log = FALSE) {
    if (!is.logical(log.arg <- log))
        stop("bad input for argument 'log'")
    rm(log)


    L = max(length(x1), length(x2), length(loc1), length(loc2),
            length(scale1), length(scale2))
    x1 = rep(x1, length.out = L); x2 = rep(x2, length.out = L);
    loc1 = rep(loc1, length.out = L); loc2 = rep(loc2, length.out = L);
    scale1 = rep(scale1, length.out = L); scale2 = rep(scale2, length.out = L);
    zedd1 = (-(x1-loc1)/scale1)
    zedd2 = (-(x2-loc2)/scale2)
    logdensity = log(2) + log(zedd1) + log(zedd2) - log(scale1) - 
                 log(scale1) - 3 * log1p(exp(zedd1) + exp(zedd2))
    if (log.arg) logdensity else exp(logdensity)
}



pbilogis4 = function(q1, q2, loc1 = 0, scale1 = 1, loc2 = 0, scale2 = 1) {
  if (!is.Numeric(q1)) stop("bad input for 'q1'")
  if (!is.Numeric(q2)) stop("bad input for 'q2'")
  if (!is.Numeric(scale1, positive = TRUE)) stop("bad input for 'scale1'")
  if (!is.Numeric(scale2, positive = TRUE)) stop("bad input for 'scale2'")


  1 / (1 + exp(-(q1-loc1)/scale1) + exp(-(q2-loc2)/scale2))
}



rbilogis4 = function(n, loc1 = 0, scale1 = 1, loc2 = 0, scale2 = 1) {
    if (!is.Numeric(n, positive = TRUE,
                    allowable.length = 1,integer.valued = TRUE))
      stop("bad input for 'n'")
    if (!is.Numeric(scale1, positive = TRUE))
      stop("bad input for 'scale1'")
    if (!is.Numeric(scale2, positive = TRUE))
      stop("bad input for 'scale2'")
    y1 = rlogis(n, location = loc1, scale = scale1)
    ezedd1 = exp(-(y1-loc1)/scale1)
    y2 = loc2 - scale2 * log(1/sqrt(runif(n) / (1 + ezedd1)^2) - 1 - ezedd1)
    cbind(y1, y2)
}



 freund61 = function(la  = "loge",
                     lap = "loge",
                     lb  = "loge",
                     lbp = "loge",
                     ea  = list(),
                     eap = list(),
                     eb  = list(),
                     ebp = list(),
                     ia = NULL, iap = NULL, ib = NULL, ibp = NULL,
                     independent = FALSE,
                     zero = NULL) {
  if (mode(la) != "character" && mode(la) != "name")
    la = as.character(substitute(la))
  if (mode(lap) != "character" && mode(lap) != "name")
    lap = as.character(substitute(lap))
  if (mode(lb) != "character" && mode(lb) != "name")
    lb = as.character(substitute(lb))
  if (mode(lbp) != "character" && mode(lbp) != "name")
    lbp = as.character(substitute(lbp))

  if (!is.list(ea )) ea  = list()
  if (!is.list(eap)) eap = list()
  if (!is.list(eb )) eb  = list()
  if (!is.list(ebp)) ebp = list()

  new("vglmff",
  blurb = c("Freund (1961) bivariate exponential distribution\n",
         "Links:    ",
         namesof("a",  la,  earg = ea ), ", ",
         namesof("ap", lap, earg = eap), ", ",
         namesof("b",  lb,  earg = eb ), ", ",
         namesof("bp", lbp, earg = ebp)),
  constraints = eval(substitute(expression({
    constraints <- cm.vgam(matrix(c(1, 1,0,0, 0,0, 1, 1), M, 2), x,
                           .independent, constraints,
                           intercept.apply = TRUE)
    constraints = cm.zero.vgam(constraints, x, .zero, M)
  }), list(.independent = independent, .zero = zero))),
  initialize = eval(substitute(expression({
    if (!is.matrix(y) || ncol(y) != 2)
        stop("the response must be a 2 column matrix") 

    predictors.names =
      c(namesof("a",  .la,  earg = .ea , short = TRUE), 
        namesof("ap", .lap, earg = .eap, short = TRUE), 
        namesof("b",  .lb,  earg = .eb , short = TRUE), 
        namesof("bp", .lbp, earg = .ebp, short = TRUE))
    extra$y1.lt.y2 = y[, 1] < y[, 2]

    if (!(arr <- sum(extra$y1.lt.y2)) || arr==n)
        stop("identifiability problem: either all y1<y2 or y2<y1")

    if (!length(etastart)) {
        sumx = sum(y[extra$y1.lt.y2, 1]); sumxp = sum(y[!extra$y1.lt.y2, 1])
        sumy = sum(y[extra$y1.lt.y2, 2]); sumyp = sum(y[!extra$y1.lt.y2, 2])
        if (FALSE) { # Noise:
            arr = min(arr + n/10, n*0.95)
            sumx = sumx * 1.1; sumxp = sumxp * 1.2;
            sumy = sumy * 1.2; sumyp = sumyp * 1.3;
        }
        ainit  = if (length(.ia))  rep(.ia, length.out = n) else
           arr / (sumx + sumyp)
        apinit = if (length(.iap)) rep(.iap,length.out = n) else
           (n-arr)/(sumxp-sumyp)
        binit  = if (length(.ib))  rep(.ib, length.out = n) else
           (n-arr)/(sumx +sumyp)
        bpinit = if (length(.ib))  rep(.ibp,length.out = n) else
           arr / (sumy - sumx)

        etastart =
          cbind(theta2eta(rep(ainit,  length.out = n), .la,  earg = .ea  ),
                theta2eta(rep(apinit, length.out = n), .lap, earg = .eap ),
                theta2eta(rep(binit,  length.out = n), .lb,  earg = .eb  ),
                theta2eta(rep(bpinit, length.out = n), .lbp, earg = .ebp ))
    }
  }), list(.la = la, .lap = lap, .lb = lb, .lbp = lbp,
           .ea = ea, .eap = eap, .eb = eb, .ebp = ebp,
           .ia = ia, .iap = iap, .ib = ib, .ibp = ibp))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    alpha  = eta2theta(eta[, 1], .la,  earg = .ea  )
    alphap = eta2theta(eta[, 2], .lap, earg = .eap )
    beta   = eta2theta(eta[, 3], .lb,  earg = .eb  )
    betap  = eta2theta(eta[, 4], .lbp, earg = .ebp )
    cbind((alphap+beta) / (alphap*(alpha+beta)),
          (alpha+betap) / (betap*(alpha+beta)))
  }, list(.la = la, .lap = lap, .lb = lb, .lbp = lbp,
          .ea = ea, .eap = eap, .eb = eb, .ebp = ebp ))),
  last = eval(substitute(expression({
    misc$link = c("a"= .la, "ap"= .lap, "b"= .lb, "bp"= .lbp)
  }), list(.la = la, .lap = lap, .lb = lb, .lbp = lbp))),
  loglikelihood = eval(substitute(
        function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
    alpha  = eta2theta(eta[, 1], .la,  earg = .ea  )
    alphap = eta2theta(eta[, 2], .lap, earg = .eap )
    beta   = eta2theta(eta[, 3], .lb,  earg = .eb  )
    betap  = eta2theta(eta[, 4], .lbp, earg = .ebp )
    if (residuals) stop("loglikelihood residuals not ",
                        "implemented yet") else {
        tmp88 = extra$y1.lt.y2
        ell1 = log(alpha[tmp88]) + log(betap[tmp88]) -
               betap[tmp88] * y[tmp88, 2] -
               (alpha+beta-betap)[tmp88] * y[tmp88, 1]
        ell2 = log(beta[!tmp88]) + log(alphap[!tmp88]) -
               alphap[!tmp88] * y[!tmp88, 1] -
               (alpha+beta-alphap)[!tmp88] * y[!tmp88, 2]
    sum(w[tmp88] * ell1) + sum(w[!tmp88] * ell2) }
  }, list(.la = la, .lap = lap, .lb = lb, .lbp = lbp,
          .ea = ea, .eap = eap, .eb = eb, .ebp = ebp ))),
  vfamily = c("freund61"),
  deriv = eval(substitute(expression({
    tmp88 = extra$y1.lt.y2
    alpha  = eta2theta(eta[, 1], .la,  earg = .ea  )
    alphap = eta2theta(eta[, 2], .lap, earg = .eap )
    beta   = eta2theta(eta[, 3], .lb,  earg = .eb  )
    betap  = eta2theta(eta[, 4], .lbp, earg = .ebp )

    dalpha.deta  = dtheta.deta(alpha,  .la,  earg = .ea  )
    dalphap.deta = dtheta.deta(alphap, .lap, earg = .eap )
    dbeta.deta   = dtheta.deta(beta,   .lb,  earg = .eb  )
    dbetap.deta  = dtheta.deta(betap,  .lbp, earg = .ebp )

    d1 = 1/alpha - y[, 1]
    d1[!tmp88] = -y[!tmp88, 2]
    d2 = 0 * alphap
    d2[!tmp88] = 1/alphap[!tmp88] - y[!tmp88, 1] + y[!tmp88, 2]
    d3 = -y[, 1]
    d3[!tmp88] = 1/beta[!tmp88] - y[!tmp88, 2]
    d4 = 1/betap - y[, 2] + y[, 1]
    d4[!tmp88] = 0

    c(w) * cbind(d1 * dalpha.deta,
                 d2 * dalphap.deta,
                 d3 * dbeta.deta,
                 d4 * dbetap.deta)
  }), list(.la = la, .lap = lap, .lb = lb, .lbp = lbp,
           .ea = ea, .eap = eap, .eb = eb, .ebp = ebp ))),
  weight = eval(substitute(expression({
    py1.lt.y2 = alpha / (alpha+beta)
    d11 = py1.lt.y2 / alpha^2
    d22 = (1-py1.lt.y2) / alphap^2
    d33 = (1-py1.lt.y2) / beta^2
    d44 = py1.lt.y2 / betap^2
    wz = matrix(0, n, M)  # diagonal
    wz[, iam(1, 1, M)] = dalpha.deta^2  * d11
    wz[, iam(2, 2, M)] = dalphap.deta^2 * d22
    wz[, iam(3, 3, M)] = dbeta.deta^2   * d33
    wz[, iam(4, 4, M)] = dbetap.deta^2  * d44
    c(w) * wz
  }), list(.la = la, .lap = lap, .lb = lb, .lbp = lbp,
           .ea = ea, .eap = eap, .eb = eb, .ebp = ebp ))))
}








 bivgamma.mckay = function(lscale = "loge",
                        lshape1 = "loge",
                        lshape2 = "loge",
                        iscale = NULL,
                        ishape1 = NULL,
                        ishape2 = NULL,
                        imethod = 1,
                        zero = 1) {
    if (mode(lscale) != "character" && mode(lscale) != "name")
        lscale = as.character(substitute(lscale))
    if (mode(lshape1) != "character" && mode(lshape1) != "name")
        lshape1 = as.character(substitute(lshape1))
    if (mode(lshape2) != "character" && mode(lshape2) != "name")
        lshape2 = as.character(substitute(lshape2))
    if (!is.null(iscale))
        if (!is.Numeric(iscale, positive = TRUE))
            stop("'iscale' must be positive or NULL")
    if (!is.null(ishape1))
        if (!is.Numeric(ishape1, positive = TRUE))
            stop("'ishape1' must be positive or NULL")
    if (!is.null(ishape2))
        if (!is.Numeric(ishape2, positive = TRUE))
            stop("'ishape2' must be positive or NULL")
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
        constraints = cm.zero.vgam(constraints, x, .zero, M)
    }), list( .zero = zero ))),
    initialize = eval(substitute(expression({
        if (!is.matrix(y) || ncol(y) != 2)
            stop("the response must be a 2 column matrix")
        if (any(y[, 1] >= y[, 2]))
            stop("the second column minus the first column must be a vector ",
                  "of positive values")
        predictors.names = c(namesof("scale",  .lscale,  short = TRUE), 
                             namesof("shape1", .lshape1, short = TRUE), 
                             namesof("shape2", .lshape2, short = TRUE))
        if (!length(etastart)) {
            momentsY = if ( .imethod == 1) {
                cbind(median(y[, 1]),  # This may not be monotonic
                      median(y[, 2])) + 0.01
            } else {
                cbind(weighted.mean(y[, 1], w),
                      weighted.mean(y[, 2], w))
            }

            mcg2.loglik = function(thetaval, y, x, w, extraargs) {
                ainit = a = thetaval
                momentsY = extraargs$momentsY
                p = (1/a) * abs(momentsY[1]) + 0.01
                q = (1/a) * abs(momentsY[2] - momentsY[1]) + 0.01
                sum(w * (-(p+q)*log(a) - lgamma(p) - lgamma(q) +
                     (p - 1)*log(y[, 1]) + (q - 1)*log(y[, 2]-y[, 1]) - y[, 2] / a ))
            }

            a.grid = if (length( .iscale )) c( .iscale ) else
               c(0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100)
            extraargs = list(momentsY = momentsY)
            ainit = getMaxMin(a.grid, objfun=mcg2.loglik,
                              y=y,  x=x, w=w, maximize = TRUE,
                              extraargs = extraargs)
            ainit = rep(if(is.Numeric( .iscale )) .iscale else ainit, length.out = n)
            pinit = (1/ainit) * abs(momentsY[1]) + 0.01
            qinit = (1/ainit) * abs(momentsY[2] - momentsY[1]) + 0.01

            pinit = rep(if(is.Numeric( .ishape1 )) .ishape1 else pinit, length.out = n)
            qinit = rep(if(is.Numeric( .ishape2 )) .ishape2 else qinit, length.out = n)

            etastart = cbind(theta2eta(ainit, .lscale),
                             theta2eta(pinit, .lshape1),
                             theta2eta(qinit, .lshape2))
        }
    }), list( .lscale=lscale, .lshape1=lshape1, .lshape2=lshape2,
              .iscale=iscale, .ishape1=ishape1, .ishape2=ishape2,
              .imethod = imethod ))),
    linkinv = eval(substitute(function(eta, extra = NULL) {
        a = eta2theta(eta[, 1], .lscale)
        p = eta2theta(eta[, 2], .lshape1)
        q = eta2theta(eta[, 3], .lshape2)
        cbind("y1"=p*a, "y2"=(p+q)*a)
    }, list( .lscale=lscale, .lshape1=lshape1, .lshape2=lshape2 ))),
    last = eval(substitute(expression({
        misc$link = c("scale"= .lscale, "shape1"= .lshape1, "shape2"= .lshape2)
        misc$ishape1 = .ishape1
        misc$ishape2 = .ishape2
        misc$iscale = .iscale
        misc$expected = TRUE
    }), list( .lscale=lscale, .lshape1=lshape1, .lshape2=lshape2,
              .iscale=iscale, .ishape1=ishape1, .ishape2=ishape2 ))),
    loglikelihood = eval(substitute(
            function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
        a = eta2theta(eta[, 1], .lscale)
        p = eta2theta(eta[, 2], .lshape1)
        q = eta2theta(eta[, 3], .lshape2)
        if (residuals) stop("loglikelihood residuals not ",
                            "implemented yet") else
        sum(w * (-(p+q)*log(a) - lgamma(p) - lgamma(q) +
                  (p - 1)*log(y[, 1]) + (q - 1)*log(y[, 2]-y[, 1]) - y[, 2] / a))
    }, list( .lscale=lscale, .lshape1=lshape1, .lshape2=lshape2 ))),
    vfamily = c("bivgamma.mckay"),
    deriv = eval(substitute(expression({
        aparam = eta2theta(eta[, 1], .lscale)
        shape1 = eta2theta(eta[, 2], .lshape1)
        shape2 = eta2theta(eta[, 3], .lshape2)
        dl.da = (-(shape1+shape2) + y[, 2] / aparam) / aparam
        dl.dshape1 = -log(aparam) - digamma(shape1) + log(y[, 1])
        dl.dshape2 = -log(aparam) - digamma(shape2) + log(y[, 2]-y[, 1])
        c(w) * cbind(dl.da      * dtheta.deta(aparam, .lscale),
                     dl.dshape1 * dtheta.deta(shape1, .lshape1),
                     dl.dshape2 * dtheta.deta(shape2, .lshape2))
    }), list( .lscale=lscale, .lshape1=lshape1, .lshape2=lshape2 ))),
    weight = eval(substitute(expression({
        d11 = (shape1+shape2) / aparam^2
        d22 = trigamma(shape1)
        d33 = trigamma(shape2)
        d12 = 1 / aparam
        d13 = 1 / aparam
        d23 = 0
        wz = matrix(0, n, dimm(M))
        wz[, iam(1, 1, M)] = dtheta.deta(aparam, .lscale)^2 * d11
        wz[, iam(2, 2, M)] = dtheta.deta(shape1, .lshape1)^2 * d22
        wz[, iam(3, 3, M)] = dtheta.deta(shape2, .lshape2)^2 * d33
        wz[, iam(1, 2, M)] = dtheta.deta(aparam, .lscale) *
                          dtheta.deta(shape1, .lshape1) * d12
        wz[, iam(1, 3, M)] = dtheta.deta(aparam, .lscale) *
                          dtheta.deta(shape2, .lshape2) * d13
        wz[, iam(2, 3, M)] = dtheta.deta(shape1, .lshape1) *
                          dtheta.deta(shape2, .lshape2) * d23
        c(w) * wz
    }), list( .lscale=lscale, .lshape1=lshape1, .lshape2=lshape2 ))))
}










rfrank = function(n, alpha) {
    if (!is.Numeric(n, positive = TRUE,
                    allowable.length = 1, integer.valued = TRUE))
      stop("bad input for argument 'n'")
    if (!is.Numeric(alpha, positive = TRUE))
      stop("bad input for argument 'alpha'")
    alpha = rep(alpha, length.out = n)
    U = runif(n)
    V = runif(n)
    T = alpha^U + (alpha - alpha^U) * V
    X = U
    index = abs(alpha - 1) < .Machine$double.eps
    Y = U
    if (any(!index))
        Y[!index] = logb(T[!index]/(T[!index]+(1-alpha[!index])*V[!index]),
                         base=alpha[!index])
    ans = matrix(c(X,Y), nrow=n, ncol = 2)
    if (any(index)) {
        ans[index, 1] = runif(sum(index)) # Uniform density for alpha == 1
        ans[index, 2] = runif(sum(index))
    }
    ans
}

pfrank = function(q1, q2, alpha) {
    if (!is.Numeric(q1)) stop("bad input for 'q1'")
    if (!is.Numeric(q2)) stop("bad input for 'q2'")
    if (!is.Numeric(alpha, positive = TRUE)) stop("bad input for 'alpha'")

    L = max(length(q1), length(q2), length(alpha))
    alpha = rep(alpha, length.out = L)
    q1 = rep(q1, length.out = L)
    q2 = rep(q2, length.out = L)

    x=q1; y=q2
    index = (x >= 1 & y <  1) | (y >= 1 & x <  1) |
            (x <= 0 | y <= 0) | (x >= 1 & y >= 1) |
            (abs(alpha - 1) < .Machine$double.eps)
    ans = as.numeric(index)
    if (any(!index))
    ans[!index] = logb(1 + ((alpha[!index])^(x[!index]) - 1)*
                  ((alpha[!index])^(y[!index]) - 1)/(alpha[!index] - 1), 
                  base=alpha[!index])
    ind2 = (abs(alpha - 1) < .Machine$double.eps)
    ans[ind2] = x[ind2] * y[ind2]
    ans[x >= 1 & y <  1] = y[x >= 1 & y < 1]   # P(Y2 < q2) = q2
    ans[y >= 1 & x <  1] = x[y >= 1 & x < 1]   # P(Y1 < q1) = q1
    ans[x <= 0 | y <= 0] = 0
    ans[x >= 1 & y >= 1] = 1
    ans
}

dfrank = function(x1, x2, alpha, log = FALSE) {
    if (!is.logical(log.arg <- log))
        stop("bad input for argument 'log'")
    rm(log)

    if (!is.Numeric(x1)) stop("bad input for 'x1'")
    if (!is.Numeric(x2)) stop("bad input for 'x2'")
    if (!is.Numeric(alpha, positive = TRUE)) stop("bad input for 'alpha'")

    L = max(length(x1), length(x2), length(alpha))
    alpha = rep(alpha, length.out = L)
    x1 = rep(x1, length.out = L)
    x2 = rep(x2, length.out = L)

    if (log.arg) {
        denom = alpha-1 + (alpha^x1  - 1) * (alpha^x2  - 1)
        denom = abs(denom)  # Needed; Genest (1987) uses this too, eqn (4.1)
        log((alpha - 1) * log(alpha)) + (x1+x2)*log(alpha) - 2 * log(denom)
    } else {
        temp = (alpha - 1) + (alpha^x1 - 1) * (alpha^x2 - 1)
        index = (abs(alpha - 1) < .Machine$double.eps)
        ans = x1
        if (any(!index))
            ans[!index] = (alpha[!index] - 1) * log(alpha[!index]) *
                (alpha[!index])^(x1[!index]+x2[!index]) / (temp[!index])^2
        ans[x1 <= 0 | x2 <= 0 | x1 >= 1 | x2 >= 1] = 0
        ans[index] = 1
        ans
    }
}




frank.control <- function(save.weight = TRUE, ...)
{
    list(save.weight=save.weight)
}



 frank = function(lapar = "loge", eapar = list(), iapar = 2, nsimEIM = 250) {
    if (mode(lapar) != "character" && mode(lapar) != "name")
      lapar = as.character(substitute(lapar))
    if (!is.Numeric(iapar, positive = TRUE))
      stop("'iapar' must be positive")

    if (!is.list(eapar)) eapar = list()
    if (length(nsimEIM) &&
       (!is.Numeric(nsimEIM, allowable.length = 1,
                    integer.valued = TRUE) ||
        nsimEIM <= 50))
      stop("'nsimEIM' should be an integer greater than 50")

    new("vglmff",
    blurb = c("Frank's bivariate distribution\n",
           "Links:    ",
           namesof("apar", lapar, earg = eapar )),
    initialize = eval(substitute(expression({
        if (!is.matrix(y) || ncol(y) != 2)
            stop("the response must be a 2 column matrix") 
        if (any(y <= 0) || any(y >= 1))
            stop("the response must have values between 0 and 1") 
        predictors.names =
          c(namesof("apar", .lapar, earg = .eapar, short = TRUE))
        if (length(dimnames(y)))
            extra$dimnamesy2 = dimnames(y)[[2]]
        if (!length(etastart)) {
            apar.init = rep(.iapar, length.out = n)
            etastart = cbind(theta2eta(apar.init, .lapar, earg = .eapar ))
        }
    }), list( .lapar = lapar, .eapar=eapar, .iapar=iapar))),
    linkinv = eval(substitute(function(eta, extra = NULL) {
        apar = eta2theta(eta, .lapar, earg = .eapar )
        fv.matrix = matrix(0.5, length(apar), 2)
        if (length(extra$dimnamesy2))
            dimnames(fv.matrix) = list(names(eta), extra$dimnamesy2)
        fv.matrix
    }, list(.lapar = lapar, .eapar=eapar ))),
    last = eval(substitute(expression({
        misc$link = c("apar"= .lapar)
        misc$earg = list("apar"= .eapar )
        misc$expected = TRUE
        misc$nsimEIM = .nsimEIM
        misc$pooled.weight = pooled.weight
    }), list(.lapar = lapar, .eapar=eapar, .nsimEIM = nsimEIM ))),
    loglikelihood = eval(substitute(
            function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
        apar = eta2theta(eta, .lapar, earg = .eapar )
        if (residuals) stop("loglikelihood residuals not ",
                            "implemented yet") else {
            sum(w * dfrank(x1=y[, 1], x2=y[, 2], alpha=apar, log = TRUE))
        }
    }, list(.lapar = lapar, .eapar=eapar ))),
    vfamily = c("frank"),
    deriv = eval(substitute(expression({
        apar = eta2theta(eta, .lapar, earg = .eapar )
        dapar.deta = dtheta.deta(apar, .lapar, earg = .eapar )

        de3 = deriv3(~ (log((apar - 1) * log(apar)) + (y1+y2)*log(apar) -
                          2 * log(apar-1 + (apar^y1  - 1) * (apar^y2  - 1))),
                        name = "apar", hessian= TRUE)

        denom = apar-1 + (apar^y[, 1]  - 1) * (apar^y[, 2]  - 1)
        tmp700 = 2*apar^(y[, 1]+y[, 2]) - apar^y[, 1] - apar^y[, 2]
        numerator = 1 + y[, 1] * apar^(y[, 1] - 1) * (apar^y[, 2]  - 1) + 
                        y[, 2] * apar^(y[, 2] - 1) * (apar^y[, 1]  - 1)
        Dl.dapar = 1/(apar - 1) + 1/(apar*log(apar)) + (y[, 1]+y[, 2])/apar -
                   2 * numerator / denom
        w * Dl.dapar * dapar.deta
    }), list(.lapar = lapar, .eapar=eapar, .nsimEIM = nsimEIM ))),
    weight = eval(substitute(expression({
    if ( is.Numeric( .nsimEIM)) {

        pooled.weight = FALSE  # For @last


        run.mean = 0
        for(ii in 1:( .nsimEIM )) {
            ysim = rfrank(n,alpha=apar)
            y1 = ysim[, 1]; y2 = ysim[, 2];
            eval.de3 = eval(de3)
            d2l.dthetas2 =  attr(eval.de3, "hessian")
            rm(ysim)
            temp3 = -d2l.dthetas2[, 1, 1]   # M = 1
            run.mean = ((ii - 1) * run.mean + temp3) / ii
        }
        wz = if (intercept.only)
            matrix(mean(run.mean), n, dimm(M)) else run.mean

        wz = wz * dapar.deta^2
        c(w) * wz
    } else {
        nump = apar^(y[, 1]+y[, 2]-2) * (2 * y[, 1] * y[, 2] +
                     y[, 1]*(y[, 1] - 1) + y[, 2]*(y[, 2] - 1)) - 
                     y[, 1]*(y[, 1] - 1) * apar^(y[, 1]-2) - 
                     y[, 2]*(y[, 2] - 1) * apar^(y[, 2]-2)
        D2l.dapar2 = 1/(apar - 1)^2 + (1+log(apar))/(apar*log(apar))^2 +
                     (y[, 1]+y[, 2])/apar^2 + 2 *
                     (nump / denom - (numerator/denom)^2)
        d2apar.deta2 = d2theta.deta2(apar, .lapar)
        wz = w * (dapar.deta^2 * D2l.dapar2 - Dl.dapar * d2apar.deta2)
        if (TRUE && intercept.only) {
            wz = cbind(wz)
            sumw = sum(w)
            for(iii in 1:ncol(wz))
                wz[,iii] = sum(wz[,iii]) / sumw
            pooled.weight = TRUE
            wz = c(w) * wz   # Put back the weights
        } else
            pooled.weight = FALSE
        wz
    }
    }), list( .lapar = lapar, .eapar=eapar, .nsimEIM = nsimEIM ))))
}



 gammahyp = function(ltheta = "loge", itheta = NULL, expected = FALSE) {
    if (mode(ltheta) != "character" && mode(ltheta) != "name")
        ltheta = as.character(substitute(ltheta))
    if (!is.logical(expected) || length(expected) != 1)
        stop("argument 'expected' must be a single logical")

    new("vglmff",
    blurb = c("Gamma hyperbola bivariate distribution\n",
           "Links:    ",
           namesof("theta", ltheta)),
    initialize = eval(substitute(expression({
        if (!is.matrix(y) || ncol(y) != 2)
            stop("the response must be a 2 column matrix") 
        if (any(y[, 1] <= 0) || any(y[, 2] <= 1))
            stop("the response has values that are out of range") 
        predictors.names = c(namesof("theta", .ltheta, short = TRUE))
        if (!length(etastart)) {
            theta.init = if (length( .itheta)) rep(.itheta, length.out = n) else {
                1 / (y[, 2] - 1 + 0.01)
            }
            etastart = cbind(theta2eta(theta.init, .ltheta))
        }
    }), list(.ltheta=ltheta, .itheta=itheta))),
    linkinv = eval(substitute(function(eta, extra = NULL) {
        theta = eta2theta(eta, .ltheta)
        cbind(theta*exp(theta), 1+1/theta)
    }, list(.ltheta=ltheta))),
    last = eval(substitute(expression({
        misc$link = c("theta"= .ltheta)
        misc$expected = .expected 
    }), list(.ltheta=ltheta, .expected=expected))),
    loglikelihood = eval(substitute(
            function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
        theta = eta2theta(eta, .ltheta)
        if (residuals) stop("loglikelihood residuals not ",
                            "implemented yet") else {
            sum(w * (-exp(-theta)*y[, 1]/theta - theta*y[, 2]))
        }
    }, list(.ltheta=ltheta))),
    vfamily = c("gammahyp"),
    deriv = eval(substitute(expression({
        theta = eta2theta(eta, .ltheta)
        Dl.dtheta = exp(-theta) * y[, 1] * (1+theta) / theta^2 - y[, 2]
        Dtheta.deta = dtheta.deta(theta, .ltheta)
        w * Dl.dtheta * Dtheta.deta
    }), list(.ltheta=ltheta))),
    weight = eval(substitute(expression({
        temp300 = 2 + theta * (2 + theta)
        if ( .expected) {
            D2l.dtheta2 = temp300 / theta^2
            wz = w * Dtheta.deta^2 * D2l.dtheta2
        } else {
            D2l.dtheta2 = temp300 * y[, 1] * exp(-theta) / theta^3
            D2theta.deta2 = d2theta.deta2(theta, .ltheta)
            wz = w * (Dtheta.deta^2 * D2l.dtheta2 - Dl.dtheta * D2theta.deta2)
        }
        wz
    }), list( .expected=expected, .ltheta=ltheta))))
}



 morgenstern = function(lapar = "rhobit", earg  = list(), iapar = NULL, tola0 = 0.01,
                        imethod = 1) {
    if (mode(lapar) != "character" && mode(lapar) != "name")
        lapar = as.character(substitute(lapar))
    if (!is.list(earg)) earg = list()

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
        if (!is.matrix(y) || ncol(y) != 2)
            stop("the response must be a 2 column matrix") 
        if (any(y < 0))
            stop("the response must have non-negative values only") 
        predictors.names = c(namesof("apar", .lapar, earg = .earg , short = TRUE))
        if (length(dimnames(y)))
            extra$dimnamesy2 = dimnames(y)[[2]]
        if (!length(etastart)) {
            ainit  = if (length(.iapar))  rep(.iapar, length.out = n) else {
                mean1 = if ( .imethod == 1) median(y[, 1]) else mean(y[, 1])
                mean2 = if ( .imethod == 1) median(y[, 2]) else mean(y[, 2])
                Finit = 0.01 + mean(y[, 1] <= mean1 & y[, 2] <= mean2)
                ((Finit+expm1(-mean1)+exp(-mean2)) / exp(-mean1-mean2) - 1)/(
                 expm1(-mean1) * expm1(-mean2))
              }
            etastart = theta2eta(rep(ainit, length.out = n), .lapar, earg = .earg )
        }
    }), list( .iapar=iapar, .lapar = lapar, .earg = earg,
              .imethod = imethod ))),
    linkinv = eval(substitute(function(eta, extra = NULL) {
        alpha = eta2theta(eta, .lapar, earg = .earg )
        fv.matrix = matrix(1, length(alpha), 2)
        if (length(extra$dimnamesy2))
            dimnames(fv.matrix) = list(names(eta), extra$dimnamesy2)
        fv.matrix
    }, list( .lapar = lapar, .earg = earg ))),
    last = eval(substitute(expression({
        misc$link = c("apar"= .lapar)
        misc$earg = list(apar = .earg)
        misc$expected = FALSE
        misc$pooled.weight = pooled.weight
    }), list( .lapar = lapar, .earg = earg ))),
    loglikelihood = eval(substitute(
            function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
        alpha  = eta2theta(eta, .lapar, earg = .earg )
        alpha[abs(alpha) < .tola0 ] = .tola0
        if (residuals) stop("loglikelihood residuals not ",
                            "implemented yet") else {
        denom = (1 + alpha - 2*alpha*(exp(-y[, 1]) + exp(-y[, 2])) +
                4*alpha*exp(-y[, 1] - y[, 2]))
        sum(w * (-y[, 1] - y[, 2] + log(denom)))
        }
    }, list( .lapar = lapar, .earg = earg, .tola0=tola0 ))),
    vfamily = c("morgenstern"),
    deriv = eval(substitute(expression({
        alpha  = eta2theta(eta, .lapar, earg = .earg )
        alpha[abs(alpha) < .tola0 ] = .tola0
        numerator = 1 - 2*(exp(-y[, 1]) + exp(-y[, 2])) + 4*exp(-y[, 1] - y[, 2])
        denom = (1 + alpha - 2*alpha*(exp(-y[, 1]) + exp(-y[, 2])) +
                4 *alpha*exp(-y[, 1] - y[, 2]))
        dl.dalpha = numerator / denom
        dalpha.deta = dtheta.deta(alpha,  .lapar, earg = .earg )
        c(w) * cbind(dl.dalpha * dalpha.deta)
    }), list( .lapar = lapar, .earg = earg, .tola0=tola0 ))),
    weight = eval(substitute(expression({
        d2l.dalpha2 = dl.dalpha^2
        d2alpha.deta2 = d2theta.deta2(alpha,  .lapar, earg = .earg )
        wz = w * (dalpha.deta^2 * d2l.dalpha2 - d2alpha.deta2 * dl.dalpha)
        if (TRUE &&
           intercept.only) {
            wz = cbind(wz)
            sumw = sum(w)
            for(iii in 1:ncol(wz))
                wz[,iii] = sum(wz[,iii]) / sumw
            pooled.weight = TRUE
            wz = c(w) * wz   # Put back the weights
        } else
            pooled.weight = FALSE
        wz
    }), list( .lapar = lapar, .earg = earg ))))
}




rfgm = function(n, alpha) {
  if (!is.Numeric(n, positive = TRUE,
                  allowable.length = 1, integer.valued = TRUE))
    stop("bad input for argument 'n'")
  if (!is.Numeric(alpha))
    stop("bad input for argument 'alpha'")
  if (any(abs(alpha) > 1))
    stop("argument 'alpha' has values out of range")

  y1 = V1 = runif(n)
  V2 = runif(n)
  temp = 2*y1 - 1
  A = alpha * temp - 1
  B = sqrt(1 - 2 * alpha * temp + (alpha*temp)^2 + 4 * alpha * V2 * temp)
  y2 = 2 * V2 / (B - A)
  matrix(c(y1, y2), nrow = n, ncol = 2)
}



dfgm = function(x1, x2, alpha, log = FALSE) {
    log.arg = log
    rm(log)
    if (!is.Numeric(alpha)) stop("bad input for 'alpha'")
    if (any(abs(alpha) > 1)) stop("'alpha' values out of range")
    if ( !is.logical( log.arg ) || length( log.arg ) != 1 )
        stop("bad input for argument 'log'")

    L = max(length(x1), length(x2), length(alpha))
    if (length(x1) != L)  x1 = rep(x1, length.out = L)
    if (length(x2) != L)  x2 = rep(x2, length.out = L)
    if (length(alpha) != L)  alpha = rep(alpha, length.out = L)
    ans = 0 * x1
    xnok = (x1 <= 0) | (x1 >= 1) | (x2 <= 0) | (x2 >= 1)
    if ( log.arg ) {
        ans[!xnok] = log1p(alpha[!xnok] * (1-2*x1[!xnok]) * (1-2*x2[!xnok]))
        ans[xnok] = log(0)
    } else {
        ans[!xnok] = 1 + alpha[!xnok] * (1-2*x1[!xnok]) * (1-2*x2[!xnok])
        ans[xnok] = 0
        if (any(ans<0))
            stop("negative values in the density (alpha out of range)")
    }
    ans
}


pfgm = function(q1, q2, alpha) {
    if (!is.Numeric(q1)) stop("bad input for 'q1'")
    if (!is.Numeric(q2)) stop("bad input for 'q2'")
    if (!is.Numeric(alpha)) stop("bad input for 'alpha'")
    if (any(abs(alpha) > 1)) stop("'alpha' values out of range")

    L = max(length(q1), length(q2), length(alpha))
    if (length(q1) != L)  q1 = rep(q1, length.out = L)
    if (length(q2) != L)  q2 = rep(q2, length.out = L)
    if (length(alpha) != L)  alpha = rep(alpha, length.out = L)

    x=q1; y=q2
    index = (x >= 1 & y<1) | (y >= 1 & x<1) | (x <= 0 | y <= 0) | (x >= 1 & y >= 1)
    ans = as.numeric(index)
    if (any(!index)) {
        ans[!index] = q1[!index] * q2[!index] * (1 + alpha[!index] *
                      (1-q1[!index])*(1-q2[!index]))
    }
    ans[x >= 1 & y<1] = y[x >= 1 & y<1]   # P(Y2 < q2) = q2
    ans[y >= 1 & x<1] = x[y >= 1 & x<1]   # P(Y1 < q1) = q1
    ans[x <= 0 | y <= 0] = 0
    ans[x >= 1 & y >= 1] = 1
    ans
}



fgm.control <- function(save.weight = TRUE, ...)
{
    list(save.weight=save.weight)
}



 fgm = function(lapar = "rhobit", earg  = list(), iapar = NULL,
                imethod = 1, nsimEIM = 200) {
    if (mode(lapar) != "character" && mode(lapar) != "name")
        lapar = as.character(substitute(lapar))
    if (!is.list(earg)) earg = list()

    if (!is.Numeric(imethod, allowable.length = 1,
                    integer.valued = TRUE, positive = TRUE) ||
       imethod > 2.5)
        stop("argument 'imethod' must be 1 or 2")
    if (!length(nsimEIM) ||
       (!is.Numeric(nsimEIM, allowable.length = 1,
                    integer.valued = TRUE) ||
        nsimEIM <= 50))
      stop("'nsimEIM' should be an integer greater than 50")
    if (length(iapar) &&
       (abs(iapar) >= 1))
      stop("'iapar' should be less than 1 in absolute value")


    new("vglmff",
    blurb = c("Farlie-Gumbel-Morgenstern distribution\n",
           "Links:    ",
           namesof("apar", lapar, earg = earg )),
    initialize = eval(substitute(expression({
        if (!is.matrix(y) || ncol(y) != 2)
            stop("the response must be a 2 column matrix") 
        if (any(y < 0) || any(y > 1))
            stop("the response must have values in the unit square")
        predictors.names = namesof("apar", .lapar, earg = .earg, short = TRUE)
        if (length(dimnames(y)))
            extra$dimnamesy2 = dimnames(y)[[2]]
        if (!length(etastart)) {
            ainit  = if (length( .iapar ))  .iapar else {
                mean1 = if ( .imethod == 1) weighted.mean(y[, 1],w) else
                        median(y[, 1])
                mean2 = if ( .imethod == 1) weighted.mean(y[, 2],w) else
                        median(y[, 2])
                Finit = weighted.mean(y[, 1] <= mean1 & y[, 2] <= mean2, w)
                (Finit / (mean1 * mean2) - 1) / ((1-mean1) * (1-mean2))
            }

            ainit = min(0.95, max(ainit, -0.95))

            etastart = theta2eta(rep(ainit, length.out = n), .lapar, earg = .earg )
        }
    }), list( .iapar=iapar, .lapar = lapar, .earg = earg,
              .imethod = imethod ))),
    linkinv = eval(substitute(function(eta, extra = NULL) {
        alpha = eta2theta(eta, .lapar, earg = .earg )
        fv.matrix = matrix(0.5, length(alpha), 2)
        if (length(extra$dimnamesy2))
            dimnames(fv.matrix) = list(names(eta), extra$dimnamesy2)
        fv.matrix
    }, list( .lapar = lapar, .earg = earg ))),
    last = eval(substitute(expression({
        misc$link = c("apar"= .lapar)
        misc$earg = list(apar = .earg)
        misc$expected = FALSE
        misc$nsimEIM = .nsimEIM
    }), list(.lapar = lapar, .earg = earg, .nsimEIM = nsimEIM ))),
    loglikelihood = eval(substitute(
            function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
        alpha = eta2theta(eta, .lapar, earg = .earg )
        if (residuals) stop("loglikelihood residuals not ",
                            "implemented yet") else {
            sum(w * dfgm(x1=y[, 1], x2=y[, 2], alpha=alpha, log = TRUE))
        }
    }, list( .lapar = lapar, .earg = earg ))),
    vfamily = c("fgm"),
    deriv = eval(substitute(expression({
        alpha  = eta2theta(eta, .lapar, earg = .earg )
        dalpha.deta = dtheta.deta(alpha, .lapar, earg = .earg )
        numerator = (1 - 2 * y[, 1])  * (1 - 2 * y[, 2])
        denom = 1 + alpha * numerator
            mytolerance = .Machine$double.eps
            bad <- (denom <= mytolerance)   # Range violation
            if (any(bad)) {
                cat("There are some range violations in @deriv\n")
                flush.console()
                denom[bad] = 2 * mytolerance
            }
        dl.dalpha = numerator / denom
        c(w) * cbind(dl.dalpha * dalpha.deta)
    }), list( .lapar = lapar, .earg = earg, .nsimEIM = nsimEIM ))),
    weight = eval(substitute(expression({
        run.var = 0
        for(ii in 1:( .nsimEIM )) {
            ysim = rfgm(n, alpha=alpha)
            numerator = (1 - 2 * ysim[, 1])  * (1 - 2 * ysim[, 2])
            denom = 1 + alpha * numerator
            dl.dalpha = numerator / denom
            rm(ysim)
            temp3 = dl.dalpha
            run.var = ((ii - 1) * run.var + temp3^2) / ii
        }
        wz = if (intercept.only)
            matrix(colMeans(cbind(run.var)),
                   n, dimm(M), byrow = TRUE) else cbind(run.var)

        wz = wz * dalpha.deta^2
        c(w) * wz
    }), list( .lapar = lapar, .earg = earg, .nsimEIM = nsimEIM ))))
}



 gumbelIbiv = function(lapar = "identity", earg  = list(),
                       iapar = NULL, imethod = 1) {
    if (mode(lapar) != "character" && mode(lapar) != "name")
        lapar = as.character(substitute(lapar))
    if (!is.list(earg)) earg = list()

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
        if (!is.matrix(y) || ncol(y) != 2)
            stop("the response must be a 2 column matrix") 
        if (any(y < 0))
            stop("the response must have non-negative values only")
        predictors.names = c(namesof("apar", .lapar, earg = .earg , short = TRUE))
        if (!length(etastart)) {
            ainit  = if (length( .iapar ))  rep( .iapar, length.out = n) else {
                mean1 = if ( .imethod == 1) median(y[, 1]) else mean(y[, 1])
                mean2 = if ( .imethod == 1) median(y[, 2]) else mean(y[, 2])
                Finit = 0.01 + mean(y[, 1] <= mean1 & y[, 2] <= mean2)
                (log(Finit+expm1(-mean1)+exp(-mean2))+mean1+mean2)/(mean1*mean2)
            }
            etastart = theta2eta(rep(ainit,  length.out = n), .lapar, earg = .earg )
        }
    }), list( .iapar=iapar, .lapar = lapar, .earg = earg,
              .imethod = imethod ))),
    linkinv = eval(substitute(function(eta, extra = NULL) {
        alpha = eta2theta(eta, .lapar, earg = .earg )
        cbind(rep(1, len=length(alpha)),
              rep(1, len=length(alpha)))
    }, list( .lapar = lapar ))),
    last = eval(substitute(expression({
        misc$link = c("apar"= .lapar)
        misc$earg = list(apar = .earg)
        misc$expected = FALSE
        misc$pooled.weight = pooled.weight
    }), list( .lapar = lapar, .earg = earg ))),
    loglikelihood = eval(substitute(
            function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
        alpha  = eta2theta(eta, .lapar, earg = .earg )
        if (residuals) stop("loglikelihood residuals not ",
                            "implemented yet") else {
            denom = (alpha*y[, 1] - 1) * (alpha*y[, 2] - 1) + alpha
            mytolerance = .Machine$double.xmin
            bad <- (denom <= mytolerance)   # Range violation
            if (any(bad)) {
                cat("There are some range violations in @deriv\n")
                flush.console()
            }
            sum(bad) * (-1.0e10) + 
            sum(w[!bad] * (-y[!bad, 1] - y[!bad, 2] +
                alpha[!bad]*y[!bad, 1]*y[!bad, 2] + log(denom[!bad])))
        }
    }, list( .lapar = lapar, .earg = earg ))),
    vfamily = c("gumbelIbiv"),
    deriv = eval(substitute(expression({
        alpha  = eta2theta(eta, .lapar, earg = .earg )
        numerator = (alpha*y[, 1] - 1)*y[, 2] + (alpha*y[, 2] - 1)*y[, 1] + 1
        denom = (alpha*y[, 1] - 1) * (alpha*y[, 2] - 1) + alpha
        denom = abs(denom)
        dl.dalpha = numerator / denom + y[, 1]*y[, 2]
        dalpha.deta = dtheta.deta(alpha,  .lapar, earg = .earg )
        c(w) * cbind(dl.dalpha * dalpha.deta)
    }), list( .lapar = lapar, .earg = earg ))),
    weight = eval(substitute(expression({
        d2l.dalpha2 = (numerator/denom)^2 - 2*y[, 1]*y[, 2] / denom
        d2alpha.deta2 = d2theta.deta2(alpha, .lapar, earg = .earg )
        wz = w * (dalpha.deta^2 * d2l.dalpha2 - d2alpha.deta2 * dl.dalpha)
        if (TRUE &&
           intercept.only) {
            wz = cbind(wz)
            sumw = sum(w)
            for(iii in 1:ncol(wz))
                wz[,iii] = sum(wz[,iii]) / sumw
            pooled.weight = TRUE
            wz = c(w) * wz   # Put back the weights
        } else
            pooled.weight = FALSE
        wz
    }), list( .lapar = lapar, .earg = earg ))))
}







pplack = function(q1, q2, oratio) {
    if (!is.Numeric(q1)) stop("bad input for 'q1'")
    if (!is.Numeric(q2)) stop("bad input for 'q2'")
    if (!is.Numeric(oratio, positive = TRUE)) stop("bad input for 'oratio'")

    L = max(length(q1), length(q2), length(oratio))
    if (length(q1) != L)  q1 = rep(q1, length.out = L)
    if (length(q2) != L)  q2 = rep(q2, length.out = L)
    if (length(oratio) != L)  oratio = rep(oratio, length.out = L)

    x=q1; y=q2
    index = (x >= 1 & y <  1) | (y >= 1 & x <  1) |
            (x <= 0 | y <= 0) | (x >= 1 & y >= 1) |
            (abs(oratio - 1) < 1.0e-6)  #  .Machine$double.eps
    ans = as.numeric(index)
    if (any(!index)) {
        temp1 = 1 + (oratio[!index]  - 1) * (q1[!index] + q2[!index])
        temp2 = temp1 - sqrt(temp1^2 - 4 * oratio[!index] *
                (oratio[!index] - 1) * q1[!index] * q2[!index])
        ans[!index] = 0.5 * temp2 / (oratio[!index] - 1)
    }

    ind2 = (abs(oratio - 1) < 1.0e-6) # .Machine$double.eps
    ans[ind2] = x[ind2] * y[ind2]
    ans[x >= 1 & y<1] = y[x >= 1 & y<1]   # P(Y2 < q2) = q2
    ans[y >= 1 & x<1] = x[y >= 1 & x<1]   # P(Y1 < q1) = q1
    ans[x <= 0 | y <= 0] = 0
    ans[x >= 1 & y >= 1] = 1
    ans
}



rplack = function(n, oratio) {
    if (!is.Numeric(n, positive = TRUE,
                    allowable.length = 1, integer.valued = TRUE))
      stop("bad input for 'n'")
    if (!is.Numeric(oratio, positive = TRUE))
      stop("bad input for 'oratio'")
    if (length(oratio) != n)  oratio = rep(oratio, length.out = n)

    y1 = U = runif(n)
    V = runif(n)
    Z = V * (1-V)
    y2 = (2*Z*(y1*oratio^2 + 1 - y1) + oratio * (1 - 2 * Z) -
          (1 - 2 * V) *
          sqrt(oratio * (oratio + 4*Z*y1*(1-y1)*(1-oratio)^2))) / (oratio +
          Z*(1-oratio)^2)
    matrix(c(y1, 0.5 * y2), nrow=n, ncol = 2)
}



dplack = function(x1, x2, oratio, log = FALSE) {
    log.arg = log
    rm(log)

    if (!is.Numeric(oratio, positive = TRUE))
      stop("bad input for 'oratio'")
    L = max(length(x1), length(x2), length(oratio))
    if (length(x1) != L)  x1 = rep(x1, length.out = L)
    if (length(x2) != L)  x2 = rep(x2, length.out = L)
    if (length(oratio) != L)  oratio = rep(oratio, length.out = L)
    if ( !is.logical( log.arg ) || length( log.arg ) != 1 )
        stop("bad input for argument 'log'")

    if ( log.arg ) {
        ans = log(oratio) + log1p((oratio - 1) *
              (x1+x2-2*x1*x2)) - 1.5 *
              log((1 + (x1+x2)*(oratio - 1))^2 - 4 * oratio * (oratio - 1)*x1*x2)
        ans[(x1 < 0) | (x1 > 1) | (x2 < 0) | (x2 > 1)] = log(0)
    } else {
        ans = oratio * ((oratio  - 1) * (x1+x2-2*x1*x2) + 1) / ((1 +
              (x1+x2)*(oratio - 1))^2 - 4 * oratio * (oratio - 1)*x1*x2)^1.5
        ans[(x1 < 0) | (x1 > 1) | (x2 < 0) | (x2 > 1)] = 0
    }
    ans
}



plackett.control <- function(save.weight = TRUE, ...)
{
    list(save.weight=save.weight)
}



 plackett = function(link = "loge", earg  = list(),
                     ioratio = NULL, imethod = 1, nsimEIM = 200) {
    if (mode(link) != "character" && mode(link) != "name")
        link = as.character(substitute(link))
    if (!is.list(earg)) earg = list()
    if (length(ioratio) && (!is.Numeric(ioratio, positive = TRUE)))
        stop("'ioratio' must be positive")
    if (!is.Numeric(imethod, allowable.length = 1,
                    integer.valued = TRUE, positive = TRUE) ||
       imethod > 2) stop("imethod must be 1 or 2")

    new("vglmff",
    blurb = c("Plackett distribution\n",
           "Links:    ",
           namesof("oratio", link, earg = earg )),
    initialize = eval(substitute(expression({
        if (!is.matrix(y) || ncol(y) != 2)
            stop("the response must be a 2 column matrix") 
        if (any(y < 0) || any(y > 1))
            stop("the response must have values in the unit square")
        predictors.names = namesof("oratio", .link, earg = .earg, short = TRUE)
        if (length(dimnames(y)))
            extra$dimnamesy2 = dimnames(y)[[2]]
        if (!length(etastart)) {
            orinit = if (length( .ioratio ))  .ioratio else {
                if ( .imethod == 2) {
                    scorp = cor(y)[1, 2]
                    if (abs(scorp) <= 0.1) 1 else
                    if (abs(scorp) <= 0.3) 3^sign(scorp) else
                    if (abs(scorp) <= 0.6) 5^sign(scorp) else
                    if (abs(scorp) <= 0.8) 20^sign(scorp) else 40^sign(scorp)
                } else {
                    y10 = weighted.mean(y[, 1], w)
                    y20 = weighted.mean(y[, 2], w)
                    (0.5 + sum(w[(y[, 1] <  y10) & (y[, 2] <  y20)])) *
                    (0.5 + sum(w[(y[, 1] >= y10) & (y[, 2] >= y20)])) / (
                    ((0.5 + sum(w[(y[, 1] <  y10) & (y[, 2] >= y20)])) *
                     (0.5 + sum(w[(y[, 1] >= y10) & (y[, 2] <  y20)]))))
                }
            }
            etastart = theta2eta(rep(orinit, length.out = n), .link, earg = .earg)
        }
    }), list( .ioratio=ioratio, .link = link, .earg = earg,
              .imethod = imethod ))),
    linkinv = eval(substitute(function(eta, extra = NULL) {
        oratio = eta2theta(eta, .link, earg = .earg )
        fv.matrix = matrix(0.5, length(oratio), 2)
        if (length(extra$dimnamesy2))
            dimnames(fv.matrix) = list(dimnames(eta)[[1]], extra$dimnamesy2)
        fv.matrix
    }, list( .link = link, .earg = earg ))),
    last = eval(substitute(expression({
        misc$link = c("oratio"= .link)
        misc$earg = list(oratio = .earg)
        misc$expected = FALSE
        misc$nsimEIM = .nsimEIM
    }), list( .link = link, .earg = earg,
              .nsimEIM = nsimEIM ))),
    loglikelihood = eval(substitute(
            function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
        oratio = eta2theta(eta, .link, earg = .earg )
        if (residuals) stop("loglikelihood residuals not ",
                            "implemented yet") else {
            sum(w * dplack(x1= y[, 1], x2= y[, 2], oratio=oratio, log = TRUE))
        }
    }, list( .link = link, .earg = earg ))),
    vfamily = c("plackett"),
    deriv = eval(substitute(expression({
        oratio  = eta2theta(eta, .link, earg = .earg )
        doratio.deta = dtheta.deta(oratio, .link, earg = .earg )
        y1 = y[, 1]
        y2 = y[, 2]
        de3 = deriv3(~ (log(oratio) + log(1+(oratio - 1) *
              (y1+y2-2*y1*y2)) - 1.5 *
              log((1 + (y1+y2)*(oratio - 1))^2 - 4 * oratio * (oratio - 1)*y1*y2)),
                        name = "oratio", hessian= FALSE)
        eval.de3 = eval(de3)
        dl.doratio =  attr(eval.de3, "gradient")
        w * dl.doratio * doratio.deta
    }), list( .link = link, .earg = earg ))),
    weight = eval(substitute(expression({
        sd3 = deriv3(~ (log(oratio) + log(1+(oratio - 1) *
              (y1sim+y2sim-2*y1sim*y2sim)) - 1.5 *
              log((1 + (y1sim+y2sim)*(oratio - 1))^2 -
              4 * oratio * (oratio - 1)*y1sim*y2sim)),
                        name = "oratio", hessian= FALSE)
        run.var = 0
        for(ii in 1:( .nsimEIM )) {
            ysim = rplack(n, oratio=oratio)
            y1sim = ysim[, 1]
            y2sim = ysim[, 1]
            eval.sd3 = eval(sd3)
            dl.doratio =  attr(eval.sd3, "gradient")
            rm(ysim, y1sim, y2sim)
            temp3 = dl.doratio
            run.var = ((ii - 1) * run.var + temp3^2) / ii
        }
        wz = if (intercept.only)
            matrix(colMeans(cbind(run.var)),
                   n, dimm(M), byrow = TRUE) else cbind(run.var)

        wz = wz * doratio.deta^2
        c(w) * wz
    }), list( .link = link, .earg = earg, .nsimEIM = nsimEIM ))))
}




damh = function(x1, x2, alpha, log = FALSE) {
    log.arg = log
    rm(log)
    if (!is.Numeric(x1)) stop("bad input for 'x1'")
    if (!is.Numeric(x2)) stop("bad input for 'x2'")
    if (!is.Numeric(alpha)) stop("bad input for 'alpha'")
    if (any(abs(alpha) > 1)) stop("'alpha' values out of range")
    L = max(length(x1), length(x2), length(alpha))
    alpha = rep(alpha, length.out = L)
    x1 = rep(x1, length.out = L)
    x2 = rep(x2, length.out = L)
    temp = 1-alpha*(1-x1)*(1-x2)
    if (log.arg) {
        ans = log1p(-alpha+2*alpha*x1*x2/temp) - 2*log(temp)
        ans[(x1 <= 0) | (x1 >= 1) | (x2 <= 0) | (x2 >= 1)] = log(0)
    } else {
        ans = (1-alpha+2*alpha*x1*x2/temp) / (temp^2)
        ans[(x1 <= 0) | (x1 >= 1) | (x2 <= 0) | (x2 >= 1)] = 0
    }
    ans
}

pamh = function(q1, q2, alpha) {
    if (!is.Numeric(q1)) stop("bad input for 'q1'")
    if (!is.Numeric(q2)) stop("bad input for 'q2'")
    if (!is.Numeric(alpha)) stop("bad input for 'alpha'")
    if (any(abs(alpha) > 1)) stop("'alpha' values out of range")

    L = max(length(q1), length(q2), length(alpha))
    if (length(q1) != L)  q1 = rep(q1, length.out = L)
    if (length(q2) != L)  q2 = rep(q2, length.out = L)
    if (length(alpha) != L)  alpha = rep(alpha, length.out = L)

    x=q1; y=q2
    index = (x >= 1 & y < 1) | (y >= 1 & x <  1) |
            (x <= 0 | y<= 0) | (x >= 1 & y >= 1)
    ans = as.numeric(index)
    if (any(!index)) {
        ans[!index] = (q1[!index]*q2[!index]) / (1 -
                      alpha[!index]*(1-q1[!index])*(1-q2[!index]))
    }
    ans[x >= 1 & y <  1] = y[x >= 1 & y < 1]   # P(Y2 < q2) = q2
    ans[y >= 1 & x <  1] = x[y >= 1 & x < 1]   # P(Y1 < q1) = q1
    ans[x <= 0 | y <= 0] = 0
    ans[x >= 1 & y >= 1] = 1
    ans
}

ramh = function(n, alpha) {
    if (!is.Numeric(n, positive = TRUE, allowable.length = 1,
                    integer.valued = TRUE))
      stop("bad input for 'n'")
    if (!is.Numeric(alpha))
      stop("bad input for 'alpha'")
    if (any(abs(alpha) > 1))
      stop("'alpha' values out of range")

    U1 = V1 = runif(n)
    V2 = runif(n)
    b = 1-V1
    A = -alpha*(2*b*V2+1)+2*alpha^2*b^2*V2+1
    B = alpha^2*(4*b^2*V2-4*b*V2+1)+alpha*(4*V2-4*b*V2-2)+1
    U2 = (2*V2*(alpha*b - 1)^2)/(A+sqrt(B))
    matrix(c(U1,U2), nrow=n, ncol = 2)
}


amh.control <- function(save.weight = TRUE, ...)
{
    list(save.weight=save.weight)
}


 amh = function(lalpha = "rhobit", ealpha = list(), ialpha = NULL,
                imethod = 1, nsimEIM = 250)
{
    if (mode(lalpha) != "character" && mode(lalpha) != "name")
      lalpha = as.character(substitute(lalpha))
    if (!is.list(ealpha)) ealpha = list()

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
        if (!is.matrix(y) || ncol(y) != 2)
            stop("the response must be a 2 column matrix")
        if (any(y < 0) || any(y > 1))
            stop("the response must have values in the unit square")
        predictors.names=c(namesof("alpha", .lalpha, earg = .ealpha, short = TRUE))
        if (length(dimnames(y)))
            extra$dimnamesy2 = dimnames(y)[[2]]
        if (!length(etastart)) {
            ainit  = if (length( .ialpha ))  .ialpha else {
                mean1 = if ( .imethod == 1) weighted.mean(y[, 1],w) else
                        median(y[, 1])
                mean2 = if ( .imethod == 1) weighted.mean(y[, 2],w) else
                        median(y[, 2])
                Finit = weighted.mean(y[, 1] <= mean1 & y[, 2] <= mean2, w)
                (1 - (mean1 * mean2 / Finit)) / ((1-mean1) * (1-mean2))
            }
            ainit = min(0.95, max(ainit, -0.95))
            etastart = theta2eta(rep(ainit, length.out = n), .lalpha, earg = .ealpha )
        }
    }), list( .lalpha = lalpha, .ealpha = ealpha, .ialpha=ialpha,
              .imethod = imethod))),
    linkinv = eval(substitute(function(eta, extra = NULL) {
        alpha = eta2theta(eta, .lalpha, earg = .ealpha )
        fv.matrix = matrix(0.5, length(alpha), 2)
        if (length(extra$dimnamesy2))
            dimnames(fv.matrix) = list(names(eta), extra$dimnamesy2)
        fv.matrix
    }, list(.lalpha = lalpha, .ealpha = ealpha ))),
    last = eval(substitute(expression({
        misc$link = c("alpha"= .lalpha)
        misc$earg = list("alpha"= .ealpha )
        misc$expected = TRUE
        misc$nsimEIM = .nsimEIM
    }), list(.lalpha = lalpha, .ealpha = ealpha, .nsimEIM = nsimEIM ))),
    loglikelihood = eval(substitute(
            function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
        alpha = eta2theta(eta, .lalpha, earg = .ealpha )
        if (residuals) stop("loglikelihood residuals not ",
                            "implemented yet") else {
            sum(w * damh(x1=y[, 1], x2=y[, 2], alpha=alpha, log = TRUE))
        }
    }, list( .lalpha = lalpha, .earg = ealpha ))),
    vfamily = c("amh"),
    deriv = eval(substitute(expression({
        alpha = eta2theta(eta, .lalpha, earg = .ealpha )
        dalpha.deta = dtheta.deta(alpha, .lalpha, earg = .ealpha )
        y1 = y[, 1]
        y2 = y[, 2]
        de3 = deriv3(~ (log(1-alpha+(2*alpha*y1*y2/(1-alpha*(1-y1)*(1-y2))))-
                        2*log(1-alpha*(1-y1)*(1-y2))) ,
                        name = "alpha", hessian= FALSE)
        eval.de3 = eval(de3)
        dl.dalpha =  attr(eval.de3, "gradient")
        w * dl.dalpha * dalpha.deta
    }), list(.lalpha = lalpha, .ealpha = ealpha ))),
    weight = eval(substitute(expression({
        sd3 = deriv3(~ (log(1-alpha+
                        (2*alpha*y1sim*y2sim/(1-alpha*(1-y1sim)*(1-y2sim))))-
                        2*log(1-alpha*(1-y1sim)*(1-y2sim))) ,
                        name = "alpha", hessian= FALSE)
        run.var = 0
        for(ii in 1:( .nsimEIM )) {
            ysim = ramh(n, alpha=alpha)
            y1sim = ysim[, 1]
            y2sim = ysim[, 1]
            eval.sd3 = eval(sd3)
            dl.alpha =  attr(eval.sd3, "gradient")
            rm(ysim, y1sim, y2sim)
            temp3 = dl.dalpha
            run.var = ((ii - 1) * run.var + temp3^2) / ii
        }
        wz = if (intercept.only)
            matrix(colMeans(cbind(run.var)),
                   n, dimm(M), byrow = TRUE) else cbind(run.var)
        wz = wz * dalpha.deta^2
        c(w) * wz
    }), list( .lalpha = lalpha, .ealpha = ealpha, .nsimEIM = nsimEIM ))))
}















dbinorm = function(x1, x2, mean1 = 0, mean2 = 0, sd1 = 1, sd2 = 1,
                   rho = 0, log = FALSE) {
  log.arg = log
  rm(log)

  temp5 = 1 - rho^2
  zedd1 = (x1 - mean1) / sd1
  zedd2 = (x2 - mean2) / sd2
  logpdf = -log(2 * pi) - log(sd1 ) - log(sd2) -
            0.5 * log1p(-rho^2) +
           -(0.5 / temp5)  * (zedd1^2 - 2 * rho * zedd1 * zedd2 + zedd2^2)
  if (log.arg) logpdf else exp(logpdf)
}




 binormal = function(lmean1 = "identity", emean1 = list(),
                     lmean2 = "identity", emean2 = list(),
                     lsd1   = "loge",     esd1   = list(),
                     lsd2   = "loge",     esd2   = list(),
                     lrho   = "rhobit",   erho   = list(),
                     imean1 = NULL,       imean2 = NULL,
                     isd1   = NULL,       isd2   = NULL,
                     irho   = NULL,       imethod = 1,
                     equalmean = FALSE,   equalsd = FALSE,
                     zero = 3:5) {
  if (mode(lmean1) != "character" && mode(lmean1) != "name")
    lmean1 = as.character(substitute(lmean1))
  if (mode(lmean2) != "character" && mode(lmean2) != "name")
    lmean2 = as.character(substitute(lmean2))
  if (mode(lsd1  ) != "character" && mode(lsd1  ) != "name")
    lsd1   = as.character(substitute(lsd1  ))
  if (mode(lsd2  ) != "character" && mode(lsd2  ) != "name")
    lsd2   = as.character(substitute(lsd2  ))
  if (mode(lrho  ) != "character" && mode(lrho  ) != "name")
    lrho   = as.character(substitute(lrho  ))

  if (!is.list(emean1)) emean1 = list()
  if (!is.list(emean2)) emean2 = list()
  if (!is.list(esd1  )) esd1   = list()
  if (!is.list(esd2  )) esd2   = list()
  if (!is.list(erho  )) erho   = list()

  trivial1 = is.logical(equalmean) && length(equalmean) == 1 && !equalmean
  trivial2 = is.logical(equalsd  ) && length(equalsd  ) == 1 && !equalsd
  if(!trivial1 && !trivial2)
    stop("only one of 'equalmean' and 'equalsd' can be assigned a value")

  if (!is.Numeric(imethod, allowable.length = 1,
                    integer.valued = TRUE, positive = TRUE) ||
     imethod > 2) stop("argument 'imethod' must be 1 or 2")

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
    constraints <- cm.vgam(temp8.m, x, .equalmean,
                           constraints, intercept.apply = TRUE)
    constraints <- cm.vgam(temp8.s, x, .equalsd,
                           constraints, intercept.apply = TRUE)
    constraints = cm.zero.vgam(constraints, x, .zero, M)
  }), list( .zero = zero,
            .equalsd   = equalsd,
            .equalmean = equalmean ))),
  initialize = eval(substitute(expression({
    if (!is.matrix(y) || ncol(y) != 2)
      stop("the response must be a 2 column matrix") 

    predictors.names = c(
      namesof("mean1", .lmean1, earg = .emean1, short = TRUE),
      namesof("mean2", .lmean2, earg = .emean2, short = TRUE),
      namesof("sd1",   .lsd1,   earg = .esd1,   short = TRUE),
      namesof("sd2",   .lsd2,   earg = .esd2,   short = TRUE),
      namesof("rho",   .lrho,   earg = .erho,   short = TRUE))

    if (length(dimnames(y)))
      extra$dimnamesy2 = dimnames(y)[[2]]

    if (!length(etastart)) {
      imean1 = rep(if (length( .imean1 )) .imean1 else
                   weighted.mean(y[, 1], w = w), length.out = n)
      imean2 = rep(if (length( .imean2 )) .imean2 else
                   weighted.mean(y[, 2], w = w), length.out = n)
      isd1   = rep(if (length( .isd1 )) .isd1 else  sd(y[, 1]), length.out = n)
      isd2   = rep(if (length( .isd2 )) .isd2 else  sd(y[, 2]), length.out = n)
      irho   = rep(if (length( .irho )) .irho else cor(y[, 1], y[, 2]),
                   length.out = n)

      if ( .imethod == 2) {
        imean1 = abs(imean1) + 0.01
        imean2 = abs(imean2) + 0.01
      }
      etastart = cbind(theta2eta(imean1, .lmean1, earg = .emean1),
                       theta2eta(imean2, .lmean2, earg = .emean2),
                       theta2eta(isd1,   .lsd1,   earg = .esd1),
                       theta2eta(isd2,   .lsd2,   earg = .esd2),
                       theta2eta(irho,   .lrho,   earg = .erho))
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
    mean1 = eta2theta(eta[, 1], .lmean1, earg = .emean1)
    mean2 = eta2theta(eta[, 2], .lmean2, earg = .emean2)
    fv.matrix = cbind(mean1, mean2)
    if (length(extra$dimnamesy2))
      dimnames(fv.matrix) = list(names(eta), extra$dimnamesy2)
    fv.matrix
  }  , list( .lmean1 = lmean1, .lmean2 = lmean2,
             .emean1 = emean1, .emean2 = emean2,
             .lsd1   = lsd1  , .lsd2   = lsd2  , .lrho = lrho,
             .esd1   = esd1  , .esd2   = esd2  , .erho = erho ))),

  last = eval(substitute(expression({
    misc$link =   c("mean1" = .lmean1,
                    "mean2" = .lmean2,
                    "sd1"   = .lsd1,
                    "sd2"   = .lsd2,
                    "rho"   = .lrho)
    misc$earg = list("mean1" = .emean1,
                     "mean2" = .emean2, 
                     "sd1"   = .esd1,
                     "sd2"   = .esd2,
                     "rho"   = .erho)
    misc$expected = TRUE
  }) , list( .lmean1 = lmean1, .lmean2 = lmean2,
             .emean1 = emean1, .emean2 = emean2,
             .lsd1   = lsd1  , .lsd2   = lsd2  , .lrho = lrho,
             .esd1   = esd1  , .esd2   = esd2  , .erho = erho ))),
  loglikelihood = eval(substitute(
          function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
    mean1 = eta2theta(eta[, 1], .lmean1, earg = .emean1)
    mean2 = eta2theta(eta[, 2], .lmean2, earg = .emean2)
    sd1   = eta2theta(eta[, 3], .lsd1  , earg = .esd1  )
    sd2   = eta2theta(eta[, 4], .lsd2  , earg = .esd2  )
    Rho   = eta2theta(eta[, 5], .lrho  , earg = .erho  )

    if (residuals) stop("loglikelihood residuals not ",
                          "implemented yet") else {
      sum(w * dbinorm(x1 = y[, 1], x2 = y[, 2],
                      mean1 = mean1, mean2 = mean2,
                      sd1 = sd1, sd2 = sd2, rho = Rho, log = TRUE))
    }
  } , list( .lmean1 = lmean1, .lmean2 = lmean2,
            .emean1 = emean1, .emean2 = emean2,
            .lsd1   = lsd1  , .lsd2   = lsd2  , .lrho = lrho,
            .esd1   = esd1  , .esd2   = esd2  , .erho = erho,
            .imethod = imethod ))),
  vfamily = c("binormal"),
  deriv = eval(substitute(expression({
    mean1 = eta2theta(eta[, 1], .lmean1, earg = .emean1)
    mean2 = eta2theta(eta[, 2], .lmean2, earg = .emean2)
    sd1   = eta2theta(eta[, 3], .lsd1  , earg = .esd1  )
    sd2   = eta2theta(eta[, 4], .lsd2  , earg = .esd2  )
    Rho   = eta2theta(eta[, 5], .lrho  , earg = .erho  )

    zedd1 = (y[, 1] - mean1) / sd1
    zedd2 = (y[, 2] - mean2) / sd2
    temp5 = 1 - Rho^2

    SigmaInv = matrix(0, n, dimm(2))
    SigmaInv[, iam(1, 1, M = 2)] = 1 / ((sd1^2) * temp5)
    SigmaInv[, iam(2, 2, M = 2)] = 1 / ((sd2^2) * temp5)
    SigmaInv[, iam(1, 2, M = 2)] = -Rho / (sd1 * sd2 * temp5)
    dl.dmeans = mux22(t(SigmaInv), y - cbind(mean1, mean2), M = 2,
                      as.matrix = TRUE)
    dl.dsd1   = -1 / sd1 + zedd1 * (zedd1 - Rho * zedd2) / (sd1 * temp5)
    dl.dsd2   = -1 / sd2 + zedd2 * (zedd2 - Rho * zedd1) / (sd2 * temp5)
    dl.drho   = -Rho * (zedd1^2 - 2 * Rho * zedd1 * zedd2 +
                        zedd2^2) / temp5^2 +
                zedd1 * zedd2 / temp5 +
                Rho / temp5

    dmean1.deta = dtheta.deta(mean1, .lmean1) 
    dmean2.deta = dtheta.deta(mean2, .lmean2) 
    dsd1.deta   = dtheta.deta(sd1  , .lsd1  ) 
    dsd2.deta   = dtheta.deta(sd2  , .lsd2  ) 
    drho.deta   = dtheta.deta(Rho  , .lrho  ) 
    dthetas.detas  = cbind(dmean1.deta,
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
    wz = matrix(0.0, n, dimm(M))
    wz[, iam(1, 1, M)] = SigmaInv[, iam(1, 1, M = 2)]
    wz[, iam(2, 2, M)] = SigmaInv[, iam(2, 2, M = 2)]
    wz[, iam(1, 2, M)] = SigmaInv[, iam(1, 2, M = 2)]
    wz[, iam(3, 3, M)] = (1 + 1 / temp5) / sd1^2
    wz[, iam(4, 4, M)] = (1 + 1 / temp5) / sd2^2
    wz[, iam(3, 4, M)] = -(Rho^2) / (temp5 * sd1 * sd2)
    wz[, iam(5, 5, M)] = 2 * (1 + 2 * Rho^2) / temp5^2 -
                         (1 + Rho^2) / temp5^2
    wz[, iam(3, 5, M)] = -Rho / (sd1 * temp5)
    wz[, iam(4, 5, M)] = -Rho / (sd2 * temp5)
    for (ilocal in 1:M)
      for (jlocal in ilocal:M)
        wz[, iam(ilocal, jlocal, M)] = wz[, iam(ilocal, jlocal, M)] *
                                       dthetas.detas[, ilocal] *
                                       dthetas.detas[, jlocal]
      c(w) * wz
  }), list( .lmean1 = lmean1, .lmean2 = lmean2,
            .emean1 = emean1, .emean2 = emean2,
            .lsd1   = lsd1  , .lsd2   = lsd2  , .lrho = lrho,
            .esd1   = esd1  , .esd2   = esd2  , .erho = erho,
            .imethod = imethod ))))
}



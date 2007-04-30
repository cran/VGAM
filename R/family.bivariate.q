# These functions are
# Copyright (C) 1998-2007 T.W. Yee, University of Auckland. All rights reserved.







bilogistic4.control <- function(save.weight=TRUE, ...)
{
    list(save.weight=save.weight)
}

bilogistic4 = function(llocation="identity",
                       lscale="loge",
                       iloc1=NULL, iscale1=NULL,
                       iloc2=NULL, iscale2=NULL,
                       method.init=1, zero=NULL) {
    if(mode(llocation) != "character" && mode(llocation) != "name")
        llocation = as.character(substitute(llocation))
    if(mode(lscale) != "character" && mode(lscale) != "name")
        lscale = as.character(substitute(lscale))
    if(!is.Numeric(method.init, allow=1, integ=TRUE, posit=TRUE) ||
       method.init > 2) stop("method.init must be 1 or 2")

    new("vglmff",
    blurb=c("Bivariate logistic distribution\n\n",
            "Link:    ",
            namesof("location1", llocation), ", ",
            namesof("scale1", lscale), ", ",
            namesof("location2", llocation), ", ",
            namesof("scale2", lscale),
            "\n", "\n",
            "Means:     location1, location2"),
    constraints=eval(substitute(expression({
        constraints = cm.zero.vgam(constraints, x, .zero, M)
    }), list(.zero=zero))),
    initialize=eval(substitute(expression({
        if(!is.matrix(y) || ncol(y) != 2)
            stop("the response must be a 2 column matrix") 
        predictors.names = c(namesof("location1", .llocation, tag= FALSE),
                             namesof("scale1", .lscale, tag= FALSE),
                             namesof("location2", .llocation, tag= FALSE),
                             namesof("scale2", .lscale, tag= FALSE))
        if(!length(etastart)) {
            if( .method.init == 1) {
                location.init1 = y[,1]
                scale.init1 = sqrt(3) * sd(y[,1]) / pi
                location.init2 = y[,2]
                scale.init2 = sqrt(3) * sd(y[,2]) / pi
            } else {
                location.init1 = median(rep(y[,1], w))
                location.init2 = median(rep(y[,2], w))
                scale.init1=sqrt(3)*sum(w*(y[,1]-location.init1)^2)/(sum(w)*pi)
                scale.init2=sqrt(3)*sum(w*(y[,2]-location.init2)^2)/(sum(w)*pi)
            }
            loc1.init = if(length(.iloc1)) rep(.iloc1, len=n) else
                             rep(location.init1, len=n)
            loc2.init = if(length(.iloc2)) rep(.iloc2, len=n) else
                             rep(location.init2, len=n)
            scale1.init = if(length(.iscale1)) rep(.iscale1, len=n) else
                             rep(1, len=n)
            scale2.init = if(length(.iscale2)) rep(.iscale2, len=n) else
                             rep(1, len=n)
            if(.llocation=="loge") location.init1 = abs(location.init1) + 0.001
            if(.llocation=="loge") location.init2 = abs(location.init2) + 0.001
            etastart = cbind(theta2eta(location.init1, .llocation),
                             theta2eta(scale1.init, .lscale),
                             theta2eta(location.init2, .llocation),
                             theta2eta(scale2.init, .lscale))
        }
    }), list(.method.init=method.init, .iloc1=iloc1, .iloc2=iloc2,
             .llocation=llocation,
             .iscale1=iscale1, .iscale2=iscale2, .lscale=lscale))),
    inverse=function(eta, extra=NULL) {
        cbind(eta[,1], eta[,2])
    },
    last=eval(substitute(expression({
        misc$link = c(location1= .llocation, scale1= .lscale,
                      location2= .llocation, scale2= .lscale)
        misc$expected = FALSE
        misc$BFGS = TRUE
    }), list(.lscale=lscale, .llocation=llocation))),
    loglikelihood=eval(substitute(
        function(mu,y,w,residuals= FALSE,eta, extra=NULL) {
        loc1 = eta2theta(eta[,1], .llocation)
        Scale1 = eta2theta(eta[,2], .lscale)
        loc2 = eta2theta(eta[,3], .llocation)
        Scale2 = eta2theta(eta[,4], .lscale)
        zedd1 = (y[,1]-loc1) / Scale1
        zedd2 = (y[,2]-loc2) / Scale2
        if(residuals) stop("loglikelihood residuals not implemented yet") else
        sum(w * (-zedd1 - zedd2 - 3 * log(1+exp(-zedd1)+exp(-zedd2)) -
                 log(Scale1) - log(Scale2)))
    }, list(.lscale=lscale, .llocation=llocation))),
    vfamily=c("bilogistic4"),
    deriv=eval(substitute(expression({
        loc1 = eta2theta(eta[,1], .llocation)
        Scale1 = eta2theta(eta[,2], .lscale)
        loc2 = eta2theta(eta[,3], .llocation)
        Scale2 = eta2theta(eta[,4], .lscale)
        zedd1 = (y[,1]-loc1) / Scale1
        zedd2 = (y[,2]-loc2) / Scale2
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
        if(iter == 1) {
            etanew = eta
        } else {
            derivold = derivnew
            etaold = etanew
            etanew = eta
        }
        derivnew = w * cbind(dl.dloc1 * dloc1.deta,
                             dl.dscale1 * dscale1.deta,
                             dl.dloc2 * dloc2.deta,
                             dl.dscale2 * dscale2.deta)
        derivnew
    }), list(.lscale=lscale, .llocation=llocation))),
    weight=eval(substitute(expression({
        if(iter == 1) {
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



rbilogis4 = function(n, loc1=0, scale1=1, loc2=0, scale2=1) {
    if(!is.Numeric(n, posit=TRUE, allow=1, integ=TRUE)) stop("bad input for n")
    if(!is.Numeric(scale1, posit=TRUE)) stop("bad input for \"scale1\"")
    if(!is.Numeric(scale2, posit=TRUE)) stop("bad input for \"scale2\"")
    y1 = rlogis(n, loc=loc1, scale=scale1)
    ezedd1 = exp(-(y1-loc1)/scale1)
    y2 = loc2 - scale2 * log(1/sqrt(runif(n) / (1 + ezedd1)^2) - 1 - ezedd1)
    cbind(y1, y2)
}

pbilogis4 = function(q1, q2, loc1=0, scale1=1, loc2=0, scale2=1) {
    if(!is.Numeric(q1)) stop("bad input for \"q1\"")
    if(!is.Numeric(q2)) stop("bad input for \"q2\"")
    if(!is.Numeric(scale1, posit=TRUE)) stop("bad input for \"scale1\"")
    if(!is.Numeric(scale2, posit=TRUE)) stop("bad input for \"scale2\"")


    1 / (1 + exp(-(q1-loc1)/scale1) + exp(-(q2-loc2)/scale2))
}

dbilogis4 = function(x1, x2, loc1=0, scale1=1, loc2=0, scale2=1) {
    if(!is.Numeric(x1)) stop("bad input for \"x1\"")
    if(!is.Numeric(x2)) stop("bad input for \"x2\"")
    if(!is.Numeric(scale1, posit=TRUE)) stop("bad input for \"scale1\"")
    if(!is.Numeric(scale2, posit=TRUE)) stop("bad input for \"scale2\"")
    ezedd1 = exp(-(x1-loc1)/scale1)
    ezedd2 = exp(-(x2-loc2)/scale2)
    2 * ezedd1 * ezedd2 / (scale1 * scale2 * (1 + ezedd1 + ezedd2)^3)
}




freund61 = function(la="loge",
                    lap="loge",
                    lb="loge",
                    lbp="loge",
                    ia=NULL, iap=NULL, ib=NULL, ibp=NULL,
                    independent=FALSE,
                    zero=NULL) {
    if(mode(la) != "character" && mode(la) != "name")
        la = as.character(substitute(la))
    if(mode(lap) != "character" && mode(lap) != "name")
        lap = as.character(substitute(lap))
    if(mode(lb) != "character" && mode(lb) != "name")
        lb = as.character(substitute(lb))
    if(mode(lbp) != "character" && mode(lbp) != "name")
        lbp = as.character(substitute(lbp))

    new("vglmff",
    blurb=c("Freund (1961) Bivariate Exponential Distribution\n",
           "Links:    ",
           namesof("a", la), ", ",
           namesof("ap", lap), ", ",
           namesof("b", lb), ", ",
           namesof("bp", lbp)),
    constraints=eval(substitute(expression({
        constraints <- cm.vgam(matrix(c(1,1,0,0, 0,0,1,1),M,2), x,
                               .independent, constraints, intercept.apply=TRUE)
        constraints = cm.zero.vgam(constraints, x, .zero, M)
    }), list(.independent=independent, .zero=zero))),
    initialize=eval(substitute(expression({
        if(!is.matrix(y) || ncol(y) != 2)
            stop("the response must be a 2 column matrix") 
        predictors.names = c(namesof("a", .la, short=TRUE), 
                             namesof("ap", .lap, short=TRUE), 
                             namesof("b", .lb, short=TRUE), 
                             namesof("bp", .lbp, short=TRUE))
        extra$y1.lt.y2 = y[,1] < y[,2]
        if(!(arr <- sum(extra$y1.lt.y2)) || arr==n)
            stop("identifiability problem: either all y1<y2 or y2<y1")
        if(!length(etastart)) {
            sumx = sum(y[extra$y1.lt.y2,1]); sumxp = sum(y[!extra$y1.lt.y2,1])
            sumy = sum(y[extra$y1.lt.y2,2]); sumyp = sum(y[!extra$y1.lt.y2,2])
            if(FALSE) { # Noise:
                arr = min(arr + n/10, n*0.95)
                sumx = sumx * 1.1; sumxp = sumxp * 1.2;
                sumy = sumy * 1.2; sumyp = sumyp * 1.3;
            }
            ainit  = if(length(.ia))  rep(.ia, len=n) else arr / (sumx + sumyp)
            apinit = if(length(.iap)) rep(.iap,len=n) else (n-arr)/(sumxp-sumyp)
            binit  = if(length(.ib))  rep(.ib, len=n) else (n-arr)/(sumx +sumyp)
            bpinit = if(length(.ib))  rep(.ibp,len=n) else arr / (sumy - sumx)
            etastart = cbind(theta2eta(rep(ainit,  len=n), .la),
                             theta2eta(rep(apinit, len=n), .lap),
                             theta2eta(rep(binit,  len=n), .lb),
                             theta2eta(rep(bpinit, len=n), .lbp))
        }
    }), list(.la=la, .lap=lap, .lb=lb, .lbp=lbp, .ia=ia, .iap=iap,
             .ib=ib, .ibp=ibp))),
    inverse=eval(substitute(function(eta, extra=NULL) {
        alpha  = eta2theta(eta[,1], .la)
        alphap = eta2theta(eta[,2], .lap)
        beta   = eta2theta(eta[,3], .lb)
        betap  = eta2theta(eta[,4], .lbp)
        cbind((alphap+beta) / (alphap*(alpha+beta)),
              (alpha+betap) / (betap*(alpha+beta)))
    }, list(.la=la, .lap=lap, .lb=lb, .lbp=lbp))),
    last=eval(substitute(expression({
        misc$link = c("a"= .la, "ap"= .lap, "b"= .lb, "bp"= .lbp)
    }), list(.la=la, .lap=lap, .lb=lb, .lbp=lbp))),
    loglikelihood= eval(substitute(
            function(mu, y, w, residuals = FALSE, eta, extra=NULL) {
        alpha  = eta2theta(eta[,1], .la)
        alphap = eta2theta(eta[,2], .lap)
        beta   = eta2theta(eta[,3], .lb)
        betap  = eta2theta(eta[,4], .lbp)
        if(residuals) stop("loglikelihood residuals not implemented yet") else {
            tmp88 = extra$y1.lt.y2
            ell1 = log(alpha[tmp88]) + log(betap[tmp88]) -
                   betap[tmp88] * y[tmp88,2] -
                   (alpha+beta-betap)[tmp88] * y[tmp88,1]
            ell2 = log(beta[!tmp88]) + log(alphap[!tmp88]) -
                   alphap[!tmp88] * y[!tmp88,1] -
                   (alpha+beta-alphap)[!tmp88] * y[!tmp88,2]
        sum(w[tmp88] * ell1) + sum(w[!tmp88] * ell2) }
    }, list(.la=la, .lap=lap, .lb=lb, .lbp=lbp))),
    vfamily=c("freund61"),
    deriv=eval(substitute(expression({
        tmp88 = extra$y1.lt.y2
        alpha  = eta2theta(eta[,1], .la)
        alphap = eta2theta(eta[,2], .lap)
        beta   = eta2theta(eta[,3], .lb)
        betap  = eta2theta(eta[,4], .lbp)
        d1 = 1/alpha - y[,1]
        d1[!tmp88] = -y[!tmp88,2]
        d2 = 0 * alphap
        d2[!tmp88] = 1/alphap[!tmp88] - y[!tmp88,1] + y[!tmp88,2]
        d3 = -y[,1]
        d3[!tmp88] = 1/beta[!tmp88] - y[!tmp88,2]
        d4 = 1/betap - y[,2] + y[,1]
        d4[!tmp88] = 0
        w * cbind(d1 * dtheta.deta(alpha,  .la),
                  d2 * dtheta.deta(alphap, .lap),
                  d3 * dtheta.deta(beta,   .lb),
                  d4 * dtheta.deta(betap,  .lbp))
    }), list(.la=la, .lap=lap, .lb=lb, .lbp=lbp))),
    weight=eval(substitute(expression({
        py1.lt.y2 = alpha / (alpha+beta)
        d11 = py1.lt.y2 / alpha^2
        d22 = (1-py1.lt.y2) / alphap^2
        d33 = (1-py1.lt.y2) / beta^2
        d44 = py1.lt.y2 / betap^2
        wz = matrix(0, n, M)  # diagonal
        wz[,iam(1,1,M)] = dtheta.deta(alpha,  .la)^2 * d11
        wz[,iam(2,2,M)] = dtheta.deta(alphap, .lap)^2 * d22
        wz[,iam(3,3,M)] = dtheta.deta(beta,   .lb)^2 * d33
        wz[,iam(4,4,M)] = dtheta.deta(betap,  .lbp)^2 * d44
        w * wz
    }), list(.la=la, .lap=lap, .lb=lb, .lbp=lbp))))
}





mckaygamma2 = function(la="loge",
                       lp="loge",
                       lq="loge",
                       ia=NULL,
                       ip=1,
                       iq=1,
                       zero=NULL) {
    if(mode(la) != "character" && mode(la) != "name")
        la = as.character(substitute(la))
    if(mode(lp) != "character" && mode(lp) != "name")
        lp = as.character(substitute(lp))
    if(mode(lq) != "character" && mode(lq) != "name")
        lq = as.character(substitute(lq))
    if(!is.Numeric(ip, positive = TRUE) || !is.Numeric(iq, positive = TRUE))
        stop("initial values for ip and iq must be positive")
    if(is.Numeric(ia) && any(ia <= 0))
        stop("ia must be positive or NULL")

    new("vglmff",
    blurb=c("McKay's Bivariate Gamma Distribution\n",
           "Links:    ",
           namesof("a", la), ", ",
           namesof("p", lp), ", ",
           namesof("q", lq)),
    constraints=eval(substitute(expression({
        constraints = cm.zero.vgam(constraints, x, .zero, M)
    }), list(.zero=zero))),
    initialize=eval(substitute(expression({
        if(!is.matrix(y) || ncol(y) != 2)
            stop("the response must be a 2 column matrix") 
        sorty1 = pmin(y[,1], y[,2])
        sorty2 = pmax(y[,1], y[,2])
        if(any(sorty2-sorty1 <= 0))
            stop("Delete those observations that are identical")
        predictors.names = c(namesof("a", .la, short=TRUE), 
                      namesof("p", .lp, short=TRUE), 
                      namesof("q", .lq, short=TRUE))
        if(!length(etastart)) {
            pinit = rep(.ip, len=n)
            qinit = rep(.iq, len=n)
            # Computing ainit from pinit and qinit is easy
            ainit = if(length(.ia)) rep(.ia, len=n) else (pinit+qinit)/(sorty1+0.1)
            etastart = cbind(theta2eta(ainit, .la),
                             theta2eta(pinit, .lp),
                             theta2eta(qinit, .lq))
        }
    }), list(.la=la, .lp=lp, .lq=lq, .ia=ia, .ip=ip, .iq=iq))),
    inverse=eval(substitute(function(eta, extra=NULL) {
        a = eta2theta(eta[,1], .la)
        p = eta2theta(eta[,2], .lp)
        q = eta2theta(eta[,3], .lq)
        cbind("pmin(y1,y2)"=(p+q)/a, "pmax(y1,y2)"=NA)
    }, list(.la=la, .lp=lp, .lq=lq))),
    last=eval(substitute(expression({
        misc$link = c("a"= .la, "p"= .lp, "q"= .lq)
    }), list(.la=la, .lp=lp, .lq=lq))),
    loglikelihood= eval(substitute(
            function(mu, y, w, residuals = FALSE, eta, extra=NULL) {
        a = eta2theta(eta[,1], .la)
        p = eta2theta(eta[,2], .lp)
        q = eta2theta(eta[,3], .lq)
        y = cbind(pmin(y[,1], y[,2]), pmax(y[,1], y[,2])) # Sort so y[,1]<y[,2]
        # Note that, after sorting, y[,1] < y[,2] is needed:
        if(residuals) stop("loglikelihood residuals not implemented yet") else
        sum(w * ((p+q)*log(a) - lgamma(p) - lgamma(q) +
                 (p-1)*log(y[,1]) + (q-1)*log(y[,2]-y[,1]) - a*y[,2] ))
    }, list(.la=la, .lp=lp, .lq=lq))),
    vfamily=c("mckaygamma2"),
    deriv=eval(substitute(expression({
        a = eta2theta(eta[,1], .la)
        p = eta2theta(eta[,2], .lp)
        q = eta2theta(eta[,3], .lq)
        sorty = y
        sorty[,1] = pmin(y[,1], y[,2])
        sorty[,2] = pmax(y[,1], y[,2])
        d1 = (p+q)/a - sorty[,2]
        d2 = log(a) - digamma(p) + log(sorty[,1])
        d3 = log(a) - digamma(q) + log(sorty[,2]-sorty[,1])
        w * cbind(d1 * dtheta.deta(a, .la),
                  d2 * dtheta.deta(p, .lp),
                  d3 * dtheta.deta(q, .lq))
    }), list(.la=la, .lp=lp, .lq=lq))),
    weight=eval(substitute(expression({
        d11 = (p+q)/a^2
        d22 = trigamma(p)
        d33 = trigamma(q)
        d12 = -1/a
        d13 = -1/a
        d23 = 0
        wz = matrix(0, n, dimm(M))
        wz[,iam(1,1,M)] = dtheta.deta(a, .la)^2 * d11
        wz[,iam(2,2,M)] = dtheta.deta(p, .lp)^2 * d22
        wz[,iam(3,3,M)] = dtheta.deta(q, .lq)^2 * d33
        wz[,iam(1,2,M)] = dtheta.deta(a, .la) * dtheta.deta(p, .lp) * d12
        wz[,iam(1,3,M)] = dtheta.deta(a, .la) * dtheta.deta(q, .lq) * d13
        wz[,iam(2,3,M)] = dtheta.deta(p, .lp) * dtheta.deta(q, .lq) * d23
        w * wz
    }), list(.la=la, .lp=lp, .lq=lq))))
}



rfrank = function(n, alpha) {
    if(!is.Numeric(n, posit=TRUE, allow=1, integ=TRUE)) stop("bad input for n")
    if(!is.Numeric(alpha, posit=TRUE)) stop("bad input for \"alpha\"")
    alpha = rep(alpha, len=n)
    U = runif(n)
    V = runif(n)
    T = alpha^U + (alpha - alpha^U) * V
    X = U
    index = abs(alpha-1) < .Machine$double.eps
    Y = U
    if(any(!index))
        Y[!index] = logb(T[!index]/(T[!index]+(1-alpha[!index])*V[!index]),
                         base=alpha[!index])
    ans = matrix(c(X,Y), nrow=n, ncol=2) # Want to suppress column names
    if(any(index)) {
        ans[index,1] = runif(sum(index)) # Uniform density for alpha==1
        ans[index,2] = runif(sum(index))
    }
    ans
}

pfrank = function(q1, q2, alpha) {
    if(!is.Numeric(q1)) stop("bad input for \"q1\"")
    if(!is.Numeric(q2)) stop("bad input for \"q2\"")
    if(!is.Numeric(alpha, posit=TRUE)) stop("bad input for \"alpha\"")

    L = max(length(q1), length(q2), length(alpha))
    alpha = rep(alpha, len=L)
    q1 = rep(q1, len=L)
    q2 = rep(q2, len=L)

    x=q1; y=q2
    index = (x>=1 & y<1) | (y>=1 & x<1) | (x<=0 | y<=0) | (x>=1 & y>=1) |
            (abs(alpha-1) < .Machine$double.eps)
    ans = as.numeric(index)
    if(any(!index))
    ans[!index] = logb(1 + ((alpha[!index])^(x[!index])-1)*
                  ((alpha[!index])^(y[!index])-1)/(alpha[!index]-1), 
                  base=alpha[!index])
    ind2 = (abs(alpha-1) < .Machine$double.eps)
    ans[ind2] = x[ind2] * y[ind2]
    ans[x>=1 & y<1] = y[x>=1 & y<1]   # P(Y2 < q2) = q2
    ans[y>=1 & x<1] = x[y>=1 & x<1]   # P(Y1 < q1) = q1
    ans[x<=0 | y<=0] = 0
    ans[x>=1 & y>=1] = 1
    ans
}

dfrank = function(x1, x2, alpha) {
    if(!is.Numeric(x1)) stop("bad input for \"x1\"")
    if(!is.Numeric(x2)) stop("bad input for \"x2\"")
    if(!is.Numeric(alpha, posit=TRUE)) stop("bad input for \"alpha\"")

    L = max(length(x1), length(x2), length(alpha))
    alpha = rep(alpha, len=L)
    x1 = rep(x1, len=L)
    x2 = rep(x2, len=L)

    temp = (alpha-1) + (alpha^x1 - 1) * (alpha^x2 - 1)
    index = (abs(alpha-1) < .Machine$double.eps)
    ans = x1
    if(any(!index))
        ans[!index] = (alpha[!index]-1) * log(alpha[!index]) *
            (alpha[!index])^(x1[!index]+x2[!index]) / (temp[!index])^2
    ans[x1<=0 | x2<=0 | x1>=1 | x2>=1] = 0
    ans[index] = 1
    ans
}




frank = function(lapar="loge", eapar=list(), iapar=2) {
    if(mode(lapar) != "character" && mode(lapar) != "name")
        lapar = as.character(substitute(lapar))
    if(!is.Numeric(iapar, positive = TRUE))
        stop("\"iapar\" must be positive")
    if(!is.list(eapar)) eapar = list()

    new("vglmff",
    blurb=c("Frank's Bivariate Distribution\n",
           "Links:    ",
           namesof("apar", lapar, earg= eapar )),
    initialize=eval(substitute(expression({
        if(!is.matrix(y) || ncol(y) != 2)
            stop("the response must be a 2 column matrix") 
        if(any(y <= 0) || any(y >= 1))
            stop("the response must have values between 0 and 1") 
        predictors.names = c(namesof("apar", .lapar, earg= .eapar, short=TRUE))
        if(!length(etastart)) {
            apar.init = rep(.iapar, len=n)
            etastart = cbind(theta2eta(apar.init, .lapar, earg= .eapar ))
        }
    }), list( .lapar=lapar, .eapar=eapar, .iapar=iapar))),
    inverse=eval(substitute(function(eta, extra=NULL) {
        apar = eta2theta(eta, .lapar, earg= .eapar )
        cbind(rep(0.5, len=length(eta)), rep(0.5, len=length(eta)))
    }, list(.lapar=lapar, .eapar=eapar ))),
    last=eval(substitute(expression({
        misc$link = c("apar"= .lapar)
        misc$earg = list("apar"= .eapar )
        misc$pooled.weight = pooled.weight
    }), list(.lapar=lapar, .eapar=eapar ))),
    loglikelihood= eval(substitute(
            function(mu, y, w, residuals = FALSE, eta, extra=NULL) {
        apar = eta2theta(eta, .lapar, earg= .eapar )
        if(residuals) stop("loglikelihood residuals not implemented yet") else {
            denom = apar-1 + (apar^y[,1] -1) * (apar^y[,2] -1)
            denom = abs(denom)  # Needed; Genest (1987) uses this too, eqn (4.1)
            sum(w * (log((apar-1) * log(apar)) + (y[,1]+y[,2])*log(apar) -
                    2 * log(denom)))
        }
    }, list(.lapar=lapar, .eapar=eapar ))),
    vfamily=c("frank"),
    deriv=eval(substitute(expression({
        apar = eta2theta(eta, .lapar, earg= .eapar )
        denom = apar-1 + (apar^y[,1] -1) * (apar^y[,2] -1)
        tmp700 = 2*apar^(y[,1]+y[,2]) - apar^y[,1] - apar^y[,2]
        numerator = 1 + y[,1] * apar^(y[,1]-1) * (apar^y[,2] -1) + 
                        y[,2] * apar^(y[,2]-1) * (apar^y[,1] -1)
        Dl.dapar = 1/(apar-1) + 1/(apar*log(apar)) + (y[,1]+y[,2])/apar -
                   2 * numerator / denom
        dapar.deta = dtheta.deta(apar, .lapar, earg= .eapar )

        w * Dl.dapar * dapar.deta
    }), list(.lapar=lapar, .eapar=eapar ))),
    weight=eval(substitute(expression({
        nump = apar^(y[,1]+y[,2]-2) * (2 * y[,1] * y[,2] +
                     y[,1]*(y[,1]-1) + y[,2]*(y[,2]-1)) - 
                     y[,1]*(y[,1]-1) * apar^(y[,1]-2) - 
                     y[,2]*(y[,2]-1) * apar^(y[,2]-2)
        D2l.dapar2 = 1/(apar-1)^2 + (1+log(apar))/(apar*log(apar))^2 +
                     (y[,1]+y[,2])/apar^2 + 2 *
                     (nump / denom - (numerator/denom)^2)
        d2apar.deta2 = d2theta.deta2(apar, .lapar)
        wz = w * (dapar.deta^2 * D2l.dapar2 - Dl.dapar * d2apar.deta2)

        if(TRUE && intercept.only) {
            wz = cbind(wz)
            sumw = sum(w)
            for(iii in 1:ncol(wz))
                wz[,iii] = sum(wz[,iii]) / sumw
            pooled.weight = TRUE
            wz = w * wz   # Put back the weights
        } else
            pooled.weight = FALSE

        wz
    }), list( .lapar=lapar, .eapar=eapar ))))
}



gammahyp = function(ltheta="loge", itheta=NULL, expected=FALSE) {
    if(mode(ltheta) != "character" && mode(ltheta) != "name")
        ltheta = as.character(substitute(ltheta))
    if(!is.logical(expected) || length(expected)!=1)
        stop("\"expected\" must be a single logical")

    new("vglmff",
    blurb=c("Gamma Hyperbola Bivariate Distribution\n",
           "Links:    ",
           namesof("theta", ltheta)),
    initialize=eval(substitute(expression({
        if(!is.matrix(y) || ncol(y) != 2)
            stop("the response must be a 2 column matrix") 
        if(any(y[,1] <= 0) || any(y[,2] <= 1))
            stop("the response has values that are out of range") 
        predictors.names = c(namesof("theta", .ltheta, short=TRUE))
        if(!length(etastart)) {
            theta.init = if(length( .itheta)) rep(.itheta, len=n) else {
                1 / (y[,2] - 1 + 0.01)
            }
            etastart = cbind(theta2eta(theta.init, .ltheta))
        }
    }), list(.ltheta=ltheta, .itheta=itheta))),
    inverse=eval(substitute(function(eta, extra=NULL) {
        theta = eta2theta(eta, .ltheta)
        cbind(theta*exp(theta), 1+1/theta)
    }, list(.ltheta=ltheta))),
    last=eval(substitute(expression({
        misc$link = c("theta"= .ltheta)
        misc$expected = .expected 
    }), list(.ltheta=ltheta, .expected=expected))),
    loglikelihood= eval(substitute(
            function(mu, y, w, residuals = FALSE, eta, extra=NULL) {
        theta = eta2theta(eta, .ltheta)
        if(residuals) stop("loglikelihood residuals not implemented yet") else {
            sum(w * (-exp(-theta)*y[,1]/theta - theta*y[,2]))
        }
    }, list(.ltheta=ltheta))),
    vfamily=c("gammahyp"),
    deriv=eval(substitute(expression({
        theta = eta2theta(eta, .ltheta)
        Dl.dtheta = exp(-theta) * y[,1] * (1+theta) / theta^2 - y[,2]
        Dtheta.deta = dtheta.deta(theta, .ltheta)
        w * Dl.dtheta * Dtheta.deta
    }), list(.ltheta=ltheta))),
    weight=eval(substitute(expression({
        temp300 = 2 + theta * (2 + theta)
        if( .expected) {
            D2l.dtheta2 = temp300 / theta^2
            wz = w * Dtheta.deta^2 * D2l.dtheta2
        } else {
            D2l.dtheta2 = temp300 * y[,1] * exp(-theta) / theta^3
            D2theta.deta2 = d2theta.deta2(theta, .ltheta)
            wz = w * (Dtheta.deta^2 * D2l.dtheta2 - Dl.dtheta * D2theta.deta2)
        }
        wz
    }), list( .expected=expected, .ltheta=ltheta))))
}



morgenstern = function(lapar="rhobit", earg=list(), iapar=NULL, tola0=0.01,
                       method.init=1) {
    if(mode(lapar) != "character" && mode(lapar) != "name")
        lapar = as.character(substitute(lapar))
    if(!is.list(earg)) earg = list()
    if(length(iapar) && (!is.Numeric(iapar, allow=1) || abs(iapar) >= 1))
        stop("'iapar' must be a single number between -1 and 1") 
    if(!is.Numeric(tola0, allow=1, posit=TRUE))
        stop("'tola0' must be a single positive number") 
    if(length(iapar) && abs(iapar) <= tola0)
        stop("'iapar' must not be between -tola0 and tola0") 
    if(!is.Numeric(method.init, allow=1, integ=TRUE, positi=TRUE) ||
       method.init > 2.5)
        stop("argument \"method.init\" must be 1 or 2")

    new("vglmff",
    blurb=c("Morgenstern's Bivariate Exponential Distribution\n",
           "Links:    ",
           namesof("apar", lapar, earg= earg )),
    initialize=eval(substitute(expression({
        if(!is.matrix(y) || ncol(y) != 2)
            stop("the response must be a 2 column matrix") 
        if(any(y < 0))
            stop("the response must have non-negative values only") 
        predictors.names = c(namesof("apar", .lapar, earg= .earg , short=TRUE))
        if(!length(etastart)) {
            ainit  = if(length(.iapar))  rep(.iapar, len=n) else {
                mean1 = if( .method.init == 1) median(y[,1]) else mean(y[,1])
                mean2 = if( .method.init == 1) median(y[,2]) else mean(y[,2])
                Finit = 0.01 + mean(y[,1] <= mean1 & y[,2] <= mean2)
                ((Finit-1+exp(-mean1)+exp(-mean2)) / exp(-mean1-mean2)  -
                 1) / ((1-exp(-mean1)) * (1-exp(-mean2)))
              }
            etastart = theta2eta(rep(ainit, len=n), .lapar, earg= .earg )
        }
    }), list( .iapar=iapar, .lapar=lapar, .earg=earg,
              .method.init=method.init ))),
    inverse=eval(substitute(function(eta, extra=NULL) {
        alpha = eta2theta(eta, .lapar, earg= .earg )
        cbind(rep(1, len=length(alpha)),
              rep(1, len=length(alpha)))
    }, list( .lapar=lapar, .earg=earg ))),
    last=eval(substitute(expression({
        misc$link = c("apar"= .lapar)
        misc$earg = list(apar = .earg)
        misc$expected = FALSE
        misc$pooled.weight = pooled.weight
    }), list( .lapar=lapar, .earg=earg ))),
    loglikelihood= eval(substitute(
            function(mu, y, w, residuals = FALSE, eta, extra=NULL) {
        alpha  = eta2theta(eta, .lapar, earg= .earg )
        alpha[abs(alpha) < .tola0 ] = .tola0
        if(residuals) stop("loglikelihood residuals not implemented yet") else {
        denom = (1 + alpha - 2*alpha*(exp(-y[,1]) + exp(-y[,2])) +
                4*alpha*exp(-y[,1] - y[,2]))
        sum(w * (-y[,1] - y[,2] + log(denom)))
        }
    }, list( .lapar=lapar, .earg=earg, .tola0=tola0 ))),
    vfamily=c("morgenstern"),
    deriv=eval(substitute(expression({
        alpha  = eta2theta(eta, .lapar, earg= .earg )
        alpha[abs(alpha) < .tola0 ] = .tola0
        numerator = 1 - 2*(exp(-y[,1]) + exp(-y[,2])) + 4*exp(-y[,1] - y[,2])
        denom = (1 + alpha - 2*alpha*(exp(-y[,1]) + exp(-y[,2])) +
                4 *alpha*exp(-y[,1] - y[,2]))
        dl.dalpha = numerator / denom
        dalpha.deta = dtheta.deta(alpha,  .lapar, earg= .earg )
        w * cbind(dl.dalpha * dalpha.deta)
    }), list( .lapar=lapar, .earg=earg, .tola0=tola0 ))),
    weight=eval(substitute(expression({
        d2l.dalpha2 = dl.dalpha^2
        d2alpha.deta2 = d2theta.deta2(alpha,  .lapar, earg= .earg )
        wz = w * (dalpha.deta^2 * d2l.dalpha2 - d2alpha.deta2 * dl.dalpha)
        if(TRUE &&
           intercept.only) {
            wz = cbind(wz)
            sumw = sum(w)
            for(iii in 1:ncol(wz))
                wz[,iii] = sum(wz[,iii]) / sumw
            pooled.weight = TRUE
            wz = w * wz   # Put back the weights
        } else
            pooled.weight = FALSE
        wz
    }), list( .lapar=lapar, .earg=earg ))))
}




dfgm = function(x1, x2, alpha) {
    if(!is.Numeric(alpha)) stop("bad input for \"alpha\"")
    if(any(alpha < -1 | alpha > 1)) stop("\"alpha\" values out of range")
    L = max(length(x1), length(x2), length(alpha))
    if(length(x1) != L)  x1 = rep(x1, len=L)
    if(length(x2) != L)  x2 = rep(x2, len=L)
    if(length(alpha) != L)  alpha = rep(alpha, len=L)
    ans = 1 + alpha * (1-2*x1) * (1-2*x2)
    ans[(x1 <= 0) | (x1 >= 1) | (x2 <= 0) | (x2 >= 1)] = 0
    if(any(ans<0))
        stop("negative values in the density (alpha out of range)") else
    ans
}



fgm = function(lapar="identity", earg=list(), iapar=NULL,
               method.init=1) { # , tola0=0.01
    if(mode(lapar) != "character" && mode(lapar) != "name")
        lapar = as.character(substitute(lapar))
    if(!is.list(earg)) earg = list()
    if(length(iapar) && !is.Numeric(iapar, allow=1))
        stop("'iapar' must be a single number")
    if(!is.Numeric(method.init, allow=1, integ=TRUE, positi=TRUE) ||
       method.init > 2.5)
        stop("argument \"method.init\" must be 1 or 2")

    new("vglmff",
    blurb=c("Farlie-Gumbel-Morgenstern Distribution\n",
           "Links:    ",
           namesof("apar", lapar, earg= earg )),
    initialize=eval(substitute(expression({
        if(!is.matrix(y) || ncol(y) != 2)
            stop("the response must be a 2 column matrix") 
        if(any(y < 0) || any(y > 1))
            stop("the response must have values in the unit square")
        predictors.names = namesof("apar", .lapar, earg= .earg, short=TRUE)
        if(!length(etastart)) {
            ainit  = if(length( .iapar ))  rep( .iapar, len=n) else {
                mean1 = if( .method.init == 1) median(y[,1]) else mean(y[,1])
                mean2 = if( .method.init == 1) median(y[,2]) else mean(y[,2])
                Finit = 0.01 + mean(y[,1] <= mean1 & y[,2] <= mean2)
                (Finit / (mean1 * mean2) - 1) / ((1-mean1) * (1-mean2))
            } 
            etastart = theta2eta(rep(ainit, len=n), .lapar, earg= .earg )
        }
    }), list( .iapar=iapar, .lapar=lapar, .earg=earg,
              .method.init=method.init ))),
    inverse=eval(substitute(function(eta, extra=NULL) {
        alpha = eta2theta(eta, .lapar, earg= .earg )
        cbind(rep(0.5, len=length(alpha)),
              rep(0.5, len=length(alpha)))
    }, list( .lapar=lapar, .earg=earg ))),
    last=eval(substitute(expression({
        misc$link = c("apar"= .lapar)
        misc$earg = list(apar = .earg)
        misc$expected = FALSE
        misc$pooled.weight = pooled.weight
    }), list( .lapar=lapar, .earg=earg ))),
    loglikelihood= eval(substitute(
            function(mu, y, w, residuals = FALSE, eta, extra=NULL) {
        alpha = eta2theta(eta, .lapar, earg= .earg )
        if(residuals) stop("loglikelihood residuals not implemented yet") else {
            denom = 1 + alpha * (1 - 2 * y[,1])  * (1 - 2 * y[,2])
            mytolerance = .Machine$double.eps
            bad <- (denom <= mytolerance)   # Range violation
            if(any(bad)) {
                cat("There are some range violations in @loglikelihood\n")
                if(exists("flush.console")) flush.console()
            }
            sum(bad) * (-1.0e10) + 
            sum(w[!bad] * log(denom[!bad]))
        }
    }, list( .lapar=lapar, .earg=earg ))),
    vfamily=c("fgm"),
    deriv=eval(substitute(expression({
        alpha  = eta2theta(eta, .lapar, earg= .earg )
        numerator = (1 - 2 * y[,1])  * (1 - 2 * y[,2])
        denom = 1 + alpha * numerator
            mytolerance = .Machine$double.eps
            bad <- (denom <= mytolerance)   # Range violation
            if(any(bad)) {
                cat("There are some range violations in @deriv\n")
                if(exists("flush.console")) flush.console()
                denom[bad] = 2 * mytolerance
            }
        dl.dalpha = numerator / denom
        dalpha.deta = dtheta.deta(alpha, .lapar, earg= .earg )
        w * cbind(dl.dalpha * dalpha.deta)
    }), list( .lapar=lapar, .earg=earg ))),
    weight=eval(substitute(expression({
        d2l.dalpha2 = dl.dalpha^2
        d2alpha.deta2 = d2theta.deta2(alpha, .lapar, earg= .earg )
        wz = w * (dalpha.deta^2 * d2l.dalpha2 - d2alpha.deta2 * dl.dalpha)
        if(TRUE &&
           intercept.only) {
            wz = cbind(wz)
            sumw = sum(w)
            for(iii in 1:ncol(wz))
                wz[,iii] = sum(wz[,iii]) / sumw
            pooled.weight = TRUE
            wz = w * wz   # Put back the weights
        } else
            pooled.weight = FALSE
        wz
    }), list( .lapar=lapar, .earg=earg ))))
}



gumbelIbiv = function(lapar="identity", earg=list(), iapar=NULL, method.init=1) {
    if(mode(lapar) != "character" && mode(lapar) != "name")
        lapar = as.character(substitute(lapar))
    if(!is.list(earg)) earg = list()
    if(length(iapar) && !is.Numeric(iapar, allow=1))
        stop("'iapar' must be a single number")
    if(!is.Numeric(method.init, allow=1, integ=TRUE, positi=TRUE) ||
       method.init > 2.5)
        stop("argument \"method.init\" must be 1 or 2")

    new("vglmff",
    blurb=c("Gumbel's Type I Bivariate Distribution\n",
           "Links:    ",
           namesof("apar", lapar, earg= earg )),
    initialize=eval(substitute(expression({
        if(!is.matrix(y) || ncol(y) != 2)
            stop("the response must be a 2 column matrix") 
        if(any(y < 0))
            stop("the response must have non-negative values only")
        predictors.names = c(namesof("apar", .lapar, earg= .earg , short=TRUE))
        if(!length(etastart)) {
            ainit  = if(length( .iapar ))  rep( .iapar, len=n) else {
                mean1 = if( .method.init == 1) median(y[,1]) else mean(y[,1])
                mean2 = if( .method.init == 1) median(y[,2]) else mean(y[,2])
                Finit = 0.01 + mean(y[,1] <= mean1 & y[,2] <= mean2)
                (log(Finit-1+exp(-mean1)+exp(-mean2))+mean1+mean2)/(mean1*mean2)
            }
            etastart = theta2eta(rep(ainit,  len=n), .lapar, earg= .earg )
        }
    }), list( .iapar=iapar, .lapar=lapar, .earg=earg,
              .method.init=method.init ))),
    inverse=eval(substitute(function(eta, extra=NULL) {
        alpha = eta2theta(eta, .lapar, earg= .earg )
        cbind(rep(1, len=length(alpha)),
              rep(1, len=length(alpha)))
    }, list( .lapar=lapar ))),
    last=eval(substitute(expression({
        misc$link = c("apar"= .lapar)
        misc$earg = list(apar = .earg)
        misc$expected = FALSE
        misc$pooled.weight = pooled.weight
    }), list( .lapar=lapar, .earg=earg ))),
    loglikelihood= eval(substitute(
            function(mu, y, w, residuals = FALSE, eta, extra=NULL) {
        alpha  = eta2theta(eta, .lapar, earg= .earg )
        if(residuals) stop("loglikelihood residuals not implemented yet") else {
            denom = (alpha*y[,1] - 1) * (alpha*y[,2] - 1) + alpha
            mytolerance = .Machine$double.xmin
            bad <- (denom <= mytolerance)   # Range violation
            if(any(bad)) {
                cat("There are some range violations in @deriv\n")
                if(exists("flush.console")) flush.console()
            }
            sum(bad) * (-1.0e10) + 
            sum(w[!bad] * (-y[!bad,1] - y[!bad,2] +
                alpha[!bad]*y[!bad,1]*y[!bad,2] + log(denom[!bad])))
        }
    }, list( .lapar=lapar, .earg=earg ))),
    vfamily=c("gumbelIbiv"),
    deriv=eval(substitute(expression({
        alpha  = eta2theta(eta, .lapar, earg= .earg )
        numerator = (alpha*y[,1] - 1)*y[,2] + (alpha*y[,2] - 1)*y[,1] + 1
        denom = (alpha*y[,1] - 1) * (alpha*y[,2] - 1) + alpha
        denom = abs(denom)
        dl.dalpha = numerator / denom + y[,1]*y[,2]
        dalpha.deta = dtheta.deta(alpha,  .lapar, earg= .earg )
        w * cbind(dl.dalpha * dalpha.deta)
    }), list( .lapar=lapar, .earg=earg ))),
    weight=eval(substitute(expression({
        d2l.dalpha2 = (numerator/denom)^2 - 2*y[,1]*y[,2] / denom
        d2alpha.deta2 = d2theta.deta2(alpha, .lapar, earg= .earg )
        wz = w * (dalpha.deta^2 * d2l.dalpha2 - d2alpha.deta2 * dl.dalpha)
        if(TRUE &&
           intercept.only) {
            wz = cbind(wz)
            sumw = sum(w)
            for(iii in 1:ncol(wz))
                wz[,iii] = sum(wz[,iii]) / sumw
            pooled.weight = TRUE
            wz = w * wz   # Put back the weights
        } else
            pooled.weight = FALSE
        wz
    }), list( .lapar=lapar, .earg=earg ))))
}



# These functions are
# Copyright (C) 1998-2007 T.W. Yee, University of Auckland. All rights reserved.
















getMaxMin = function(vov, objfun, y, x, w, extraargs=NULL, maximize=TRUE,
                     abs.arg=FALSE) {
    if(!is.vector(vov)) stop("vov must be a vector")
    objvals = vov
    for(ii in 1:length(vov))
        objvals[ii] = objfun(vov[ii], y=y, x=x, w=w, extraargs=extraargs)
    try.this = if(abs.arg) {
                   if(maximize) vov[abs(objvals) == max(abs(objvals))] else
                   vov[abs(objvals) == min(abs(objvals))]
               } else {
                   if(maximize) vov[objvals == max(objvals)] else
                   vov[objvals == min(objvals)]
               }
    if(!length(try.this)) stop("something has gone wrong!")
    if(length(try.this) == 1) try.this else sample(try.this, size=1)
}



mccullagh89 = function(ltheta="rhobit", lnu="logoff",
                       itheta=NULL, inu=NULL,
                       etheta=list(),
                       enu=if(lnu == "logoff") list(offset=0.5) else list(),
                       zero=NULL)
{
    if(mode(ltheta) != "character" && mode(ltheta) != "name")
        ltheta = as.character(substitute(ltheta))
    if(mode(lnu) != "character" && mode(lnu) != "name")
        lnu = as.character(substitute(lnu))
    if(length(zero) && !is.Numeric(zero, integer=TRUE, posit=TRUE))
        stop("bad input for argument \"zero\"")
    if(!is.list(etheta)) etheta = list()
    if(!is.list(enu)) enu = list()

    new("vglmff",
    blurb=c("McCullagh (1989)'s distribution \n",
    "f(y) = (1-2*theta*y+theta^2)^(-nu) * [1 - y^2]^(nu-1/2) /\n",
            "       Beta[nu+1/2, 1/2], ",
            "  -1 < y < 1, -1 < theta < 1, nu > -1/2\n",
            "Links:     ",
            namesof("theta", ltheta, earg=etheta), ", ",
            namesof("nu", lnu, earg=enu),
            "\n",
            "\n",
            "Mean:     nu*theta/(1+nu)"),
    constraints=eval(substitute(expression({
        constraints = cm.zero.vgam(constraints, x, .zero, M)
    }), list( .zero=zero ))),
    initialize=eval(substitute(expression({
        y = as.numeric(y)
        if(any(y <= -1 | y >= 1))
            stop("all y values must be in (-1,1)")
        predictors.names= c(namesof("theta", .ltheta, earg= .etheta,tag=FALSE),
                            namesof("nu",    .lnu,    earg= .enu,   tag=FALSE))
        if(!length(etastart)) {
            theta.init = if(length(.itheta)) rep(.itheta, length=n) else {
                mccullagh89.aux = function(thetaval, y, x, w, extraargs)
                mean((y-thetaval)*(thetaval^2-1)/(1-2*thetaval*y+thetaval^2))
                theta.grid = seq(-0.9, 0.9, by=0.05)
                try.this = getMaxMin(theta.grid, objfun=mccullagh89.aux,
                                     y=y,  x=x, w=w, maximize=FALSE,
                                     abs.arg=TRUE)
                try.this = rep(try.this, len=n)
                try.this
            }
            tmp = y / (theta.init-y)
            tmp[tmp < -0.4] = -0.4
            tmp[tmp > 10.0] = 10.0
            nu.init = rep(if(length(.inu)) .inu else tmp, length=n)
            nu.init[!is.finite(nu.init)] = 0.4
            etastart = cbind(theta2eta(theta.init, .ltheta, earg=.etheta ),
                             theta2eta(nu.init, .lnu, earg= .enu ))
        }
    }), list( .ltheta=ltheta, .lnu=lnu, .inu=inu, .itheta=itheta,
              .etheta = etheta, .enu=enu ))),
    inverse=eval(substitute(function(eta, extra=NULL) {
        Theta = eta2theta(eta[,1], .ltheta, earg= .etheta )
        nu = eta2theta(eta[,2], .lnu, earg= .enu )
        nu*Theta/(1+nu)
    }, list( .ltheta=ltheta, .lnu=lnu,
             .etheta = etheta, .enu=enu ))),
    last=eval(substitute(expression({
        misc$link = c("theta"= .ltheta, "nu"= .lnu)
        misc$earg =  list(theta = .etheta, nu= .enu )
    }), list( .ltheta=ltheta, .lnu=lnu, .etheta = etheta, .enu=enu ))),
    loglikelihood=eval(substitute(
        function(mu,y,w,residuals= FALSE,eta, extra=NULL) {
        Theta = eta2theta(eta[,1], .ltheta, earg= .etheta )
        nu = eta2theta(eta[,2], .lnu, earg= .enu )
        if(residuals) stop("loglikelihood residuals not implemented yet") else
        sum(w * ((nu-0.5)*log1p(-y^2) - nu * log1p(-2*Theta*y + Theta^2) -
                lbeta(nu+0.5,0.5 )))
    }, list( .ltheta=ltheta, .lnu=lnu, .etheta = etheta, .enu=enu ))),
    vfamily=c("mccullagh89"),
    deriv=eval(substitute(expression({
        Theta = eta2theta(eta[,1], .ltheta, earg= .etheta )
        nu = eta2theta(eta[,2], .lnu, earg= .enu )
        dTheta.deta = dtheta.deta(Theta, .ltheta, earg= .etheta )
        dnu.deta = dtheta.deta(nu, .lnu, earg= .enu )
        dl.dTheta = 2 * nu * (y-Theta) / (1 -2*Theta*y + Theta^2)
        dl.dnu = log1p(-y^2) - log1p(-2*Theta*y + Theta^2) -
                 digamma(nu+0.5) + digamma(nu+1)
        w * cbind(dl.dTheta * dTheta.deta, dl.dnu * dnu.deta)
    }), list( .ltheta=ltheta, .lnu=lnu, .etheta = etheta, .enu=enu ))),
    weight=eval(substitute(expression({
        d2l.dTheta2 = (2 * nu^2 / (1+nu)) / (1-Theta^2)
        d2l.dnu2 = trigamma(nu+0.5) - trigamma(nu+1)
        wz = matrix(as.numeric(NA), n, M)  #diagonal matrix
        wz[,iam(1,1,M)] = d2l.dTheta2 * dTheta.deta^2
        wz[,iam(2,2,M)] = d2l.dnu2 * dnu.deta^2
        w * wz
    }), list( .ltheta=ltheta, .lnu=lnu ))))
}





hzeta = function(link="loglog", earg=list(), init.alpha=NULL)
{
    if(length(init.alpha) && !is.Numeric(init.alpha, positive=TRUE))
        stop("'init.alpha' must be > 0")

    if(mode(link) != "character" && mode(link) != "name")
        link = as.character(substitute(link))
    if(!is.list(earg)) earg = list()

    new("vglmff",
    blurb=c(
    "Haight's Zeta distribution f(y) = (2y-1)^(-alpha) - (2y+1)^(-alpha),\n",
            "    alpha>0, y=1,2,..\n\n",
            "Link:    ",
            namesof("alpha", link, earg=earg), "\n\n",
            "Mean:     (1-2^(-alpha)) * zeta(alpha) if alpha>1",
            "\n",
            "Variance: (1-2^(1-alpha)) * zeta(alpha-1) - mean^2 if alpha>2"),
    initialize=eval(substitute(expression({
        y = as.numeric(y)
        if(ncol(cbind(y)) != 1)
            stop("response must be a vector or a one-column matrix")
        if(any(y < 1))
            stop("all y values must be in 1,2,3,...")
        predictors.names = namesof("alpha", .link, earg= .earg, tag=FALSE)
        if(!length(etastart)) {
            ainit = if(length( .init.alpha)) .init.alpha else {
                if((meany <- mean(y)) < 1.5) 3.0 else
                if(meany < 2.5) 1.4 else 1.1 
            }
            ainit = rep(ainit, length=n) 
            etastart = theta2eta(ainit, .link, earg= .earg )
        }
    }), list( .link=link, .earg=earg, .init.alpha=init.alpha ))),
    inverse=eval(substitute(function(eta, extra=NULL) {
        alpha = eta2theta(eta, .link, earg= .earg )
        mu = (1-2^(-alpha)) * zeta(alpha)
        mu[alpha <= 1] = Inf
        mu
    }, list( .link=link, .earg=earg ))),
    last=eval(substitute(expression({
        misc$d3 = d3    # because save.weights=F
        misc$link = c(alpha= .link)
        misc$earg = list(alpha= .earg)
        misc$pooled.weight = pooled.weight
    }), list( .link=link, .earg=earg ))),
    loglikelihood=eval(substitute(
        function(mu,y,w,residuals= FALSE,eta, extra=NULL) {
        alpha = eta2theta(eta, .link, earg= .earg )
        if(residuals) stop("loglikelihood residuals not implemented yet") else
        sum(w * log((2*y-1)^(-alpha) - (2*y+1)^(-alpha )))
    }, list( .link=link, .earg=earg ))),
    vfamily=c("hzeta"),
    deriv=eval(substitute(expression({
        if(iter==1) {
            d3 = deriv3(~ w * log((2*y-1)^(-alpha) - (2*y+1)^(-alpha)),
                        "alpha", hessian= TRUE)
        }

        alpha = eta2theta(eta, .link, earg= .earg ) 
        eval.d3 = eval(d3)
        dl.dalpha =  attr(eval.d3, "gradient")
        dalpha.deta = dtheta.deta(alpha, .link, earg= .earg )
        dl.dalpha * dalpha.deta
    }), list( .link=link, .earg=earg ))),
    weight=eval(substitute(expression({
        d2l.dalpha2 =  as.vector(attr(eval.d3, "hessian"))
        wz = -dalpha.deta^2 * d2l.dalpha2  -
              dl.dalpha * d2theta.deta2(alpha, .link, earg= .earg )

        if(FALSE && intercept.only) {
            sumw = sum(w)
            for(i in 1:ncol(wz))
                wz[,i] = sum(wz[,i]) / sumw
            pooled.weight = TRUE
            wz = w * wz   # Put back the weights
        } else
            pooled.weight = FALSE
       c(wz)
    }), list( .link=link, .earg=earg ))))
}



dhzeta = function(x, alpha) 
{
    if(!is.Numeric(alpha, posit=TRUE))
        stop("'alpha' must be numeric and have positive values")
    nn = max(length(x), length(alpha))
    x = rep(x, len=nn)
    alpha = rep(alpha, len=nn)
    ox = !is.finite(x)
    zero = ox | round(x) != x | x < 1
    ans = 0 * x
    ans[!zero] = (2*x[!zero]-1)^(-alpha[!zero]) - (2*x[!zero]+1)^(-alpha[!zero])
    ans
}


phzeta = function(q, alpha) 
{
    if(!is.Numeric(alpha, posit=TRUE))
        stop("'alpha' must be numeric and have positive values")
    nn = max(length(q), length(alpha))
    q = rep(q, len=nn)
    alpha = rep(alpha, len=nn)
    oq = !is.finite(q)
    zero = oq | q < 1
    q = floor(q)
    ans = 0 * q
    ans[!zero] = 1 - (2*q[!zero]+1)^(-alpha[!zero])
    ans
}


qhzeta = function(p, alpha) 
{
    if(!is.Numeric(alpha, posit=TRUE))
        stop("'alpha' must be numeric and have positive values")
    if(!is.Numeric(p, posit=TRUE) || any(p >= 1))
        stop("argument \"p\" must have values inside the interval (0,1)")
    nn = max(length(p), length(alpha))
    p = rep(p, len=nn)
    alpha = rep(alpha, len=nn)
    ans = (((1 - p)^(-1/alpha) - 1) / 2) # p is in (0,1)
    floor(ans+1)
}

rhzeta = function(n, alpha) 
{
    if(!is.Numeric(alpha, posit=TRUE))
        stop("'alpha' must be numeric and have positive values")
    if(!is.Numeric(n, posit=TRUE, integ=TRUE, allow=1))
        stop("argument \"n\" must be a positive integer")
    ans = ((runif(n)^(-1/alpha) - 1) / 2)
    floor(ans+1)
}


dirmultinomial = function(lphi="logit", ephi = list(),
                          iphi = 0.10, parallel= FALSE, zero="M")
{

    if(mode(lphi) != "character" && mode(lphi) != "name")
        lphi = as.character(substitute(lphi))
    if(length(zero) && 
       !(is.Numeric(zero, integer=TRUE, posit=TRUE) || is.character(zero )))
        stop("bad input for argument \"zero\"")
    if(!is.Numeric(iphi, positive=TRUE) || max(iphi) >= 1.0)
        stop("bad input for argument \"iphi\"")
    if(!is.list(ephi)) ephi = list()

    new("vglmff",
    blurb=c("Dirichlet-multinomial distribution\n\n",
            "Links:    ",
            "log(prob[1]/prob[M]), ..., log(prob[M-1]/prob[M]), ",
            namesof("phi", lphi, earg=ephi), "\n", "\n",
            "Mean:     shape_j / sum_j(shape_j)"),
    constraints=eval(substitute(expression({
        .ZERO = .zero
        if(is.character(.ZERO)) .ZERO = eval(parse(text = .ZERO))
        .PARALLEL = .parallel
        if(is.logical(.PARALLEL) && .PARALLEL) {
            mycmatrix = if(length(.ZERO))
                stop("can only handle parallel=TRUE when zero=NULL") else
                cbind(rbind(matrix(1,M-1,1), 0), rbind(matrix(0,M-1,1), 1))
        } else
            mycmatrix = if(M==1) diag(1) else diag(M)
        constraints=cm.vgam(mycmatrix, x, .PARALLEL, constraints, int=TRUE)
        constraints = cm.zero.vgam(constraints, x, .ZERO, M)
    }), list( .parallel=parallel, .zero=zero ))),
    initialize=eval(substitute(expression({
        delete.zero.colns <- TRUE
        eval(process.categorical.data.vgam)

        y = as.matrix(y)
        ycount = as.matrix(y * w)
        M = ncol(y)
        if(max(abs(ycount - round(ycount ))) > 1.0e-6)
            warning("there appears to be non-integer responses")
        if(min(ycount) < 0)
            stop("all values of the response (matrix) must be non-negative")
        predictors.names =
            c(paste("log(prob[,",1:(M-1),"]/prob[,",M,"])", sep=""),
              namesof("phi", .lphi, short=TRUE))
        extra$n2 = w  # aka omega, must be integer # as.vector(apply(y, 1, sum))
        if(!length(etastart)) {
            prob.init = apply(ycount, 2, sum)
            prob.init = prob.init / sum(prob.init)
            prob.init = matrix(prob.init, n, M, byrow=TRUE)
            phi.init = rep( .iphi, len=n)
            etastart = cbind(log(prob.init[,-M]/prob.init[,M]),
                             theta2eta(phi.init, .lphi, earg= .ephi ))
        }
    }), list( .lphi=lphi, .ephi=ephi, .iphi=iphi ))),
    inverse=eval(substitute(function(eta, extra=NULL) {
        M = if(is.matrix(eta)) ncol(eta) else 1
        temp = cbind(exp(eta[,-M]), 1)
        temp / as.vector(temp %*% rep(1, M))
    }, list( .ephi=ephi, .lphi=lphi ))),
    last=eval(substitute(expression({
        misc$link = c(rep("noLinkFunction", length=M-1), .lphi)
        names(misc$link) = c(paste("prob", 1:(M-1), sep=""), "phi")
        misc$earg = vector("list", M)
        names(misc$earg) = names(misc$link)
        for(ii in 1:(M-1)) misc$earg[[ii]] = list()
        misc$earg[[M]] = .ephi
        misc$expected = TRUE
        if(intercept.only) {
            misc$shape=probs[1,]*(1/phi[1]-1) # phi & probs computed in @deriv
        }
    }), list( .ephi=ephi, .lphi=lphi ))),
    loglikelihood=eval(substitute(
        function(mu,y,w,residuals= FALSE,eta, extra=NULL) {
        M = if(is.matrix(eta)) ncol(eta) else 1
        probs = cbind(exp(eta[,-M]), 1)
        probs = probs / as.vector(probs %*% rep(1, M))
        phi = eta2theta(eta[,M], .lphi, earg= .ephi )
        n = length(phi)
        ycount = as.matrix(y * w)
        if(residuals) stop("loglikelihood residuals not implemented yet") else {
            ans = rep(0.0, len=n)
            omega = extra$n2
            for(jay in 1:M) {
                maxyj = max(ycount[,jay])
                loopOveri = n < maxyj
                if(loopOveri) {
                    for(iii in 1:n) {
                        rrr = 1:ycount[iii,jay] # a vector
                        if(ycount[iii,jay] > 0)
                        ans[iii] = ans[iii] + sum(log((1-phi[iii]) *
                                   probs[iii,jay] + (rrr-1)*phi[iii]))

                    }
                } else {
                    for(rrr in 1:maxyj) {
                        index = (rrr <= ycount[,jay]) & (ycount[,jay] > 0)
                        if(any(index))
                            ans[index] = ans[index] + log((1-phi[index]) *
                                         probs[index,jay] + (rrr-1)*phi[index])
                    }
                }
            } # end of jay loop

            maxomega = max(omega)
            loopOveri = n < maxomega
            if(loopOveri) {
                for(iii in 1:n) {
                    rrr = 1:omega[iii]
                    ans[iii]= ans[iii] - sum(log1p(-phi[iii] + (rrr-1)*phi[iii]))
                }
            } else {
                for(rrr in 1:maxomega) {
                    ind8 = rrr <= omega
                    ans[ind8] = ans[ind8] - log1p(-phi[ind8] + (rrr-1)*phi[ind8])
                }
            }
            sum(ans)
        }
    }, list( .ephi=ephi, .lphi=lphi ))),
    vfamily=c("dirmultinomial "),
    deriv=eval(substitute(expression({
        probs = cbind(exp(eta[,-M]), 1)
        probs = probs / as.vector(probs %*% rep(1, M))
        phi = eta2theta(eta[,M], .lphi, earg= .ephi )
        dl.dprobs = matrix(0.0, n, M-1)
        dl.dphi = rep(0.0, len=n)
        omega = extra$n2
        ycount = as.matrix(y * w)
        for(jay in 1:M) {
            maxyj = max(ycount[,jay])
            loopOveri = n < maxyj
            if(loopOveri) {
                for(iii in 1:n) {
                    rrr = 1:ycount[iii,jay]
                    if(ycount[iii,jay] > 0) {
                        PHI = phi[iii]
                        dl.dphi[iii]=dl.dphi[iii] + sum((rrr-1-probs[iii,jay]) /
                                       ((1-PHI)*probs[iii,jay] + (rrr-1)*PHI))

                        tmp9 = (1-PHI) / ((1-PHI)*probs[iii,jay] + (rrr-1)*PHI)
                        if(jay < M) {
                            dl.dprobs[iii,jay] = dl.dprobs[iii,jay] + sum(tmp9)
                        } else {
                            for(jay2 in 1:(M-1))
                               dl.dprobs[iii,jay2]=dl.dprobs[iii,jay2]-sum(tmp9)
                        }
                    }
                }
            } else {
                for(rrr in 1:maxyj) {
                    index = (rrr <= ycount[,jay]) & (ycount[,jay] > 0)
                    PHI = phi[index]
                    dl.dphi[index] = dl.dphi[index] + (rrr-1-probs[index,jay]) /
                        ((1-PHI)*probs[index,jay] + (rrr-1)*PHI)
                    tmp9 = (1-PHI) / ((1-PHI)*probs[index,jay] + (rrr-1)*PHI)
                    if(jay < M) {
                        dl.dprobs[index,jay] = dl.dprobs[index,jay] + tmp9
                    } else {
                        for(jay2 in 1:(M-1))
                            dl.dprobs[index,jay2] = dl.dprobs[index,jay2] - tmp9
                    }
                }
            }
        } # end of jay loop
        maxomega = max(omega)
        loopOveri = n < maxomega
        if(loopOveri) {
            for(iii in 1:n) {
                rrr = 1:omega[iii]
                dl.dphi[iii]=dl.dphi[iii] - sum((rrr-2)/(1 + (rrr-2)*phi[iii]))
            }
        } else {
            for(rrr in 1:maxomega) {
                index = rrr <= omega
                dl.dphi[index]=dl.dphi[index] - (rrr-2)/(1 + (rrr-2)*phi[index])
            }
        }
        dprobs.deta = probs[,-M] * (1 - probs[,-M])    # n x (M-1)
        dphi.deta = dtheta.deta(phi, .lphi, earg= .ephi )
        ans = cbind(dl.dprobs * dprobs.deta, dl.dphi * dphi.deta)
        ans
    }), list( .ephi=ephi, .lphi=lphi ))),
    weight=eval(substitute(expression({
        wz = matrix(0, n, dimm(M))
        loopOveri = n < maxomega
        if(loopOveri) {
            for(iii in 1:n) {
                rrr = 1:omega[iii]  # A vector
                PHI = phi[iii]
                pYiM.ge.rrr = 1 - pbetabin.ab(q=rrr-1, size=omega[iii],
                    shape1=probs[iii,M]*(1/PHI-1),
                    shape2=(1-probs[iii,M])*(1/PHI-1))  # A vector
                denomM = ((1-PHI)*probs[iii,M] + (rrr-1)*PHI)^2  # A vector
                wz[iii,iam(M,M,M)] = wz[iii,iam(M,M,M)] +
                        sum(probs[iii,M]^2 * pYiM.ge.rrr / denomM) -
                        sum(1 / (1 + (rrr-2)*PHI)^2)
                for(jay in 1:(M-1)) {
                    denomj = ((1-PHI)*probs[iii,jay] + (rrr-1)*PHI)^2
                    pYij.ge.rrr = 1 - pbetabin.ab(q=rrr-1, size=omega[iii],
                        shape1=probs[iii,jay]*(1/PHI-1),
                        shape2=(1-probs[iii,jay])*(1/PHI-1))
                    wz[iii,iam(jay,jay,M)] = wz[iii,iam(jay,jay,M)] + 
                        sum(pYij.ge.rrr / denomj) + 
                        sum(pYiM.ge.rrr / denomM)
                    for(kay in jay:(M-1)) if(kay > jay) {
                        wz[iii,iam(jay,kay,M)] = wz[iii,iam(jay,kay,M)] + 
                            sum(pYiM.ge.rrr / denomM)
                    }
                    wz[iii,iam(jay,M,M)] = wz[iii,iam(jay,M,M)] +
                            sum(probs[iii,jay] * pYij.ge.rrr / denomj) -
                            sum(probs[iii,M]   * pYiM.ge.rrr / denomM)
                    wz[iii,iam(M,M,M)] = wz[iii,iam(M,M,M)] +
                            sum(probs[iii,jay]^2 * pYij.ge.rrr / denomj)
                } # end of jay loop
            } # end of iii loop
        } else {
            for(rrr in 1:maxomega) {
                ind5 = rrr <= omega
                PHI = phi[ind5]
                pYiM.ge.rrr = 1 - pbetabin.ab(q=rrr-1, size=omega[ind5],
                    shape1=probs[ind5,M]*(1/PHI-1),
                    shape2=(1-probs[ind5,M])*(1/PHI-1))
                denomM = ((1-PHI)*probs[ind5,M] + (rrr-1)*PHI)^2
                wz[ind5,iam(M,M,M)] = wz[ind5,iam(M,M,M)] +
                        probs[ind5,M]^2 * pYiM.ge.rrr / denomM -
                        1 / (1 + (rrr-2)*PHI)^2
                for(jay in 1:(M-1)) {
                    denomj = ((1-PHI)*probs[ind5,jay] + (rrr-1)*PHI)^2
                    pYij.ge.rrr = 1 - pbetabin.ab(q=rrr-1, size=omega[ind5],
                        shape1=probs[ind5,jay]*(1/PHI-1),
                        shape2=(1-probs[ind5,jay])*(1/PHI-1))
                    wz[ind5,iam(jay,jay,M)] = wz[ind5,iam(jay,jay,M)] + 
                        pYij.ge.rrr / denomj + pYiM.ge.rrr / denomM 
                    for(kay in jay:(M-1)) if(kay > jay) {
                        wz[ind5,iam(jay,kay,M)] = wz[ind5,iam(jay,kay,M)] + 
                            pYiM.ge.rrr / denomM 
                    }
                    wz[ind5,iam(jay,M,M)] = wz[ind5,iam(jay,M,M)] +
                            probs[ind5,jay] * pYij.ge.rrr / denomj -
                            probs[ind5,M]   * pYiM.ge.rrr / denomM
                    wz[ind5,iam(M,M,M)] = wz[ind5,iam(M,M,M)] +
                            probs[ind5,jay]^2 * pYij.ge.rrr / denomj
                } # end of jay loop
            } # end of rrr loop
        }

        for(jay in 1:(M-1))
            for(kay in jay:(M-1))
                wz[,iam(jay,kay,M)] = wz[,iam(jay,kay,M)] * (1-phi)^2
        for(jay in 1:(M-1))
            wz[,iam(jay,M,M)] = wz[,iam(jay,M,M)] * (phi-1) / phi
        wz[,iam(M,M,M)] = wz[,iam(M,M,M)] / phi^2

        d1Thetas.deta = cbind(dprobs.deta, dphi.deta)
        index = iam(NA, NA, M, both = TRUE, diag = TRUE)
        wz = wz * d1Thetas.deta[,index$row] * d1Thetas.deta[,index$col]
        wz
    }), list( .ephi=ephi, .lphi=lphi ))))
}


dirmul.old = function(link="loge", earg=list(), init.alpha = 0.01,
                      parallel= FALSE, zero=NULL)
{

    if(mode(link) != "character" && mode(link) != "name")
        link = as.character(substitute(link))
    if(length(zero) && !is.Numeric(zero, integer=TRUE, posit=TRUE))
        stop("bad input for argument \"zero\"")
    if(!is.Numeric(init.alpha, posit=TRUE))
        stop("'init.alpha' must contain positive values only")
    if(!is.list(earg)) earg = list()

    new("vglmff",
    blurb=c("Dirichlet-Multinomial distribution\n\n",
            "Links:     ",
            namesof("shape1", link, earg=earg), ", ..., ",
            namesof("shapeM", link, earg=earg), "\n\n",
            "Posterior mean:    (n_j + shape_j)/(2*sum(n_j) + sum(shape_j))\n"),
    constraints=eval(substitute(expression({
        constraints = cm.vgam(matrix(1,M,1), x, .parallel, constraints, int= TRUE)
        constraints = cm.zero.vgam(constraints, x, .zero, M)
    }), list( .parallel=parallel, .zero=zero ))),
    initialize=eval(substitute(expression({
        y = as.matrix(y)
        M = ncol(y)
        if(any(y != round(y )))
            stop("all y values must be integer-valued")

        predictors.names = namesof(paste("shape", 1:M, sep=""), .link,
            earg=.earg, short=TRUE)
        extra$n2 = as.vector(apply(y, 1, sum))  # Nb. don't multiply by 2
        extra$y  = y
        if(!length(etastart)) {
            yy = if(is.numeric(.init.alpha)) 
                matrix(.init.alpha, n, M, byrow= TRUE) else
                matrix(runif(n*M), n, M)
            etastart = theta2eta(yy, .link, earg=.earg)
        }
    }), list( .link=link, .earg=earg, .init.alpha=init.alpha ))),
    inverse=eval(substitute(function(eta, extra=NULL) {
        shape = eta2theta(eta, .link, earg=.earg)
        M = if(is.matrix(eta)) ncol(eta) else 1
        sumshape = as.vector(shape %*% rep(1, len=M))
        (extra$y + shape) / (extra$n2 + sumshape)
    }, list( .link=link, .earg=earg ))),
    last=eval(substitute(expression({
        misc$link = rep(.link, length=M)
        names(misc$link) = paste("shape", 1:M, sep="")
        misc$earg = vector("list", M)
        names(misc$earg) = names(misc$link)
        for(ii in 1:M) misc$earg[[ii]] = .earg
        misc$pooled.weight = pooled.weight
    }), list( .link=link, .earg=earg ))),
    loglikelihood=eval(substitute(
        function(mu,y,w,residuals= FALSE,eta, extra=NULL) {
        shape = eta2theta(eta, .link, earg=.earg)
        M = if(is.matrix(eta)) ncol(eta) else 1
        sumshape = as.vector(shape %*% rep(1, len=M))
        if(residuals) stop("loglikelihood residuals not implemented yet") else
        sum(w*(lgamma(sumshape) - lgamma(extra$n2 + sumshape ))) +
            sum(w * (lgamma(y + shape) - lgamma(shape )))
    }, list( .link=link, .earg=earg ))),
    vfamily=c("dirmul.old"),
    deriv=eval(substitute(expression({
        shape = eta2theta(eta, .link, earg=.earg)
        sumshape = as.vector(shape %*% rep(1, len=M))
        dl.dsh = digamma(sumshape) - digamma(extra$n2 + sumshape) +
                 digamma(y + shape) - digamma(shape)
        dsh.deta = dtheta.deta(shape, .link, earg=.earg)
        w * dl.dsh * dsh.deta
    }), list( .link=link, .earg=earg ))),
    weight=eval(substitute(expression({
        index = iam(NA, NA, M, both = TRUE, diag = TRUE)
        wz = matrix(trigamma(sumshape)-trigamma(extra$n2 + sumshape),
                    nrow=n, ncol=dimm(M))
        wz[,1:M] = wz[,1:M] + trigamma(y + shape) - trigamma(shape)
        wz = -wz * dsh.deta[, index$row] * dsh.deta[, index$col]


        if(TRUE && intercept.only) {
            sumw = sum(w)
            for(i in 1:ncol(wz))
                wz[,i] = sum(wz[,i]) / sumw
            pooled.weight = TRUE
            wz = w * wz   # Put back the weights
        } else
            pooled.weight = FALSE

        wz
    }), list( .link=link, .earg=earg ))))
}




rdiric = function(n, shape, dimension=NULL) {
    if(!is.numeric(dimension))
        dimension = length(shape)
    shape = rep(shape, len=dimension)

    ans = if(is.R()) rgamma(n*dimension, rep(shape, rep(n, dimension ))) else
               rgamma(n*dimension, rep(shape, each=n)) 
    dim(ans) = c(n, dimension) 


    ans = ans / apply(ans, 1, sum)
    ans 
}


dirichlet = function(link="loge", earg=list(), zero=NULL)
{
    if(mode(link) != "character" && mode(link) != "name")
        link = as.character(substitute(link))
    if(length(zero) && !is.Numeric(zero, integer=TRUE, posit=TRUE))
        stop("bad input for argument \"zero\"")
    if(!is.list(earg)) earg = list()

    new("vglmff",
    blurb=c("Dirichlet distribution\n\n",
            "Links:     ",
            namesof("shapej", link, earg=earg), "\n\n",
            "Mean:     shape_j/(1 + sum(shape_j)), j=1,..,ncol(y)"),
    constraints=eval(substitute(expression({
        constraints = cm.zero.vgam(constraints, x, .zero, M)
    }), list( .zero=zero ))),
    initialize=eval(substitute(expression({
        y = as.matrix(y)
        M = ncol(y)
        if(any(y <= 0) || any(y>=1))
            stop("all y values must be > 0 and < 1")
        predictors.names = namesof(paste("shape", 1:M, sep=""), .link,
           earg=.earg, short=TRUE)
        if(!length(etastart)) {
            yy = matrix(t(y) %*% rep(1/nrow(y), nrow(y)), nrow(y), M, byrow= TRUE)
            etastart = theta2eta(yy, .link, earg= .earg )
        }
    }), list( .link=link, .earg=earg ))),
    inverse=eval(substitute(function(eta, extra=NULL) {
        shape = eta2theta(eta, .link, earg= .earg )
        M = if(is.matrix(eta)) ncol(eta) else 1
        sumshape = as.vector(shape %*% rep(1, len=M))  # apply(shape, 1, sum)
        shape / sumshape
    }, list( .link=link, .earg=earg ))),
    last=eval(substitute(expression({
        misc$link = c(shape= .link)
        temp.names = paste("shape", 1:M, sep="")
        misc$link = rep( .link, len=M)
        names(misc$link) = temp.names
        misc$earg = vector("list", M)
        names(misc$earg) = names(misc$link)
        for(ii in 1:M) misc$earg[[ii]] = .earg
    }), list( .link=link, .earg=earg ))),
    loglikelihood=eval(substitute(
        function(mu,y,w,residuals= FALSE,eta, extra=NULL) {
        shape = eta2theta(eta, .link, earg= .earg )
        M = if(is.matrix(eta)) ncol(eta) else 1
        sumshape = as.vector(shape %*% rep(1, len=M))
        if(residuals) stop("loglikelihood residuals not implemented yet") else
        sum(w * (lgamma(sumshape) - lgamma(shape) + (shape-1)*log(y )))
    }, list( .link=link, .earg=earg ))),
    vfamily=c("dirichlet"),
    deriv=eval(substitute(expression({
        shape = eta2theta(eta, .link, earg= .earg )
        sumshape = as.vector(shape %*% rep(1, len=M))
        dl.dsh = digamma(sumshape) - digamma(shape) + log(y)
        dsh.deta = dtheta.deta(shape, .link, earg= .earg )
        w * dl.dsh * dsh.deta
    }), list( .link=link, .earg=earg ))),
    weight=expression({
        index = iam(NA, NA, M, both = TRUE, diag = TRUE)
        wz = matrix(trigamma(sumshape), nrow=n, ncol=dimm(M))
        wz[,1:M] = wz[,1:M] - trigamma(shape)
        wz = -w * wz * dsh.deta[, index$row] * dsh.deta[, index$col]
        wz
    }))
}



zeta = function(x, deriv=0) {


    deriv.arg = deriv
    if(!is.Numeric(deriv.arg, allow=1, integer=TRUE, positi=TRUE) && deriv.arg!=0)
        stop("'deriv' must be a single non-negative integer")
    if(!(deriv.arg==0 || deriv.arg==1 || deriv.arg==2))
        stop("'deriv' must be 0, 1, or 2")



    if(deriv.arg > 0)
        return(zeta.derivative(x, deriv=deriv))



    if(any(special <- Re(x) <= 1)) {
        ans <- x
        ans[special] <- Inf   # For Re(x)==1

        special3 <- Re(x) < 1
        ans[special3] <- NA # For 0 < Re(x) < 1

        special4 <- (0 < Re(x)) & (Re(x) < 1) & (Im(x) == 0)
        ans[special4] <- zeta.derivative(x[special4], deriv=deriv)


        special2 <- Re(x) < 0
        if(any(special2)) {
            x2 = x[special2]
            cx = 1-x2
            ans[special2] = 2^(x2) * pi^(x2-1) * sin(pi*x2/2) * gamma(cx) * Recall(cx)
        }

        if(any(!special)) {
            ans[!special] <- Recall(x[!special])
        }
        return(ans)
    }

    a=12; k=8  # Markman paper 
    B = c(1/6, -1/30,1/42,-1/30,5/66,-691/2730,7/6,-3617/510)
    ans = 0
    for(i in 1:(a-1))
       ans = ans + 1.0 / i^x
    ans = ans + 1.0 / ((x-1.0)* a^(x-1.0)) + 1.0 / (2.0 * a^x)

    term = (x/2) / a^(x+1)
    ans = ans + term * B[1]

    for(m in 2:k) {
        term = term * (x+2*m-2) * (x+2*m-3) / (a*a* 2*m * (2*m-1))
        ans = ans + term * B[m]
    }
    ans
}


zeta.derivative = function(x, deriv=0) 
{


    deriv.arg = deriv
    if(!is.Numeric(deriv.arg, allow=1, integer=TRUE, positi=TRUE) && deriv.arg!=0)
        stop("'deriv' must be a single non-negative integer")
    if(!(deriv.arg==0 || deriv.arg==1 || deriv.arg==2))
        stop("'deriv' must be 0, 1, or 2")

    if(any(Im(x) != 0))
        stop("Sorry, currently can only handle x real, not complex")
    if(any(x < 0))
        stop("Sorry, currently cannot handle x < 0")

    ok = is.finite(x) & x > 0 & x != 1   # Handles NAs
    ans = rep(as.numeric(NA), length(x))
    nn = sum(ok)  # Effective length (excludes x < 0 and x = 1 values)
    if(nn)
        ans[ok] = dotFortran(name="vzetawr", as.double(x[ok]), ans=double(nn),
                  as.integer(deriv.arg), as.integer(nn))$ans



    if(deriv==0)
        ans[is.finite(x) & abs(x) < 1.0e-12] = -0.5 

    ans
}


dzeta = function(x, p) 
{
    if(!is.Numeric(p, allow=1, posit=TRUE) || p <= 1)
        stop("'p' must be numeric and > 1")
    p = rep(p, len=length(x))

    ox = !is.finite(x)
    zero = ox | round(x) != x | x < 1
    if(any(zero)) warning("non-integer x and/or x < 1 or NAs")
    ans = rep(0.0, len=length(x))
    if(any(!zero))
        ans[!zero] = x[!zero]^(-p[!zero]) / zeta(p[!zero]) 
    if(any(ox)) ans[ox] = NA
    ans
}

zetaff = function(link="loge", earg=list(), init.p=NULL)
{

    if(length(init.p) && !is.Numeric(init.p, positi=TRUE))
        stop("argument \"init.p\" must be > 0")
    if(mode(link) != "character" && mode(link) != "name")
        link = as.character(substitute(link))
    if(!is.list(earg)) earg = list()

    new("vglmff",
    blurb=c("Zeta distribution f(y) = 1/(y^(p+1) zeta(p+1)), p>0, y=1,2,..\n\n",
            "Link:    ",
            namesof("p", link, earg=earg), "\n\n",
            "Mean:     zeta(p) / zeta(p+1), provided p>1\n",
            "Variance: zeta(p-1) / zeta(p+1) - mean^2, provided p>2"),
    initialize=eval(substitute(expression({
        y = as.numeric(y)
        if(any(y < 1))
            stop("all y values must be in 1,2,3,...")
        if(any(y != round(y )))
            warning("y should be integer-valued")
        if(ncol(cbind(y)) != 1)
            stop("response must be a vector or a one-column matrix")

        predictors.names = namesof("p", .link, earg=.earg, tag=FALSE) 

        if(!length(etastart)) {
            zetaff.Loglikfun = function(pp, y, x, w, extraargs) {
                sum(w * (-(pp+1) * log(y) - log(zeta(pp+1))))
            }
            p.grid = seq(0.1, 3.0, len=19)
            pp.init = if(length( .init.p )) .init.p else
                      getMaxMin(p.grid, objfun=zetaff.Loglikfun, y=y,  x=x, w=w)
            pp.init = rep(pp.init, length=length(y))
            if( .link == "loglog") pp.init[pp.init <= 1] = 1.2
            etastart = theta2eta(pp.init, .link, earg=.earg)
        }
    }), list( .link=link, .earg=earg, .init.p=init.p ))),
    inverse=eval(substitute(function(eta, extra=NULL) {
        ans = pp = eta2theta(eta, .link, earg=.earg)
        ans[pp>1] = zeta(pp[pp>1]) / zeta(pp[pp>1]+1)
        ans[pp<=1] = NA
        ans
    }, list( .link=link, .earg=earg ))),
    last=eval(substitute(expression({
        misc$link = c(pp= .link)
        misc$earg = list(pp = .earg)
    }), list( .link=link, .earg=earg ))),
    loglikelihood=eval(substitute(
        function(mu,y,w,residuals= FALSE,eta, extra=NULL) {
        pp = eta2theta(eta, .link, earg=.earg)
        if(residuals) stop("loglikelihood residuals not implemented yet") else
        sum(w * (-(pp+1) * log(y) - log(zeta(pp+1 ))))
    }, list( .link=link, .earg=earg ))),
    vfamily=c("zeta"),
    deriv=eval(substitute(expression({
        pp = eta2theta(eta, .link, earg=.earg)
        fred1 = zeta(pp+1)
        fred2 = zeta(pp+1, deriv=1)
        dl.dpp = -log(y) - fred2 / fred1
        dpp.deta = dtheta.deta(pp, .link, earg=.earg)
        w * dl.dpp * dpp.deta
    }), list( .link=link, .earg=earg ))),
    weight=expression({
        ed2l.dpp2 = zeta(pp+1, deriv=2) / fred1 - (fred2/fred1)^2
        wz = w * dpp.deta^2 * ed2l.dpp2
        wz
    }))
}



gharmonic = function(n, s=1, lognexponent=0) {
    if(!is.Numeric(n, integ=TRUE, posit=TRUE))
        stop("bad input for argument \"n\"")
    if(!is.Numeric(s, posit=TRUE))
        stop("bad input for argument \"s\"")
    if(!is.Numeric(lognexponent, allow=1))
        stop("bad input for argument \"lognexponent\"")
    if(length(n) == 1 && length(s) == 1) {
        if(lognexponent != 0) sum(log(1:n)^lognexponent * (1:n)^(-s)) else
            sum((1:n)^(-s))
    } else {
        LEN = max(length(n), length(s))
        n = rep(n, len=LEN)
        ans = s = rep(s, len=LEN)
        if(lognexponent != 0) {
            for(i in 1:LEN)
                ans[i] = sum(log(1:n[i])^lognexponent * (1:n[i])^(-s[i]))
        } else
            for(i in 1:LEN)
                ans[i] = sum((1:n[i])^(-s[i]))
        ans
    }
}

dzipf = function(x, N, s)
{
    if(!is.Numeric(x))
        stop("bad input for argument \"x\"")
    if(!is.Numeric(N, integ=TRUE, posit=TRUE))
        stop("bad input for argument \"N\"")
    if(!is.Numeric(s, posit=TRUE))
        stop("bad input for argument \"s\"")
    nn = max(length(x), length(N), length(s))
    x = rep(x, len=nn); N = rep(N, len=nn); s = rep(s, len=nn);
    ox = !is.finite(x)
    zero = ox | round(x) != x | x < 1 | x > N
    ans = 0 * x
    if(any(!zero))
        ans[!zero] = x[!zero]^(-s[!zero]) / gharmonic(N[!zero], s[!zero])
    ans
}



pzipf = function(q, N, s) {
    if(!is.Numeric(q))
        stop("bad input for argument \"q\"")
    if(!is.Numeric(N, integ=TRUE, posit=TRUE))
        stop("bad input for argument \"N\"")
    if(!is.Numeric(s, posit=TRUE))
        stop("bad input for argument \"s\"")

    nn = max(length(q), length(N), length(s))
    q = rep(q, len=nn); N = rep(N, len=nn); s = rep(s, len=nn);
    oq = !is.finite(q)
    zeroOR1 = oq | q < 1 | q >= N
    floorq = floor(q)
    ans = 0 * floorq
    ans[oq | q >= N] = 1
    if(any(!zeroOR1))
        ans[!zeroOR1] = gharmonic(floorq[!zeroOR1], s[!zeroOR1]) /
                        gharmonic(N[!zeroOR1], s[!zeroOR1])
    ans
}


zipf = function(N=NULL, link="loge", earg=list(), init.s=NULL)
{
    if(length(N) &&
      (!is.Numeric(N, positi=TRUE, integ=TRUE, allow=1) || N <= 1))
        stop("bad input for argument \"N\"")
    enteredN = length(N)
    if(length(init.s) && !is.Numeric(init.s, positi=TRUE))
        stop("argument \"init.s\" must be > 0")

    if(mode(link) != "character" && mode(link) != "name")
        link = as.character(substitute(link))
    if(!is.list(earg)) earg = list()

    new("vglmff",
    blurb=c("Zipf distribution f(y;s) = y^(-s) / sum((1:N)^(-s)),",
            " s>0, y=1,2,...,N", ifelse(enteredN, paste("=",N,sep=""), ""),
            "\n\n",
            "Link:    ",
            namesof("s", link, earg=earg),
            "\n\n",
            "Mean:    gharmonic(N,s-1) / gharmonic(N,s)"),
    initialize=eval(substitute(expression({
        if(ncol(cbind(y)) != 1)
            stop("response must be a vector or a one-column matrix")
        y = as.numeric(y)
        if(any(y != round(y )))
            stop("y must be integer-valued")
        predictors.names = namesof("s", .link, earg= .earg, tag=FALSE) 
        NN = .N
        if(!is.Numeric(NN, allow=1, posit=TRUE, integ=TRUE))
            NN = max(y)
        if(max(y) > NN)
            stop("maximum of the response is greater than argument \"N\"")
        if(any(y < 1))
            stop(paste("all response values must be in 1,2,3,...,N=",NN,sep=""))
        extra$N = NN
        if(!length(etastart)) {
            llfun = function(ss, y, N, w) {
                sum(w * ((-ss) * log(y) - log(gharmonic(N, ss))))
            }
            ss.init = if(length( .init.s )) .init.s else
                getInitVals(gvals=seq(0.1, 3.0, len=19), llfun=llfun,
                            y=y, N=extra$N, w=w)
            ss.init = rep(ss.init, length=length(y))
            if( .link == "loglog") ss.init[ss.init <= 1] = 1.2
            etastart = theta2eta(ss.init, .link, earg= .earg)
        }
    }), list( .link=link, .earg=earg, .init.s=init.s, .N=N ))),
    inverse=eval(substitute(function(eta, extra=NULL) {
        ss = eta2theta(eta, .link, earg= .earg)
        gharmonic(extra$N, s=ss-1) / gharmonic(extra$N, s=ss)
    }, list( .link=link, .earg=earg ))),
    last=eval(substitute(expression({
        misc$expected = FALSE
        misc$link = c(s= .link)
        misc$earg = list(s= .earg )
        misc$N = extra$N
    }), list( .link=link, .earg=earg ))),
    loglikelihood=eval(substitute(
        function(mu,y,w,residuals= FALSE,eta, extra=NULL) {
        ss = eta2theta(eta, .link, earg= .earg)
        if(residuals) stop("loglikelihood residuals not implemented yet") else
        sum(w * ((-ss) * log(y) - log(gharmonic(extra$N, ss ))))
    }, list( .link=link, .earg=earg ))),
    vfamily=c("zipf"),
    deriv=eval(substitute(expression({
        ss = eta2theta(eta, .link, earg= .earg)
        fred1 = gharmonic(extra$N, ss)
        fred2 = gharmonic(extra$N, ss, lognexp=1)
        dl.dss = -log(y) + fred2 / fred1
        dss.deta = dtheta.deta(ss, .link, earg= .earg)
        d2ss.deta2 = d2theta.deta2(ss, .link, earg= .earg)
        w * dl.dss * dss.deta
    }), list( .link=link, .earg=earg ))),
    weight=expression({
        d2l.dss = gharmonic(extra$N, ss, lognexp=2) / fred1 - (fred2/fred1)^2
        wz = w * (dss.deta^2 * d2l.dss - d2ss.deta2 * dl.dss)
        wz
    }))
}








cauchy1 = function(scale.arg=1, llocation="identity",
                   elocation=list(),
                   ilocation=NULL, method.init=1)
{
    if(mode(llocation) != "character" && mode(llocation) != "name")
        llocation = as.character(substitute(llocation))
    if(!is.Numeric(scale.arg, posit=TRUE)) stop("bad input for \"scale.arg\"")
    if(!is.Numeric(method.init, allow=1, integ=TRUE, posit=TRUE) ||
       method.init > 3)
        stop("'method.init' must be 1 or 2 or 3")
    if(!is.list(elocation)) elocation = list()

    new("vglmff",
    blurb=c("One parameter Cauchy distribution (location unknown, scale known)\n\n",
            "Link:    ",
            namesof("location", llocation, earg=elocation), "\n\n",
            "Mean:     NA\n",
            "Variance: NA"),
    initialize=eval(substitute(expression({
        predictors.names = namesof("location", .llocation,
            earg=.elocation, tag=FALSE)
        if(ncol(cbind(y)) != 1)
            stop("response must be a vector or a one-column matrix")

        if(!length(etastart)) {
            loc.init = if(length(.ilocation)) .ilocation else {
                if( .method.init == 2) median(rep(y,w)) else 
                if( .method.init == 3) y else {
                    cauchy1.Loglikfun = function(loc, y, x, w, extraargs) {
                         scal = extraargs
                         sum(w * (-log1p(((y-loc)/scal)^2) - log(scal)))
                     }
                     loc.grid = quantile(y, probs=seq(0.1, 0.9, by=0.05))
                     try.this = getMaxMin(loc.grid, objfun=cauchy1.Loglikfun,
                                          y=y,  x=x, w=w, extraargs= .scale.arg)
                    try.this = rep(try.this, len=n)
                    try.this
                }
            }
            loc.init = rep(loc.init, len=n)
            if(.llocation == "loge") loc.init = abs(loc.init)+0.01
            etastart = theta2eta(loc.init, .llocation, earg=.elocation)
        }
    }), list( .scale.arg=scale.arg, .ilocation=ilocation,
              .elocation=elocation, .llocation=llocation,
              .method.init=method.init ))),
    inverse=function(eta, extra=NULL) {
        rep(as.numeric(NA), length(eta)) 
    },
    last=eval(substitute(expression({
        misc$link = c("location"= .llocation)
        misc$earg = list(location= .elocation )
        misc$scale.arg = .scale.arg 
    }), list( .scale.arg=scale.arg, .elocation=elocation,
             .llocation=llocation ))),
    loglikelihood=eval(substitute(
        function(mu,y,w,residuals= FALSE,eta, extra=NULL) {
        location = eta2theta(eta, .llocation, earg=.elocation)
        temp = (y-location)/ .scale.arg
        if(residuals) stop("loglikelihood residuals not implemented yet") else
        sum(w * (-log1p(temp^2) - log(pi) - log(.scale.arg )))
    }, list( .scale.arg=scale.arg, .elocation=elocation,
             .llocation=llocation ))),
    vfamily=c("cauchy1"),
    deriv=eval(substitute(expression({
        location = eta2theta(eta, .llocation, earg=.elocation)
        temp = (y-location)/.scale.arg
        dl.dlocation = 2 * temp / ((1 + temp^2) * .scale.arg)
        dlocation.deta = dtheta.deta(location, .llocation, earg=.elocation)
        w * dl.dlocation * dlocation.deta
    }), list( .scale.arg=scale.arg, .elocation=elocation,
              .llocation=llocation ))),
    weight=eval(substitute(expression({
        wz = w * dlocation.deta^2 / (.scale.arg^2 * 2)
        wz
    }), list( .scale.arg=scale.arg, .elocation=elocation,
              .llocation=llocation ))))
}






logistic1 = function(llocation="identity",
                     elocation=list(),
                     scale.arg=1, method.init=1)
{
    if(mode(llocation) != "character" && mode(llocation) != "name")
        llocation = as.character(substitute(llocation))
    if(!is.Numeric(scale.arg, allow=1, posit=TRUE))
        stop("'scale.arg' must be a single positive number")
    if(!is.Numeric(method.init, allow=1, integ=TRUE, posit=TRUE) ||
       method.init > 2)
        stop("'method.init' must be 1 or 2")
    if(!is.list(elocation)) elocation = list()

    new("vglmff",
    blurb=c("One-parameter logistic distribution (location unknown, scale known)\n\n",
            "Link:    ",
            namesof("location", llocation, earg=elocation), "\n\n",
            "Mean:     location", "\n",
            "Variance: (pi*scale)^2 / 3"),
    initialize=eval(substitute(expression({
        predictors.names = namesof("location", .llocation, 
            earg= .elocation, tag=FALSE)
        if(!length(etastart)) {
            location.init = if( .method.init == 1) y else median(rep(y, w))
            location.init = rep(location.init, len=n)
            if(.llocation == "loge") location.init = abs(location.init) + 0.001
            etastart = theta2eta(location.init, .llocation, earg= .elocation)
        }
    }), list( .method.init=method.init, .llocation=llocation,
              .elocation=elocation ))),
    inverse=eval(substitute(function(eta, extra=NULL) {
        eta2theta(eta, .llocation, earg= .elocation)
    }, list( .llocation=llocation,
             .elocation=elocation ))),
    last=eval(substitute(expression({
        misc$link = c(location= .llocation)
        misc$earg = list(location= .elocation )
        misc$scale.arg = .scale.arg 
    }), list( .llocation=llocation, 
              .elocation=elocation, .scale.arg=scale.arg ))),
    loglikelihood=eval(substitute(
        function(mu,y,w,residuals= FALSE,eta, extra=NULL) {
        location = eta2theta(eta, .llocation, earg= .elocation)
        zedd = (y-location)/.scale.arg
        if(residuals) stop("loglikelihood residuals not implemented yet") else
        sum(w * (-zedd - 2 * log1p(exp(-zedd)) - log(.scale.arg )))
    }, list( .llocation=llocation,
             .elocation=elocation, .scale.arg=scale.arg ))),
    vfamily=c("logistic1"),
    deriv=eval(substitute(expression({
        location = eta2theta(eta, .llocation, earg= .elocation)
        ezedd = exp(-(y-location)/.scale.arg)
        dl.dlocation = (1 - ezedd) / ((1 + ezedd) * .scale.arg)
        dlocation.deta = dtheta.deta(location, .llocation, earg= .elocation)
        w * dl.dlocation * dlocation.deta
    }), list( .llocation=llocation,
              .elocation=elocation, .scale.arg=scale.arg ))),
    weight=eval(substitute(expression({
        wz = w * dlocation.deta^2 / (.scale.arg^2 * 3) 
        wz
    }), list( .scale.arg=scale.arg ))))
}




erlang = function(shape.arg, link="loge", earg=list(), method.init=1)
{

    if(!is.Numeric(shape.arg, allow=1, integer=TRUE, positi=TRUE))
        stop("\"shape\" must be a positive integer")
    if(!is.Numeric(method.init, allow=1, integer=TRUE, positi=TRUE) ||
       method.init > 2)
        stop("\"method.init\" must be 1 or 2")

    if(mode(link) != "character" && mode(link) != "name")
        link = as.character(substitute(link))
    if(!is.list(earg)) earg = list()

    new("vglmff",
    blurb=c("Erlang distribution\n\n",
            "Link:    ", namesof("scale", link, earg=earg), "\n", "\n",
            "Mean:     shape * scale", "\n",
            "Variance: shape * scale^2"),
    initialize=eval(substitute(expression({
        if(ncol(y <- as.matrix(y)) > 1)
            stop("erlang cannot handle matrix responses yet")
        if(any(y < 0))
            stop("all y values must be >= 0")

        predictors.names = namesof("scale", .link, earg=.earg, tag=FALSE)

        if(!length(etastart)) {
            if(.method.init==1) 
                sc.init = y / .shape.arg
            if(.method.init==2) {
                sc.init = median(y) / .shape.arg
                sc.init = rep(sc.init, length=n) 
            }
            etastart = theta2eta(sc.init, .link, earg=.earg)
        }
    }), list( .link=link, .earg=earg,
              .shape.arg=shape.arg, .method.init=method.init ))),
    inverse=eval(substitute(function(eta, extra=NULL) {
        sc = eta2theta(eta, .link, earg=.earg)
        .shape.arg * sc 
    }, list( .link=link, .earg=earg, .shape.arg=shape.arg ))),
    last=eval(substitute(expression({
        misc$link = c(scale= .link)
        misc$earg = list(scale= .earg )
        misc$shape.arg = .shape.arg 
    }), list( .link=link, .earg=earg, .shape.arg=shape.arg ))),
    loglikelihood=eval(substitute(
        function(mu,y,w,residuals= FALSE,eta, extra=NULL) {
        sc = eta2theta(eta, .link, earg=.earg)
        if(residuals) stop("loglikelihood residuals not implemented yet") else
        sum(w * ((.shape.arg - 1) * log(y) - y / sc - .shape.arg * log(sc) -
                 lgamma( .shape.arg )))
    }, list( .link=link, .earg=earg, .shape.arg=shape.arg ))),
    vfamily=c("erlang"),
    deriv=eval(substitute(expression({
        sc = eta2theta(eta, .link, earg=.earg)
        dl.dsc = (y / sc - .shape.arg) / sc
        dsc.deta = dtheta.deta(sc, .link, earg=.earg)
        w * dl.dsc * dsc.deta
    }), list( .link=link, .earg=earg, .shape.arg=shape.arg ))),
    weight=eval(substitute(expression({
        ed2l.dsc2 = .shape.arg / sc^2 # Use the expected info matrix
        wz = w * dsc.deta^2 * ed2l.dsc2
        wz
    }), list( .earg=earg, .shape.arg=shape.arg ))))
}



borel.tanner = function(shape.arg, link="logit", earg=list())
{

    if(!is.Numeric(shape.arg, allow=1, integ=TRUE))
        stop("bad input for argument \"shape.arg\"")

    if(mode(link) != "character" && mode(link) != "name")
        link = as.character(substitute(link))
    if(!is.list(earg)) earg = list()

    new("vglmff",
    blurb=c("Borel-Tanner distribution\n\n",
            "Link:    ",
            namesof("a", link, earg=earg), "\n\n",
            "Mean:     n/(1-a)",
            "\n",
            "Variance: n*a / (1-a)^3"),
    initialize=eval(substitute(expression({
        y = as.numeric(y)
        if(any(y < .shape.arg))
            stop("all y values must be >= n")
        if(max(abs(y - round(y )))>0.00001)
            stop("response must be integer-valued")
        if(ncol(cbind(y)) != 1)
            stop("response must be a vector or a one-column matrix")

        predictors.names = namesof("a", .link, earg=.earg, tag=FALSE)


        if(!length(etastart)) {
            a.init = .shape.arg / y 
            etastart = theta2eta(a.init, .link, earg=.earg)
        }
    }), list( .link=link, .earg=earg, .shape.arg=shape.arg ))),
    inverse=eval(substitute(function(eta, extra=NULL) {
        a = eta2theta(eta, .link, earg=.earg)
        .shape.arg / (1 - a)
    }, list( .link=link, .earg=earg, .shape.arg=shape.arg ))),
    last=eval(substitute(expression({
        misc$link = c(a= .link)
        misc$earg = list(a= .earg )
        misc$shape.arg = .shape.arg 
    }), list( .link=link, .earg=earg, .shape.arg=shape.arg ))),
    loglikelihood=eval(substitute(
        function(mu,y,w,residuals= FALSE,eta, extra=NULL) {
        a = eta2theta(eta, .link, earg=.earg)
        if(residuals) stop("loglikelihood residuals not implemented yet") else
        sum(w * ((y- .shape.arg) * log(a) - a * y))
    }, list( .link=link, .earg=earg, .shape.arg=shape.arg ))),
    vfamily=c("borel.tanner"),
    deriv=eval(substitute(expression({
        a = eta2theta(eta, .link, earg=.earg)
        dl.da = (y- .shape.arg)/a - y 
        da.deta = dtheta.deta(a, .link, earg=.earg)
        w * dl.da * da.deta
    }), list( .link=link, .earg=earg, .shape.arg=shape.arg ))),
    weight=eval(substitute(expression({
        ed2l.da2 = .shape.arg/(a*(1-a))   # Use the expected info matrix
        wz = w * da.deta^2 * ed2l.da2
        wz
    }), list( .shape.arg=shape.arg ))))
}


dsnorm = function(x, location=0, scale=1, shape=0) {
    if(!is.Numeric(scale, posit=TRUE))
        stop("bad input for argument \"scale\"")
    zedd = (x - location) / scale
    2 * dnorm(zedd) * pnorm(shape * zedd) / scale
}



rsnorm = function(n, location=0, scale=1, shape=0) {
    if(!is.Numeric(n, posit=TRUE, integ=TRUE, allow=1))
        stop("bad input for argument \"n\"")
    if(!is.Numeric(scale, posit=TRUE))
        stop("bad input for argument \"scale\"")
    if(!is.Numeric(shape)) stop("bad input for argument \"shape\"")
    rho = shape / sqrt(1 + shape^2)
    u0 = rnorm(n)
    v = rnorm(n)
    u1 = rho*u0 + sqrt(1 - rho^2) * v
    location + scale * ifelse(u0 >= 0, u1, -u1)
}


skewnormal1 = function(lshape="identity", earg = list(), ishape=NULL)
{
    if(mode(lshape) != "character" && mode(lshape) != "name")
        lshape = as.character(substitute(lshape))
    if(!is.list(earg)) earg = list()

    new("vglmff",
    blurb=c("1-parameter Skew-normal distribution\n\n",
            "Link:     ",
            namesof("shape", lshape, earg=earg), "\n",
            "Mean:     shape * sqrt(2 / (pi * (1+shape^2 )))\n",
            "Variance: 1-mu^2"),
    initialize=eval(substitute(expression({
        y = cbind(y)
        if(ncol(y) != 1)
            stop("response must be a vector or a one-column matrix")
        predictors.names = namesof("shape", .lshape, earg=.earg, tag=FALSE)
        if(!length(etastart)) {
            init.shape = if(length( .ishape)) rep( .ishape, len=n) else {
                temp = y
                index = abs(y) < sqrt(2/pi)-0.01
                temp[!index] = y[!index]
                temp[index] = sign(y[index])/sqrt(2/(pi*y[index]*y[index])-1)
                temp
            }
            etastart = matrix(init.shape, n, ncol(y))
        }
    }), list( .lshape=lshape, .earg=earg, .ishape=ishape ))), 
    inverse=eval(substitute(function(eta, extra=NULL) {
        alpha = eta2theta(eta, .lshape, earg=.earg)
        alpha * sqrt(2/(pi * (1+alpha^2 )))
    }, list( .earg=earg, .lshape=lshape ))),
    last=eval(substitute(expression({
        misc$link = c(shape= .lshape) 
        misc$earg = list(shape= .earg )
    }), list( .earg=earg, .lshape=lshape ))),
    link=eval(substitute(function(mu, extra=NULL) {
        alpha = mu / sqrt(2/pi - mu^2)
        theta2eta(alpha, .lshape, earg=.earg)
    }, list( .earg=earg, .lshape=lshape ))),
    loglikelihood=eval(substitute(
         function(mu, y, w, residuals = FALSE, eta, extra=NULL) {
            alpha = eta2theta(eta, .lshape, earg=.earg)
        if(residuals) stop("loglikelihood residuals not implemented yet") else
            sum(w * (pnorm(y*alpha, log=TRUE )))
    }, list( .earg=earg, .lshape=lshape ))), 
    vfamily=c("skewnormal1"),
    deriv=eval(substitute(expression({
        alpha = eta2theta(eta, .lshape, earg=.earg)
        zedd = y*alpha
        tmp76 = pnorm(zedd)
        tmp86 = dnorm(zedd)
        dl.dshape = tmp86 * y / tmp76
        dshape.deta = dtheta.deta(alpha, .lshape, earg=.earg)
        w * dl.dshape * dshape.deta
    }), list( .earg=earg, .lshape=lshape ))),
    weight=eval(substitute(expression({
        d2shape.deta2 = d2theta.deta2(alpha, .lshape, earg=.earg)
        d2l.dshape = -y*y * tmp86 * (tmp76 * zedd + tmp86) / tmp76^2
        wz = -(dshape.deta^2) * d2l.dshape - d2shape.deta2 * dl.dshape
        wz = w * wz
        wz[wz < .Machine$double.eps] = .Machine$double.eps
        wz
    }), list( .earg=earg, .lshape=lshape ))))
}


betaff = function(link="loge", earg=list(),
                  i1=NULL, i2=NULL, trim=0.05,
                  A=0, B=1, zero=NULL)
{
    if(mode(link) != "character" && mode(link) != "name")
        link = as.character(substitute(link))
    if(length(zero) && !is.Numeric(zero, integer=TRUE, posit=TRUE))
        stop("bad input for argument \"zero\"")

    if(!is.Numeric(A, allow=1) || !is.Numeric(B, allow=1) || A >= B)
        stop("A must be < B, and both must be of length one")
    stdbeta = (A==0 && B==1)  # stdbeta==T iff standard beta distribution
    if(!is.list(earg)) earg = list()

    new("vglmff",
    blurb=c("Two-parameter Beta distribution\n",
            if(stdbeta)
            paste("y^(shape1-1) * (1-y)^(shape2-1) / B(shape1,shape2),",
            "0<=y<=1, shape1>0, shape2>0\n\n")
            else
            paste("(y-",A,")^(shape1-1) * (",B,
            "-y)^(shape2-1) / [B(shape1,shape2) * (",
            B, "-", A, ")^(shape1+shape2-1)], ",
             A,"<=y<=",B," shape1>0, shape2>0\n\n", sep=""),
            "Links:    ",
            namesof("shape1", link, earg=earg),  ", ",
            namesof("shape2", link, earg=earg)),
    constraints=eval(substitute(expression({
        constraints = cm.zero.vgam(constraints, x, .zero, M)
    }), list( .zero=zero ))),
    initialize=eval(substitute(expression({
        if(min(y) <= .A || max(y) >= .B)
            stop("data not within (A, B)")
        if(ncol(cbind(y)) != 1)
            stop("response must be a vector or a one-column matrix")
        predictors.names = c(namesof("shape1", .link, earg= .earg, short=TRUE),
                             namesof("shape2", .link, earg= .earg, short=TRUE))
        if(is.numeric( .i1 ) && is.numeric( .i2 )) {
            vec = c(.i1, .i2)
            vec = c(theta2eta(vec[1], .link, earg= .earg ),
                    theta2eta(vec[2], .link, earg= .earg ))
            etastart = matrix(vec, n, 2, byrow= TRUE)
        }
        if(!length(etastart)) {
            mu1d = mean(y, trim=.trim)
            uu = (mu1d-.A) / (.B - .A) 
            DD = (.B - .A)^2 
            pinit = uu^2 * (1-uu)*DD/var(y) - uu   # But var(y) is not robust
            qinit = pinit * (1-uu) / uu
            etastart = matrix(theta2eta(c(pinit,qinit), .link, earg= .earg ),
                              n, 2, byrow=TRUE)
        }
    }), list( .link=link, .i1=i1, .i2=i2, .trim=trim, .A=A, .B=B,
              .earg=earg ))),
    inverse=eval(substitute(function(eta, extra=NULL) {
        shapes = eta2theta(eta, .link, earg= .earg )
        .A + (.B-.A) * shapes[,1] / (shapes[,1] + shapes[,2])
    }, list( .link=link, .A=A, .B=B, .earg=earg ))),
    last=eval(substitute(expression({
        misc$link = c(shape1 = .link, shape2 = .link)
        misc$limits = c(.A, .B)
        misc$earg = list(shape1= .earg, shape2= .earg)
    }), list( .link=link, .A=A, .B=B, .earg=earg ))),
    loglikelihood=eval(substitute(
         function(mu, y, w, residuals= FALSE, eta, extra=NULL){
        shapes = eta2theta(eta, .link, earg= .earg )
        temp = if(is.R()) lbeta(shapes[,1], shapes[,2]) else
               lgamma(shapes[,1]) + lgamma(shapes[,2]) -
               lgamma(shapes[,1]+shapes[,2])
        if(residuals) stop("loglikelihood residuals not implemented yet") else
        sum(w * ((shapes[,1]-1) * log(y- .A) + (shapes[,2]-1) * log( .B -y) -
                 temp - (shapes[,1]+shapes[,2]-1) * log( .B - .A )))
    }, list( .link=link, .A=A, .B=B, .earg=earg ))),
    vfamily="betaff",
    deriv=eval(substitute(expression({
        shapes = eta2theta(eta, .link, earg= .earg )
        dshapes.deta = dtheta.deta(shapes, .link, earg= .earg )
        dl.dshapes = cbind(log(y-.A), log(.B-y)) - digamma(shapes) +
                     digamma(shapes[,1] + shapes[,2]) - log(.B - .A)
        w * dl.dshapes * dshapes.deta
    }), list( .link=link, .A=A, .B=B, .earg=earg ))),
    weight=expression({
        temp2 = trigamma(shapes[,1]+shapes[,2])
        d2l.dshape12 = temp2 - trigamma(shapes[,1])
        d2l.dshape22 = temp2 - trigamma(shapes[,2])
        d2l.dshape1shape2 = temp2

        wz = matrix(as.numeric(NA), n, dimm(M))   #3=dimm(M)
        wz[,iam(1,1,M)] = d2l.dshape12 * dshapes.deta[,1]^2
        wz[,iam(2,2,M)] = d2l.dshape22 * dshapes.deta[,2]^2
        wz[,iam(1,2,M)] = d2l.dshape1shape2 * dshapes.deta[,1] * dshapes.deta[,2]

        -w * wz
    }))
}



beta4 = function(link="loge", earg=list(),
                 i1=2.3, i2=2.4, iA=NULL, iB=NULL)
{



    if(mode(link) != "character" && mode(link) != "name")
        link = as.character(substitute(link))
    if(!is.list(earg)) earg = list()

    new("vglmff",
    blurb=c("Four-parameter Beta distribution\n",
            "(y-A)^(shape1-1) * (B-y)^(shape2-1), A < y < B \n\n",
            "Links:    ",
            namesof("shape1", link, earg=earg),  ", ",
            namesof("shape2", link, earg=earg), ", ",
            " A, B"),
    initialize=eval(substitute(expression({
        if(!is.vector(y) || (is.matrix(y) && ncol(y) != 1))
            stop("y must be a vector or a one-column matrix")

        if(length(.iA) && any(y < .iA))
            stop("initial A value out of range")
        if(length(.iB) && any(y > .iB))
            stop("initial B value out of range")

        predictors.names = c(
            namesof("shape1", .link, earg=.earg, short=TRUE),
            namesof("shape2", .link, earg=.earg, short=TRUE), "A", "B")
        my.range = diff(range(y))
        if(!length(etastart)) {
            etastart = cbind(shape1= rep(.i1, len=length(y)),
                             shape2= .i2,
                             A = if(length(.iA)) .iA else min(y)-my.range/70,
                             B = if(length(.iB)) .iB else max(y)+my.range/70)
        }
    }), list( .i1=i1, .i2=i2, .iA=iA, .iB=iB, .link=link, .earg=earg ))), 
    inverse=eval(substitute(function(eta, extra=NULL) {
        shapes = eta2theta(eta[,1:2], .link, earg=.earg)
        .A = eta[,3]
        .B = eta[,4]
        .A + (.B-.A) * shapes[,1] / (shapes[,1] + shapes[,2])
    }, list( .link=link, .earg=earg ))),
    last=eval(substitute(expression({
        misc$link = c(shape1 = .link, shape2 = .link, 
                      A="identity", B="identity")
        misc$earg = list(shape1 = .earg, shape2 = .earg, 
                         A=list(), B=list())
    }), list( .link=link, .earg=earg ))), 
    loglikelihood=eval(substitute(
        function(mu,y,w,residuals= FALSE,eta,extra=NULL) {
        shapes = eta2theta(eta[,1:2], .link, earg=.earg)
        .A = eta[,3]
        .B = eta[,4]
        temp = if(is.R()) lbeta(shapes[,1], shapes[,2]) else
               lgamma(shapes[,1]) + lgamma(shapes[,2]) -
               lgamma(shapes[,1]+shapes[,2])
        if(residuals) stop("loglikelihood residuals not implemented yet") else
        sum(w * ((shapes[,1]-1)*log(y-.A) + (shapes[,2]-1)*log(.B-y) - temp -
            (shapes[,1]+shapes[,2]-1)*log(.B-.A )))
    }, list( .link=link, .earg=earg ))), 
    vfamily="beta4",
    deriv=eval(substitute(expression({
        shapes = eta2theta(eta[,1:2], .link, earg=.earg)
        .A = eta[,3]
        .B = eta[,4]
        dshapes.deta = dtheta.deta(shapes, .link, earg=.earg)
        rr1 = (.B - .A)
        temp3 = (shapes[,1] + shapes[,2] - 1)
        temp1 = temp3 / rr1
        dl.dshapes = cbind(log(y-.A), log(.B-y)) - digamma(shapes) +
                     digamma(shapes[,1] + shapes[,2]) - log(.B - .A)
        dl.dA = -(shapes[,1]-1) / (y- .A)  + temp1
        dl.dB =  (shapes[,2]-1) / (.B - y) - temp1
        w * cbind(dl.dshapes * dshapes.deta, dl.dA, dl.dB)
    }), list( .link=link, .earg=earg ))), 
    weight=expression({

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


        -w * wz
    }))
}



simple.exponential = function()
{
    new("vglmff",
    blurb=c("Simple Exponential distribution\n",
            "Link:    log(rate)\n"),
    deviance= function(mu, y, w, residuals = FALSE, eta, extra=NULL) {
        devy = -log(y) - 1
        devmu = -log(mu) - y/mu
        devi = 2 * (devy - devmu)
        if(residuals) sign(y - mu) * sqrt(abs(devi) * w) else sum(w * devi)
    },
    initialize=expression({
        predictors.names = "log(rate)"
        mustart = y + (y == 0) / 8
    }),
    inverse=function(eta, extra=NULL)
        exp(-eta),
    link=function(mu, extra=NULL)
        -log(mu),
    vfamily="simple.exponential",
    deriv=expression({
        rate = 1 / mu
        dl.drate = mu - y
        drate.deta = dtheta.deta(rate, "loge")
        w * dl.drate * drate.deta
    }),
    weight=expression({
        ed2l.drate2 = -1 / rate^2
        wz = -w * drate.deta^2 * ed2l.drate2
        wz
    }))
}


exponential = function(link="loge", earg=list(), location=0, expected=TRUE)
{
    if(!is.Numeric(location, allow=1))
        stop("bad input for argument \"location\"")

    if(mode(link) != "character" && mode(link) != "name")
        link = as.character(substitute(link))
    if(!is.list(earg)) earg = list()

    new("vglmff",
    blurb=c("Exponential distribution\n\n",
            "Link:     ", namesof("rate", link, tag= TRUE), "\n",
            "Mean:     ", "mu =", location, "+ 1 / ",
            namesof("rate", link, tag= TRUE, earg=earg), "\n",
            "Variance: ",
            if(location==0) "Exponential: mu^2" else
            paste("(mu-", location, ")^2", sep="")),
    initialize=eval(substitute(expression({
        if(ncol(cbind(y)) != 1)
            stop("response must be a vector or a one-column matrix")
        extra$loc = .location   # This is passed into, e.g., link, deriv etc.
        if(any(y <= extra$loc))
            stop(paste("all responses must be greater than",extra$loc))
        predictors.names = namesof("rate", .link, tag=FALSE)
        mu = y + (y == extra$loc) / 8
        if(!length(etastart))
            etastart = theta2eta(1/(mu-extra$loc), .link, earg=.earg)
    }), list( .location=location, .link=link, .earg=earg ))),
    inverse=eval(substitute(function(eta, extra=NULL)
        extra$loc + 1 / eta2theta(eta, .link, earg=.earg),
    list( .link=link, .earg=earg ))),
    last=eval(substitute(expression({
        misc$location = extra$loc
        misc$link = c(rate = .link)
        misc$earg = list(rate = .earg)
    }), list( .link=link, .earg=earg ))),
    link=eval(substitute(function(mu, extra=NULL) 
        theta2eta(1/(mu-extra$loc), .link, earg=.earg),
    list( .link=link, .earg=earg ))),
    deviance=eval(substitute(
        function(mu,y,w,residuals= FALSE,eta,extra=NULL) {
        devy = -log(y - .location) - 1
        devmu = -log(mu - .location) - (y - .location)/(mu - .location)
        devi = 2 * (devy - devmu)
        if(residuals)
            sign(y - mu) * sqrt(abs(devi) * w) else 
            sum(w * devi)
    }, list( .location=location, .earg=earg ))),
    vfamily=c("exponential"),
    deriv=eval(substitute(expression({
        rate = 1 / (mu - extra$loc)
        dl.drate = mu - y
        drate.deta = dtheta.deta(rate, .link, earg=.earg)
        w * dl.drate * drate.deta
    }), list( .link=link, .earg=earg ))),
    weight=eval(substitute(expression({
        d2l.drate2 = - ((mu-extra$loc)^2)
        wz = -(drate.deta^2) * d2l.drate2
        if(! .expected) {
            # Use the observed info matrix rather than the expected
            d2rate.deta2 = d2theta.deta2(rate, .link, earg=.earg)
            wz = wz - dl.drate * d2rate.deta2
        }
        w * wz
    }), list( .link=link, .expected=expected, .earg=earg ))))
}




gamma1 = function(link="loge", earg=list())
{
    if(mode(link) != "character" && mode(link) != "name")
        link = as.character(substitute(link))
    if(!is.list(earg)) earg = list()

    new("vglmff",
    blurb=c("1-parameter Gamma distribution\n",
            "Link:     ",
            namesof("shape", link, earg=earg, tag= TRUE), "\n", 
            "Mean:       mu (=shape)\n",
            "Variance:   mu (=shape)"),
    initialize=eval(substitute(expression({
        if(any(y <= 0))
            stop("all responses must be positive")
        M = if(is.matrix(y)) ncol(y) else 1
        temp.names = if(M == 1) "shape" else paste("shape", 1:M, sep="")
        predictors.names = namesof(temp.names, .link, earg=.earg, short=TRUE)
        if(!length(etastart))
            etastart = cbind(theta2eta(y+1/8, .link, earg=.earg ))
    }), list( .link=link, .earg=earg ))), 
    inverse=eval(substitute(function(eta, extra=NULL)
        eta2theta(eta, .link, earg=.earg)),
    list( .link=link, .earg=earg )),
    last=eval(substitute(expression({
        temp.names = if(M == 1) "shape" else paste("shape", 1:M, sep="")
        misc$link = rep( .link, length=M)
        names(misc$link) = temp.names
        misc$earg = vector("list", M)
        names(misc$earg) = names(misc$link)
        for(ii in 1:M) misc$earg[[ii]] = .earg
        misc$expected = TRUE
    }), list( .link=link, .earg=earg ))),
    link=eval(substitute(function(mu, extra=NULL)
        theta2eta(mu, .link, earg=.earg)),
    list( .link=link, .earg=earg )),
    loglikelihood= function(mu, y, w, residuals = FALSE, eta, extra=NULL)
        if(residuals) stop("loglikelihood residuals not implemented yet") else
        sum(w * ((mu-1)*log(y) - y - lgamma(mu))),
    vfamily=c("gamma1"),
    deriv=eval(substitute(expression({
        shape = mu
        dl.dshape = log(y) - digamma(shape)
        dshape.deta = dtheta.deta(shape, .link, earg=.earg)
        w * dl.dshape * dshape.deta
    }), list( .link=link, .earg=earg ))),
    weight=expression({
        d2l.dshape = -trigamma(shape)
        wz = -(dshape.deta^2) * d2l.dshape
        w * wz
    }))
}


gamma2.ab = function(lrate="loge", lshape="loge",
                     erate=list(), eshape=list(),
                     irate=NULL, ishape=NULL, expected=TRUE, zero=2)
{
    if(mode(lrate) != "character" && mode(lrate) != "name")
        lrate = as.character(substitute(lrate))
    if(mode(lshape) != "character" && mode(lshape) != "name")
        lshape = as.character(substitute(lshape))
    if(length( irate) && !is.Numeric(irate, posit=TRUE))
        stop("bad input for argument \"irate\"")
    if(length( ishape) && !is.Numeric(ishape, posit=TRUE))
        stop("bad input for argument \"ishape\"")
    if(length(zero) && !is.Numeric(zero, integer=TRUE, posit=TRUE))
        stop("bad input for argument \"zero\"")
    if(!is.logical(expected) || length(expected) != 1)
        stop("bad input for argument \"expected\"")
    if(!is.list(erate)) erate = list()
    if(!is.list(eshape)) eshape = list()

    new("vglmff",
    blurb=c("2-parameter Gamma distribution\n",
            "Links:    ",
            namesof("rate", lrate, earg=erate), ", ", 
            namesof("shape", lshape, earg=eshape), "\n",
            "Mean:     mu = shape/rate\n",
            "Variance: (mu^2)/shape = shape/rate^2"),
    constraints=eval(substitute(expression({
        constraints = cm.zero.vgam(constraints, x, .zero, M)
    }), list( .zero=zero ))),
    initialize=eval(substitute(expression({
        # Error check
        if(ncol(cbind(y)) != 1)
            stop("response must be a vector or a one-column matrix")
        if(any(y <= 0))
            stop("all responses must be positive")
        predictors.names = c(namesof("rate",  .lrate, earg=.erate,  tag=FALSE),
                             namesof("shape", .lshape, earg=.eshape, tag=FALSE))
        if(!length(etastart)) {
            mymu = y + 0.167 * (y == 0)
            junk = lsfit(x, y, wt = w, intercept = FALSE)
            var.y.est = sum(w * junk$resid^2) / (nrow(x) - length(junk$coef))
            init.shape =  if(length( .ishape)) .ishape else mymu^2 / var.y.est
            init.rate =  if(length( .irate)) .irate else init.shape / mymu
            init.rate = rep(init.rate, len=n)
            init.shape = rep(init.shape, len=n)
            if( .lshape == "loglog")
                init.shape[init.shape <= 1] = 3.1 #Hopefully value is big enough
            etastart = cbind(theta2eta(init.rate, .lrate, earg=.erate),
                             theta2eta(init.shape, .lshape, earg=.eshape))
        }
    }), list( .lrate=lrate, .lshape=lshape, .irate=irate, .ishape=ishape,
              .erate=erate, .eshape=eshape ))),
    inverse=eval(substitute(function(eta, extra=NULL) {
        eta2theta(eta[,2], .lshape, earg=.eshape) / eta2theta(eta[,1], .lrate,
        earg=.erate)
    }, list( .lrate=lrate, .lshape=lshape,
             .erate=erate, .eshape=eshape ))),
    last=eval(substitute(expression({
        misc$link = c(rate= .lrate, shape= .lshape)
        misc$earg = list(rate= .erate, shape= .eshape)
    }), list( .lrate=lrate, .lshape=lshape,
              .erate=erate, .eshape=eshape ))),
    loglikelihood=eval(substitute(
        function(mu, y, w, residuals= FALSE,eta, extra=NULL) {
        rate = eta2theta(eta[,1], .lrate, earg=.erate)
        shape = eta2theta(eta[,2], .lshape, earg=.eshape)
        if(residuals) stop("loglikelihood residuals not implemented yet") else
        sum(w*(-rate * y + (shape-1)*log(y) + shape*log(rate) - lgamma(shape )))
    }, list( .lrate=lrate, .lshape=lshape,
             .erate=erate, .eshape=eshape ))),
    vfamily=c("gamma2.ab"),
    deriv=eval(substitute(expression({
        rate = eta2theta(eta[,1], .lrate, earg=.erate)
        shape = eta2theta(eta[,2], .lshape, earg=.eshape)
        dl.drate = mu - y
        dl.dshape = log(y*rate) - digamma(shape)
        dratedeta = dtheta.deta(rate, .lrate, earg=.erate)
        dshape.deta = dtheta.deta(shape, .lshape, earg=.eshape)
        w * cbind(dl.drate * dratedeta, dl.dshape * dshape.deta)
    }), list( .lrate=lrate, .lshape=lshape,
              .erate=erate, .eshape=eshape ))),
    weight=eval(substitute(expression({
        d2l.dshape2 = -trigamma(shape)
        d2l.drate2 = -shape/(rate^2)
        d2l.drateshape = 1/rate
        wz = matrix(as.numeric(NA), n, dimm(M))  #3=dimm(M)
        wz[,iam(1,1,M)] = -d2l.drate2 * dratedeta^2
        wz[,iam(2,2,M)] = -d2l.dshape2 * dshape.deta^2
        wz[,iam(1,2,M)] = -d2l.drateshape * dratedeta * dshape.deta
        if(! .expected) {
            d2ratedeta2 = d2theta.deta2(rate, .lrate, earg=.erate)
            d2shapedeta2 = d2theta.deta2(shape, .lshape, earg=.eshape)
            wz[,iam(1,1,M)] = wz[,iam(1,1,M)] - dl.drate * d2ratedeta2
            wz[,iam(2,2,M)] = wz[,iam(2,2,M)] - dl.dshape * d2shapedeta2
        }
        w * wz
    }), list( .lrate=lrate, .lshape=lshape,
              .erate=erate, .eshape=eshape, .expected=expected ))))
}



gamma2 = function(lmu="loge", lshape="loge",
                  emu=list(), eshape=list(),
                  method.init=1,
                  deviance.arg=FALSE, ishape=NULL, zero=-2)
{
    if(mode(lmu) != "character" && mode(lmu) != "name")
        lmu = as.character(substitute(lmu))
    if(mode(lshape) != "character" && mode(lshape) != "name")
        lshape = as.character(substitute(lshape))
    if(length(zero) && !is.Numeric(zero, integer=TRUE))
        stop("bad input for argument \"zero\"")
    if(length( ishape) && !is.Numeric(ishape, posit=TRUE))
        stop("bad input for argument \"ishape\"")
    if(!is.Numeric(method.init, allow=1, integ=TRUE, posit=TRUE) ||
       method.init > 2)
        stop("'method.init' must be 1 or 2")
    if(!is.list(emu)) emu = list()
    if(!is.list(eshape)) eshape = list()

    ans = 
    new("vglmff",
    blurb=c("2-parameter Gamma distribution",
            " (McCullagh and Nelder 1989 parameterization)\n",
            "Links:    ",
            namesof("mu", lmu, earg=emu), ", ", 
            namesof("shape", lshape, earg=eshape), "\n",
            "Mean:     mu\n",
            "Variance: (mu^2)/shape"),
    constraints=eval(substitute(expression({
        temp752 = .zero
        if(length(temp752) && all(temp752 == -2))
            temp752 = 2*(1:ncol(y))
        constraints = cm.zero.vgam(constraints, x, temp752, M)
    }), list( .zero=zero ))),
    initialize=eval(substitute(expression({
        if(is.R()) assign("CQO.FastAlgorithm",
            ( .lmu == "loge" && .lshape == "loge"), envir = VGAMenv) else
            CQO.FastAlgorithm <<- ( .lmu == "loge" && .lshape == "loge")
        if(any(function.name == c("cqo","cao")) &&
           is.Numeric( .zero, allow=1) && .zero != -2)
            stop("argument zero=-2 is required")

        y = as.matrix(y)
        M = 2 * ncol(y)
        NOS = ncoly = ncol(y)  # Number of species
        temp1.names = if(NOS==1) "mu" else paste("mu", 1:NOS, sep="")
        temp2.names = if(NOS==1) "shape" else paste("shape", 1:NOS, sep="")
        predictors.names =
            c(namesof(temp1.names, .lmu, earg=.emu, tag=FALSE),
              namesof(temp2.names, .lshape, earg=.eshape, tag=FALSE))
        predictors.names = predictors.names[interleave.VGAM(M, M=2)]


        # Error check
        if(any(y <= 0))
            stop("all responses must be positive") # see @loglikelihood
        if(!length(etastart)) {
            init.shape = matrix(1.0, n, NOS)
            mymu = y # + 0.167 * (y == 0)  # method.init == 1 (the default)
            if( .method.init == 2) {
                for(ii in 1:ncol(y)) {
                    mymu[,ii] = weighted.mean(y[,ii], w=w)
                }
            }
            for(spp in 1:NOS) {
                junk = lsfit(x, y[,spp], wt = w, intercept = FALSE)
                var.y.est = sum(w * junk$resid^2) / (n - length(junk$coef))
                init.shape[,spp] = if(length( .ishape)) .ishape else
                    mymu[,spp]^2 / var.y.est
                if( .lshape == "loglog") init.shape[init.shape[,spp] <=
                             1,spp] = 3.1 # Hopefully value is big enough
            }
            etastart = cbind(theta2eta(mymu, .lmu, earg=.emu ),
                             theta2eta(init.shape, .lshape, earg=.eshape ))
            etastart = etastart[,interleave.VGAM(M, M=2),drop=FALSE]
        }
    }), list( .lmu=lmu, .lshape=lshape, .ishape=ishape, .zero=zero,
              .emu=emu, .eshape=eshape,
              .method.init=method.init ))),
    inverse=eval(substitute(function(eta, extra=NULL) {
        NOS = ncol(eta) / 2
        eta2theta(eta[,2*(1:NOS)-1,drop=FALSE], .lmu, earg=.emu )
    }, list( .lmu=lmu, .emu=emu ))),
    last=eval(substitute(expression({
       if(is.R()) {
            if(exists("CQO.FastAlgorithm", envir = VGAMenv))
                rm("CQO.FastAlgorithm", envir = VGAMenv)
        } else {
            while(exists("CQO.FastAlgorithm"))
                remove("CQO.FastAlgorithm")
        }
        tmp34 = c(rep( .lmu, length=NOS), rep( .lshape, length=NOS))
        names(tmp34) = c(if(NOS==1) "mu" else paste("mu", 1:NOS, sep=""), 
                         if(NOS==1) "shape" else paste("shape", 1:NOS, sep=""))
        tmp34 = tmp34[interleave.VGAM(M, M=2)]
        misc$link = tmp34 # Already named
        misc$earg = vector("list", M)
        names(misc$earg) = names(misc$link)
        for(ii in 1:NOS) {
            misc$earg[[2*ii-1]] = .emu
            misc$earg[[2*ii  ]] = .eshape
        }
        misc$expected = TRUE
    }), list( .lmu=lmu, .lshape=lshape,
              .emu=emu, .eshape=eshape ))),
    link=eval(substitute(function(mu, extra=NULL) {
        temp = theta2eta(mu, .lmu, earg=.emu )
        temp = cbind(temp, NA * temp)
        temp[,interleave.VGAM(ncol(temp), M=2),drop=FALSE]
    }, list( .lmu=lmu, .emu=emu ))),
    loglikelihood=eval(substitute(
        function(mu, y, w, residuals= FALSE,eta, extra=NULL) {
        NOS = ncol(eta) / 2
        mymu = mu  # eta2theta(eta[,2*(1:NOS)-1], .lmu, earg=.emu )
        shapemat = eta2theta(eta[,2*(1:NOS),drop=FALSE], .lshape, earg=.eshape )
        if(residuals) stop("loglikelihood residuals not implemented yet") else
        sum(w*((shapemat - 1) * log(y) + shapemat *
                (log(shapemat) - y / mymu - log(mymu)) - lgamma(shapemat )))
    }, list( .lmu=lmu, .lshape=lshape,
             .emu=emu, .eshape=eshape ))),
    vfamily=c("gamma2"),
    deriv=eval(substitute(expression({
        NOS = ncol(eta) / 2
        mymu = eta2theta(eta[,2*(1:NOS)-1], .lmu, earg=.emu )
        shape = eta2theta(eta[,2*(1:NOS)], .lshape, earg=.eshape )
        dl.dmu = shape * (y / mymu - 1) / mymu
        dl.dshape = log(y) + log(shape) - log(mymu) + 1 - digamma(shape) -
                    y / mymu
        dmu.deta = dtheta.deta(mymu, .lmu, earg=.emu )
        dshape.deta = dtheta.deta(shape, .lshape, earg=.eshape )
        myderiv = w * cbind(dl.dmu * dmu.deta, dl.dshape * dshape.deta)
        myderiv[,interleave.VGAM(M, M=2)]
    }), list( .lmu=lmu, .lshape=lshape,
              .emu=emu, .eshape=eshape ))),
    weight=eval(substitute(expression({
        ed2l.dmu2 = shape / (mymu^2)
        ed2l.dshape2 = trigamma(shape) - 1 / shape
        wz = matrix(as.numeric(NA), n, M)  #2=M; diagonal!
        NOS = M / 2 
        wz[,2*(1:NOS)-1] = ed2l.dmu2 * dmu.deta^2
        wz[,2*(1:NOS)] = ed2l.dshape2 * dshape.deta^2
        w * wz
    }), list( .lmu=lmu ))))

    if(deviance.arg) ans@deviance=eval(substitute(
        function(mu,y,w,residuals= FALSE,eta,extra=NULL) {
        NOS = ncol(eta) / 2
        temp300 =  eta[,2*(1:NOS),drop=FALSE]
        if( .lshape == "loge") {
            bigval = 28
            temp300[temp300 >  bigval] =  bigval
            temp300[temp300 < -bigval] = -bigval
        } else stop("can only handle the 'loge' link")
        shape =  eta2theta(temp300, .lshape, earg=.eshape )
        devi = -2 * (log(y/mu) - y/mu + 1)
        if(residuals) {
           warning("not 100% sure about these deviance residuals!")
           sign(y - mu) * sqrt(abs(devi) * w)
        } else
           sum(w * devi)
    }, list( .lshape=lshape )))
    ans
}


geometric =function(link="logit", earg=list(), expected=TRUE)
{
    if(!is.logical(expected) || length(expected) != 1)
        stop("bad input for argument \"expected\"")
    if(mode(link) != "character" && mode(link) != "name")
        link = as.character(substitute(link))
    if(!is.list(earg)) earg = list()

    new("vglmff",
    blurb=c("Geometric distribution (P[Y=y] = prob*(1-prob)^y, y=0,1,2,...)\n",
            "Link:     ",
            namesof("prob", link, earg=earg), "\n",
            "Mean:     mu = (1-prob)/prob\n",
            "Variance: mu*(1+mu) = (1-prob)/prob^2"),
    initialize=eval(substitute(expression({
        if(ncol(cbind(y)) != 1)
            stop("response must be a vector or a 1-column matrix")
        if(any(y < 0)) stop("all responses must be >= 0")
        if(any(y!=round(y ))) stop("response should be integer-valued")
        predictors.names = namesof("prob", .link, earg=.earg, tag=FALSE)
        if(!length(etastart)) {
            prob.init = 1 / (1 + y + 1/16)
            etastart = theta2eta(prob.init, .link, earg= .earg)
        }
    }), list( .link=link, .earg=earg ))),
    inverse=eval(substitute(function(eta, extra=NULL) {
        prob = eta2theta(eta, .link, earg= .earg)
        (1-prob)/prob 
    }, list( .link=link, .earg=earg ))),
    last=eval(substitute(expression({
        misc$link = c(prob= .link)
        misc$earg = list(prob= .earg )
        misc$expected = .expected
    }), list( .link=link, .earg=earg, .expected=expected ))),
    loglikelihood=eval(substitute(
        function(mu, y, w, residuals= FALSE,eta, extra=NULL) {
        prob = eta2theta(eta, .link, earg= .earg)
        if(residuals) stop("loglikelihood residuals not implemented yet") else {
            if(is.R()) sum(w * dgeom(x=y, prob=prob, log=TRUE)) else
            sum(w*(y * log1p(-prob) + log(prob )))
        }
    }, list( .link=link, .earg=earg ))),
    vfamily=c("geometric"),
    deriv=eval(substitute(expression({
        prob = eta2theta(eta, .link, earg= .earg)
        dl.dprob = -y/(1-prob) + 1/prob 
        dprobdeta = dtheta.deta(prob, .link, earg= .earg)
        w * cbind(dl.dprob * dprobdeta)
    }), list( .link=link, .earg=earg, .expected=expected ))),
    weight=eval(substitute(expression({
        ed2l.dprob2 = if( .expected ) 1 / (prob^2 * (1-prob)) else
            y / (1-prob)^2 + 1 / prob^2
        wz = ed2l.dprob2 * dprobdeta^2
        if( !( .expected )) wz = wz - dl.dprob * d2theta.deta2(prob, .link, earg= .earg)
        w * wz
    }), list( .link=link, .earg=earg, .expected=expected ))))
}


dbetageom = function(x, shape1, shape2, log=FALSE) {
    if(!is.Numeric(x)) stop("bad input for argument \"x\"")
    if(!is.Numeric(shape1, pos=TRUE)) stop("bad input for argument \"shape1\"")
    if(!is.Numeric(shape2, pos=TRUE)) stop("bad input for argument \"shape2\"")
    N = max(length(x), length(shape1), length(shape2))
    x = rep(x, len=N); shape1 = rep(shape1, len=N); shape2 = rep(shape2, len=N)
    if(!is.R()) {
        beta  = function(x,y) gamma(x) * gamma(y) / gamma(x+y)
        lbeta = function(x,y) lgamma(x) + lgamma(y) - lgamma(x+y)
    }
    ans = if(log) lbeta(1+shape1, shape2+abs(x)) - lbeta(shape1, shape2) else
          beta(1+shape1, shape2+abs(x)) / beta(shape1, shape2)
    ifelse(x == round(x) & x >= 0, ans, if(log) -Inf else 0)
}


pbetageom = function(q, shape1, shape2, log.p=FALSE) {
    if(!is.Numeric(q)) stop("bad input for argument \"q\"")
    if(!is.Numeric(shape1, pos=TRUE)) stop("bad input for argument \"shape1\"")
    if(!is.Numeric(shape2, pos=TRUE)) stop("bad input for argument \"shape2\"")
    N = max(length(q), length(shape1), length(shape2))
    q = rep(q, len=N); shape1 = rep(shape1, len=N); shape2 = rep(shape2, len=N)
    ans = q * 0  # Retains names(q)
    if(max(abs(shape1-shape1[1])) < 1.0e-08 &&
       max(abs(shape2-shape2[1])) < 1.0e-08) {
        qstar = floor(q)
        temp = if(max(qstar) >= 0) dbetageom(x=0:max(qstar), 
               shape1=shape1[1], shape2=shape2[1]) else 0*qstar
        unq = unique(qstar)
        for(i in unq) {
            index = qstar == i
            ans[index] = if(i >= 0) sum(temp[1:(1+i)]) else 0
        }
    } else
    for(i in 1:N) {
        qstar = floor(q[i])
        ans[i] = if(qstar >= 0) sum(dbetageom(x=0:qstar, 
                 shape1=shape1[i], shape2=shape2[i])) else 0
    }
    if(log.p) log(ans) else ans
}

rbetageom = function(n, shape1, shape2) {
    if(!is.Numeric(n, integ=TRUE,allow=1)) stop("bad input for argument \"n\"")
    if(!is.Numeric(shape1, pos=TRUE)) stop("bad input for argument \"shape1\"")
    if(!is.Numeric(shape2, pos=TRUE)) stop("bad input for argument \"shape2\"")
    rgeom(n=n, prob = rbeta(n=n, shape1=shape1, shape2=shape2))
}




tobit = function(Lower = 0, Upper = Inf, lmu="identity",
                 lsd="loge", emu=list(), esd=list(), imethod=1, zero=2)
{
    if(mode(lmu) != "character" && mode(lmu) != "name")
        lmu = as.character(substitute(lmu))
    if(mode(lsd) != "character" && mode(lsd) != "name")
        lsd = as.character(substitute(lsd))
    if(!is.Numeric(imethod, allow=1, integer=TRUE, positi=TRUE) || imethod > 2)
        stop("imethod must be 1 or 2")
    if(length(Lower) != 1 || length(Upper) != 1 ||
       !is.numeric(Lower) || !is.numeric(Upper) || Lower >= Upper)
    stop("Lower and Upper must have length 1 and be numeric with Lower < Upper")
    if(length(zero) && !is.Numeric(zero, integer=TRUE, posit=TRUE))
        stop("bad input for argument \"zero\"")
    if(!is.list(emu)) emu = list()
    if(!is.list(esd)) esd = list()

    new("vglmff",
    blurb=c("Tobit model\n\n",
            "Links:    ", namesof("mu", lmu, earg=emu, tag= TRUE), "; ",
                          namesof("sd", lsd, earg=esd, tag= TRUE), "\n",
            "Conditional variance: sd^2"),
    constraints=eval(substitute(expression({
        constraints = cm.zero.vgam(constraints, x, .zero, M)
    }), list( .zero=zero ))),
    initialize=eval(substitute(expression({
        y = cbind(y)
        if(ncol(y)!=1)stop("the response must be a vector or a 1-column matrix")
        extra$censoredL = (y <= .Lower)
        extra$censoredU = (y >= .Upper)
        if(min(y) < .Lower) {
             warning(paste("replacing response values less than the value ",
                           .Lower, "by", .Lower))
             y[y < .Lower] = .Lower
        }
        if(max(y) > .Upper) {
             warning(paste("replacing response values greater than the value", 
                            .Upper, "by", .Upper))
             y[y > .Upper] = .Upper
        }
        predictors.names = c(namesof("mu", .lmu, earg=.emu, tag=FALSE),
                             namesof("sd", .lsd, earg=.esd, tag=FALSE))
        if(!length(etastart)) {
            anyc = extra$censoredL | extra$censoredU
            i11 = if( .imethod == 1) anyc else FALSE  # can be all data
            junk=if(is.R()) lm.wfit(x=cbind(x[!i11,]),y=y[!i11],w=w[!i11]) else
                   lm.wfit(x=cbind(x[!i11,]), y=y[!i11], w=w[!i11],method="qr")
            sd.y.est = sqrt( sum(w[!i11] * junk$resid^2) / junk$df.residual )
            etastart = cbind(mu=y, rep(theta2eta(sd.y.est, .lsd, earg= .esd), length=n))
            if(any(anyc)) etastart[anyc,1] = x[anyc,,drop=FALSE] %*% junk$coeff
        }
   }), list( .Lower=Lower, .Upper=Upper, .lmu=lmu, .lsd=lsd,
             .emu=emu, .esd=esd, .imethod=imethod ))),
    inverse=eval(substitute( function(eta, extra=NULL) {
        eta2theta(eta[,1], .lmu, earg= .emu)
    }, list( .lmu=lmu, .emu=emu ))),
    last=eval(substitute(expression({
        misc$link = c("mu"= .lmu, "sd"= .lsd)
        misc$earg = list("mu"= .emu, "sd"= .esd)
        misc$expected = TRUE
        misc$Lower = .Lower
        misc$Upper = .Upper
    }), list( .lmu=lmu, .lsd=lsd,
              .emu=emu, .esd=esd,
              .Lower=Lower, .Upper=Upper ))),
    loglikelihood=eval(substitute(
        function(mu,y,w,residuals= FALSE,eta, extra=NULL) {
        cenL = extra$censoredL
        cenU = extra$censoredU
        cen0 = !cenL & !cenU   # uncensored obsns
        mum = eta2theta(eta[,1], .lmu, earg= .emu)
        sd = eta2theta(eta[,2], .lsd, earg= .esd)
        ell1 = -log(sd[cen0]) - 0.5 * ((y[cen0] - mum[cen0])/sd[cen0])^2
        ell2 = log1p(-pnorm((mum[cenL] - .Lower)/sd[cenL]))
        ell3 = log1p(-pnorm(( .Upper -  mum[cenU])/sd[cenU]))
        if(residuals) stop("loglikelihood residuals not implemented yet") else
        sum(w[cen0] * ell1) + sum(w[cenL] * ell2) + sum(w[cenU] * ell3)
    }, list( .lmu=lmu, .lsd=lsd,
             .emu=emu, .esd=esd,
             .Lower=Lower, .Upper=Upper ))),
    vfamily=c("tobit"),
    deriv=eval(substitute(expression({
        cenL = extra$censoredL
        cenU = extra$censoredU
        cen0 = !cenL & !cenU   # uncensored obsns
        mum = eta2theta(eta[,1], .lmu, earg= .emu)
        sd = eta2theta(eta[,2], .lsd, earg= .esd)
        dl.dmu = (y-mum) / sd^2
        dl.dsd = (((y-mum)/sd)^2 - 1) / sd
        dmu.deta = dtheta.deta(mum, .lmu, earg= .emu)
        dsd.deta = dtheta.deta(sd, .lsd, earg= .esd)
        if(any(cenL)) {
            mumL = mum - .Lower
            temp21L = mumL[cenL] / sd[cenL]
            PhiL = pnorm(temp21L)
            phiL = dnorm(temp21L)
            fred21 = phiL / (1 - PhiL)
            dl.dmu[cenL] = -fred21 / sd[cenL]
            dl.dsd[cenL] = mumL[cenL] * fred21 / sd[cenL]^2
            rm(fred21)
        }
        if(any(cenU)) {
            mumU = .Upper - mum
            temp21U = mumU[cenU] / sd[cenU]
            PhiU = pnorm(temp21U)
            phiU = dnorm(temp21U)
            fred21 = phiU / (1 - PhiU)
            dl.dmu[cenU] = fred21 / sd[cenU]   # Negated
            dl.dsd[cenU] = mumU[cenU] * fred21 / sd[cenU]^2
            rm(fred21)
        }
        w * cbind(dl.dmu * dmu.deta, dl.dsd * dsd.deta)
    }), list( .lmu=lmu, .lsd=lsd,
              .emu=emu, .esd=esd,
              .Lower=Lower, .Upper=Upper ))),
    weight=eval(substitute(expression({
        A1 = 1 - pnorm((mum - .Lower) / sd)   # Lower
        A3 = 1 - pnorm(( .Upper - mum) / sd)  # Upper
        A2 = 1 - A1 - A3                      # Middle; uncensored
        wz = matrix(0, n, 3)
        wz[,iam(1,1,M)] = A2 * 1 / sd^2  # ed2l.dmu2
        wz[,iam(2,2,M)] = A2 * 2 / sd^2  # ed2l.dsd2
        mumL = mum - .Lower
        temp21L = mumL / sd
        PhiL = pnorm(temp21L)
        phiL = dnorm(temp21L)
        temp31L = ((1-PhiL) * sd)^2 
        wz.cenL11 = phiL * (phiL - (1-PhiL)*temp21L) / temp31L
        wz.cenL22 = mumL * phiL * ((1-PhiL) * (2 - temp21L^2) +
                    mumL * phiL / sd) / (sd * temp31L)
        wz.cenL12 = phiL * ((1-PhiL)*(temp21L^2 - 1) - temp21L*phiL) / temp31L
        wz.cenL11[!is.finite(wz.cenL11)] = 0
        wz.cenL22[!is.finite(wz.cenL22)] = 0
        wz.cenL12[!is.finite(wz.cenL12)] = 0
        wz[,iam(1,1,M)] = wz[,iam(1,1,M)] + A1 * wz.cenL11
        wz[,iam(2,2,M)] = wz[,iam(2,2,M)] + A1 * wz.cenL22
        wz[,iam(1,2,M)] = A1 * wz.cenL12
        mumU = .Upper - mum    # often Inf
        temp21U = mumU / sd    # often Inf
        PhiU = pnorm(temp21U)  # often 1
        phiU = dnorm(temp21U)  # often 0
        temp31U = ((1-PhiU) * sd)^2  # often 0
        tmp8 = (1-PhiU)*temp21U
        wzcenU11 = phiU * (phiU - tmp8) / temp31U
        tmp9 = (1-PhiU) * (2 - temp21U^2)
        wzcenU22 = mumU * phiU * (tmp9 + mumU * phiU / sd) / (sd * temp31U)
        wzcenU12 = -phiU * ((1-PhiU)*(temp21U^2 - 1) - temp21U*phiU) / temp31U
        wzcenU11[!is.finite(wzcenU11)] = 0  # Needed when .Upper==Inf
        wzcenU22[!is.finite(wzcenU22)] = 0  # Needed when .Upper==Inf
        wzcenU12[!is.finite(wzcenU12)] = 0  # Needed when .Upper==Inf
        wz[,iam(1,1,M)] = wz[,iam(1,1,M)] + A3 * wzcenU11
        wz[,iam(2,2,M)] = wz[,iam(2,2,M)] + A3 * wzcenU22
        wz[,iam(1,2,M)] = wz[,iam(1,2,M)] + A3 * wzcenU12
        wz[,iam(1,1,M)] = w * wz[,iam(1,1,M)] * dmu.deta^2
        wz[,iam(2,2,M)] = w * wz[,iam(2,2,M)] * dsd.deta^2
        wz[,iam(1,2,M)] = w * wz[,iam(1,2,M)] * dmu.deta * dsd.deta
        wz
    }), list( .lmu=lmu, .Lower=Lower, .Upper=Upper, .lsd=lsd ))))
}



normal1 = function(lmean="identity", lsd="loge",
                   emean=list(), esd=list(), zero=NULL)
{

    if(mode(lmean) != "character" && mode(lmean) != "name")
        lmean = as.character(substitute(lmean))
    if(mode(lsd) != "character" && mode(lsd) != "name")
        lsd = as.character(substitute(lsd))
    if(length(zero) && !is.Numeric(zero, integer=TRUE, posit=TRUE))
        stop("bad input for argument \"zero\"")
    if(!is.list(emean)) emean = list()
    if(!is.list(esd)) esd = list()

    new("vglmff",
    blurb=c("Univariate Normal distribution\n\n",
            "Links:    ",
            namesof("mean", lmean, earg=emean, tag= TRUE), "; ",
            namesof("sd", lsd, earg=esd, tag= TRUE),
            "\n",
            "Variance: sd^2"),
    constraints=eval(substitute(expression({
        constraints = cm.zero.vgam(constraints, x, .zero, M)
    }), list( .zero=zero ))),
    initialize=eval(substitute(expression({
        predictors.names = c(namesof("mean", .lmean, earg=.emean, tag=FALSE),
                             namesof("sd",   .lsd, earg=.esd, tag=FALSE))
        if(ncol(y <- cbind(y)) != 1)
            stop("response must be a vector or a one-column matrix")
        if(!length(etastart)) {
            junk = if(is.R()) lm.wfit(x=x, y=y, w=w) else 
                              lm.wfit(x=x, y=y, w=w, method="qr")
            sd.y.est = sqrt( sum(w * junk$resid^2) / junk$df.residual )
            mean.init = if( .lmean == "loge") pmax(1/1024, y) else y
            etastart = cbind(theta2eta(mean.init, .lmean, earg= .emean),
                             theta2eta(sd.y.est, .lsd, earg= .esd))
        }
    }), list( .lmean=lmean, .lsd=lsd, .emean=emean, .esd=esd ))),
    inverse=eval(substitute(function(eta, extra=NULL) {
        eta2theta(eta[,1], .lmean, earg= .emean)
    }, list( .lmean=lmean, .emean=emean, .esd=esd ))),
    last=eval(substitute(expression({
        misc$link = c("mu"= .lmean, "sd"= .lsd)
        misc$earg = list("mu"= .emean, "sd"= .esd)
        misc$expected = TRUE
    }), list( .lmean=lmean, .lsd=lsd, .emean=emean, .esd=esd ))),
    loglikelihood=eval(substitute(
        function(mu,y,w,residuals= FALSE,eta, extra=NULL) {
        sd = eta2theta(eta[,2], .lsd, earg= .esd)
        if(residuals) stop("loglikelihood residuals not implemented yet") else {
        if(is.R())
            sum(w*dnorm(y, m=mu, sd=sd, log=TRUE)) else
            sum(w * (-log(sd*sqrt(2*pi)) - 0.5 * ((y - mu)/sd)^2))
        }
    }, list( .lsd=lsd, .emean=emean, .esd=esd ))),
    vfamily=c("normal1"),
    deriv=eval(substitute(expression({
        mymu = eta2theta(eta[,1], .lmean, earg= .emean)
        sd = eta2theta(eta[,2], .lsd, earg= .esd)
        dl.dmu = (y-mymu) / sd^2
        dl.dsd = -1/sd + (y-mymu)^2 / sd^3
        dmu.deta = dtheta.deta(mymu, .lmean, earg= .emean)
        dsd.deta = dtheta.deta(sd, .lsd, earg= .esd)
        cbind(w * dl.dmu * dmu.deta, w * dl.dsd * dsd.deta)
    }), list( .lmean=lmean, .lsd=lsd, .emean=emean, .esd=esd ))),
    weight=expression({
        wz = matrix(as.numeric(NA), n, 2)  # diagonal matrix; y is one-column too
        ed2l.dmu2 = -1 / sd^2
        ed2l.dsd2 = -2 / sd^2    # zz; replace 2 by 0.5 ??
        wz[,iam(1,1,M)] = -w * ed2l.dmu2 * dmu.deta^2
        wz[,iam(2,2,M)] = -w * ed2l.dsd2 * dsd.deta^2
        wz
    }))
}





lognormal = function(lmeanlog="identity", lsdlog="loge",
                     emeanlog=list(), esdlog=list(), zero=NULL)
{
    if(mode(lmeanlog) != "character" && mode(lmeanlog) != "name")
        lmeanlog = as.character(substitute(lmeanlog))
    if(mode(lsdlog) != "character" && mode(lsdlog) != "name")
        lsdlog = as.character(substitute(lsdlog))
    if(length(zero) && (!is.Numeric(zero, integer=TRUE, posit=TRUE) ||
                        zero > 2)) stop("bad input for argument argument \"zero\"")
    if(!is.list(emeanlog)) emeanlog = list()
    if(!is.list(esdlog)) esdlog = list()

    new("vglmff",
    blurb=c("Two-parameter (univariate) lognormal distribution\n\n",
            "Links:    ", namesof("meanlog", lmeanlog, earg=emeanlog, tag= TRUE), ", ",
                          namesof("sdlog", lsdlog, earg=esdlog, tag= TRUE)),
    constraints=eval(substitute(expression({
        constraints = cm.zero.vgam(constraints, x, .zero, M)
    }), list( .zero=zero ))),
    initialize=eval(substitute(expression({
        if(ncol(cbind(y)) != 1)
            stop("response must be a vector or a one-column matrix")
        if(min(y) <= 0) stop("response must be positive")
        predictors.names = c(namesof("meanlog", .lmeanlog, earg=.emeanlog, tag=FALSE),
                             namesof("sdlog", .lsdlog, earg=.esdlog, tag=FALSE))
        if(!length(etastart)) {
            junk = if(is.R()) lm.wfit(x=x, y=log(y), w=w) else 
                              lm.wfit(x=x, y=log(y), w=w, method="qr")
            sdlog.y.est = sqrt( sum(w * junk$resid^2) / junk$df.residual )
            etastart = cbind(
            meanlog= rep(theta2eta(log(median(y)), .lmeanlog, earg= .emeanlog), length=n),
            sdlog= rep(theta2eta(sdlog.y.est, .lsdlog, earg= .esdlog), length=n))
        }
    }), list( .lmeanlog = lmeanlog, .lsdlog=lsdlog,
              .emeanlog = emeanlog, .esdlog=esdlog ))),
    inverse=eval(substitute(function(eta, extra=NULL) {
        mulog = eta2theta(eta[,1], .lmeanlog, earg= .emeanlog)
        sdlog = eta2theta(eta[,2], .lsdlog, earg= .esdlog)
        exp(mulog + 0.5 * sdlog^2)
    }, list( .lmeanlog = lmeanlog, .lsdlog=lsdlog,
              .emeanlog = emeanlog, .esdlog=esdlog ))),
    last=eval(substitute(expression({
        misc$link = c("meanlog"= .lmeanlog, "sdlog"= .lsdlog)
        misc$earg = list("meanlog"= .emeanlog, "sdlog"= .esdlog)
        misc$expected = TRUE
    }), list( .lmeanlog = lmeanlog, .lsdlog=lsdlog,
              .emeanlog = emeanlog, .esdlog=esdlog ))),
    loglikelihood=eval(substitute(
        function(mu,y,w,residuals= FALSE,eta, extra=NULL) {
        mulog = eta2theta(eta[,1], .lmeanlog, earg= .emeanlog)
        sdlog = eta2theta(eta[,2], .lsdlog, earg= .esdlog)
        if(residuals) stop("loglikelihood residuals not implemented yet") else {
        if(is.R())
            sum(w*dlnorm(y, meanlog=mulog, sdlog=sdlog, log=TRUE)) else
            sum(w * (-log(y*sdlog*sqrt(2*pi)) - 0.5 * ((log(y) - mulog)/sdlog)^2))
        }
    }, list( .lmeanlog = lmeanlog, .lsdlog=lsdlog,
              .emeanlog = emeanlog, .esdlog=esdlog ))),
    vfamily=c("lognormal"),
    deriv=eval(substitute(expression({
        mulog = eta2theta(eta[,1], .lmeanlog, earg= .emeanlog)
        sdlog = eta2theta(eta[,2], .lsdlog, earg= .esdlog)
        dl.dmulog = (log(y)-mulog) / sdlog^2
        dl.dsdlog = -1/sdlog + (log(y)-mulog)^2 / sdlog^3
        dl.dlambda = (1 + (log(y)-mulog) / sdlog^2) / y
        dmulog.deta = dtheta.deta(mulog, .lmeanlog, earg= .emeanlog)
        dsdlog.deta = dtheta.deta(sdlog, .lsdlog, earg= .esdlog)
        w * cbind(dl.dmulog * dmulog.deta, 
                  dl.dsdlog * dsdlog.deta)
    }), list( .lmeanlog = lmeanlog, .lsdlog=lsdlog,
              .emeanlog = emeanlog, .esdlog=esdlog ))),
    weight=expression({
        wz = matrix(as.numeric(NA), n, 2)  # Diagonal!
        ed2l.dmulog2 = 1 / sdlog^2
        ed2l.dsdlog2 = 2 * ed2l.dmulog2
        wz[,iam(1,1,M)] = ed2l.dmulog2 * dmulog.deta^2
        wz[,iam(2,2,M)] = ed2l.dsdlog2 * dsdlog.deta^2
        wz = w * wz
        wz
    }))
}

if(!is.R()) {

    qlognormal = function(p, meanlog=0, sdlog=1, lambda=0)
        lambda + exp(qnorm(p=p, mean=meanlog, sd=sdlog))

    dlognormal = function(x, meanlog=0, sdlog=1, lambda=0)
        (x > lambda) *
        dnorm(x=log(abs(x-lambda)), mean=meanlog, sd=sdlog) / (x-lambda)

    rlognormal = function(n, meanlog=0, sdlog=1, lambda=0)
        lambda + exp(rnorm(n, mean=meanlog, sd=sdlog))

    plognormal = function(q, meanlog=0, sdlog=1, lambda=0)
        (q>lambda) * pnorm(q=log(abs(q-lambda)), mean=meanlog, sd=sdlog)
}




lognormal3 = function(lmeanlog="identity", lsdlog="loge",
                      emeanlog=list(), esdlog=list(),
                      powers.try = (-3):3,
                      delta=NULL, zero=NULL)
{
    if(length(delta) && !is.Numeric(delta, positive=TRUE))
        stop("bad input for argument argument \"delta\"")
    if(mode(lmeanlog) != "character" && mode(lmeanlog) != "name")
        lmeanlog = as.character(substitute(lmeanlog))
    if(mode(lsdlog) != "character" && mode(lsdlog) != "name")
        lsdlog = as.character(substitute(lsdlog))
    if(length(zero) && (!is.Numeric(zero, integer=TRUE, posit=TRUE) ||
                        zero > 3))
        stop("bad input for argument argument \"zero\"")
    if(!is.list(emeanlog)) emeanlog = list()
    if(!is.list(esdlog)) esdlog = list()

    new("vglmff",
    blurb=c("Three-parameter (univariate) lognormal distribution\n\n",
            "Links:    ",
            namesof("meanlog", lmeanlog, earg=emeanlog, tag= TRUE),
            "; ", namesof("sdlog", lsdlog, earg=esdlog, tag= TRUE),
            "; ", namesof("lambda", "identity", earg=list(), tag= TRUE)),
    constraints=eval(substitute(expression({
        constraints = cm.zero.vgam(constraints, x, .zero, M)
    }), list( .zero=zero ))),
    initialize=eval(substitute(expression({
        if(ncol(cbind(y)) != 1)
            stop("response must be a vector or a one-column matrix")
        predictors.names = 
        c(namesof("meanlog", .lmeanlog, earg=.emeanlog, tag=FALSE), 
          namesof("sdlog",   .lsdlog, earg=.esdlog, tag=FALSE), "lambda")

        if(!length(etastart)) {
            miny = min(y)
            if(length( .delta)) {
                lambda.init = rep(miny- .delta, length=n)
            } else {
                pvalue.vec = NULL
                powers.try = .powers.try
                for(delta in 10^powers.try) {
                    pvalue.vec = c(pvalue.vec,
                                   shapiro.test(sample(log(y-miny+delta),
                                   size=min(5000, length(y ))))$p.value) 
                }
                index.lambda=(1:length(powers.try))[pvalue.vec==max(pvalue.vec)]
                lambda.init = miny - 10^powers.try[index.lambda]
            }
            junk = if(is.R()) lm.wfit(x=x, y=log(y-lambda.init), w=w) else 
                              lm.wfit(x=x, y=log(y-lambda.init), w=w, method="qr")
            sdlog.y.est = sqrt( sum(w * junk$resid^2) / junk$df.residual )
            etastart = cbind(mu=log(median(y - lambda.init)),
            sdlog=rep(theta2eta(sdlog.y.est, .lsdlog, earg= .esdlog), length=n),
            lambda = lambda.init)
        }
    }), list( .lmeanlog=lmeanlog, .lsdlog=lsdlog,
              .emeanlog = emeanlog, .esdlog=esdlog,
              .delta = delta, .powers.try=powers.try ))),
    inverse=eval(substitute(function(eta, extra=NULL) {
        mymu = eta2theta(eta[,1], .lmeanlog, earg= .emeanlog)
        sdlog = eta2theta(eta[,2], .lsdlog, earg= .esdlog)
        lambda = eta2theta(eta[,3], "identity", earg=list())
        lambda + exp(mymu + 0.5 * sdlog^2)
    }, list( .lmeanlog=lmeanlog, .lsdlog=lsdlog,
              .emeanlog = emeanlog, .esdlog=esdlog ))),
    last=eval(substitute(expression({
        misc$link = c("meanlog"= .lmeanlog,"sdlog"= .lsdlog,"lambda"="identity")
        misc$earg = list("meanlog"= .emeanlog, "sdlog"= .esdlog,
                         "lambda"=list())
        misc$expected = TRUE
    }), list( .lmeanlog=lmeanlog, .lsdlog=lsdlog,
              .emeanlog = emeanlog, .esdlog=esdlog ))),
    loglikelihood=eval(substitute(
        function(mu,y,w,residuals= FALSE,eta, extra=NULL) {
        mymu = eta2theta(eta[,1], .lmeanlog, earg= .emeanlog)
        sdlog = eta2theta(eta[,2], .lsdlog, earg= .esdlog)
        lambda = eta2theta(eta[,3], "identity", earg=list())
        if(any(y < lambda))
            cat("warning: bad y\n")
        if(residuals) stop("loglikelihood residuals not implemented yet") else {
        if(is.R())
            sum(w*dlnorm(y-lambda, meanlog=mymu, sdlog=sdlog, log=TRUE)) else
            sum(w * (-log((y-lambda)*sdlog*sqrt(2*pi)) -
                     0.5 * ((log(y-lambda) - mymu)/sdlog)^2))
        }
    }, list( .lmeanlog=lmeanlog, .lsdlog=lsdlog,
              .emeanlog = emeanlog, .esdlog=esdlog ))),
    vfamily=c("lognormal3"),
    deriv=eval(substitute(expression({
        mymu = eta2theta(eta[,1], .lmeanlog, earg= .emeanlog)
        sdlog = eta2theta(eta[,2], .lsdlog, earg= .esdlog)
        lambda = eta2theta(eta[,3], "identity", earg=list())
        if(any(y < lambda))
            cat("warning: bad y\n")
        dl.dmymu = (log(y-lambda)-mymu) / sdlog^2
        dl.dsdlog = -1/sdlog + (log(y-lambda)-mymu)^2 / sdlog^3
        dl.dlambda = (1 + (log(y-lambda)-mymu) / sdlog^2) / (y-lambda)
        dmymu.deta = dtheta.deta(mymu, .lmeanlog, earg= .emeanlog)
        dsdlog.deta = dtheta.deta(sdlog, .lsdlog, earg= .esdlog)
        dlambda.deta = dtheta.deta(lambda, "identity", earg=list())
        w * cbind(dl.dmymu * dmymu.deta, 
                  dl.dsdlog * dsdlog.deta, 
                  dl.dlambda * dlambda.deta)
    }), list( .lmeanlog=lmeanlog, .lsdlog=lsdlog,
              .emeanlog = emeanlog, .esdlog=esdlog ))),
    weight=expression({
        wz = matrix(0, n, dimm(M))
        ed2l.dmymu2 = 1 / sdlog^2
        ed2l.dsdlog = 2 / sdlog^2
        temp9 = exp(-mymu+sdlog^2 / 2)
        ed2l.dlambda2 = exp(2*(-mymu+sdlog^2)) * (1+sdlog^2) / sdlog^2
        wz[,iam(1,1,M)] = ed2l.dmymu2 * dmymu.deta^2
        wz[,iam(2,2,M)] = ed2l.dsdlog * dsdlog.deta^2
        wz[,iam(3,3,M)] = ed2l.dlambda2 * dlambda.deta^2
        wz[,iam(1,3,M)] = temp9 * dmymu.deta * dlambda.deta / sdlog^2
        wz[,iam(2,3,M)] = -2 * temp9 / sdlog * dsdlog.deta * dlambda.deta
        wz = w * wz
        wz
    }))
}








interleave.VGAM = function(L, M) c(matrix(1:L, nrow=M, byrow=TRUE))

negbinomial = function(lmu = "loge", lk = "loge",
                       emu =list(), ek=list(),
                       ik = NULL, cutoff = 0.995, Maxiter=5000,
                       deviance.arg=FALSE, method.init=1, shrinkage.init=0.95,
                       zero = -2)
{




    if(length(ik) && !is.Numeric(ik, posit=TRUE))
        stop("bad input for argument \"ik\"")
    if(!is.Numeric(cutoff, allow=1) || cutoff<0.8 || cutoff>=1)
        stop("range error in the argument \"cutoff\"")
    if(!is.Numeric(Maxiter, integ=TRUE, allow=1) || Maxiter < 100)
        stop("bad input for argument \"Maxiter\"")
    if(!is.Numeric(method.init, allow=1, integ=TRUE, posit=TRUE) ||
       method.init > 2) stop("argument \"method.init\" must be 1 or 2")
    if(!is.Numeric(shrinkage.init, allow=1) || shrinkage.init < 0 ||
       shrinkage.init > 1) stop("bad input for argument \"shrinkage.init\"")

    if(mode(lmu) != "character" && mode(lmu) != "name")
        lmu = as.character(substitute(lmu))
    if(mode(lk) != "character" && mode(lk) != "name")
        lk = as.character(substitute(lk))
    if(!is.list(emu)) emu = list()
    if(!is.list(ek)) ek = list()

    ans = 
    new("vglmff",
    blurb=c("Negative-binomial distribution\n\n",
            "Links:    ",
            namesof("mu", lmu, earg=emu), ", ",
            namesof("k", lk, earg=ek), "\n",
            "Mean:     mu\n",
            "Variance: mu * (1 + mu/k)"),
    constraints=eval(substitute(expression({
        temp752 = .zero
        if(length(temp752) && all(temp752 == -2))
            temp752 = 2*(1:ncol(y))
        constraints = cm.zero.vgam(constraints, x, temp752, M)
    }), list( .zero=zero ))),
    initialize=eval(substitute(expression({
        if(is.R())
            assign("CQO.FastAlgorithm", ( .lmu == "loge") &&
                   ( .lk == "loge"), envir = VGAMenv) else
            CQO.FastAlgorithm <<- ( .lmu == "loge") && ( .lk == "loge")
        if(any(function.name == c("cqo","cao")) &&
           is.Numeric( .zero, allow=1) && .zero != -2)
            stop("argument zero=-2 is required")

        y = as.matrix(y) 
        M = 2 * ncol(y) 
        NOS = ncoly = ncol(y)  # Number of species
        predictors.names = c(namesof(if(NOS==1) "mu" else
            paste("mu", 1:NOS, sep=""), .lmu, earg=.emu, tag=FALSE),
            namesof(if(NOS==1) "k" else paste("k", 1:NOS, sep=""), .lk, earg=.ek,
            tag=FALSE))
        predictors.names = predictors.names[interleave.VGAM(M, M=2)]
        if(!length(etastart)) {
            mu.init = y
            for(iii in 1:ncol(y)) {
                use.this = if( .method.init == 2) {
                    weighted.mean(y[,iii], w) + 1/16
                } else {
                    median(y[,iii]) + 1/16
                }
                mu.init[,iii] = (1- .sinit) * (y[,iii]+1/16) + .sinit * use.this
            }


            if( is.Numeric( .k.init )) {
                kay.init = matrix( .k.init, nr=n, nc=NOS, byrow=TRUE)
            } else {
                negbinomial.Loglikfun = function(kmat, y, x, w, extraargs) {
                     mu = extraargs
                     sum(w * (y * log(mu/(mu+kmat)) +
                              kmat*log(kmat/(mu+kmat)) + lgamma(y+kmat) -
                              lgamma(kmat) - lgamma(y+1)))
                }
                k.grid = 2^((-3):6)
                kay.init = matrix(0, nr=n, nc=NOS)
                for(spp. in 1:NOS) {
                    kay.init[,spp.] = getMaxMin(k.grid,
                                      objfun=negbinomial.Loglikfun,
                                      y=y[,spp.], x=x, w=w,
                                      extraargs= mu.init[,spp.])
                }
            }
            etastart = cbind(theta2eta(mu.init, .lmu, earg= .emu),
                             theta2eta(kay.init, .lk, earg= .ek))
            etastart = etastart[,interleave.VGAM(M, M=2),drop=FALSE]
        }
    }), list( .lmu=lmu, .lk=lk, .k.init=ik, .zero=zero,
              .emu=emu, .ek=ek,
              .sinit=shrinkage.init,
              .method.init=method.init ))),
    inverse=eval(substitute(function(eta, extra=NULL) {
        NOS = ncol(eta) / 2
        eta2theta(eta[,2*(1:NOS)-1,drop=FALSE], .lmu, earg= .emu)
    }, list( .lmu=lmu, .emu=emu, .ek=ek ))),
    last=eval(substitute(expression({
        if(is.R()) {
            if(exists("CQO.FastAlgorithm", envir = VGAMenv))
                rm("CQO.FastAlgorithm", envir = VGAMenv)
        } else {
            while(exists("CQO.FastAlgorithm"))
                remove("CQO.FastAlgorithm")
        }
        temp0303 = c(rep( .lmu, length=NOS), rep( .lk, length=NOS))
        names(temp0303) = c(if(NOS==1) "mu" else paste("mu", 1:NOS, sep=""), 
                            if(NOS==1) "k" else paste("k", 1:NOS, sep=""))
        temp0303 = temp0303[interleave.VGAM(M, M=2)]
        misc$link = temp0303 # Already named
        misc$earg = vector("list", M)
        names(misc$earg) = names(misc$link)
        for(ii in 1:NOS) {
            misc$earg[[2*ii-1]] = .emu
            misc$earg[[2*ii  ]] = .ek
        }
        misc$cutoff = .cutoff 
        misc$method.init = .method.init 
    }), list( .lmu=lmu, .lk=lk, .cutoff=cutoff,
              .emu=emu, .ek=ek,
              .method.init=method.init ))),
    link=eval(substitute(function(mu, extra=NULL) {
        temp = theta2eta(mu, .lmu, earg= .emu)
        temp = cbind(temp, NA * temp)
        temp[,interleave.VGAM(ncol(temp), M=2),drop=FALSE]
    }, list( .lmu=lmu, .emu=emu, .ek=ek ))),
    loglikelihood=eval(substitute(
        function(mu,y,w,residuals= FALSE,eta, extra=NULL) {
        NOS = ncol(eta) / 2
        temp300 = eta[,2*(1:NOS),drop=FALSE]
        if( .lk == "loge") {
            bigval = 28
            temp300 = ifelse(temp300 >  bigval,  bigval, temp300)
            temp300 = ifelse(temp300 < -bigval, -bigval, temp300)
        }
        kmat = eta2theta(temp300, .lk, earg= .ek)
        ans = 
        sum(w * (y * log(mu/(mu+kmat)) + kmat*log(kmat/(mu+kmat)) +
                 lgamma(y+kmat) - lgamma(kmat) - lgamma(y+1 )))
        if(residuals) stop("loglikelihood residuals not implemented yet") else
        ans
    }, list( .lk=lk, .emu=emu, .ek=ek ))),
    vfamily=c("negbinomial"),
    deriv=eval(substitute(expression({
        NOS = ncol(eta) / 2
        M = ncol(eta)
        temp3 = eta[,2*(1:NOS),drop=FALSE]
        bigval = 28
        temp3 = ifelse(temp3 >  bigval,  bigval, temp3)
        temp3 = ifelse(temp3 < -bigval, -bigval, temp3)
        kmat = eta2theta(temp3, .lk, earg= .ek)
        dl.dmu = y/mu - (y+kmat)/(kmat+mu)
        dl.dk = digamma(y+kmat) - digamma(kmat) - (y+kmat)/(mu+kmat) + 1 +
                log(kmat/(kmat+mu))
        dmu.deta = dtheta.deta(mu, .lmu, earg= .emu)
        dk.deta = dtheta.deta(kmat, .lk, earg= .ek)
        myderiv = w * cbind(dl.dmu * dmu.deta, dl.dk * dk.deta)
        myderiv[,interleave.VGAM(M, M=2)]
    }), list( .lmu=lmu, .lk=lk, .emu=emu, .ek=ek ))),
    weight=eval(substitute(expression({
        wz = matrix(as.numeric(NA), n, M)  # wz is 'diagonal' 
        ed2l.dmu2 = 1/mu - 1/(mu+kmat)
        fred2 = dotFortran(name="enbin9",
                      ans=double(n*NOS), as.double(kmat),
                      as.double(mu), as.double( .cutoff ),
                      as.integer(n), ok=as.integer(1), as.integer(NOS),
                      sumpdf=double(1), macheps=as.double(.Machine$double.eps),
                      as.integer( .Maxiter ))
        if(fred2$ok != 1)
            stop("error in Fortran subroutine exnbin9")
        dim(fred2$ans) = c(n, NOS)
        ed2l.dk2 = -fred2$ans - 1/kmat + 1/(kmat+mu)
        wz[,2*(1:NOS)-1] = dmu.deta^2 * ed2l.dmu2
        wz[,2*(1:NOS)] = dk.deta^2 * ed2l.dk2
        w * wz
    }), list( .cutoff=cutoff, .Maxiter=Maxiter ))))

    if(deviance.arg) ans@deviance=eval(substitute(
        function(mu,y,w,residuals= FALSE,eta,extra=NULL) {
        NOS = ncol(eta) / 2
        temp300 =  eta[,2*(1:NOS),drop=FALSE]
        if( .lk == "loge") {
            bigval = 28
            temp300[temp300 >  bigval] =  bigval
            temp300[temp300 < -bigval] = -bigval
        } else stop("can only handle the 'loge' link")
        k =  eta2theta(temp300, .lk, earg= .ek)
        devi = 2 * (y*log(ifelse(y<1, 1, y)/mu) + (y+k)*log((mu+k)/(k+y )))
        if(residuals)
           sign(y - mu) * sqrt(abs(devi) * w) else
           sum(w * devi)
    }, list( .lk=lk, .emu=emu, .ek=ek,)))
    ans
}


negbin.ab = function(link.alpha ="loge", link.k ="loge",
                     ealpha=list(), ek=list(),
                     k.init=1,
                     zero=2,
                     cutoff=0.995)
{



    if(!is.Numeric(cutoff, allow=1) || cutoff<0.8 || cutoff>=1)
        stop("range error in the argument cutoff")
    if(length(zero) && !is.Numeric(zero, integer=TRUE, posit=TRUE))
        stop("bad input for argument \"zero\"")

    if(mode(link.alpha) != "character" && mode(link.alpha) != "name")
        link.alpha = as.character(substitute(link.alpha))
    if(mode(link.k) != "character" && mode(link.k) != "name")
        link.k = as.character(substitute(link.k))
    if(!is.list(ealpha)) ealpha = list()
    if(!is.list(ek)) ek = list()

    new("vglmff",
    blurb=c("Negative-binomial distribution\n\n",
            "Links:    ",
            namesof("alpha", link.alpha, earg=ealpha), ", ",
            namesof("k", link.k, earg=ek),
            "\n",
            "Mean:     alpha * k",
            "\n",
            "Variance: alpha * k * (1 + alpha)"),
    constraints=eval(substitute(expression({
        constraints = cm.zero.vgam(constraints, x, .zero, M)
    }), list( .zero=zero ))),
    initialize=eval(substitute(expression({
        if(ncol(cbind(y)) != 1)
            stop("response must be a vector or a one-column matrix")
        predictors.names =
        c(namesof("alpha", .link.alpha, earg=.ealpha, tag=FALSE),
          namesof("k", .link.k, earg=.ek, tag=FALSE))

        if(!length(etastart)) {
            etastart = cbind(
            theta2eta((y + 1/8) / .k.init, .link.alpha, earg= .ealpha),
            theta2eta( .k.init, .link.k, earg= .ek))
        }
    }), list( .link.alpha=link.alpha, .link.k=link.k, .k.init=k.init,
              .ealpha=ealpha, .ek=ek ))),
    inverse=eval(substitute(function(eta, extra=NULL) {
        alpha = eta2theta(eta[,1], .link.alpha, earg= .ealpha)
        k = eta2theta(eta[,2], .link.k, earg= .ek)
        alpha * k 
    }, list( .link.alpha=link.alpha, .link.k=link.k,
              .ealpha=ealpha, .ek=ek ))),
    last=eval(substitute(expression({
        misc$link = c(alpha= .link.alpha, k= .link.k)
        misc$earg = list(alpha= .ealpha, k= .ek )
    }), list( .link.alpha=link.alpha, .link.k=link.k,
              .ealpha=ealpha, .ek=ek ))),
    loglikelihood=eval(substitute(
            function(mu,y,w,residuals= FALSE,eta, extra=NULL) {
        alpha = eta2theta(eta[,1], .link.alpha, earg= .ealpha)
        k = eta2theta(eta[,2], .link.k, earg= .ek)
        if(residuals) stop("loglikelihood residuals not implemented yet") else
        sum(w * (y * log(alpha) - (y+k)*log(alpha+1) + lgamma(y+k) -
                 lgamma(k) - lgamma(y+1 )))
    }, list( .link.alpha=link.alpha, .link.k=link.k,
              .ealpha=ealpha, .ek=ek ))),
    vfamily=c("negbin.ab"),
    deriv=eval(substitute(expression({
        alpha = eta2theta(eta[,1], .link.alpha, earg= .ealpha)
        k = eta2theta(eta[,2], .link.k, earg= .ek)
        dl.dalpha = (y/alpha - k)/(1+alpha)
        dl.dk = digamma(y+k) - digamma(k) - log1p(alpha)
        dalpha.deta = dtheta.deta(alpha, .link.alpha, earg= .ealpha)
        dk.deta = dtheta.deta(k, .link.k, earg= .ek)
        cbind(w * dl.dalpha * dalpha.deta, w * dl.dk * dk.deta)
    }), list( .link.alpha=link.alpha, .link.k=link.k,
              .ealpha=ealpha, .ek=ek ))),
    weight=eval(substitute(expression({
        wz = matrix(as.numeric(NA), n, dimm(M))    # 3==dimm(M)
        ed2l.dalpha2 = k/(alpha*(1+alpha))
        ed2l.dalphak =  1/(1+alpha)   # Not -1/(1+alpha)

        fred = dotFortran(name="enbin8",
                      ans=double(n),
                      as.double(k),
                      as.double(1/(1+alpha)),
                      as.double( .cutoff ),
                      as.integer(n), ok=as.integer(1), as.integer(1),
                      sumpdf=double(1), macheps=as.double(.Machine$double.eps))
        if(fred$ok != 1)
            stop("error in Fortran subroutine enbin8")
        ed2l.dk2 = -fred$ans

        wz[,iam(1,1,M)] = dalpha.deta^2 * ed2l.dalpha2
        wz[,iam(2,2,M)] = dk.deta^2 * ed2l.dk2
        wz[,iam(1,2,M)] = dk.deta * dalpha.deta * ed2l.dalphak 

        w * wz
    }), list( .cutoff=cutoff,
              .ealpha=ealpha, .ek=ek ))))
}



if(FALSE)
nbmud = function(lmu = c("loge","identity","reciprocal"),
                 k.init = 1,
                 zero = -2,
                 cutoff = 0.995,
                 deviance.arg=FALSE)
{
    ans = negbinomial(link.mu = lmu[1],
                    link.k = "reciprocal",
                    k.init = k.init,
                    zero = zero,
                    cutoff = cutoff,
                    deviance.arg=deviance.arg)
    ans@vfamily = "nbmud"
    ans
}




if(FALSE)
neg.binomial = function(link.p="logit", link.k="loge",
                        ep=list(), ek=list(),
                        zero=2,
                        ik=NULL,
                        cutoff=0.995)
{



    if(!is.Numeric(cutoff, allow=1) || cutoff<0.8 || cutoff>=1)
        stop("range error in the argument cutoff")

    if(mode(link.p) != "character" && mode(link.p) != "name")
        link.p = as.character(substitute(link.p))
    if(link.p=="canonical")
        link.p = "logc"
    if(mode(link.k) != "character" && mode(link.k) != "name")
        link.k = as.character(substitute(link.k))
    if(!is.list(ep)) ep = list()
    if(!is.list(ek)) ek = list()

    new("vglmff",
    blurb=c("Negative-binomial distribution\n\n",
            "Links:    ",
            namesof("p", link.p, earg=ep), ", ",
            namesof("k", link.k, earg=ek), "; mu=k*(1-p)/p",
            "\n",
            "Variance: mu(1 + mu/k)"),
    constraints=eval(substitute(expression({
        constraints = cm.zero.vgam(constraints, x, .zero, M)
    }), list( .zero=zero ))),
    initialize=eval(substitute(expression({
        if(ncol(cbind(y)) != 1)
            stop("response must be a vector or a one-column matrix")
        y = as.numeric(y)
        if(any(y < 0))
            stop("response must be non-negative")
        if(max(abs(y - round(y )))>0.00001)
            stop("response must be integer-valued")

        predictors.names =
        c(namesof("p", .link.p, earg=.ep, tag=FALSE),
          namesof("k", .link.k, earg=.ek, tag=FALSE))



        junk = lm(y ~ x - 1, weight=w, smart= FALSE) # singular.ok = FALSE,
        var.y.est = summary(junk)$sigma^2
        mu.adj = fitted(junk)

        if(FALSE) { 
            mu = rep(weighted.mean(y, w=w), len=length(y))
            mu = rep(median(rep(y+0.167, times=w)), len=length(y))

            k = mean(rep(mu^2 / (var.y.est - mu), w), trim=0.05)
            k = rep(k, length(mu))
        } else {
            mu = mu.adj
            mu[mu <= 0] = min(mu[mu > 0])
            k = mu.adj^2 / (var.y.est - mu.adj)
            k[k <= 0] = quantile(k[k>0], prob=0.02)
        }

        if(length( .ik )) {
            mu = median(rep(y, times=w))
            k =  rep( .ik , len=length(y))
        }
    
        if(!length(etastart)) {
            prob = k / (k + mu)
            etastart = cbind(theta2eta(prob, .link.p, earg= .ep),
                              theta2eta(k, .link.k, earg= .ek))
        }
    }), list( .link.p=link.p, .link.k=link.k, .ik=ik,
              .ep=ep, .ek=ek ))),
    inverse=eval(substitute(function(eta, extra=NULL) {
        prob = eta2theta(eta[,1], .link.p, earg= .ep)
        k = eta2theta(eta[,2], .link.k, earg= .ek)
        k * (1 - prob) / prob
    }, list( .link.p=link.p, .link.k=link.k,
              .ep=ep, .ek=ek ))),
    last=eval(substitute(expression({
        misc$link = c(p= .link.p, k= .link.k )
        misc$earg = list(p= .ep, k= .ek )
    }), list( .link.p=link.p, .link.k=link.k,
              .ep=ep, .ek=ek ))),
    loglikelihood=eval(substitute(
        function(mu,y,w,residuals= FALSE,eta, extra=NULL) {
        prob = eta2theta(eta[,1], .link.p, earg= .ep)
        k = eta2theta(eta[,2], .link.k, earg= .ek)
        if(residuals) stop("loglikelihood residuals not implemented yet") else
        sum(w * (y * log1p(-prob) + k * log(prob) + lgamma(y+k) -
                 lgamma(k) - lgamma(y+1 )))
    }, list( .link.p=link.p, .link.k=link.k,
              .ep=ep, .ek=ek ))),
    vfamily=c("neg.binomial"),
    deriv=eval(substitute(expression({
        prob = eta2theta(eta[,1], .link.p, earg= .ep)
        k = eta2theta(eta[,2], .link.k, earg= .ek)
        dl.dp = k/prob - y/(1-prob)
        dl.dk = log(prob) + digamma(y+k) - digamma(k)
        dp.deta = dtheta.deta(prob, .link.p, earg= .ep)
        dk.deta = dtheta.deta(k, .link.k, earg= .ek)
        w * cbind(dl.dp * dp.deta, dl.dk * dk.deta)
    }), list( .link.p=link.p, .link.k=link.k,
              .ep=ep, .ek=ek ))),
    weight=eval(substitute(expression({
        wz = matrix(as.numeric(NA), n, dimm(M))    # 3==dimm(M)
        d2l.dpk = 1/prob


        ed2l.dp2 = -k/(prob^2 * (1-prob))      # "e" for expected value
        fred = dotFortran(name="exnbin",
                      ans=double(n),
                      as.double(k),
                      as.double(prob),
                      as.double( .cutoff ),
                      as.integer(n), ok=as.integer(1), as.integer(1),
                      sumpdf=double(1))
        if(fred$ok != 1)
            stop("error in Fortran subroutine exnbin")

        ed2l.dk2 = fred$ans

        wz[,iam(1,1,M)] = dp.deta^2 * ed2l.dp2
        wz[,iam(1,2,M)] = d2l.dpk * dp.deta * dk.deta #ed2l.dpk=d2l.dpk
        wz[,iam(2,2,M)] = dk.deta^2 * ed2l.dk2
        wz = -w * wz
        wz
    }), list( .cutoff=cutoff,
              .ep=ep, .ek=ek ))))
}



if(FALSE)
neg.binomial.k = function(k, link="logit", earg=list(), expected=TRUE, ...)
{

    if(!is.Numeric(k, allow=1, posit=TRUE))
        stop("bad input for argument argument \"k\"")
    if(mode(link) != "character" && mode(link) != "name")
        link = as.character(substitute(link))
    if(!is.list(earg)) earg = list()

    new("vglmff",
    blurb=c("Negative-binomial distribution with k known and p unknown\n",
            "(k=", k, ") ", 
            if(k==1) "Geometric\n\n" else "\n\n",
            "Links:    ",
            namesof("p", link, earg=earg), "; p=",k,"/(",k,"+mu)",
            "\n",
            "Variance: ",
            if(k==1) "Geometric: mu(1+mu)" else 
                   paste("mu(1 + mu/",k,")", sep="")),
    deviance=eval(substitute(
        function(mu,y,w,residuals= FALSE,eta,extra=NULL) {
        prob = .k / ( .k + mu)
        devy = .k * log( .k / ( .k + y))
        nz = y != 0
        devy[nz] = devy[nz] + y[nz] * log(y[nz] / ( .k + y[nz]))
        devmu = y * log1p(-prob) + .k * log(prob)
        devi = 2 * (devy - devmu)
        if(residuals)
           sign(y - mu) * sqrt(abs(devi) * w) else
           sum(w * devi)
    }, list( .link=link, .earg=earg, .k=k ))),
    initialize=eval(substitute(expression({
        predictors.names = namesof("p", .link, earg=.ep, tag=FALSE)
        mu = y + 0.167 * (y == 0)

        if(!length(etastart)) {
            prob = .k / ( .k + mu)
            etastart = theta2eta(prob, .link, earg= .earg)
        }
    }), list( .link=link, .earg=earg, .k=k ))),
    inverse=eval(substitute(function(eta, extra=NULL) {
        prob = eta2theta(eta, .link, earg= .earg)
        .k * (1 - prob) / prob
    }, list( .link=link, .earg=earg, .k=k ))),
    last=eval(substitute(expression({
        misc$link = c(p = .link)
        misc$k = .k
    }), list( .link=link, .earg=earg, .k=k ))),
    loglikelihood=eval(substitute(
        function(mu,y,w,residuals= FALSE,eta, extra=NULL) {
        prob = eta2theta(eta, .link, earg= .earg)
        if(residuals) stop("loglikelihood residuals not implemented yet") else
        sum(w * (y * log1p(-prob) + .k * log(prob) + lgamma(y+ .k) -
                 lgamma( .k ) - lgamma(y+1 )))
    }, list( .link=link, .earg=earg, .k=k ))),
    vfamily=c("neg.binomial.k"),
    deriv=eval(substitute(expression({
        prob = .k / ( .k + mu)
        dp.deta = dtheta.deta(prob, .link, earg= .earg)
        w * ( .k/prob - y/(1-prob)) * dp.deta
    }), list( .link=link, .earg=earg, .k=k ))),
    weight=eval(substitute(expression({
        wz = dp.deta^2 * (y/(1 - prob)^2 + .k/prob^2) 
        if(! .expected) {
            d2pdeta2 = d2theta.deta2(prob, .link, earg= .earg)
            wz = wz - d2pdeta2 * ( .k/prob - y/(1-prob))
        }
        w * wz
    }), list( .link=link, .earg=earg, .k=k, .expected=expected ))))
}





simple.poisson = function()
{
    new("vglmff",
    blurb=c("Poisson distribution\n\n",
            "Link:     log(lambda)",
            "\n",
            "Variance: lambda"),
    deviance= function(mu, y, w, residuals = FALSE, eta, extra=NULL) {
        nz = y > 0
        devi =  - (y - mu)
        devi[nz] = devi[nz] + y[nz] * log(y[nz]/mu[nz])
        if(residuals) sign(y - mu) * sqrt(2 * abs(devi) * w) else
            2 * sum(w * devi)
    },
    initialize=expression({
        if(ncol(cbind(y)) != 1)
            stop("response must be a vector or a one-column matrix")
        predictors.names = "log(lambda)"
        mu = y + 0.167 * (y == 0)
        if(!length(etastart))
            etastart = log(mu)
    }), 
    inverse=function(eta, extra=NULL)
        exp(eta),
    last=expression({
        misc$link = c(lambda = "loge")
        misc$earg = list(lambda = list())
    }),
    link=function(mu, extra=NULL)
        log(mu),
    vfamily="simple.poisson",
    deriv=expression({
        lambda = mu
        dl.dlambda = -1 + y/lambda
        dlambda.deta = dtheta.deta(theta=lambda, link="loge", earg= list())
        w * dl.dlambda * dlambda.deta
    }),
    weight=expression({
        d2l.dlambda2 = 1 / lambda
        w * d2l.dlambda2 * dlambda.deta^2
    }))
}



studentt =  function(link.df="loglog", earg=list())
{

    if(mode(link.df) != "character" && mode(link.df) != "name")
        link.df = as.character(substitute(link.df))
    if(!is.list(earg)) earg = list()

    new("vglmff",
    blurb=c("Student t-distribution\n\n",
            "Link:     ",
            namesof("df", link.df, earg=earg),
            "\n",
            "Variance: df/(df-2) if df > 2\n"),
    initialize=eval(substitute(expression({
        if(ncol(cbind(y)) != 1)
            stop("response must be a vector or a one-column matrix")
        predictors.names = namesof("df", .link.df, earg=.earg, tag=FALSE)
        if(!length(etastart)) {
            init.df = (2*var(y)/(var(y)-1))
            if(is.na(init.df) || init.df<1)
                init.df = 4
            etastart = rep(theta2eta(init.df, .link.df, earg= .earg),
                           len=length(y))
        }
    }), list( .link.df=link.df, .earg=earg ))), 
    inverse=eval(substitute(function(eta, extra=NULL) {
        df =  eta2theta(eta, .link.df, earg= .earg)
        ifelse(df > 1, 0, NA)
    }, list( .link.df=link.df, .earg=earg ))),
    last=eval(substitute(expression({
        misc$link = c(df = .plink )
        misc$earg = list(df = .earg )
    }), list( .plink=link.df, .earg=earg ))),
    link=eval(substitute(function(mu, extra=NULL) {
        alpha = mu / sqrt(2/pi - mu^2)
        theta2eta(alpha, .plink, earg= .earg)
    }, list( .plink=link.df, .earg=earg ))),
    loglikelihood=eval(substitute(function(mu,  y,  w,  residuals = FALSE,  eta, 
        extra=NULL) {
        df =  eta2theta(eta, .link.df, earg= .earg)
        temp1 =  y^2 / df
        if(residuals) stop("loglikelihood residuals not implemented yet") else {
        if(is.R()) sum(w * dt(x=y, df=df, log=TRUE)) else
            sum(w * (-log(pi*df)/2 - (df+1)*log1p(temp1)/2 +
                    lgamma((df+1)/2) - lgamma(df/2)))
        }
    }, list( .link.df=link.df, .earg=earg ))), 
    vfamily=c("studentt"),
    deriv=eval(substitute(expression({
        df = eta2theta(eta, .link.df, earg= .earg)
        temp = 1/df
        temp1 = y^2 * temp
        dl.ddf = 0.5*(-temp -log1p(temp1) +(df+1)*y^2/(df^2 * (1+temp1)) +
                 digamma((df+1)/2)-digamma(df/2))
        ddf.deta =  dtheta.deta(theta=df, .link.df, earg= .earg)
        w * dl.ddf * ddf.deta
    }), list( .link.df=link.df, .earg=earg ))),
    weight=eval(substitute(expression({
        temp2 =  (df+1)/2
        d2df.deta2 = d2theta.deta2(theta=df,  .link.df, earg= .earg)
        negative =  -trigamma(df/2)/4 -
                     0.5*y^2*( (1+temp)/(df+y^2) + temp^2 )/(df+y^2)
        positive = 0.5*temp^2 +trigamma(temp2)/4 + 0.5*y^2*temp/(df+y^2)
        d2l.ddf2 =  positive + negative 
        wz =  -ddf.deta^2 * d2l.ddf2 - dl.ddf * d2df.deta2
        wz * w
    }), list( .link.df=link.df, .earg=earg ))))
}


 
chisq = function(link = "loge", earg=list())
{
    if(mode(link) != "character" && mode(link) != "name")
        link = as.character(substitute(link))
    if(!is.list(earg)) earg = list()

    new("vglmff",
    blurb=c("Chi-squared distribution\n\n",
            "Link:     ",
            namesof("df", link, earg=earg)),
    inverse =eval(substitute(function(eta,extra=NULL) {
        eta2theta(eta, .link, earg= .earg)
    }, list( .link = link, .earg=earg ))),
    initialize =eval(substitute(expression({
        if(ncol(cbind(y)) != 1)
            stop("response must be a vector or a one-column matrix")
        predictors.names = namesof("df", .link, earg=.earg, tag = FALSE)
        mu = y + 0.167 * (y == 0)
    }), list( .link = link, .earg=earg ))),
    last =eval(substitute(expression({
        misc$link = c(df = .link)
        misc$earg = list(df = .earg )
    }), list( .link = link, .earg=earg ))),
    link=eval(substitute(function(mu, extra = NULL) {
        theta2eta(mu, .link, earg= .earg)
    }, list( .link = link, .earg=earg ))),
    loglikelihood =eval(substitute(
        function(mu,y,w,residuals= FALSE,eta,extra=NULL) {
        df = eta2theta(eta, .link, earg= .earg)
        if(residuals) stop("loglikelihood residuals not implemented yet") else
        sum(w * (df*log(0.5)/2 + (df/2 - 1)*log(y) - y/2 - 
            lgamma(df/2 )))
    }, list( .link = link, .earg=earg ))),
    vfamily="chisq",
    deriv=eval(substitute(expression({
        df = eta2theta(eta, .link, earg= .earg)
        dl.dv = (log(y/2) - digamma(df/2)) / 2 
        dv.deta = dtheta.deta(df, .link, earg= .earg)
        w * dl.dv * dv.deta
    }), list( .link = link, .earg=earg ))),
    weight =eval(substitute(expression({
        ed2l.dv2 = -trigamma(df/2) / 4
        wz = -ed2l.dv2 * dv.deta^2
        wz * w
    }), list( .link = link, .earg=earg ))))
}







simplex = function(lmu="logit", lsigma="loge",
                   emu=list(), esigma=list(), imu=NULL, isigma=NULL)
{

    if(mode(lmu) != "character" && mode(lmu) != "name")
        lmu = as.character(substitute(lmu))
    if(mode(lsigma) != "character" && mode(lsigma) != "name")
        lsigma = as.character(substitute(lsigma))
    if(!is.list(emu)) emu = list()
    if(!is.list(esigma)) esigma = list()

    new("vglmff",
    blurb=c("Univariate Simplex distribution \n",
            "f(y) = [2*pi*sigma^2*(y*(1-y))^3]^(-0.5) * \n",
            "       exp[-0.5*(y-mu)^2 / (y*(1-y)*mu^2*(1-mu)^2)/sigma^2], ",
            "  0 < y < 1,\n",
            "Links:     ",
            namesof("mu", lmu, earg=emu), ", ",
            namesof("sigma", lsigma, earg=esigma), "\n\n",
            "Mean:     mu\n",
            "Variance: sigma^2"),
    initialize=eval(substitute(expression({
        y = as.numeric(y)
        if(any(y <= 0 | y >= 1))
            stop("all y values must be in (0,1)")

        predictors.names = c(namesof("mu", .lmu, earg=.emu, tag=FALSE),
                             namesof("sigma", .lsigma, earg=.esigma, tag=FALSE))

        if(!length(etastart)) {
            mu.init = rep(if(length( .imu)) .imu else
                           median(y), length=n)
            sigma.init = rep(if(length( .isigma)) .isigma else
                           sqrt(var(y)), length=n)
            etastart = cbind(theta2eta(mu.init, .lmu, earg= .emu),
                             theta2eta(sigma.init, .lsigma, earg= .esigma))
        }
    }), list( .lmu=lmu, .lsigma=lsigma,
              .emu=emu, .esigma=esigma,
              .imu=imu, .isigma=isigma ))),
    inverse=eval(substitute(function(eta, extra=NULL) {
        eta2theta(eta[,1], .lmu, earg= .emu)
    }, list( .lmu=lmu,
              .emu=emu, .esigma=esigma ))),
    last=eval(substitute(expression({
        misc$d3 = d3    # because save.weights=F
        misc$link = c(mu= .lmu, sigma= .lsigma)
        misc$earg = list(mu= .emu, sigma= .esigma)
        misc$pooled.weight = pooled.weight
    }), list( .lmu=lmu, .lsigma=lsigma,
              .emu=emu, .esigma=esigma ))),
    loglikelihood=eval(substitute(function(mu,y,w,residuals= FALSE,eta, extra=NULL) {
        sigma = eta2theta(eta[,2], .lsigma, earg= .esigma)
        if(residuals) stop("loglikelihood residuals not implemented yet") else
        sum(w * (-0.5*log(2*pi*sigma^2*(y*(1-y))^3) -
                          (0.5/sigma^2)*(y-mu)^2 / (y*(1-y)*mu^2*(1-mu)^2 )))
    }, list( .lsigma=lsigma,
             .emu=emu, .esigma=esigma ))),
    vfamily=c("simplex1"),
    deriv=eval(substitute(expression({
        if(iter==1) {
            d3 = deriv3(~ w * (-0.5*log(2*pi*sigma^2*(y*(1-y))^3) -
                          (0.5/sigma^2)*(y-mu)^2 / (y*(1-y)*mu^2*(1-mu)^2)),
                        c("mu", "sigma"), hessian= TRUE)
        }

        sigma = eta2theta(eta[,2], .lsigma, earg= .esigma)

        eval.d3 = eval(d3)
        dl.dthetas =  attr(eval.d3, "gradient")

        dmu.deta = dtheta.deta(mu, .lmu, earg= .emu)
        dsigma.deta = dtheta.deta(sigma, .lsigma, earg= .esigma)
        dtheta.detas = cbind(dmu.deta, dsigma.deta)

        dl.dthetas * dtheta.detas
    }), list( .lmu=lmu, .lsigma=lsigma,
             .emu=emu, .esigma=esigma ))),
    weight=eval(substitute(expression({
        d2l.dthetas2 =  attr(eval.d3, "hessian")

        wz = matrix(as.numeric(NA), n, dimm(M))  #3=dimm(M)
        wz[,iam(1,1,M)] = -d2l.dthetas2[,1,1] * dtheta.detas[,1]^2
        wz[,iam(2,2,M)] = -d2l.dthetas2[,2,2] * dtheta.detas[,2]^2
        wz[,iam(1,2,M)] = -d2l.dthetas2[,1,2] * dtheta.detas[,1] *
                                                dtheta.detas[,2]
        if(!.expected) {
            d2mudeta2 = d2theta.deta2(mu, .lmu, earg= .emu)
            d2sigmadeta2 = d2theta.deta2(sigma, .lsigma, earg= .esigma)
            wz[,iam(1,1,M)] = wz[,iam(1,1,M)] - dl.dthetas[,1] * d2mudeta2
            wz[,iam(2,2,M)] = wz[,iam(2,2,M)] - dl.dthetas[,2] * d2sigmadeta2
        }

        if(intercept.only) {
            sumw = sum(w)
            for(i in 1:ncol(wz))
                wz[,i] = sum(wz[,i]) / sumw
            pooled.weight = TRUE
            wz = w * wz   # Put back the weights
        } else
            pooled.weight = FALSE

        wz
    }), list( .lmu=lmu, .lsigma=lsigma, .expected=FALSE,
              .emu=emu, .esigma=esigma ))))
}



rig = function(lmu="identity", llambda="loge",
               emu=list(), elambda=list(), imu=NULL, ilambda=1)
{

    if(mode(lmu) != "character" && mode(lmu) != "name")
        lmu = as.character(substitute(lmu))
    if(mode(llambda) != "character" && mode(llambda) != "name")
        llambda = as.character(substitute(llambda))
    if(!is.Numeric(ilambda, posit=TRUE))
        stop("bad input for \"ilambda\"")
    if(!is.list(emu)) emu = list()
    if(!is.list(elambda)) elambda = list()

    new("vglmff",
    blurb=c("Reciprocal inverse Gaussian distribution \n",
            "f(y) = [lambda/(2*pi*y)]^(0.5) * \n",
            "       exp[-0.5*(lambda/y) * (y-mu)^2], ",
            "  0 < y,\n",
            "Links:     ",
            namesof("mu", lmu, earg=emu), ", ",
            namesof("lambda", llambda, earg=elambda), "\n\n",
            "Mean:     mu"),
    initialize=eval(substitute(expression({
        if(ncol(cbind(y)) != 1)
            stop("response must be a vector or a one-column matrix")
        y = as.numeric(y)
        if(any(y <= 0))
            stop("all y values must be > 0")
        predictors.names = 
        c(namesof("mu", .lmu, earg=.emu, tag=FALSE),
          namesof("lambda", .llambda, earg=.elambda, tag=FALSE))
        if(!length(etastart)) {
            mu.init = rep(if(length( .imu)) .imu else
                           median(y), length=n)
            lambda.init = rep(if(length( .ilambda )) .ilambda else
                           sqrt(var(y)), length=n)
            etastart = cbind(theta2eta(mu.init, .lmu, earg= .emu),
                             theta2eta(lambda.init, .llambda, earg= .elambda))
        }
    }), list( .lmu=lmu, .llambda=llambda,
              .emu=emu, .elambda=elambda,
              .imu=imu, .ilambda=ilambda ))),
    inverse=eval(substitute(function(eta, extra=NULL) {
        eta2theta(eta[,1], .lmu, earg= .emu)
    }, list( .lmu=lmu,
             .emu=emu, .elambda=elambda ))),
    last=eval(substitute(expression({
        misc$d3 = d3    # because save.weights=FALSE
        misc$link = c(mu= .lmu, lambda= .llambda)
        misc$earg = list(mu= .emu, lambda= .elambda)
        misc$pooled.weight = pooled.weight
    }), list( .lmu=lmu, .llambda=llambda,
              .emu=emu, .elambda=elambda ))),
    loglikelihood=eval(substitute(
                  function(mu,y,w,residuals= FALSE,eta, extra=NULL) {
        lambda = eta2theta(eta[,2], .llambda, earg= .elambda)
        if(residuals) stop("loglikelihood residuals not implemented yet") else
        sum(w * (-0.5*log(y) + 0.5*log(lambda) - (0.5*lambda/y) * (y-mu)^2))
    }, list( .llambda=llambda,
             .emu=emu, .elambda=elambda ))),
    vfamily=c("rig"),
    deriv=eval(substitute(expression({
        if(iter==1) {
            d3 = deriv3(~ w * 
                 (-0.5*log(y) + 0.5*log(lambda) - (0.5*lambda/y) * (y-mu)^2),
                        c("mu", "lambda"), hessian= TRUE)
        }

        lambda = eta2theta(eta[,2], .llambda, earg= .elambda)

        eval.d3 = eval(d3)
        dl.dthetas =  attr(eval.d3, "gradient")

        dmu.deta = dtheta.deta(mu, .lmu, earg= .emu)
        dlambda.deta = dtheta.deta(lambda, .llambda, earg= .elambda)
        dtheta.detas = cbind(dmu.deta, dlambda.deta)

        dl.dthetas * dtheta.detas
    }), list( .lmu=lmu, .llambda=llambda,
              .emu=emu, .elambda=elambda ))),
    weight=eval(substitute(expression({
        d2l.dthetas2 =  attr(eval.d3, "hessian")

        wz = matrix(as.numeric(NA), n, dimm(M))  #3=dimm(M)
        wz[,iam(1,1,M)] = -d2l.dthetas2[,1,1] * dtheta.detas[,1]^2
        wz[,iam(2,2,M)] = -d2l.dthetas2[,2,2] * dtheta.detas[,2]^2
        wz[,iam(1,2,M)] = -d2l.dthetas2[,1,2] * dtheta.detas[,1] *
                                                 dtheta.detas[,2]
        if(!.expected) {
            d2mudeta2 = d2theta.deta2(mu, .lmu, earg= .emu)
            d2lambda = d2theta.deta2(lambda, .llambda, earg= .elambda)
            wz[,iam(1,1,M)] = wz[,iam(1,1,M)] - dl.dthetas[,1] * d2mudeta2
            wz[,iam(2,2,M)] = wz[,iam(2,2,M)] - dl.dthetas[,2] * d2lambda
        }

        if(intercept.only) {
            sumw = sum(w)
            for(i in 1:ncol(wz))
                wz[,i] = sum(wz[,i]) / sumw
            pooled.weight = TRUE
            wz = w * wz   # Put back the weights
        } else
            pooled.weight = FALSE

        wz
    }), list( .lmu=lmu, .llambda=llambda, .expected=FALSE,
              .emu=emu, .elambda=elambda ))))
}



hypersecant = function(link.theta="elogit",
    earg=if(link.theta=="elogit") list(min=-pi/2, max=pi/2) else list(),
    init.theta=NULL)
{

    if(mode(link.theta) != "character" && mode(link.theta) != "name")
        link.theta = as.character(substitute(link.theta))
    if(!is.list(earg)) earg = list()

    new("vglmff",
    blurb=c("Hyperbolic Secant distribution \n",
            "f(y) = exp(theta*y + log(cos(theta ))) / (2*cosh(pi*y/2))\n",
            "  for all y,\n",
            "Link:     ",
            namesof("theta", link.theta, earg=earg), "\n\n",
            "Mean:     tan(theta)"),
    initialize=eval(substitute(expression({
        if(ncol(cbind(y)) != 1)
            stop("response must be a vector or a one-column matrix")
        predictors.names = namesof("theta", .link.theta, earg=.earg, tag=FALSE)
        if(!length(etastart)) {
            theta.init = rep(if(length( .init.theta)) .init.theta else
                             median(y), length=n)
            etastart = theta2eta(theta.init, .link.theta, earg= .earg)
        }
    }), list( .link.theta=link.theta, .earg=earg,
              .init.theta=init.theta ))),
    inverse=eval(substitute(function(eta, extra=NULL) {
        theta = eta2theta(eta, .link.theta, earg= .earg)
        tan(theta)
    }, list( .link.theta=link.theta, .earg=earg ))),
    last=eval(substitute(expression({
        misc$link = c(theta= .link.theta )
        misc$earg = list(theta= .earg )
        misc$expected = TRUE
    }), list( .link.theta=link.theta, .earg=earg ))),
    loglikelihood=eval(substitute(function(mu,y,w,residuals= FALSE,eta, extra=NULL) {
        theta = eta2theta(eta, .link.theta, earg= .earg)
        if(residuals) stop("loglikelihood residuals not implemented yet") else
        sum(w * (theta*y + log(cos(theta)) - log(cosh(pi*y/2 ))))
    }, list( .link.theta=link.theta, .earg=earg ))),
    vfamily=c("hypersecant"),
    deriv=eval(substitute(expression({
        theta = eta2theta(eta, .link.theta, earg= .earg)
        dl.dthetas =  y - tan(theta)
        dparam.deta = dtheta.deta(theta, .link.theta, earg= .earg)
        w * dl.dthetas * dparam.deta
    }), list( .link.theta=link.theta, .earg=earg ))),
    weight=expression({
        d2l.dthetas2 =  1 / cos(theta)^2
        wz = w * d2l.dthetas2 * dparam.deta^2
        wz
    }))
}



hypersecant.1 = function(link.theta="elogit",
    earg=if(link.theta=="elogit") list(min=-pi/2, max=pi/2) else list(),
    init.theta=NULL)
{

    if(mode(link.theta) != "character" && mode(link.theta) != "name")
        link.theta = as.character(substitute(link.theta))
    if(!is.list(earg)) earg = list()

    new("vglmff",
    blurb=c("Hyperbolic Secant distribution \n",
            "f(y) = (cos(theta)/pi) * y^(-0.5+theta/pi) * \n",
            "       (1-y)^(-0.5-theta/pi), ",
            "  0 < y < 1,\n",
            "Link:     ",
            namesof("theta", link.theta, earg=earg), "\n\n",
            "Mean:     0.5 + theta/pi", "\n",
            "Variance: (pi^2 - 4*theta^2) / (8*pi^2)"),
    initialize=eval(substitute(expression({
        if(ncol(cbind(y)) != 1)
            stop("response must be a vector or a one-column matrix")
        y = as.numeric(y)
        if(any(y <= 0 | y >= 1))
            stop("all y values must be in (0,1)")
        predictors.names = namesof("theta", .link.theta, earg=.earg, tag=FALSE)
        if(!length(etastart)) {
            theta.init = rep(if(length( .init.theta)) .init.theta else
                           median(y), length=n)

            etastart = theta2eta(theta.init, .link.theta, earg= .earg)
        }
    }), list( .link.theta=link.theta, .earg=earg,
              .init.theta=init.theta ))),
    inverse=eval(substitute(function(eta, extra=NULL) {
        theta = eta2theta(eta, .link.theta, earg= .earg)
        0.5 + theta/pi
    }, list( .link.theta=link.theta, .earg=earg ))),
    last=eval(substitute(expression({
        misc$link = c(theta= .link.theta)
        misc$earg = list(theta= .earg )
        misc$expected = TRUE
    }), list( .link.theta=link.theta, .earg=earg ))),
    loglikelihood=eval(substitute(function(mu,y,w,residuals= FALSE,eta, extra=NULL) {
        theta = eta2theta(eta, .link.theta, earg= .earg)
        if(residuals) stop("loglikelihood residuals not implemented yet") else
        sum(w * (log(cos(theta)) + (-0.5+theta/pi)*log(y) +
                (-0.5-theta/pi)*log1p(-y )))
    }, list( .link.theta=link.theta, .earg=earg ))),
    vfamily=c("hypersecant.1"),
    deriv=eval(substitute(expression({
        theta = eta2theta(eta, .link.theta, earg= .earg)
        dl.dthetas =  -tan(theta) + log(y/(1-y)) / pi 
        dparam.deta = dtheta.deta(theta, .link.theta, earg= .earg)
        w * dl.dthetas * dparam.deta
    }), list( .link.theta=link.theta, .earg=earg ))),
    weight=expression({
        d2l.dthetas2 =  1 / cos(theta)^2
        wz = w * d2l.dthetas2 * dparam.deta^2
        wz
    }))
}



leipnik = function(lmu="logit", llambda="loge",
                   emu=list(), elambda=list(), imu=NULL, ilambda=NULL)
{


    if(mode(lmu) != "character" && mode(lmu) != "name")
        lmu = as.character(substitute(lmu))
    if(mode(llambda) != "character" && mode(llambda) != "name")
        llambda = as.character(substitute(llambda))
    if(is.Numeric(ilambda) && any(ilambda <= -1))
        stop("ilambda must be > -1")
    if(!is.list(emu)) emu = list()
    if(!is.list(elambda)) elambda = list()

    new("vglmff",
    blurb=c("Leipnik's distribution \n",
    "f(y) = (y(1-y))^(-1/2) * [1 + (y-mu)^2 / (y*(1-y))]^(-lambda/2) /\n",
            "       Beta[(lambda+1)/2, 1/2], ",
            "  0 < y < 1,  lambda > -1\n",
            "Links:     ",
            namesof("mu", lmu, earg=emu), ", ",
            namesof("lambda", llambda, earg=elambda), "\n\n",
            "Mean:     mu\n",
            "Variance: mu*(1-mu)"),
    initialize=eval(substitute(expression({
        if(ncol(cbind(y)) != 1)
            stop("response must be a vector or a one-column matrix")
        y = as.numeric(y)
        if(any(y <= 0 | y >= 1))
            stop("all y values must be in (0,1)")
        predictors.names =
        c(namesof("mu", .lmu, earg=.emu, tag=FALSE),
          namesof("lambda", .llambda, earg=.elambda, tag=FALSE))
        if(!length(etastart)) {
            mu.init = rep(if(length( .imu)) .imu else
                          (y), length=n)
            lambda.init = rep(if(length( .ilambda)) .ilambda else
                           1/var(y), length=n)
            etastart = cbind(theta2eta(mu.init, .lmu, earg= .emu),
                             theta2eta(lambda.init, .llambda, earg= .elambda))
        }
    }), list( .lmu=lmu, .llambda=llambda, .imu=imu, .ilambda=ilambda,
              .emu=emu, .elambda=elambda ))),
    inverse=eval(substitute(function(eta, extra=NULL) {
        eta2theta(eta[,1], .lmu, earg= .emu)
    }, list( .lmu=lmu,
             .emu=emu, .elambda=elambda ))),
    last=eval(substitute(expression({
        if(!is.R()) 
            misc$d3 = d3    # because save.weights=FALSE
        misc$link = c(mu= .lmu, lambda= .llambda)
        misc$earg = list(mu= .emu, lambda= .elambda)
        misc$pooled.weight = pooled.weight
        misc$expected = FALSE
    }), list( .lmu=lmu, .llambda=llambda,
              .emu=emu, .elambda=elambda ))),
    loglikelihood=eval(substitute(
        function(mu,y,w,residuals= FALSE,eta, extra=NULL) {
        lambda = eta2theta(eta[,2], .llambda, earg= .elambda)
        if(residuals) stop("loglikelihood residuals not implemented yet") else
        sum(w * (-0.5*log(y*(1-y)) - 0.5 * lambda *
                log1p((y-mu)^2 / (y*(1-y ))) - lgamma((lambda+1)/2) +
                lgamma(1+ lambda/2 )))
    }, list( .llambda=llambda,
             .emu=emu, .elambda=elambda ))),
    vfamily=c("leipnik"),
    deriv=eval(substitute(expression({
        lambda = eta2theta(eta[,2], .llambda, earg= .elambda)
        if(is.R()) {
            dl.dthetas = w * cbind(dl.dmu=lambda*(y-mu)/(y*(1-y)+(y-mu)^2),
                                   dl.dlambda=-0.5*log1p((y-mu)^2 / (y*(1-y))) -
                                   0.5*digamma((lambda+1)/2) +
                                   0.5*digamma(1+lambda/2))
        } else {
            if(iter==1) {
                d3 = dfun(~ w * (-0.5*log(y*(1-y)) - 0.5 * lambda * log(1 +
                    (y-mu)^2 / (y*(1-y ))) - lgamma((lambda+1)/2) +
                    lgamma(1+ lambda/2)),
                    c("mu", "lambda"), hessian= TRUE)
            }
            eval.d3 = eval(d3)
            dl.dthetas =  attr(eval.d3, "gradient")
        }
        dmu.deta = dtheta.deta(mu, .lmu, earg= .emu)
        dlambda.deta = dtheta.deta(lambda, .llambda, earg= .elambda)
        dtheta.detas = cbind(dmu.deta, dlambda.deta)
        dl.dthetas * dtheta.detas
    }), list( .lmu=lmu, .llambda=llambda,
              .emu=emu, .elambda=elambda ))),
    weight=eval(substitute(expression({
        if(is.R()) {
            denominator = y*(1-y) + (y-mu)^2
            d2l.dthetas2 =  array(NA, c(n,2,2))
            d2l.dthetas2[,1,1] = w * lambda*(-y*(1-y)+(y-mu)^2)/denominator^2
            d2l.dthetas2[,1,2] = 
            d2l.dthetas2[,2,1] = w * (y-mu) / denominator
            d2l.dthetas2[,2,2] = w * (-0.25*trigamma((lambda+1)/2) +
                                       0.25*trigamma(1+lambda/2))
        } else {
            d2l.dthetas2 =  attr(eval.d3, "hessian")
        }

        wz = matrix(as.numeric(NA), n, dimm(M))  #3=dimm(M)
        wz[,iam(1,1,M)] = -d2l.dthetas2[,1,1] * dtheta.detas[,1]^2
        wz[,iam(2,2,M)] = -d2l.dthetas2[,2,2] * dtheta.detas[,2]^2
        wz[,iam(1,2,M)] = -d2l.dthetas2[,1,2] * dtheta.detas[,1] *
                                                dtheta.detas[,2]
        if(!.expected) {
            d2mudeta2 = d2theta.deta2(mu, .lmu, earg= .emu)
            d2lambda = d2theta.deta2(lambda, .llambda, earg= .elambda)
            wz[,iam(1,1,M)] = wz[,iam(1,1,M)] - dl.dthetas[,1] * d2mudeta2
            wz[,iam(2,2,M)] = wz[,iam(2,2,M)] - dl.dthetas[,2] * d2lambda
        }

        if(intercept.only) {
            sumw = sum(w)
            for(i in 1:ncol(wz))
                wz[,i] = sum(wz[,i]) / sumw
            pooled.weight = TRUE
            wz = w * wz   # Put back the weights
        } else
            pooled.weight = FALSE

        wz
    }), list( .lmu=lmu, .llambda=llambda, .expected=FALSE,
              .emu=emu, .elambda=elambda ))))
}





invbinomial = function(lrho="logit", llambda="loge",
                       erho=list(), elambda=list(),
                       irho=0.75, 
                       ilambda=NULL,
                       zero=NULL)
{

    if(mode(lrho) != "character" && mode(lrho) != "name")
        lrho = as.character(substitute(lrho))
    if(mode(llambda) != "character" && mode(llambda) != "name")
        llambda = as.character(substitute(llambda))
    if(length(zero) && !is.Numeric(zero, integer=TRUE, posit=TRUE))
        stop("bad input for argument \"zero\"")
    if(!is.list(erho)) erho = list()
    if(!is.list(elambda)) elambda = list()

    new("vglmff",
    blurb=c("Inverse binomial distribution\n\n",
            "Links:    ",
            namesof("rho", lrho, earg=erho), ", ", 
            namesof("lambda", llambda, earg=elambda), "\n", 
            "Mean:     lambda*(1-rho)/(2*rho-1)\n"),
    constraints=eval(substitute(expression({
        constraints = cm.zero.vgam(constraints, x, .zero, M)
    }), list( .zero=zero ))),
    initialize=eval(substitute(expression({
        if(ncol(cbind(y)) != 1)
            stop("response must be a vector or a one-column matrix")
        predictors.names =
        c(namesof("rho", .lrho, earg=.erho, tag=FALSE),
          namesof("lambda", .llambda, earg=.elambda, tag=FALSE))
        if(!length(etastart)) {
            rho = rep(if(length( .irho)) .irho else 0.75, length=n)
            lambda = rep(if(length( .ilambda)) .ilambda else 1, length=n)
            etastart = cbind(theta2eta(rho, .lrho, earg= .erho),
                             theta2eta(lambda, .llambda, earg= .elambda))
        }
    }), list( .llambda=llambda, .lrho=lrho,
              .elambda=elambda, .erho=erho,
              .ilambda=ilambda, .irho=irho ))),
    inverse=eval(substitute(function(eta, extra=NULL) {
        rho = eta2theta(eta[,1], .lrho, earg= .erho)
        lambda = eta2theta(eta[,2], .llambda, earg= .elambda)
        lambda*(1-rho)/(2*rho-1)
    }, list( .llambda=llambda, .lrho=lrho,
             .elambda=elambda, .erho=erho ))),
    last=eval(substitute(expression({
        misc$link = c(rho= .lrho, lambda= .llambda)
        misc$earg = list(rho= .erho, lambda= .elambda)
        misc$pooled.weight = pooled.weight
    }), list( .llambda=llambda, .lrho=lrho,
              .elambda=elambda, .erho=erho ))),
    loglikelihood=eval(substitute(
        function(mu, y, w, residuals= FALSE,eta, extra=NULL) {
        rho = eta2theta(eta[,1], .lrho, earg= .erho)
        lambda = eta2theta(eta[,2], .llambda, earg= .elambda)
        if(residuals) stop("loglikelihood residuals not implemented yet") else
        sum(w*(log(lambda) - lgamma(2*y+lambda) - lgamma(y+1) -
               lgamma(y+lambda+1) + y*log(rho*(1-rho)) + lambda*log(rho )))
    }, list( .llambda=llambda, .lrho=lrho,
             .elambda=elambda, .erho=erho ))),
    vfamily=c("invbinomial"),
    deriv=eval(substitute(expression({
        rho = eta2theta(eta[,1], .lrho, earg= .erho)
        lambda = eta2theta(eta[,2], .llambda, earg= .elambda)
        dl.drho = y * (1-2*rho)/(rho*(1-rho)) + lambda /rho
        dl.dlambda = 1/lambda - digamma(2*y+lambda) - digamma(y+lambda+1) +
                      log(rho)
        drho.deta = dtheta.deta(rho, .lrho, earg= .erho)
        dlambda.deta = dtheta.deta(lambda, .llambda, earg= .elambda)
        w * cbind( dl.drho * drho.deta, dl.dlambda * dlambda.deta )
    }), list( .llambda=llambda, .lrho=lrho,
              .elambda=elambda, .erho=erho ))),
    weight=eval(substitute(expression({
        d2l.drho2 = y * (-1+2*rho-2*rho^2) / (rho*(1-rho))^2 - lambda/rho^2
        d2l.dlambda2 = -1/(lambda^2) - trigamma(2*y+lambda) -
                        trigamma(y+lambda+1)
        d2l.dlambdarho = 1/rho
        wz = matrix(as.numeric(NA), n, dimm(M))  #3=dimm(M)
        wz[,iam(2,2,M)] = -d2l.dlambda2 * dlambda.deta^2
        wz[,iam(1,1,M)] = -d2l.drho2 * drho.deta^2
        wz[,iam(1,2,M)] = -d2l.dlambdarho * dlambda.deta * drho.deta

        d2rhodeta2 = d2theta.deta2(rho, .lrho, earg= .erho)
        d2lambda.deta2 = d2theta.deta2(lambda, .llambda, earg= .elambda)
        wz[,iam(1,1,M)] = wz[,iam(1,1,M)] - dl.dlambda * d2lambda.deta2
        wz[,iam(2,2,M)] = wz[,iam(2,2,M)] - dl.drho * d2rhodeta2
        wz = w * wz

        if(intercept.only) {
            sumw = sum(w)
            for(i in 1:ncol(wz))
                wz[,i] = sum(wz[,i]) / sumw
            pooled.weight = TRUE
            wz = w * wz   # Put back the weights
        } else
            pooled.weight = FALSE

        wz
    }), list( .llambda=llambda, .lrho=lrho,
              .elambda=elambda, .erho=erho ))))
}



genpoisson = function(llambda="logit", ltheta="loge",
                      elambda=list(), etheta=list(),
                      ilambda=0.5, itheta=NULL, zero=NULL)
{

    if(mode(llambda) != "character" && mode(llambda) != "name")
        llambda = as.character(substitute(llambda))
    if(mode(ltheta) != "character" && mode(ltheta) != "name")
        ltheta = as.character(substitute(ltheta))
    if(length(zero) && !is.Numeric(zero, integer=TRUE, posit=TRUE))
        stop("bad input for argument \"zero\"")
    if(!is.list(elambda)) elambda = list()
    if(!is.list(etheta)) etheta = list()

    new("vglmff",
    blurb=c("Generalized Poisson distribution\n\n",
            "Links:    ",
            namesof("lambda", llambda, earg=elambda), ", ", 
            namesof("theta", ltheta, earg=etheta), "\n", 
            "Mean:     theta / (1-lambda)\n",
            "Variance: theta / (1-lambda)^3"),
    constraints=eval(substitute(expression({
        constraints = cm.zero.vgam(constraints, x, .zero, M)
    }), list( .zero=zero ))),
    initialize=eval(substitute(expression({
        if(ncol(cbind(y)) != 1)
            stop("response must be a vector or a one-column matrix")
        predictors.names =
           c(namesof("lambda", .llambda, earg=.elambda, tag=FALSE),
             namesof("theta", .ltheta, earg=.etheta, tag=FALSE))
        if(!length(etastart)) {
            lambda = rep(if(length( .ilambda)) .ilambda else
                       0.5, length=n)
            theta = rep(if(length( .itheta)) .itheta else
                          median(y) * (1-lambda), length=n)
            etastart = cbind(theta2eta(lambda, .llambda, earg= .elambda),
                             theta2eta(theta, .ltheta, earg= .etheta))
        }
    }), list( .ltheta=ltheta, .llambda=llambda,
              .etheta=etheta, .elambda=elambda,
              .itheta=itheta, .ilambda=ilambda )) ),
    inverse=eval(substitute(function(eta, extra=NULL) {
        lambda = eta2theta(eta[,1], .llambda, earg= .elambda)
        theta = eta2theta(eta[,2], .ltheta, earg= .etheta)
        theta/(1-lambda)
    }, list( .ltheta=ltheta, .llambda=llambda,
             .etheta=etheta, .elambda=elambda ))),
    last=eval(substitute(expression({
        misc$link = c(lambda=.llambda, theta=.ltheta)
        misc$earg = list(lambda=.elambda, theta=.etheta)
        misc$pooled.weight = pooled.weight
    }), list( .ltheta=ltheta, .llambda=llambda,
              .etheta=etheta, .elambda=elambda ))),
    loglikelihood=eval(substitute(
        function(mu, y, w, residuals= FALSE,eta, extra=NULL) {
        lambda = eta2theta(eta[,1], .llambda, earg= .elambda)
        theta = eta2theta(eta[,2], .ltheta, earg= .etheta)
        index = y == 0 
        if(residuals) stop("loglikelihood residuals not implemented yet") else
        sum(w[index]*(-theta[index])) + 
        sum(w[!index]*(-y[!index]*lambda[!index]-theta[!index]+
            (y[!index]-1)*log(theta[!index]+y[!index]*lambda[!index]) +
            log(theta[!index] )))
    }, list( .ltheta=ltheta, .llambda=llambda,
             .etheta=etheta, .elambda=elambda ))),
    vfamily=c("genpoisson"),
    deriv=eval(substitute(expression({
        lambda = eta2theta(eta[,1], .llambda, earg= .elambda)
        theta = eta2theta(eta[,2], .ltheta, earg= .etheta)
        dl.dlambda = -y + y*(y-1)/(theta+y*lambda)
        dl.dtheta = -1 + (y-1)/(theta+y*lambda) + 1/theta
        dTHETA.deta = dtheta.deta(theta, .ltheta, earg= .etheta)
        dlambda.deta = dtheta.deta(lambda, .llambda, earg= .elambda)
        w * cbind( dl.dlambda * dlambda.deta, dl.dtheta * dTHETA.deta )
    }), list( .ltheta=ltheta, .llambda=llambda,
              .etheta=etheta, .elambda=elambda ))),
    weight=eval(substitute(expression({
        d2l.dlambda2 = -y^2 * (y-1) / (theta+y*lambda)^2
        d2l.dtheta2 = -(y-1)/(theta+y*lambda)^2 - 1 / theta^2
        d2l.dthetalambda =  -y * (y-1) / (theta+y*lambda)^2 
        wz = matrix(as.numeric(NA), n, dimm(M))  #3=dimm(M)
        wz[,iam(1,1,M)] = -d2l.dlambda2 * dlambda.deta^2
        wz[,iam(2,2,M)] = -d2l.dtheta2 * dTHETA.deta^2
        wz[,iam(1,2,M)] = -d2l.dthetalambda * dTHETA.deta * dlambda.deta

        d2THETA.deta2 = d2theta.deta2(theta, .ltheta, earg= .etheta)
        d2lambdadeta2 = d2theta.deta2(lambda, .llambda, earg= .elambda)
        wz[,iam(1,1,M)] = wz[,iam(1,1,M)] - dl.dlambda * d2lambdadeta2
        wz[,iam(2,2,M)] = wz[,iam(2,2,M)] - dl.dtheta * d2THETA.deta2
        wz = w * wz

        if(intercept.only) {
            sumw = sum(w)
            for(i in 1:ncol(wz))
                wz[,i] = sum(wz[,i]) / sumw
            pooled.weight = TRUE
            wz = w * wz   # Put back the weights
        } else
            pooled.weight = FALSE


        wz
    }), list( .ltheta=ltheta, .llambda=llambda,
              .etheta=etheta, .elambda=elambda ))))
}



lgammaff = function(link="loge", earg=list(), init.k=NULL)
{
    if(mode(link) != "character" && mode(link) != "name")
        link = as.character(substitute(link))
    if(!is.list(earg)) earg = list()

    new("vglmff",
    blurb=c("Log-gamma distribution f(y) = exp(ky - e^y)/gamma(k)), k>0\n\n",
            "Link:    ",
            namesof("k", link, earg=earg), "\n", "\n",
            "Mean:    digamma(k)", "\n"),
    initialize=eval(substitute(expression({
        if(ncol(cbind(y)) != 1)
            stop("response must be a vector or a one-column matrix")
        predictors.names = namesof("k", .link, earg=.earg, tag=FALSE) 
        if(!length(etastart)) {
            k.init = if(length( .init.k)) rep( .init.k, len=length(y)) else {
                medy = median(y)
                if(medy < 2) 5 else if(medy < 4) 20 else exp(0.7 * medy)
            }
            etastart = theta2eta(k.init, .link, earg= .earg)
        }
    }), list( .link=link, .earg=earg, .init.k=init.k ))),
    inverse=eval(substitute(function(eta, extra=NULL) {
        k = eta2theta(eta, .link, earg= .earg)
        digamma(k)
    }, list( .link=link, .earg=earg ))),
    last=eval(substitute(expression({
        misc$link = c(k= .link )
        misc$earg = list(k= .earg )
    }), list( .link=link, .earg=earg ))),
    loglikelihood=eval(substitute(
        function(mu,y,w,residuals= FALSE,eta, extra=NULL) {
        k = eta2theta(eta, .link, earg= .earg)
        if(residuals) stop("loglikelihood residuals not implemented yet") else
        sum(w * (k * y - exp(y) - lgamma(k )))
    }, list( .link=link, .earg=earg ))),
    vfamily=c("lgammaff"),
    deriv=eval(substitute(expression({
        k = eta2theta(eta, .link, earg= .earg) 
        dl.dk = y - digamma(k)
        dk.deta = dtheta.deta(k, .link, earg= .earg)
        w * dl.dk * dk.deta
    }), list( .link=link, .earg=earg ))),
    weight=eval(substitute(expression({
        ed2l.dk2 = trigamma(k)
        wz = w * dk.deta^2 * ed2l.dk2
        wz
    }), list( .link=link, .earg=earg ))))
}


dlgamma = function(x, location=0, scale=1, k=1) {
    if(!is.Numeric(scale, posit=TRUE)) stop("bad input for argument \"scale\"")
    if(!is.Numeric(k, posit=TRUE)) stop("bad input for argument \"k\"")
    z = (x-location) / scale
    exp(k * z - exp(z)) / (scale * gamma(k))
}
plgamma = function(q, location=0, scale=1, k=1) {
    if(!is.Numeric(scale, posit=TRUE)) stop("bad input for argument \"scale\"")
    if(!is.Numeric(k, posit=TRUE)) stop("bad input for argument \"k\"")
    z = (q-location)/scale
    pgamma(exp(z), k)
}
qlgamma = function(p, location=0, scale=1, k=1) {
    if(!is.Numeric(scale, posit=TRUE)) stop("bad input for argument \"scale\"")
    if(!is.Numeric(k, posit=TRUE)) stop("bad input for argument \"k\"")
    q = qgamma(p, k)
    location + scale * log(q)
}
rlgamma = function(n, location=0, scale=1, k=1) {
    if(!is.Numeric(n, posit=TRUE, integ=TRUE, allow=1)) 
        stop("bad input for argument \"n\"")
    if(!is.Numeric(scale, posit=TRUE)) stop("bad input for argument \"scale\"")
    if(!is.Numeric(k, posit=TRUE)) stop("bad input for argument \"k\"")
    y = rgamma(n, k)
    location + scale * log(y)
}


lgamma3ff = function(llocation="identity", lscale="loge", lshape="loge",
                     elocation=list(), escale=list(), eshape=list(),
                     ilocation=NULL, iscale=NULL, ishape=1, zero=NULL)
{
    if(mode(llocation) != "character" && mode(llocation) != "name")
        llocation = as.character(substitute(llocation))
    if(mode(lscale) != "character" && mode(lscale) != "name")
        lscale = as.character(substitute(lscale))
    if(mode(lshape) != "character" && mode(lshape) != "name")
        lshape = as.character(substitute(lshape))
    if(length(zero) && !is.Numeric(zero, integer=TRUE, posit=TRUE))
        stop("bad input for argument \"zero\"")
    if(!is.list(elocation)) elocation = list()
    if(!is.list(escale)) escale = list()
    if(!is.list(eshape)) eshape = list()

    new("vglmff",
    blurb=c("Log-gamma distribution",
            " f(y) = exp(k(y-a)/b - e^((y-a)/b))/(b*gamma(k)), ",
            "location=a, scale=b>0, shape=k>0\n\n",
            "Links:    ",
            namesof("location", llocation, earg=elocation), ", ",
            namesof("scale", lscale, earg=escale), ", ",
            namesof("shape", lshape, earg=eshape), "\n\n",
            "Mean:     a + b*digamma(k)", "\n"),
    constraints=eval(substitute(expression({
        constraints = cm.zero.vgam(constraints, x, .zero, M)
    }), list( .zero=zero ))),
    initialize=eval(substitute(expression({
        if(ncol(cbind(y)) != 1)
            stop("response must be a vector or a one-column matrix")
        predictors.names =
        c(namesof("location", .llocation, earg=.elocation, tag=FALSE),
          namesof("scale", .lscale, earg=.escale, tag=FALSE),
          namesof("shape", .lshape, earg=.eshape, tag=FALSE))
        if(!length(etastart)) {
            k.init = if(length( .ishape)) rep( .ishape, len=length(y)) else {
                rep(exp(median(y)), len=length(y))
            }
            scale.init = if(length( .iscale)) rep( .iscale, len=length(y)) else {
                rep(sqrt(var(y) / trigamma(k.init)), len=length(y))
            }
            loc.init = if(length( .iloc)) rep( .iloc, len=length(y)) else {
                rep(median(y) - scale.init * digamma(k.init), len=length(y))
            }
            etastart = cbind(theta2eta(loc.init, .llocation, earg= .elocation),
                             theta2eta(scale.init, .lscale, earg= .escale),
                             theta2eta(k.init, .lshape, earg= .eshape))
        }
    }), list( .llocation=llocation, .lscale=lscale, .lshape=lshape,
              .elocation=elocation, .escale=escale, .eshape=eshape,
              .iloc=ilocation, .iscale=iscale, .ishape=ishape ))),
    inverse=eval(substitute(function(eta, extra=NULL) {
        eta2theta(eta[,1], .llocation, earg= .elocation) +
        eta2theta(eta[,2], .lscale, earg= .escale) *
        digamma(eta2theta(eta[,3], .lshape, earg= .eshape))
    }, list( .llocation=llocation, .lscale=lscale, .lshape=lshape,
             .elocation=elocation, .escale=escale, .eshape=eshape ))),
    last=eval(substitute(expression({
        misc$link = c(location= .llocation, scale= .lscale, shape= .lshape)
        misc$earg = list(location= .elocation, scale= .escale, shape= .eshape)
    }), list( .llocation=llocation, .lscale=lscale, .lshape=lshape,
              .elocation=elocation, .escale=escale, .eshape=eshape ))),
    loglikelihood=eval(substitute(
        function(mu,y,w,residuals= FALSE,eta, extra=NULL) {
        a = eta2theta(eta[,1], .llocation, earg= .elocation)
        b = eta2theta(eta[,2], .lscale, earg= .escale)
        k = eta2theta(eta[,3], .lshape, earg= .eshape)
        zedd = (y-a)/b
        if(residuals) stop("loglikelihood residuals not implemented yet") else
        sum(w * (k * zedd - exp(zedd) - lgamma(k) - log(b )))
    }, list( .llocation=llocation, .lscale=lscale, .lshape=lshape,
             .elocation=elocation, .escale=escale, .eshape=eshape ))),
    vfamily=c("lgamma3ff"),
    deriv=eval(substitute(expression({
        a = eta2theta(eta[,1], .llocation, earg= .elocation)
        b = eta2theta(eta[,2], .lscale, earg= .escale)
        k = eta2theta(eta[,3], .lshape, earg= .eshape)
        zedd = (y-a)/b
        dl.da = (exp(zedd) - k) / b
        dl.db = (zedd * (exp(zedd) - k) - 1) / b
        dl.dk = zedd - digamma(k)
        da.deta = dtheta.deta(a, .llocation, earg= .elocation)
        db.deta = dtheta.deta(b, .lscale, earg= .escale)
        dk.deta = dtheta.deta(k, .lshape, earg= .eshape)
        w * cbind(dl.da * da.deta, dl.db * db.deta, dl.dk * dk.deta)
    }), list( .llocation=llocation, .lscale=lscale, .lshape=lshape,
              .elocation=elocation, .escale=escale, .eshape=eshape ))),
    weight=eval(substitute(expression({
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
        wz = w * wz
        wz
    }), list( .llocation=llocation, .lscale=lscale, .lshape=lshape,
              .elocation=elocation, .escale=escale, .eshape=eshape ))))
}

prentice74 = function(llocation="identity", lscale="loge", lshape="identity",
                      elocation=list(), escale=list(), eshape=list(),
                      ilocation=NULL, iscale=NULL, ishape=NULL, zero=NULL)
{
    if(mode(llocation) != "character" && mode(llocation) != "name")
        llocation = as.character(substitute(llocation))
    if(mode(lscale) != "character" && mode(lscale) != "name")
        lscale = as.character(substitute(lscale))
    if(mode(lshape) != "character" && mode(lshape) != "name")
        lshape = as.character(substitute(lshape))
    if(length(zero) && !is.Numeric(zero, integer=TRUE, posit=TRUE))
        stop("bad input for argument \"zero\"")
    if(!is.list(elocation)) elocation = list()
    if(!is.list(escale)) escale = list()
    if(!is.list(eshape)) eshape = list()

    new("vglmff",
    blurb=c("Log-gamma distribution (Prentice, 1974)",
            " f(y) = |q| * exp(w/q^2 - e^w) / (b*gamma(1/q^2)) ,\n",
            "w=(y-a)*q/b + digamma(1/q^2), location=a, scale=b>0, shape=q\n\n",
            "Links:    ",
            namesof("location", llocation, earg=elocation), ", ",
            namesof("scale", lscale, earg=escale), ", ",
            namesof("shape", lshape, earg=eshape), "\n", "\n",
            "Mean:     a", "\n"),
    constraints=eval(substitute(expression({
        constraints = cm.zero.vgam(constraints, x, .zero, M)
    }), list( .zero=zero ))),
    initialize=eval(substitute(expression({
        if(ncol(cbind(y)) != 1)
            stop("response must be a vector or a one-column matrix")
        predictors.names =
        c(namesof("location", .llocation, earg=.elocation, tag=FALSE),
          namesof("scale", .lscale, earg=.escale, tag=FALSE),
          namesof("shape", .lshape, earg=.eshape, tag=FALSE))
        if(!length(etastart)) {
            sdy = sqrt(var(y))
            k.init = if(length( .ishape)) rep( .ishape, len=length(y)) else {
                skewness = mean((y-mean(y))^3) / sdy^3 # <0 Left Skewed
                rep(-skewness, len=length(y))
            }
            scale.init = if(length( .iscale)) rep( .iscale, len=length(y)) else {
                rep(sdy, len=length(y))
            }
            loc.init = if(length( .iloc)) rep( .iloc, len=length(y)) else {
                rep(median(y), len=length(y))
            }
            etastart = cbind(theta2eta(loc.init, .llocation, earg= .elocation),
                             theta2eta(scale.init, .lscale, earg= .escale),
                             theta2eta(k.init, .lshape, earg= .eshape))
        }
    }), list( .llocation=llocation, .lscale=lscale, .lshape=lshape,
              .elocation=elocation, .escale=escale, .eshape=eshape,
              .iloc=ilocation, .iscale=iscale, .ishape=ishape ))),
    inverse=eval(substitute(function(eta, extra=NULL) {
        eta2theta(eta[,1], .llocation, earg= .elocation)
    }, list( .llocation=llocation, .lscale=lscale, .lshape=lshape,
             .elocation=elocation, .escale=escale, .eshape=eshape ))),
    last=eval(substitute(expression({
        misc$link = c(location= .llocation, scale= .lscale, shape= .lshape)
        misc$earg = list(location= .elocation, scale= .escale, shape= .eshape)
    }), list( .llocation=llocation, .lscale=lscale, .lshape=lshape,
              .elocation=elocation, .escale=escale, .eshape=eshape ))),
    loglikelihood=eval(substitute(
        function(mu,y,w,residuals= FALSE,eta, extra=NULL) {
        a = eta2theta(eta[,1], .llocation, earg= .elocation)
        b = eta2theta(eta[,2], .lscale, earg= .escale)
        k = eta2theta(eta[,3], .lshape, earg= .eshape)
        tmp55 = k^(-2)
        doubw = (y-a)*k/b + digamma(tmp55)
        if(residuals) stop("loglikelihood residuals not implemented yet") else
        sum(w*(log(abs(k)) -log(b) -lgamma(tmp55) + doubw*tmp55 -exp(doubw )))
    }, list( .llocation=llocation, .lscale=lscale, .lshape=lshape,
             .elocation=elocation, .escale=escale, .eshape=eshape ))),
    vfamily=c("prentice74"),
    deriv=eval(substitute(expression({
        a = eta2theta(eta[,1], .llocation, earg= .elocation)
        b = eta2theta(eta[,2], .lscale, earg= .escale)
        k = eta2theta(eta[,3], .lshape, earg= .eshape)
        tmp55 = k^(-2)
        mustar = digamma(tmp55)
        doubw = (y-a)*k/b + mustar
        sigmastar2 = trigamma(tmp55)
        dl.da = k*(exp(doubw) - tmp55) / b
        dl.db = ((doubw - mustar) * (exp(doubw) - tmp55) - 1) / b
        dl.dk = 1/k - 2 * (doubw - mustar) / k^3 - (exp(doubw) - tmp55) *
                ((doubw - mustar) / k - 2 * sigmastar2 / k^3)
        da.deta = dtheta.deta(a, .llocation, earg= .elocation)
        db.deta = dtheta.deta(b, .lscale, earg= .escale)
        dk.deta = dtheta.deta(k, .lshape, earg= .eshape)
        w * cbind(dl.da * da.deta, dl.db * db.deta, dl.dk * dk.deta)
    }), list( .llocation=llocation, .lscale=lscale, .lshape=lshape,
              .elocation=elocation, .escale=escale, .eshape=eshape ))),
    weight=eval(substitute(expression({
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
        wz = w * wz
        wz
    }), list( .llocation=llocation, .lscale=lscale, .lshape=lshape,
              .elocation=elocation, .escale=escale, .eshape=eshape ))))
}



dggamma = function(x, scale=1, d=1, k=1) {
    if(!is.Numeric(scale, posit=TRUE)) stop("bad input for argument \"scale\"")
    if(!is.Numeric(d, posit=TRUE)) stop("bad input for argument \"d\"")
    if(!is.Numeric(k, posit=TRUE)) stop("bad input for argument \"k\"")
    N = max(length(x), length(scale), length(d), length(k))
    x = rep(x, len=N); scale = rep(scale, len=N);
    d = rep(d, len=N); k = rep(k, len=N); 
    ans = rep(0.0, len=N)
    ind = x > 0
    if(any(ind)) {
        z = (x[ind]/scale[ind])^d[ind]
        ans[ind] = d[ind] * scale[ind]^(-d[ind]*k[ind]) *
                   x[ind]^(d[ind]*k[ind]-1) * exp(-z) / gamma(k[ind])
    }
    ans
}
pggamma = function(q, scale=1, d=1, k=1) {
    if(!is.Numeric(scale, posit=TRUE)) stop("bad input for argument \"scale\"")
    if(!is.Numeric(d, posit=TRUE)) stop("bad input for argument \"d\"")
    if(!is.Numeric(k, posit=TRUE)) stop("bad input for argument \"k\"")
    z = (q/scale)^d
    pgamma(z, k)
}
qggamma = function(p, scale=1, d=1, k=1) {
    if(!is.Numeric(scale, posit=TRUE)) stop("bad input for argument \"scale\"")
    if(!is.Numeric(d, posit=TRUE)) stop("bad input for argument \"d\"")
    if(!is.Numeric(k, posit=TRUE)) stop("bad input for argument \"k\"")
    q = qgamma(p, k)
    scale * q^(1/d)
}
rggamma = function(n, scale=1, d=1, k=1) {
    if(!is.Numeric(n, posit=TRUE, integ=TRUE, allow=1)) 
        stop("bad input for \"n\"")
    if(!is.Numeric(scale, posit=TRUE)) stop("bad input for \"scale\"")
    if(!is.Numeric(d, posit=TRUE)) stop("bad input for \"d\"")
    if(!is.Numeric(k, posit=TRUE)) stop("bad input for \"k\"")
    y = rgamma(n, k)
    scale * y^(1/d)
}

ggamma = function(lscale="loge", ld="loge", lk="loge",
                  escale=list(), ed=list(), ek=list(),
                  iscale=NULL, id=NULL, ik=NULL, zero=NULL)
{
    if(mode(lscale) != "character" && mode(lscale) != "name")
        lscale = as.character(substitute(lscale))
    if(mode(ld) != "character" && mode(ld) != "name")
        ld = as.character(substitute(ld))
    if(mode(lk) != "character" && mode(lk) != "name")
        lk = as.character(substitute(lk))
    if(length(zero) && !is.Numeric(zero, integer=TRUE, posit=TRUE))
        stop("bad input for argument \"zero\"")
    if(!is.list(escale)) escale = list()
    if(!is.list(ed)) ed = list()
    if(!is.list(ek)) ek = list()

    new("vglmff",
    blurb=c("Generalized gamma distribution",
            " f(y) = d * b^(-d*k) * y^(d*k-1) * exp(-(y/b)^d) /  gamma(k),\n",
            "scale=b>0, d>0, k>0, y>0\n\n",
            "Links:    ",
            namesof("scale", lscale, earg=escale), ", ",
            namesof("d", ld, earg=ed), ", ",
            namesof("k", lk, earg=ek), "\n", "\n",
            "Mean:     b*k", "\n"),
    constraints=eval(substitute(expression({
        constraints = cm.zero.vgam(constraints, x, .zero, M)
    }), list( .zero=zero ))),
    initialize=eval(substitute(expression({
        if(ncol(cbind(y)) != 1)
            stop("response must be a vector or a one-column matrix")
        if(any(y <= 0)) stop("response must be have positive values only")
        predictors.names = 
            c(namesof("scale", .lscale, earg=.escale, tag=FALSE),
              namesof("d", .ld, earg=.ed, tag=FALSE),
              namesof("k", .lk, earg=.ek, tag=FALSE))
        if(!length(etastart)) {
            b.init = if(length( .iscale)) rep( .iscale, len=length(y)) else {
                rep(mean(y^2) / mean(y), len=length(y))
            }
            k.init = if(length( .ik)) rep( .ik, len=length(y)) else {
                rep(mean(y) / b.init, len=length(y))
            }
            d.init = if(length( .id)) rep( .id, len=length(y)) else {
                rep(digamma(k.init) / mean(log(y/b.init)), len=length(y))
            }
            etastart = cbind(theta2eta(b.init, .lscale, earg= .escale),
                             theta2eta(d.init, .ld, earg= .ed),
                             theta2eta(k.init, .lk, earg= .ek))
        }
    }), list( .lscale=lscale, .ld=ld, .lk=lk,
              .escale=escale, .ed=ed, .ek=ek,
              .iscale=iscale, .id=id, .ik=ik ))),
    inverse=eval(substitute(function(eta, extra=NULL) {
        b = eta2theta(eta[,1], .lscale, earg= .escale)
        k = eta2theta(eta[,3], .lk, earg= .ek)
        b * k
    }, list( .ld=ld, .lscale=lscale, .lk=lk,
             .escale=escale, .ed=ed, .ek=ek ))),
    last=eval(substitute(expression({
        misc$link = c(scale= .lscale, d= .ld, k= .lk)
        misc$earg = list(scale= .escale, d= .ed, k= .ek)
    }), list( .lscale=lscale, .ld=ld, .lk=lk,
              .escale=escale, .ed=ed, .ek=ek ))),
    loglikelihood=eval(substitute(
        function(mu,y,w,residuals= FALSE,eta, extra=NULL) {
        b = eta2theta(eta[,1], .lscale, earg= .escale)
        d = eta2theta(eta[,2], .ld, earg= .ed)
        k = eta2theta(eta[,3], .lk, earg= .ek)
        if(residuals) stop("loglikelihood residuals not implemented yet") else
        sum(w*(log(d) - lgamma(k) + (d*k-1) * log(y) - d*k*log(b) - (y/b)^d))
    }, list( .lscale=lscale, .ld=ld, .lk=lk,
             .escale=escale, .ed=ed, .ek=ek ))),
    vfamily=c("ggamma"),
    deriv=eval(substitute(expression({
        b = eta2theta(eta[,1], .lscale, earg= .escale)
        d = eta2theta(eta[,2], .ld, earg= .ed)
        k = eta2theta(eta[,3], .lk, earg= .ek)
        tmp22 = (y/b)^d
        tmp33 = log(y/b)
        dl.db = d * (tmp22 - k) / b
        dl.dd = 1/d + tmp33 * (k - tmp22)
        dl.dk = d * tmp33 - digamma(k)
        db.deta = dtheta.deta(b, .lscale, earg= .escale)
        dd.deta = dtheta.deta(d, .ld, earg= .ed)
        dk.deta = dtheta.deta(k, .lk, earg= .ek)
        w * cbind(dl.db * db.deta, dl.dd * dd.deta, dl.dk * dk.deta)
    }), list( .lscale=lscale, .ld=ld, .lk=lk,
              .escale=escale, .ed=ed, .ek=ek ))),
    weight=eval(substitute(expression({
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
        wz = w * wz
        wz
    }), list( .lscale=lscale, .ld=ld, .lk=lk,
              .escale=escale, .ed=ed, .ek=ek ))))
}


dlog = function(x, prob) {
    if(!is.Numeric(x)) stop("bad input for argument \"x\"")
    if(!is.Numeric(prob, posit=TRUE) || max(prob) >= 1)
        stop("bad input for argument \"prob\"")
    N = max(length(x), length(prob))
    x = rep(x, len=N); prob = rep(prob, len=N); 
    ox = !is.finite(x)
    zero = round(x) != x | x < 1
    ans = rep(0.0, len=length(x))
    if(any(!zero))
        ans[!zero] = -(prob[!zero]^(x[!zero])) / (x[!zero] * log1p(-prob[!zero]))
    if(any(ox))
        ans[ox] = NA
    ans
}

plog = function(q, prob, log.p=FALSE) {
    if(!is.Numeric(q)) stop("bad input for argument \"q\"")
    if(!is.Numeric(prob, posit=TRUE) || max(prob) >= 1)
        stop("bad input for argument \"prob\"")
    N = max(length(q), length(prob))
    q = rep(q, len=N); prob = rep(prob, len=N);
    ans = q * 0  # Retains names(q)
    if(max(abs(prob-prob[1])) < 1.0e-08) {
        qstar = floor(q)
        temp = if(max(qstar) >= 1) dlog(x=1:max(qstar), 
               prob=prob[1]) else 0*qstar
        unq = unique(qstar)
        for(i in unq) {
            index = qstar == i
            ans[index] = if(i >= 1) sum(temp[1:i]) else 0
        }
    } else
    for(i in 1:N) {
        qstar = floor(q[i])
        ans[i] = if(qstar >= 1) sum(dlog(x=1:qstar, prob=prob[i])) else 0
    }
    if(log.p) log(ans) else ans
}

rlog = function(n, prob, Smallno=1.0e-6) {
    if(!is.Numeric(n, posit=TRUE, integ=TRUE))
        stop("bad input for argument \"n\"")
    if(!is.Numeric(prob, allow=1, posit=TRUE) || max(prob) >= 1)
        stop("bad input for argument \"prob\"")
    if(!is.Numeric(Smallno, posit=TRUE, allow=1) || Smallno > 0.01 ||
       Smallno < 2 * .Machine$double.eps)
        stop("bad input for argument \"Smallno\"")
    ans = rep(0.0, len=n)

    ptr1 = 1; ptr2 = 0
    a = -1 / log1p(-prob)
    mean = a*prob/(1-prob)    # E(Y)
    sigma = sqrt(a*prob*(1-a*prob)) / (1-prob)   # sd(Y)
    ymax = dlog(x=1, prob)
    while(ptr2 < n) {
        Lower = 0.5 # A continuity correction is used = 1 - 0.5.
        Upper = mean + 5 * sigma
        while(plog(q=Upper, prob) < 1-Smallno)
            Upper = Upper + sigma
        Upper = Upper + 0.5
        x = round(runif(2*n, min=Lower, max=Upper))
        index = runif(2*n, max=ymax) < dlog(x,prob)
        sindex = sum(index)
        if(sindex) {
            ptr2 = min(n, ptr1 + sindex - 1)
            ans[ptr1:ptr2] = (x[index])[1:(1+ptr2-ptr1)]
            ptr1 = ptr2 + 1
        }
    }
    ans
}


logff = function(link="logit", earg=list(), init.c=NULL)
{
    if(length(init.c) &&
       (!is.Numeric(init.c, posit=TRUE) || max(init.c) >= 1))
        stop("init.c must be in (0,1)")
    if(mode(link) != "character" && mode(link) != "name")
        link = as.character(substitute(link))
    if(!is.list(earg)) earg = list()

    new("vglmff",
    blurb=c("Logarithmic distribution f(y) = a * c^y / y, y=1,2,3,...,\n",
            "            0 < c < 1, a = -1 / log(1-c)  \n\n",
            "Link:    ", namesof("c", link, earg=earg), "\n", "\n",
            "Mean:    a * c / (1 - c)", "\n"),
    initialize=eval(substitute(expression({
        if(ncol(cbind(y)) != 1)
            stop("response must be a vector or a one-column matrix")
        predictors.names = namesof("c", .link, earg=.earg, tag=FALSE) 
        if(!length(etastart)) {
            llfun = function(cc, y, w) {
                a = -1 / log1p(-cc)
                sum(w * (log(a) + y * log(cc) - log(y)))
            }
            c.init = if(length( .init.c )) .init.c else
                getInitVals(gvals=seq(0.05, 0.95, len=9), llfun=llfun, y=y, w=w)
            c.init = rep(c.init, length=length(y))
            etastart = theta2eta(c.init, .link, earg= .earg)
        }
    }), list( .link=link, .earg=earg, .init.c=init.c ))),
    inverse=eval(substitute(function(eta, extra=NULL) {
        cc = eta2theta(eta, .link, earg= .earg)
        a = -1 / log1p(-cc)
        a * cc / (1-cc)
    }, list( .link=link, .earg=earg ))),
    last=eval(substitute(expression({
        misc$link = c(c= .link)
        misc$earg = list(c= .earg)
    }), list( .link=link, .earg=earg ))),
    loglikelihood=eval(substitute(
        function(mu,y,w,residuals= FALSE,eta, extra=NULL) {
        cc = eta2theta(eta, .link, earg= .earg)
        a = -1 / log1p(-cc)
        if(residuals) stop("loglikelihood residuals not implemented yet") else
        sum(w * (log(a) + y * log(cc) - log(y )))
    }, list( .link=link, .earg=earg ))),
    vfamily=c("logff"),
    deriv=eval(substitute(expression({
        cc = eta2theta(eta, .link, earg= .earg)
        a = -1 / log1p(-cc)
        dl.dc = 1 / ((1-cc) * log1p(-cc)) + y / cc
        dc.deta = dtheta.deta(cc, .link, earg= .earg)
        w * dl.dc * dc.deta
    }), list( .link=link, .earg=earg ))),
    weight=eval(substitute(expression({
        ed2l.dc2 = a * (1 - a * cc) / (cc * (1-cc)^2)
        wz = w * dc.deta^2 * ed2l.dc2
        wz
    }), list( .link=link, .earg=earg ))))
}


levy = function(delta=NULL, link.gamma="loge",
                earg=list(), idelta=NULL, igamma=NULL)
{



    delta.known = is.Numeric(delta, allow=1)
    if(mode(link.gamma) != "character" && mode(link.gamma) != "name")
        link.gamma = as.character(substitute(link.gamma))
    if(!is.list(earg)) earg = list()

    new("vglmff",
    blurb=c("Levy distribution f(y) = sqrt(gamma/(2*pi)) * ",
            "(y-delta)^(-3/2) * \n",
            "          exp(-gamma / (2*(y-delta ))),\n",
            "          delta < y, gamma > 0",
            if(delta.known) paste(", delta = ", delta, ",", sep=""),
            "\n\n",
            if(delta.known) "Link:    " else "Links:   ",
            namesof("gamma", link.gamma, earg=earg),
            if(! delta.known) 
                c(", ", namesof("delta", "identity", earg=list())),
            "\n\n",
            "Mean:    NA", 
            "\n"),
    initialize=eval(substitute(expression({
        if(ncol(cbind(y)) != 1)
            stop("response must be a vector or a one-column matrix")
        predictors.names = 
            c(namesof("gamma", .link.gamma, earg=.earg, tag=FALSE),
              if( .delta.known) NULL else 
              namesof("delta", "identity", earg=list(), tag=FALSE))

        if(!length(etastart)) {
            delta.init = if( .delta.known) {
                           if(min(y,na.rm= TRUE) <= .delta)
                               stop("delta must be < min(y)")
                           .delta 
                         } else {
                           if(length( .idelta)) .idelta else
                               min(y,na.rm= TRUE) - 1.0e-4 *
                               diff(range(y,na.rm= TRUE))
                         }
            gamma.init = if(length( .igamma)) .igamma else
                         median(y - delta.init) # = 1/median(1/(y-delta.init))
            gamma.init = rep(gamma.init, length=length(y))
            etastart = cbind(theta2eta(gamma.init, .link.gamma, earg= .earg),
                             if( .delta.known) NULL else delta.init)
                             
        }
    }), list( .link.gamma=link.gamma, .earg=earg,
             .delta.known=delta.known,
             .delta=delta,
             .idelta=idelta,
             .igamma=igamma ))),
    inverse=eval(substitute(function(eta, extra=NULL) {
        eta = as.matrix(eta)
        mygamma = eta2theta(eta[,1], .link.gamma, earg= .earg)
        delta = if( .delta.known) .delta else eta[,2]


        NA * mygamma
    }, list( .link.gamma=link.gamma, .earg=earg,
             .delta.known=delta.known,
             .delta=delta ))),
    last=eval(substitute(expression({
        misc$link = if( .delta.known) NULL else c(delta="identity")
        misc$link = c(gamma = .link.gamma, misc$link)
        misc$earg = if( .delta.known) list(gamma = .earg) else
                    list(gamma = .earg, delta=list())
        if( .delta.known)
            misc$delta = .delta
    }), list( .link.gamma=link.gamma, .earg=earg,
             .delta.known=delta.known,
             .delta=delta ))),
    loglikelihood=eval(substitute(
        function(mu,y,w,residuals= FALSE,eta, extra=NULL) {
        eta = as.matrix(eta)
        mygamma = eta2theta(eta[,1], .link.gamma, earg= .earg)
        delta = if( .delta.known) .delta else eta[,2]
        if(residuals) stop("loglikelihood residuals not implemented yet") else
        sum(w * 0.5 * (log(mygamma) -3*log(y-delta) - mygamma / (y-delta )))
    }, list( .link.gamma=link.gamma, .earg = earg,
             .delta.known=delta.known,
             .delta=delta ))),
    vfamily=c("levy"),
    deriv=eval(substitute(expression({
        eta = as.matrix(eta)
        mygamma = eta2theta(eta[,1], .link.gamma, earg= .earg)
        delta = if( .delta.known) .delta else eta[,2]
        if(! .delta.known)
            dl.ddelta  = (3 - mygamma / (y-delta)) / (2 * (y-delta))
        dl.dgamma = 0.5 * (1 / mygamma - 1 / (y-delta))
        dgamma.deta = dtheta.deta(mygamma, .link.gamma, earg= .earg)
        w * cbind(dl.dgamma * dgamma.deta, 
                  if( .delta.known) NULL else dl.ddelta)
    }), list( .link.gamma=link.gamma, .earg=earg,
             .delta.known=delta.known,
             .delta=delta ))),
    weight=eval(substitute(expression({
        wz = matrix(as.numeric(NA), n, dimm(M))   # M = if(delta is known) 1 else 2
        wz[,iam(1,1,M)] = 1 * dgamma.deta^2 
        if(! .delta.known) {
            wz[,iam(1,2,M)] =  3 * dgamma.deta
            wz[,iam(2,2,M)] =  21
        }
        wz = w * wz / (2 * mygamma^2) 
        wz
    }), list( .link.gamma=link.gamma, .earg=earg,
             .delta.known=delta.known,
             .delta=delta ))))
}


        

if(FALSE) 
stoppa = function(y0,
                  link.alpha="loge",
                  link.theta="loge", ealpha=list(), etheta=list(),
                  ialpha=NULL,
                  itheta=1.0,
                  zero=NULL)
{
    if(!is.Numeric(y0, allo=1) || y0 <= 0)
        stop("y0 must be a positive value")

    if(mode(link.alpha) != "character" && mode(link.alpha) != "name")
        link.alpha = as.character(substitute(link.alpha))
    if(mode(link.theta) != "character" && mode(link.theta) != "name")
        link.theta = as.character(substitute(link.theta))
    if(length(zero) && !is.Numeric(zero, integer=TRUE, posit=TRUE))
        stop("bad input for argument \"zero\"")
    if(!is.list(ealpha)) ealpha = list()
    if(!is.list(etheta)) etheta = list()

    new("vglmff",
    blurb=c("Stoppa distribution\n\n",
            "Links:    ",
            namesof("alpha", link.alpha, earg=ealpha), ", ", 
            namesof("theta", link.theta, earg=etheta), "\n", 
            if(is.R()) "Mean:     theta*y0*beta(1-1/alpha, theta)" else
                       "Mean:     theta*y0*beta(1-1/alpha, theta)"),
    constraints=eval(substitute(expression({
        constraints = cm.zero.vgam(constraints, x, .zero, M)
    }), list( .zero=zero ))),
    initialize=eval(substitute(expression({
        predictors.names = 
        c(namesof("alpha", .link.alpha, earg=.ealpha, tag=FALSE),
          namesof("theta", .link.theta, earg=.etheta, tag=FALSE))

        y0 = .y0 
        if(min(y) < y0) stop("y0 must lie in the interval (0, min(y))")
        if(!length( .ialpha) || !length( .itheta)) {
            qvec = c(.25, .5, .75)   # Arbitrary; could be made an argument
            init.theta = if(length( .itheta)) .itheta else 1
            xvec = log1p(-qvec^(1/init.theta))
            fit0 = lsfit(x=xvec, y=log(quantile(y, qvec))-log(y0), intercept=FALSE)
        }

        extra$y0 = y0
        if(!length(etastart)) {
            alpha = rep(if(length( .ialpha)) .ialpha else -1/fit0$coef[1], length=n)
            theta = rep(if(length( .itheta)) .itheta else 1.0, length=n)
            etastart = cbind(theta2eta(alpha, .link.alpha, earg= .ealpha),
                             theta2eta(theta, .link.theta, earg= .etheta))
        }
    }), list( .link.theta=link.theta, .link.alpha=link.alpha,
            .y0=y0,
            .itheta=itheta, .ialpha=ialpha ))),
    inverse=eval(substitute(function(eta, extra=NULL) {
        alpha = eta2theta(eta[,1], .link.alpha, earg= .ealpha)
        theta = eta2theta(eta[,2], .link.theta, earg= .etheta)
        theta * extra$y0 * beta(1-1/alpha, theta)
    }, list( .link.theta=link.theta, .link.alpha=link.alpha ))),
    last=eval(substitute(expression({
        misc$link = c(alpha= .link.alpha, theta= .link.theta)
    }), list( .link.theta=link.theta, .link.alpha=link.alpha ))),
    loglikelihood=eval(substitute(
            function(mu, y, w, residuals= FALSE,eta, extra=NULL) {
        alpha = eta2theta(eta[,1], .link.alpha, earg= .ealpha)
        theta = eta2theta(eta[,2], .link.theta, earg= .etheta)
        if(residuals) stop("loglikelihood residuals not implemented yet") else
        sum(w*(log(theta*alpha) + alpha*log(extra$y0) -(alpha+1)*log(y)+
               (theta-1) * log1p(-(y/extra$y0)^(-alpha))))
    }, list( .link.theta=link.theta, .link.alpha=link.alpha ))),
    vfamily=c("stoppa"),
    deriv=eval(substitute(expression({
        alpha = eta2theta(eta[,1], .link.alpha, earg= .ealpha)
        theta = eta2theta(eta[,2], .link.theta, earg= .etheta)
        temp8  = (y / extra$y0)^(-alpha)
        temp8a = log(temp8)
        temp8b = log1p(-temp8)
        dl.dalpha = 1/alpha - log(y/extra$y0) + (theta-1) * temp8 *
                    log(y / extra$y0) / (1-temp8)
        dl.dtheta = 1/theta + temp8b
        dalpha.deta = dtheta.deta(alpha, .link.alpha, earg= .ealpha)
        dTHETA.deta = dtheta.deta(theta, .link.theta, earg= .etheta)
        w * cbind( dl.dalpha * dalpha.deta, dl.dtheta * dTHETA.deta )
    }), list( .link.theta=link.theta, .link.alpha=link.alpha ))),
    weight=eval(substitute(expression({
        ed2l.dalpha = 1/alpha^2 + theta * (2 * log(extra$y0) * (digamma(2)-
                      digamma(theta+4)) -
                      (trigamma(1)+trigamma(theta+3)) / alpha^3) /
                      (alpha * (theta+1) * (theta+2) / n) # zz / sum(w)
        ed2l.dtheta = 1 / theta^2
        ed2l.dalphatheta = (digamma(2)-digamma(theta+2)) / (alpha*(theta+1))
        wz = matrix(as.numeric(NA), n, dimm(M))  #3=dimm(M)
        wz[,iam(1,1,M)] = ed2l.dalpha * dalpha.deta^2
        wz[,iam(2,2,M)] = ed2l.dtheta * dTHETA.deta^2
        wz[,iam(1,2,M)] = ed2l.dalpha * dTHETA.deta * dalpha.deta
        wz = w * wz
        wz
    }), list( .link.theta=link.theta, .link.alpha=link.alpha ))) )
}



dlino = function(x, shape1, shape2, lambda=1) {
    if(!is.Numeric(x)) stop("bad input for argument \"x\"")
    if(!is.Numeric(shape1, posit=TRUE)) 
        stop("bad input for argument \"shape1\"")
    if(!is.Numeric(shape2, posit=TRUE)) 
        stop("bad input for argument \"shape2\"")
    if(!is.Numeric(lambda, posit=TRUE)) 
        stop("bad input for argument \"lambda\"")
    dbeta(x=x, shape1=shape1, shape2=shape2) * lambda^shape1 /
        (1 - (1-lambda)*x)^(shape1+shape2)
}

plino = function(q, shape1, shape2, lambda=1) {
    if(!is.Numeric(q)) stop("bad input for \"q\"")
    if(!is.Numeric(shape1, posit=TRUE)) 
        stop("bad input for argument \"shape1\"")
    if(!is.Numeric(shape2, posit=TRUE)) 
        stop("bad input for argument \"shape2\"")
    if(!is.Numeric(lambda, posit=TRUE)) 
        stop("bad input for argument \"lambda\"")
    pbeta(q=lambda*q/(1 - (1-lambda)*q), shape1=shape1, shape2=shape2)
}

qlino = function(p, shape1, shape2, lambda=1) {
    if(!is.Numeric(p, posit=TRUE) || any(p >= 1)) 
        stop("bad input for argument \"p\"")
    if(!is.Numeric(shape1, posit=TRUE)) 
        stop("bad input for argument \"shape1\"")
    if(!is.Numeric(lambda, posit=TRUE)) 
        stop("bad input for argument \"lambda\"")
    Y = qbeta(p=p, shape1=shape1, shape2=shape2)
    Y / (lambda + (1-lambda)*Y)
}


rlino = function(n, shape1, shape2, lambda=1) {
    if(!is.Numeric(n, posit=TRUE, integ=TRUE, allow=1)) 
        stop("bad input for argument \"n\"")
    if(!is.Numeric(shape1, posit=TRUE)) 
        stop("bad input for argument \"shape1\"")
    if(!is.Numeric(shape2, posit=TRUE)) 
        stop("bad input for argument \"shape2\"")
    if(!is.Numeric(lambda, posit=TRUE)) 
        stop("bad input for argument \"lambda\"")
    Y = rbeta(n=n, shape1=shape1, shape2=shape2)
    Y / (lambda + (1-lambda)*Y)
}



lino = function(lshape1="loge",
                lshape2="loge",
                llambda="loge",
                eshape1=list(), eshape2=list(), elambda=list(),
                ishape1=NULL, ishape2=NULL, ilambda=1, zero=NULL)
{
    if(mode(lshape1) != "character" && mode(lshape1) != "name")
        lshape1 = as.character(substitute(lshape1))
    if(mode(lshape2) != "character" && mode(lshape2) != "name")
        lshape2 = as.character(substitute(lshape2))
    if(mode(llambda) != "character" && mode(llambda) != "name")
        llambda = as.character(substitute(llambda))
    if(length(zero) && !is.Numeric(zero, integer=TRUE, posit=TRUE))
        stop("bad input for argument \"zero\"")
    if(!is.Numeric(ilambda, positive=TRUE))
        stop("bad input for argument \"ilambda\"")
    if(!is.list(eshape1)) eshape1 = list()
    if(!is.list(eshape2)) eshape2 = list()
    if(!is.list(elambda)) elambda = list()

    new("vglmff",
    blurb=c("Generalized Beta distribution (Libby and Novick, 1982)\n\n",
            "Links:    ",
            namesof("shape1", lshape1, earg=eshape1), ", ", 
            namesof("shape2", lshape2, earg=eshape2), ", ", 
            namesof("lambda", llambda, earg=elambda), "\n", 
            "Mean:     something complicated"),
    constraints=eval(substitute(expression({
        constraints = cm.zero.vgam(constraints, x, .zero, M)
    }), list( .zero=zero ))),
    initialize=eval(substitute(expression({
        predictors.names = 
        c(namesof("shape1", .lshape1, earg=.eshape1, tag=FALSE),
          namesof("shape2", .lshape2, earg=.eshape2, tag=FALSE),
          namesof("lambda", .llambda, earg=.elambda, tag=FALSE))
        if(min(y) <= 0 || max(y) >= 1)
            stop("values of the response must be between 0 and 1 (0,1)")
        if(ncol(cbind(y)) != 1)
            stop("response must be a vector or a one-column matrix")
        if(!length(etastart)) {
            lambda.init = rep(if(length( .ilambda )) .ilambda else 1, length=n)
            sh1.init = if(length( .ishape1 )) rep( .ishape1, length=n) else NULL
            sh2.init = if(length( .ishape2 )) rep( .ishape2, length=n) else NULL
            txY.init = lambda.init * y / (1+lambda.init*y - y)
            mean1 = mean(txY.init)
            mean2 = mean(1/txY.init)
            if(!is.Numeric(sh1.init))
                sh1.init = rep((mean2 - 1) / (mean2 - 1/mean1), length=n)
            if(!is.Numeric(sh2.init))
                sh2.init = rep(sh1.init * (1-mean1) / mean1, length=n)
            etastart = cbind(theta2eta(sh1.init, .lshape1, earg= .eshape1),
                             theta2eta(sh2.init, .lshape2, earg= .eshape2),
                             theta2eta(lambda.init, .llambda, earg= .elambda))
        }
    }), list( .lshape1=lshape1, .lshape2=lshape2, .llambda=llambda,
              .eshape1=eshape1, .eshape2=eshape2, .elambda=elambda,
              .ishape1=ishape1, .ishape2=ishape2, .ilambda=ilambda ))),
    inverse=eval(substitute(function(eta, extra=NULL) {
        sh1 = eta2theta(eta[,1], .lshape1, earg= .eshape1)
        sh2 = eta2theta(eta[,2], .lshape2, earg= .eshape2)
        lambda = eta2theta(eta[,3], .llambda, earg= .elambda)
        rep(as.numeric(NA), length=nrow(eta))
    }, list( .lshape1=lshape1, .lshape2=lshape2, .llambda=llambda,
             .eshape1=eshape1, .eshape2=eshape2, .elambda=elambda ))),
    last=eval(substitute(expression({
        misc$link = c(shape1 = .lshape1, shape2 = .lshape2, lambda = .llambda)
        misc$earg =list(shape1 = .eshape1, shape2 = .eshape2, lambda = .elambda)
    }), list( .lshape1=lshape1, .lshape2=lshape2, .llambda=llambda,
              .eshape1=eshape1, .eshape2=eshape2, .elambda=elambda ))),
    loglikelihood=eval(substitute(
            function(mu, y, w, residuals= FALSE,eta, extra=NULL) {
        sh1 = eta2theta(eta[,1], .lshape1, earg= .eshape1)
        sh2 = eta2theta(eta[,2], .lshape2, earg= .eshape2)
        lambda = eta2theta(eta[,3], .llambda, earg= .elambda)
        if(!is.R()) lbeta = function(a,b) lgamma(a) + lgamma(b) - lgamma(a+b)
        if(residuals) stop("loglikelihood residuals not implemented yet") else
        sum(w*(sh1*log(lambda) + (sh1-1)*log(y) + (sh2-1)*log1p(-y) -
               lbeta(sh1,sh2) -(sh1+sh2)*log1p(-(1-lambda)*y)) )
    }, list( .lshape1=lshape1, .lshape2=lshape2, .llambda=llambda,
             .eshape1=eshape1, .eshape2=eshape2, .elambda=elambda ))),
    vfamily=c("lino"),
    deriv=eval(substitute(expression({
        sh1 = eta2theta(eta[,1], .lshape1, earg= .eshape1)
        sh2 = eta2theta(eta[,2], .lshape2, earg= .eshape2)
        lambda = eta2theta(eta[,3], .llambda, earg= .elambda)
        temp1 = log1p(-(1-lambda) * y)
        temp2 = digamma(sh1+sh2)
        dl.dsh1 = log(lambda) + log(y) - digamma(sh1) + temp2 - temp1
        dl.dsh2 = log1p(-y) - digamma(sh2) + temp2 - temp1
        dl.dlambda = sh1/lambda - (sh1+sh2) * y / (1 - (1-lambda) * y)
        dsh1.deta = dtheta.deta(sh1, .lshape1, earg= .eshape1)
        dsh2.deta = dtheta.deta(sh2, .lshape2, earg= .eshape2)
        dlambda.deta = dtheta.deta(lambda, .llambda, earg= .elambda)
        w * cbind( dl.dsh1 * dsh1.deta,
                   dl.dsh2 * dsh2.deta,
                   dl.dlambda * dlambda.deta)
    }), list( .lshape1=lshape1, .lshape2=lshape2, .llambda=llambda,
              .eshape1=eshape1, .eshape2=eshape2, .elambda=elambda ))),
    weight=eval(substitute(expression({
        if(!is.R()) beta = function(a,b) (gamma(a) / gamma(a+b)) * gamma(b)
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
        wz = w * wz
        wz
    }), list( .lshape1=lshape1, .lshape2=lshape2, .llambda=llambda,
              .eshape1=eshape1, .eshape2=eshape2, .elambda=elambda ))))
}


genbetaII= function(link.a="loge",
                    link.scale="loge",
                    link.p="loge",
                    link.q="loge",
                    earg.a=list(), earg.scale=list(), earg.p=list(), earg.q=list(),
                    init.a=NULL,
                    init.scale=NULL,
                    init.p=1.0,
                    init.q=1.0,
                    zero=NULL)
{

    if(mode(link.a) != "character" && mode(link.a) != "name")
        link.a = as.character(substitute(link.a))
    if(mode(link.scale) != "character" && mode(link.scale) != "name")
        link.scale = as.character(substitute(link.scale))
    if(mode(link.p) != "character" && mode(link.p) != "name")
        link.p = as.character(substitute(link.p))
    if(mode(link.q) != "character" && mode(link.q) != "name")
        link.q = as.character(substitute(link.q))
    if(length(zero) && !is.Numeric(zero, integer=TRUE, posit=TRUE))
        stop("bad input for argument \"zero\"")
    if(!is.list(earg.a)) earg.a = list()
    if(!is.list(earg.scale)) earg.scale = list()
    if(!is.list(earg.p)) earg.p = list()
    if(!is.list(earg.q)) earg.q = list()

    new("vglmff",
    blurb=c("Generalized Beta II distribution\n\n",
            "Links:    ",
            namesof("a", link.a, earg=earg.a), ", ", 
            namesof("scale", link.scale, earg=earg.scale), ", ", 
            namesof("p", link.p, earg=earg.p), ", ", 
            namesof("q", link.q, earg=earg.q), "\n", 
            "Mean:     scale*gamma(p + 1/a)*gamma(q - 1/a)/(gamma(p)*gamma(q))"),
    constraints=eval(substitute(expression({
        constraints = cm.zero.vgam(constraints, x, .zero, M)
    }), list( .zero=zero ))),
    initialize=eval(substitute(expression({
        predictors.names = 
        c(namesof("a", .link.a, earg=.earg.a, tag=FALSE),
          namesof("scale", .link.scale, earg=.earg.scale, tag=FALSE),
          namesof("p", .link.p, earg=.earg.p, tag=FALSE),
          namesof("q", .link.q, earg=.earg.q, tag=FALSE))

        if(!length(.init.a) || !length(.init.scale)) {
            qvec = c(.25, .5, .75)   # Arbitrary; could be made an argument
            init.q = if(length(.init.q)) .init.q else 1
            xvec = log( (1-qvec)^(-1/ init.q ) - 1 )
            fit0 = lsfit(x=xvec, y=log(quantile(y, qvec )))
        }

        if(!length(etastart)) {
            aa = rep(if(length(.init.a)) .init.a else 1/fit0$coef[2], length=n)
            scale = rep(if(length(.init.scale)) .init.scale else
                        exp(fit0$coef[1]), length=n)
            qq = rep(if(length(.init.q)) .init.q else 1.0, length=n)
            parg = rep(if(length(.init.p)) .init.p else 1.0, length=n)
            etastart = cbind(theta2eta(aa, .link.a, earg= .earg.a),
                             theta2eta(scale, .link.scale, earg= .earg.scale),
                             theta2eta(parg, .link.p, earg= .earg.p),
                             theta2eta(qq, .link.q, earg= .earg.q))
        }
    }), list( .link.a=link.a, .link.scale=link.scale,
              .link.p=link.p, .link.q=link.q,
              .earg.a=earg.a, .earg.scale=earg.scale, 
              .earg.p=earg.p, .earg.q=earg.q,
              .init.a=init.a, .init.scale=init.scale, 
              .init.p=init.p, .init.q=init.q ))),
    inverse=eval(substitute(function(eta, extra=NULL) {
        aa = eta2theta(eta[,1], .link.a, earg= .earg.a)
        scale = eta2theta(eta[,2], .link.scale, earg= .earg.scale)
        parg = eta2theta(eta[,3], .link.p, earg= .earg.p)
        qq = eta2theta(eta[,4], .link.q, earg= .earg.q)
        scale*gamma(parg + 1/aa)*gamma(qq-1/aa)/(gamma(parg)*gamma(qq))
    }, list( .link.a=link.a, .link.scale=link.scale,
              .earg.a=earg.a, .earg.scale=earg.scale, 
              .earg.p=earg.p, .earg.q=earg.q,
             .link.p=link.p, .link.q=link.q ))),
    last=eval(substitute(expression({
        misc$link = c(a= .link.a, scale= .link.scale,
                      p= .link.p, q= .link.q)
        misc$earg = list(a= .earg.a, scale= .earg.scale,
                      p= .earg.p, q= .earg.q)
    }), list( .link.a=link.a, .link.scale=link.scale,
              .earg.a=earg.a, .earg.scale=earg.scale, 
              .earg.p=earg.p, .earg.q=earg.q,
              .link.p=link.p, .link.q=link.q ))),
    loglikelihood=eval(substitute(
            function(mu, y, w, residuals= FALSE,eta, extra=NULL) {
        aa = eta2theta(eta[,1], .link.a, earg= .earg.a)
        scale = eta2theta(eta[,2], .link.scale, earg= .earg.scale)
        parg = eta2theta(eta[,3], .link.p, earg= .earg.p)
        qq = eta2theta(eta[,4], .link.q, earg= .earg.q)
        if(residuals) stop("loglikelihood residuals not implemented yet") else
        sum(w*(log(aa) + (aa*parg-1)*log(y) - aa*parg*log(scale) +
   (if(is.R()) -lbeta(parg, qq) else lgamma(parg+qq)-lgamma(parg)-lgamma(qq))-
            (parg+qq)*log1p((y/scale)^aa)))
    }, list( .link.a=link.a, .link.scale=link.scale,
              .earg.a=earg.a, .earg.scale=earg.scale, 
              .earg.p=earg.p, .earg.q=earg.q,
            .link.p=link.p, .link.q=link.q ))),
    vfamily=c("genbetaII"),
    deriv=eval(substitute(expression({
        aa = eta2theta(eta[,1], .link.a, earg= .earg.a)
        scale = eta2theta(eta[,2], .link.scale, earg= .earg.scale)
        parg = eta2theta(eta[,3], .link.p, earg= .earg.p)
        qq = eta2theta(eta[,4], .link.q, earg= .earg.q)

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
        da.deta = dtheta.deta(aa, .link.a, earg= .earg.a)
        dscale.deta = dtheta.deta(scale, .link.scale, earg= .earg.scale)
        dp.deta = dtheta.deta(parg, .link.p, earg= .earg.p)
        dq.deta = dtheta.deta(qq, .link.q, earg= .earg.q)
        w * cbind( dl.da * da.deta, dl.dscale * dscale.deta,
                   dl.dp * dp.deta, dl.dq * dq.deta )
    }), list( .link.a=link.a, .link.scale=link.scale,
              .earg.a=earg.a, .earg.scale=earg.scale, 
              .earg.p=earg.p, .earg.q=earg.q,
              .link.p=link.p, .link.q=link.q ))),
    weight=eval(substitute(expression({
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
        wz = w * wz
        wz
    }), list( .link.a=link.a, .link.scale=link.scale,
              .earg.a=earg.a, .earg.scale=earg.scale, 
              .earg.p=earg.p, .earg.q=earg.q,
              .link.p=link.p, .link.q=link.q ))))
}


rsinmad = function(n, a, scale, q.arg)
    qsinmad(runif(n), a, scale, q.arg)

rlomax = function(n, scale, q.arg)
    rsinmad(n, a=1, scale, q.arg)

rfisk = function(n, a, scale)
    rsinmad(n, a, scale, q.arg=1)

rparalogistic = function(n, a, scale)
    rsinmad(n, a, scale, a)

rdagum = function(n, a, scale, p.arg)
    qdagum(runif(n), a, scale, p.arg)

rinvlomax = function(n, scale, p.arg)
    rdagum(n, a=1, scale, p.arg)

rinvparalogistic = function(n, a, scale)
    rdagum(n, a, scale, a)




qsinmad = function(p, a, scale, q.arg) {
    bad = (p < 0) | (p > 1)
    ans = NA * p
    a = rep(a, len=length(p))[!bad]
    scale = rep(scale, len=length(p))[!bad]
    q = rep(q.arg, len=length(p))[!bad]
    xx = p[!bad]
    ans[!bad] = scale* ((1 - xx)^(-1/q) - 1)^(1/a)
    ans
}

qlomax = function(p, scale, q.arg)
    qsinmad(p, a=1, scale, q.arg)

qfisk = function(p, a, scale)
    qsinmad(p, a, scale, q.arg=1)

qparalogistic = function(p, a, scale)
    qsinmad(p, a, scale, a)

qdagum = function(p, a, scale, p.arg) {
    bad = (p < 0) | (p > 1)
    ans = NA * p
    a = rep(a, len=length(p))[!bad]
    scale = rep(scale, len=length(p))[!bad]
    p.arg = rep(p.arg, len=length(p))[!bad]
    xx = p[!bad]
    ans[!bad] = scale* (xx^(-1/p.arg) - 1)^(-1/a)
    ans
}

qinvlomax = function(p, scale, p.arg)
    qdagum(p, a=1, scale, p.arg)

qinvparalogistic = function(p, a, scale)
    qdagum(p, a, scale, a)






psinmad = function(q, a, scale, q.arg) {
    zero = q <= 0
    a = rep(a, len=length(q))[!zero]
    scale = rep(scale, len=length(q))[!zero]
    q.arg = rep(q.arg, len=length(q))[!zero]
    ans = 0 * q
    xx = q[!zero]
    ans[!zero] = 1 - (1 + (xx/scale)^a)^(-q.arg)
    ans
}

plomax = function(q, scale, q.arg)
    psinmad(q, a=1, scale, q.arg)

pfisk = function(q, a, scale)
    psinmad(q, a, scale, q.arg=1)

pparalogistic = function(q, a, scale)
    psinmad(q, a, scale, a)



pdagum = function(q, a, scale, p.arg) {
    zero = q <= 0
    a = rep(a, len=length(q))[!zero]
    scale = rep(scale, len=length(q))[!zero]
    p = rep(p.arg, len=length(q))[!zero]
    ans = 0 * q
    xx = q[!zero]
    ans[!zero] = (1 + (xx/scale)^(-a))^(-p)
    ans
}

pinvlomax = function(q, scale, p.arg)
    pdagum(q, a=1, scale, p.arg)

pinvparalogistic = function(q, a, scale)
    pdagum(q, a, scale, a)



dsinmad = function(x, a, scale, q.arg) {
    zero = x <= 0
    a = rep(a, len=length(x))[!zero]
    scale = rep(scale, len=length(x))[!zero]
    q = rep(q.arg, len=length(x))[!zero]
    ans = 0 * x
    xx = x[!zero]
    ans[!zero] = a * q * xx^(a-1) / (scale^a * (1 + (xx/scale)^a)^(1+q))
    ans
}

dlomax = function(x, scale, q.arg)
    dsinmad(x, a=1, scale, q.arg)

dfisk = function(x, a, scale)
    dsinmad(x, a, scale, q.arg=1)

dparalogistic = function(x, a, scale)
    dsinmad(x, a, scale, a)



ddagum = function(x, a, scale, p.arg) {
    zero = x <= 0
    a = rep(a, len=length(x))[!zero]
    scale = rep(scale, len=length(x))[!zero]
    p = rep(p.arg, len=length(x))[!zero]
    ans = 0 * x
    xx = x[!zero]
    ans[!zero] = a * p * xx^(a*p-1) / (scale^(a*p) * (1 + (xx/scale)^a)^(1+p))
    ans
}

dinvlomax = function(x, scale, p.arg)
    ddagum(x, a=1, scale, p.arg)

dinvparalogistic = function(x, a, scale)
    ddagum(x, a, scale, a)



sinmad = function(link.a="loge",
                  link.scale="loge",
                  link.q="loge",
                  earg.a=list(), earg.scale=list(), earg.q=list(),
                  init.a=NULL, 
                  init.scale=NULL,
                  init.q=1.0, 
                  zero=NULL)
{

    if(mode(link.a) != "character" && mode(link.a) != "name")
        link.a = as.character(substitute(link.a))
    if(mode(link.scale) != "character" && mode(link.scale) != "name")
        link.scale = as.character(substitute(link.scale))
    if(mode(link.q) != "character" && mode(link.q) != "name")
        link.q = as.character(substitute(link.q))
    if(length(zero) && !is.Numeric(zero, integer=TRUE, posit=TRUE))
        stop("bad input for argument \"zero\"")
    if(!is.list(earg.a)) earg.a = list()
    if(!is.list(earg.scale)) earg.scale = list()
    if(!is.list(earg.q)) earg.q = list()

    new("vglmff",
    blurb=c("Singh-Maddala distribution\n\n",
            "Links:    ",
            namesof("a", link.a, earg=earg.a), ", ", 
            namesof("scale", link.scale, earg=earg.scale), ", ", 
            namesof("q", link.q, earg=earg.q), "\n", 
            "Mean:     scale*gamma(1 + 1/a)*gamma(q - 1/a)/gamma(q)"),
    constraints=eval(substitute(expression({
        constraints = cm.zero.vgam(constraints, x, .zero, M)
    }), list( .zero=zero ))),
    initialize=eval(substitute(expression({
        if(ncol(cbind(y)) != 1)
            stop("response must be a vector or a one-column matrix")
        predictors.names = 
            c(namesof("a", .link.a, earg=.earg.a, tag=FALSE),
              namesof("scale", .link.scale, earg=.earg.scale, tag=FALSE),
              namesof("q", .link.q, earg=.earg.q, tag=FALSE))
        parg = 1

        if(!length(.init.a) || !length(.init.scale)) {
            qvec = c(.25, .5, .75)   # Arbitrary; could be made an argument
            init.q = if(length(.init.q)) .init.q else 1
            xvec = log( (1-qvec)^(-1/ init.q ) - 1 )
            fit0 = lsfit(x=xvec, y=log(quantile(y, qvec )))
        }

        if(!length(etastart)) {
            aa = rep(if(length(.init.a)) .init.a else 1/fit0$coef[2], length=n)
            scale = rep(if(length(.init.scale)) .init.scale else
                        exp(fit0$coef[1]), length=n)
            qq = rep(if(length(.init.q)) .init.q else 1.0, length=n)
            etastart = cbind(theta2eta(aa, .link.a, earg= .earg.a),
                             theta2eta(scale, .link.scale, earg= .earg.scale),
                             theta2eta(qq, .link.q, earg= .earg.q))
        }
    }), list( .link.a=link.a, .link.scale=link.scale,
              .link.q=link.q,
              .earg.a=earg.a, .earg.scale=earg.scale, 
              .earg.q=earg.q,
              .init.a=init.a, .init.scale=init.scale, 
              .init.q=init.q ))),
    inverse=eval(substitute(function(eta, extra=NULL) {
        aa = eta2theta(eta[,1], .link.a, earg= .earg.a)
        scale = eta2theta(eta[,2], .link.scale, earg= .earg.scale)
        qq = eta2theta(eta[,3], .link.q, earg= .earg.q)
        scale*gamma(1 + 1/aa)*gamma(qq-1/aa)/(gamma(qq))
    }, list( .link.a=link.a, .link.scale=link.scale,
              .earg.a=earg.a, .earg.scale=earg.scale, 
              .earg.q=earg.q,
             .link.q=link.q ))),
    last=eval(substitute(expression({
        misc$link = c(a= .link.a, scale= .link.scale, q= .link.q)
        misc$earg = list(a= .earg.a, scale= .earg.scale, q= .earg.q)
    }), list( .link.a=link.a, .link.scale=link.scale,
              .earg.a=earg.a, .earg.scale=earg.scale, 
              .earg.q=earg.q,
              .link.q=link.q ))),
    loglikelihood=eval(substitute(
            function(mu, y, w, residuals= FALSE,eta, extra=NULL) {
        aa = eta2theta(eta[,1], .link.a, earg= .earg.a)
        scale = eta2theta(eta[,2], .link.scale, earg= .earg.scale)
        parg = 1
        qq = eta2theta(eta[,3], .link.q, earg= .earg)
        if(residuals) stop("loglikelihood residuals not implemented yet") else
        sum(w*(log(aa) + (aa*parg-1)*log(y) - aa*parg*log(scale) +
   (if(is.R()) -lbeta(parg, qq) else lgamma(parg+qq)-lgamma(parg)-lgamma(qq))-
            (parg+qq)*log1p((y/scale)^aa)))
    }, list( .link.a=link.a, .link.scale=link.scale,
              .earg.a=earg.a, .earg.scale=earg.scale, 
              .earg.q=earg.q,
            .link.q=link.q ))),
    vfamily=c("sinmad"),
    deriv=eval(substitute(expression({
        aa = eta2theta(eta[,1], .link.a, earg= .earg.a)
        scale = eta2theta(eta[,2], .link.scale, earg= .earg.scale)
        parg = 1
        qq = eta2theta(eta[,3], .link.q, earg= .earg.q)

        temp1 = log(y/scale)
        temp2 = (y/scale)^aa
        temp3a = digamma(parg)
        temp3b = digamma(qq)

        dl.da = 1/aa + parg * temp1 - (parg+qq) * temp1 / (1+1/temp2)
        dl.dscale = (aa/scale) * (-parg + (parg+qq) / (1+1/temp2))
        dl.dq = digamma(parg + qq) - temp3b - log1p(temp2)
        da.deta = dtheta.deta(aa, .link.a, earg= .earg.a)
        dscale.deta = dtheta.deta(scale, .link.scale, earg= .earg.scale)
        dq.deta = dtheta.deta(qq, .link.q, earg= .earg.q)
        w * cbind( dl.da * da.deta, dl.dscale * dscale.deta,
                   dl.dq * dq.deta )
    }), list( .link.a=link.a, .link.scale=link.scale,
              .earg.a=earg.a, .earg.scale=earg.scale, 
              .earg.q=earg.q,
              .link.q=link.q ))),
    weight=eval(substitute(expression({
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
        wz = w * wz
        wz
    }), list( .link.a=link.a, .link.scale=link.scale,
              .earg.a=earg.a, .earg.scale=earg.scale, 
              .earg.q=earg.q,
              .link.q=link.q ))))
}


 dagum = function(link.a="loge",
                  link.scale="loge",
                  link.p="loge",
                  earg.a=list(), earg.scale=list(), earg.p=list(),
                  init.a=NULL, 
                  init.scale=NULL,
                  init.p=1.0, 
                  zero=NULL)
{

    if(mode(link.a) != "character" && mode(link.a) != "name")
        link.a = as.character(substitute(link.a))
    if(mode(link.scale) != "character" && mode(link.scale) != "name")
        link.scale = as.character(substitute(link.scale))
    if(mode(link.p) != "character" && mode(link.p) != "name")
        link.p = as.character(substitute(link.p))
    if(length(zero) && !is.Numeric(zero, integer=TRUE, posit=TRUE))
        stop("bad input for argument \"zero\"")
    if(!is.list(earg.a)) earg.a = list()
    if(!is.list(earg.scale)) earg.scale = list()
    if(!is.list(earg.p)) earg.p = list()

    new("vglmff",
    blurb=c("Dagum distribution\n\n",
            "Links:    ",
            namesof("a", link.a, earg=earg.a), ", ", 
            namesof("scale", link.scale, earg=earg.scale), ", ", 
            namesof("p", link.p, earg=earg.p), "\n", 
            "Mean:     scale*gamma(p + 1/a)*gamma(1 - 1/a)/gamma(p)"),
    constraints=eval(substitute(expression({
        constraints = cm.zero.vgam(constraints, x, .zero, M)
    }), list( .zero=zero ))),
    initialize=eval(substitute(expression({
        if(ncol(cbind(y)) != 1)
            stop("response must be a vector or a one-column matrix")
        predictors.names =
        c(namesof("a", .link.a, earg=.earg.a, tag=FALSE),
          namesof("scale", .link.scale, earg=.earg.scale, tag=FALSE),
          namesof("p", .link.p, earg=.earg.p, tag=FALSE))

        if(!length(.init.a) || !length(.init.scale)) {
            qvec = c(.25, .5, .75)   # Arbitrary; could be made an argument
            init.p = if(length(.init.p)) .init.p else 1
            xvec = log( qvec^(-1/ init.p ) - 1 )
            fit0 = lsfit(x=xvec, y=log(quantile(y, qvec )))
        }

        if(!length(etastart)) {
            parg = rep(if(length(.init.p)) .init.p else 1.0, length=n)
            aa = rep(if(length(.init.a)) .init.a else -1/fit0$coef[2], length=n)
            scale = rep(if(length(.init.scale)) .init.scale else
                        exp(fit0$coef[1]), length=n)
            etastart = cbind(theta2eta(aa, .link.a, earg= .earg.a),
                             theta2eta(scale, .link.scale, earg= .earg.scale),
                             theta2eta(parg, .link.p, earg= .earg.p))
        }
    }), list( .link.a=link.a, .link.scale=link.scale,
              .link.p=link.p,
              .earg.a=earg.a, .earg.scale=earg.scale, 
              .earg.p=earg.p,
              .init.a=init.a, .init.scale=init.scale, 
              .init.p=init.p ))),
    inverse=eval(substitute(function(eta, extra=NULL) {
        aa = eta2theta(eta[,1], .link.a, earg= .earg.a)
        scale = eta2theta(eta[,2], .link.scale, earg= .earg.scale)
        parg = eta2theta(eta[,3], .link.p, earg= .earg.p)
        qq = 1
        scale*gamma(parg + 1/aa)*gamma(qq-1/aa)/(gamma(parg)*gamma(qq))
    }, list( .link.a=link.a, .link.scale=link.scale,
              .earg.a=earg.a, .earg.scale=earg.scale, 
              .earg.p=earg.p,
             .link.p=link.p ))),
    last=eval(substitute(expression({
        misc$link = c(a= .link.a, scale= .link.scale, p= .link.p )
        misc$earg = list(a= .earg.a, scale= .earg.scale, p= .earg.p)
    }), list( .link.a=link.a, .link.scale=link.scale,
              .earg.a=earg.a, .earg.scale=earg.scale, 
              .earg.p=earg.p,
              .link.p=link.p ))),
    loglikelihood=eval(substitute(
            function(mu, y, w, residuals= FALSE,eta, extra=NULL) {
        aa = eta2theta(eta[,1], .link.a, earg= .earg.a)
        scale = eta2theta(eta[,2], .link.scale, earg= .earg.scale)
        parg = eta2theta(eta[,3], .link.p, earg= .earg.p)
        qq = 1
        if(residuals) stop("loglikelihood residuals not implemented yet") else
        sum(w*(log(aa) + (aa*parg-1)*log(y) - aa*parg*log(scale) +
   (if(is.R()) -lbeta(parg, qq) else lgamma(parg+qq)-lgamma(parg)-lgamma(qq))-
            (parg+qq)*log1p((y/scale)^aa)))
    }, list( .link.a=link.a, .link.scale=link.scale,
              .earg.a=earg.a, .earg.scale=earg.scale, 
              .earg.p=earg.p,
            .link.p=link.p ))),
    vfamily=c("dagum"),
    deriv=eval(substitute(expression({
        aa = eta2theta(eta[,1], .link.a, earg= .earg.a)
        scale = eta2theta(eta[,2], .link.scale, earg= .earg.scale)
        parg = eta2theta(eta[,3], .link.p, earg= .earg.p)
        qq = 1

        temp1 = log(y/scale)
        temp2 = (y/scale)^aa
        temp3a = digamma(parg)
        temp3b = digamma(qq)

        dl.da = 1/aa + parg * temp1 - (parg+qq) * temp1 / (1+1/temp2)
        dl.dscale = (aa/scale) * (-parg + (parg+qq) / (1+1/temp2))
        dl.dp = aa * temp1 + digamma(parg + qq) - temp3a - log1p(temp2)
        da.deta = dtheta.deta(aa, .link.a, earg= .earg.a)
        dscale.deta = dtheta.deta(scale, .link.scale, earg= .earg.scale)
        dp.deta = dtheta.deta(parg, .link.p, earg= .earg.p)
        w * cbind( dl.da * da.deta, dl.dscale * dscale.deta,
                   dl.dp * dp.deta )
    }), list( .link.a=link.a, .link.scale=link.scale,
              .earg.a=earg.a, .earg.scale=earg.scale, 
              .earg.p=earg.p,
              .link.p=link.p ))),
    weight=eval(substitute(expression({
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
        wz = w * wz
        wz
    }), list( .link.a=link.a, .link.scale=link.scale,
              .earg.a=earg.a, .earg.scale=earg.scale, 
              .earg.p=earg.p,
              .link.p=link.p ))))
}



betaII= function(link.scale="loge",
                 link.p="loge",
                 link.q="loge",
                  earg.scale=list(), earg.p=list(), earg.q=list(),
                 init.scale=NULL,
                 init.p=1.0, 
                 init.q=1.0, 
                 zero=NULL)
{

    if(mode(link.scale) != "character" && mode(link.scale) != "name")
        link.scale = as.character(substitute(link.scale))
    if(mode(link.p) != "character" && mode(link.p) != "name")
        link.p = as.character(substitute(link.p))
    if(mode(link.q) != "character" && mode(link.q) != "name")
        link.q = as.character(substitute(link.q))
    if(length(zero) && !is.Numeric(zero, integer=TRUE, posit=TRUE))
        stop("bad input for argument \"zero\"")
    if(!is.list(earg.scale)) earg.scale = list()
    if(!is.list(earg.p)) earg.p = list()
    if(!is.list(earg.q)) earg.q = list()

    new("vglmff",
    blurb=c("Beta II distribution\n\n",
            "Links:    ",
            namesof("scale", link.scale, earg=earg.scale), ", ", 
            namesof("p", link.p, earg=earg.p), ", ", 
            namesof("q", link.q, earg=earg.q), "\n", 
            "Mean:     scale*gamma(p + 1)*gamma(q - 1)/(gamma(p)*gamma(q))"),
    constraints=eval(substitute(expression({
        constraints = cm.zero.vgam(constraints, x, .zero, M)
    }), list( .zero=zero ))),
    initialize=eval(substitute(expression({
        if(ncol(cbind(y)) != 1)
            stop("response must be a vector or a one-column matrix")
        predictors.names = 
        c(namesof("scale", .link.scale, earg=.earg.scale, tag=FALSE),
          namesof("p", .link.p, earg=.earg.p, tag=FALSE),
          namesof("q", .link.q, earg=.earg.q, tag=FALSE))

        if(!length(.init.scale)) {
            qvec = c(.25, .5, .75)   # Arbitrary; could be made an argument
            init.q = if(length(.init.q)) .init.q else 1
            xvec = log( (1-qvec)^(-1/ init.q ) - 1 )
            fit0 = lsfit(x=xvec, y=log(quantile(y, qvec )))
        }

        if(!length(etastart)) {
            scale = rep(if(length(.init.scale)) .init.scale else
                        exp(fit0$coef[1]), length=n)
            qq = rep(if(length(.init.q)) .init.q else 1.0, length=n)
            parg = rep(if(length(.init.p)) .init.p else 1.0, length=n)
            etastart = cbind(theta2eta(scale, .link.scale, earg= .earg.scale),
                             theta2eta(parg, .link.p, earg= .earg.p),
                             theta2eta(qq, .link.q, earg= .earg.q))
        }
    }), list( .link.scale=link.scale,
              .link.p=link.p, .link.q=link.q,
              .earg.scale=earg.scale, 
              .earg.p=earg.p, .earg.q=earg.q,
              .init.scale=init.scale, 
              .init.p=init.p, .init.q=init.q ))),
    inverse=eval(substitute(function(eta, extra=NULL) {
        aa = 1
        scale = eta2theta(eta[,1], .link.scale, earg= .earg.scale)
        parg = eta2theta(eta[,2], .link.p, earg= .earg.p)
        qq = eta2theta(eta[,3], .link.q, earg= .earg.q)
        scale*gamma(parg + 1/aa)*gamma(qq-1/aa)/(gamma(parg)*gamma(qq))
    }, list( .link.scale=link.scale,
              .earg.scale=earg.scale, 
              .earg.p=earg.p, .earg.q=earg.q,
             .link.p=link.p, .link.q=link.q ))),
    last=eval(substitute(expression({
        misc$link = c(scale= .link.scale, p= .link.p, q= .link.q)
        misc$earg = list(scale= .earg.scale, p= .earg.p, q= .earg.q)
    }), list( .link.scale=link.scale,
              .earg.scale=earg.scale, 
              .earg.p=earg.p, .earg.q=earg.q,
              .link.p=link.p, .link.q=link.q ))),
    loglikelihood=eval(substitute(
            function(mu, y, w, residuals= FALSE,eta, extra=NULL) {
        aa = 1
        scale = eta2theta(eta[,1], .link.scale, earg= .earg.scale)
        parg = eta2theta(eta[,2], .link.p, earg= .earg.p)
        qq = eta2theta(eta[,3], .link.q, earg= .earg.q)
        if(residuals) stop("loglikelihood residuals not implemented yet") else
        sum(w*(log(aa) + (aa*parg-1)*log(y) - aa*parg*log(scale) +
   (if(is.R()) -lbeta(parg, qq) else lgamma(parg+qq)-lgamma(parg)-lgamma(qq))-
            (parg+qq)*log1p((y/scale)^aa)))
    }, list( .link.scale=link.scale,
              .earg.scale=earg.scale, 
              .earg.p=earg.p, .earg.q=earg.q,
             .link.p=link.p, .link.q=link.q ))),
    vfamily=c("betaII"),
    deriv=eval(substitute(expression({
        aa = 1
        scale = eta2theta(eta[,1], .link.scale, earg= .earg.scale)
        parg = eta2theta(eta[,2], .link.p, earg= .earg.p)
        qq = eta2theta(eta[,3], .link.q, earg= .earg.q)

        temp1 = log(y/scale)
        temp2 = (y/scale)^aa
        temp3 = digamma(parg + qq)
        temp3a = digamma(parg)
        temp3b = digamma(qq)
        temp4 = log1p(temp2)

        dl.dscale = (aa/scale) * (-parg + (parg+qq) / (1+1/temp2))
        dl.dp = aa * temp1 + temp3 - temp3a - temp4
        dl.dq = temp3 - temp3b - temp4
        dscale.deta = dtheta.deta(scale, .link.scale, earg= .earg.scale)
        dp.deta = dtheta.deta(parg, .link.p, earg= .earg.p)
        dq.deta = dtheta.deta(qq, .link.q, earg= .earg.q)
        w * cbind( dl.dscale * dscale.deta,
                   dl.dp * dp.deta, dl.dq * dq.deta )
    }), list( .link.scale=link.scale,
              .earg.scale=earg.scale, 
              .earg.p=earg.p, .earg.q=earg.q,
              .link.p=link.p, .link.q=link.q ))),
    weight=eval(substitute(expression({
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
        wz = w * wz
        wz
    }), list( .link.scale=link.scale,
              .earg.scale=earg.scale, 
              .earg.p=earg.p, .earg.q=earg.q,
              .link.p=link.p, .link.q=link.q ))))
}



lomax = function(link.scale="loge",
                 link.q="loge",
                 earg.scale=list(), earg.q=list(),
                 init.scale=NULL,
                 init.q=1.0, 
                 zero=NULL)
{

    if(mode(link.scale) != "character" && mode(link.scale) != "name")
        link.scale = as.character(substitute(link.scale))
    if(mode(link.q) != "character" && mode(link.q) != "name")
        link.q = as.character(substitute(link.q))
    if(length(zero) && !is.Numeric(zero, integer=TRUE, posit=TRUE))
        stop("bad input for argument \"zero\"")
    if(!is.list(earg.scale)) earg.scale = list()
    if(!is.list(earg.q)) earg.q = list()

    new("vglmff",
    blurb=c("Lomax distribution\n\n",
            "Links:    ",
            namesof("scale", link.scale, earg=earg.scale), ", ", 
            namesof("q", link.q, earg=earg.q), "\n", 
            "Mean:     scale/(q-1)"),
    constraints=eval(substitute(expression({
        constraints = cm.zero.vgam(constraints, x, .zero, M)
    }), list( .zero=zero ))),
    initialize=eval(substitute(expression({
        if(ncol(cbind(y)) != 1)
            stop("response must be a vector or a one-column matrix")
        predictors.names =
        c(namesof("scale", .link.scale, earg=.earg.scale, tag=FALSE),
          namesof("q", .link.q, earg=.earg.q, tag=FALSE))
        aa = parg = 1

        if(!length(.init.scale)) {
            qvec = c(.25, .5, .75)   # Arbitrary; could be made an argument
            init.q = if(length(.init.q)) .init.q else 1
            xvec = log( (1-qvec)^(-1/ init.q ) - 1 )
            fit0 = lsfit(x=xvec, y=log(quantile(y, qvec )))
        }

        if(!length(etastart)) {
            qq = rep(if(length(.init.q)) .init.q else 1.0, length=n)
            scale = rep(if(length(.init.scale)) .init.scale else
                        exp(fit0$coef[1]), length=n)
            etastart = cbind(theta2eta(scale, .link.scale, earg= .earg.scale),
                             theta2eta(qq, .link.q, earg= .earg.q))
        }
    }), list( .link.scale=link.scale,
              .link.q=link.q,
              .earg.scale=earg.scale, 
              .earg.q=earg.q,
              .init.scale=init.scale, 
              .init.q=init.q ))),
    inverse=eval(substitute(function(eta, extra=NULL) {
        scale = eta2theta(eta[,1], .link.scale, earg= .earg.scale)
        qq = eta2theta(eta[,2], .link.q, earg= .earg.q)
        scale/(qq-1)
    }, list( .link.scale=link.scale,
              .earg.scale=earg.scale, 
              .earg.q=earg.q,
             .link.q=link.q ))),
    last=eval(substitute(expression({
        misc$link = c(scale= .link.scale, q= .link.q)
        misc$earg = list(scale= .earg.scale, q= .earg.q)
    }), list( .link.scale=link.scale,
              .earg.scale=earg.scale, 
              .earg.q=earg.q,
              .link.q=link.q ))),
    loglikelihood=eval(substitute(
            function(mu, y, w, residuals= FALSE,eta, extra=NULL) {
        aa = 1
        scale = eta2theta(eta[,1], .link.scale, earg= .earg.scale)
        parg = 1
        qq = eta2theta(eta[,2], .link.q, earg= .earg.q)
        if(residuals) stop("loglikelihood residuals not implemented yet") else
        sum(w*(log(aa) + (aa*parg-1)*log(y) - aa*parg*log(scale) +
   (if(is.R()) -lbeta(parg, qq) else lgamma(parg+qq)-lgamma(parg)-lgamma(qq))-
            (parg+qq)*log1p((y/scale)^aa)))
    }, list( .link.scale=link.scale,
              .earg.scale=earg.scale, 
              .earg.q=earg.q,
            .link.q=link.q ))),
    vfamily=c("lomax"),
    deriv=eval(substitute(expression({
        aa = 1
        scale = eta2theta(eta[,1], .link.scale, earg= .earg.scale)
        parg = 1
        qq = eta2theta(eta[,2], .link.q, earg= .earg.q)
        temp2 = (y/scale)^aa

        dl.dscale = (aa/scale) * (-parg + (parg+qq) / (1+1/temp2))
        dl.dq = digamma(parg + qq) - digamma(qq) - log1p(temp2)
        dscale.deta = dtheta.deta(scale, .link.scale, earg= .earg.scale)
        dq.deta = dtheta.deta(qq, .link.q, earg= .earg.q)
        w * cbind( dl.dscale * dscale.deta,
                   dl.dq * dq.deta )
    }), list( .link.scale=link.scale,
              .earg.scale=earg.scale, 
              .earg.q=earg.q,
              .link.q=link.q ))),
    weight=eval(substitute(expression({
        ed2l.dscale = aa^2 * parg * qq / (scale^2 * (1+parg+qq))
        ed2l.dq = 1/qq^2 
        ed2l.dscaleq = -aa * parg / (scale*(parg+qq))
        wz = matrix(as.numeric(NA), n, dimm(M))  #M==2 means 3=dimm(M)
        wz[,iam(1,1,M)] = ed2l.dscale * dscale.deta^2
        wz[,iam(2,2,M)] = ed2l.dq * dq.deta^2
        wz[,iam(1,2,M)] = ed2l.dscaleq * dscale.deta * dq.deta
        wz = w * wz
        wz
    }), list( .link.scale=link.scale,
              .earg.scale=earg.scale, 
              .earg.q=earg.q,
              .link.q=link.q ))))
}


 fisk = function(link.a="loge",
                 link.scale="loge",
                 earg.a=list(), earg.scale=list(),
                 init.a=NULL, 
                 init.scale=NULL,
                 zero=NULL)
{

    if(mode(link.a) != "character" && mode(link.a) != "name")
        link.a = as.character(substitute(link.a))
    if(mode(link.scale) != "character" && mode(link.scale) != "name")
        link.scale = as.character(substitute(link.scale))
    if(length(zero) && !is.Numeric(zero, integer=TRUE, posit=TRUE))
        stop("bad input for argument \"zero\"")
    if(!is.list(earg.a)) earg.a = list()
    if(!is.list(earg.scale)) earg.scale = list()

    new("vglmff",
    blurb=c("Fisk distribution\n\n",
            "Links:    ",
            namesof("a", link.a, earg=earg.a), ", ", 
            namesof("scale", link.scale, earg=earg.scale), "\n", 
                "Mean:     scale * gamma(1 + 1/a) * gamma(1 - 1/a)"),
    constraints=eval(substitute(expression({
        constraints = cm.zero.vgam(constraints, x, .zero, M)
    }), list( .zero=zero ))),
    initialize=eval(substitute(expression({
        predictors.names =
        c(namesof("a", .link.a, earg=.earg.a, tag=FALSE),
          namesof("scale", .link.scale, earg=.earg.scale, tag=FALSE))
        qq = parg = 1

        if(!length(.init.scale)) {
            qvec = c(.25, .5, .75)   # Arbitrary; could be made an argument
            xvec = log( 1/qvec - 1 )
            fit0 = lsfit(x=xvec, y=log(quantile(y, qvec )))
        }

        if(!length(etastart)) {
            aa = rep(if(length(.init.a)) .init.a else -1/fit0$coef[2], length=n)
            scale = rep(if(length(.init.scale)) .init.scale else
                        exp(fit0$coef[1]), length=n)
            etastart = cbind(theta2eta(aa, .link.a, earg= .earg.a),
                             theta2eta(scale, .link.scale, earg= .earg.scale))
        }
    }), list( .link.a=link.a, .link.scale=link.scale,
              .earg.a=earg.a, .earg.scale=earg.scale, 
              .init.a=init.a, .init.scale=init.scale
            ))),
    inverse=eval(substitute(function(eta, extra=NULL) {
        aa = eta2theta(eta[,1], .link.a, earg= .earg.a)
        scale = eta2theta(eta[,2], .link.scale, earg= .earg.scale)
        qq = 1
        scale*gamma(1 + 1/aa)*gamma(1-1/aa)
    }, list( .link.a=link.a, .link.scale=link.scale,
              .earg.a=earg.a, .earg.scale=earg.scale
           ))),
    last=eval(substitute(expression({
        misc$link = c(a= .link.a, scale= .link.scale)
        misc$earg = list(a= .earg.a, scale= .earg.scale)
    }), list( .link.a=link.a, .link.scale=link.scale,
              .earg.a=earg.a, .earg.scale=earg.scale
            ))),
    loglikelihood=eval(substitute(
            function(mu, y, w, residuals= FALSE,eta, extra=NULL) {
        aa = eta2theta(eta[,1], .link.a, earg= .earg)
        scale = eta2theta(eta[,2], .link.scale, earg= .earg.scale)
        parg = qq = 1
        if(residuals) stop("loglikelihood residuals not implemented yet") else
        sum(w*(log(aa) + (aa*parg-1)*log(y) - aa*parg*log(scale) +
   (if(is.R()) -lbeta(parg, qq) else lgamma(parg+qq)-lgamma(parg)-lgamma(qq))-
            (parg+qq)*log1p((y/scale)^aa)))
    }, list( .link.a=link.a, .link.scale=link.scale,
              .earg.a=earg.a, .earg.scale=earg.scale ))),
    vfamily=c("fisk"),
    deriv=eval(substitute(expression({
        aa = eta2theta(eta[,1], .link.a, earg= .earg.a)
        scale = eta2theta(eta[,2], .link.scale, earg= .earg.scale)
        parg = qq = 1

        temp1 = log(y/scale)
        temp2 = (y/scale)^aa
        temp3a = digamma(parg)
        temp3b = digamma(qq)

        dl.da = 1/aa + parg * temp1 - (parg+qq) * temp1 / (1+1/temp2)
        dl.dscale = (aa/scale) * (-parg + (parg+qq) / (1+1/temp2))
        da.deta = dtheta.deta(aa, .link.a, earg= .earg.a)
        dscale.deta = dtheta.deta(scale, .link.scale, earg= .earg.scale)
        w * cbind( dl.da * da.deta, dl.dscale * dscale.deta )
    }), list( .link.a=link.a, .link.scale=link.scale,
              .earg.a=earg.a, .earg.scale=earg.scale
              ))),
    weight=eval(substitute(expression({
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
        wz = w * wz
        wz
    }), list( .link.a=link.a, .link.scale=link.scale,
              .earg.a=earg.a, .earg.scale=earg.scale ))))
}


invlomax = function(link.scale="loge",
                    link.p="loge",
                    earg.scale=list(), earg.p=list(),
                    init.scale=NULL,
                    init.p=1.0, 
                    zero=NULL)
{

    if(mode(link.scale) != "character" && mode(link.scale) != "name")
        link.scale = as.character(substitute(link.scale))
    if(mode(link.p) != "character" && mode(link.p) != "name")
        link.p = as.character(substitute(link.p))
    if(length(zero) && !is.Numeric(zero, integer=TRUE, posit=TRUE))
        stop("bad input for argument \"zero\"")
    if(!is.list(earg.scale)) earg.scale = list()
    if(!is.list(earg.p)) earg.p = list()

    new("vglmff",
    blurb=c("Inverse Lomax distribution\n\n",
            "Links:    ",
            namesof("scale", link.scale, earg=earg.scale), ", ", 
            namesof("p", link.p, earg=earg.p), "\n", 
            "Mean:     does not exist"),
    constraints=eval(substitute(expression({
        constraints = cm.zero.vgam(constraints, x, .zero, M)
    }), list( .zero=zero ))),
    initialize=eval(substitute(expression({
        if(ncol(cbind(y)) != 1)
            stop("response must be a vector or a one-column matrix")
        predictors.names =
        c(namesof("scale", .link.scale, earg=.earg.scale, tag=FALSE),
          namesof("p", .link.p, earg=.earg.p, tag=FALSE))
        qq = aa = 1

        if(!length(.init.scale)) {
            qvec = c(.25, .5, .75)   # Arbitrary; could be made an argument
            init.p = if(length(.init.p)) .init.p else 1
            xvec = log( qvec^(-1/ init.p ) - 1 )
            fit0 = lsfit(x=xvec, y=log(quantile(y, qvec )))
        }
        if(!length(etastart)) {
            scale = rep(if(length(.init.scale)) .init.scale else
                        exp(fit0$coef[1]), length=n)
            parg = rep(if(length(.init.p)) .init.p else 1.0, length=n)
            etastart = cbind(theta2eta(scale, .link.scale, earg= .earg.scale),
                             theta2eta(parg, .link.p, earg= .earg.p))
        }
    }), list( .link.scale=link.scale,
              .link.p=link.p,
              .earg.scale=earg.scale, 
              .earg.p=earg.p,
              .init.scale=init.scale, 
              .init.p=init.p ))),
    inverse=eval(substitute(function(eta, extra=NULL) {
        rep(as.numeric(NA), len=nrow(eta))
    }, list( .link.scale=link.scale,
              .earg.scale=earg.scale, 
              .earg.p=earg.p,
             .link.p=link.p ))),
    last=eval(substitute(expression({
        misc$link = c(scale= .link.scale, p= .link.p )
        misc$earg = list(scale= .earg.scale, p= .earg.p )
    }), list( .link.scale=link.scale,
              .earg.scale=earg.scale, 
              .earg.p=earg.p,
              .link.p=link.p ))),
    loglikelihood=eval(substitute(
            function(mu, y, w, residuals= FALSE,eta, extra=NULL) {
        aa = qq = 1
        scale = eta2theta(eta[,1], .link.scale, earg= .earg.scale)
        parg = eta2theta(eta[,2], .link.p, earg= .earg.p)
        if(residuals) stop("loglikelihood residuals not implemented yet") else
        sum(w*(log(aa) + (aa*parg-1)*log(y) - aa*parg*log(scale) +
   (if(is.R()) -lbeta(parg, qq) else lgamma(parg+qq)-lgamma(parg)-lgamma(qq))-
            (parg+qq)*log1p((y/scale)^aa)))
    }, list( .link.scale=link.scale,
              .earg.scale=earg.scale, 
              .earg.p=earg.p,
            .link.p=link.p ))),
    vfamily=c("invlomax"),
    deriv=eval(substitute(expression({
        aa = qq = 1 
        scale = eta2theta(eta[,1], .link.scale, earg= .earg.scale)
        parg = eta2theta(eta[,2], .link.p, earg= .earg.p)

        temp1 = log(y/scale)
        temp2 = (y/scale)^aa

        dl.dscale = (aa/scale) * (-parg + (parg+qq) / (1+1/temp2))
        dl.dp = aa * temp1 + digamma(parg + qq) - digamma(parg) - log1p(temp2)
        dscale.deta = dtheta.deta(scale, .link.scale, earg= .earg.scale)
        dp.deta = dtheta.deta(parg, .link.p, earg= .earg.p)
        w * cbind( dl.dscale * dscale.deta,
                   dl.dp * dp.deta )
    }), list( .link.scale=link.scale,
              .earg.scale=earg.scale, 
              .earg.p=earg.p,
              .link.p=link.p ))),
    weight=eval(substitute(expression({
        ed2l.dscale = aa^2 * parg * qq / (scale^2 * (1+parg+qq))
        ed2l.dp = 1/parg^2 
        ed2l.dscalep =  aa * qq   / (scale*(parg+qq))
        wz = matrix(as.numeric(NA), n, dimm(M))  #M==2 means 3=dimm(M)
        wz[,iam(1,1,M)] = ed2l.dscale * dscale.deta^2
        wz[,iam(2,2,M)] = ed2l.dp * dp.deta^2
        wz[,iam(1,2,M)] = ed2l.dscalep * dscale.deta * dp.deta
        wz = w * wz
        wz
    }), list( .link.scale=link.scale,
              .earg.scale=earg.scale, 
              .earg.p=earg.p,
              .link.p=link.p ))))
}


paralogistic = function(link.a="loge",
                  link.scale="loge",
                    earg.a=list(), earg.scale=list(), 
                  init.a=1.0,
                  init.scale=NULL,
                  zero=NULL)
{

    if(mode(link.a) != "character" && mode(link.a) != "name")
        link.a = as.character(substitute(link.a))
    if(mode(link.scale) != "character" && mode(link.scale) != "name")
        link.scale = as.character(substitute(link.scale))
    if(length(zero) && !is.Numeric(zero, integer=TRUE, posit=TRUE))
        stop("bad input for argument \"zero\"")
    if(!is.list(earg.a)) earg.a = list()
    if(!is.list(earg.scale)) earg.scale = list()

    new("vglmff",
    blurb=c("Paralogistic distribution\n\n",
            "Links:    ",
            namesof("a", link.a, earg=earg.a), ", ", 
            namesof("scale", link.scale, earg=earg.scale), "\n", 
            "Mean:     scale*gamma(1 + 1/a)*gamma(a - 1/a)/gamma(a)"),
    constraints=eval(substitute(expression({
        constraints = cm.zero.vgam(constraints, x, .zero, M)
    }), list( .zero=zero ))),
    initialize=eval(substitute(expression({
        if(ncol(cbind(y)) != 1)
            stop("response must be a vector or a one-column matrix")
        predictors.names =
        c(namesof("a", .link.a, earg=.earg.a, tag=FALSE),
          namesof("scale", .link.scale, earg=.earg.scale, tag=FALSE))
        parg = 1

        if(!length(.init.a) || !length(.init.scale)) {
            qvec = c(.25, .5, .75)   # Arbitrary; could be made an argument
            init.a = if(length(.init.a)) .init.a else 1
            xvec = log( (1-qvec)^(-1/ init.a ) - 1 )
            fit0 = lsfit(x=xvec, y=log(quantile(y, qvec )))
        }

        if(!length(etastart)) {
            aa = rep(if(length(.init.a)) .init.a else 1/fit0$coef[2], length=n)
            scale = rep(if(length(.init.scale)) .init.scale else
                    exp(fit0$coef[1]), length=n)
            etastart = cbind(theta2eta(aa, .link.a, earg= .earg.a),
                             theta2eta(scale, .link.scale, earg= .earg.scale))
        }
    }), list( .link.a=link.a, .link.scale=link.scale,
              .earg.a=earg.a, .earg.scale=earg.scale, 
              .init.a=init.a, .init.scale=init.scale
              ))),
    inverse=eval(substitute(function(eta, extra=NULL) {
        aa = eta2theta(eta[,1], .link.a, earg= .earg.a)
        scale = eta2theta(eta[,2], .link.scale, earg= .earg.scale)
        qq = aa
        scale*gamma(1 + 1/aa)*gamma(qq-1/aa)/(gamma(qq))
    }, list( .link.a=link.a, .link.scale=link.scale,
              .earg.a=earg.a, .earg.scale=earg.scale
             ))),
    last=eval(substitute(expression({
        misc$link = c(a= .link.a, scale= .link.scale)
        misc$earg = list(a= .earg.a, scale= .earg.scale )
    }), list( .link.a=link.a, .link.scale=link.scale,
              .earg.a=earg.a, .earg.scale=earg.scale
              ))),
    loglikelihood=eval(substitute(
            function(mu, y, w, residuals= FALSE,eta, extra=NULL) {
        aa = eta2theta(eta[,1], .link.a, earg= .earg.a)
        scale = eta2theta(eta[,2], .link.scale, earg= .earg.scale)
        parg = 1
        qq = aa
        if(residuals) stop("loglikelihood residuals not implemented yet") else
        sum(w*(log(aa) + (aa*parg-1)*log(y) - aa*parg*log(scale) +
   (if(is.R()) -lbeta(parg, qq) else lgamma(parg+qq)-lgamma(parg)-lgamma(qq))-
            (parg+qq)*log1p((y/scale)^aa)))
    }, list( .link.a=link.a, .link.scale=link.scale,
              .earg.a=earg.a, .earg.scale=earg.scale
            ))),
    vfamily=c("paralogistic"),
    deriv=eval(substitute(expression({
        aa = eta2theta(eta[,1], .link.a, earg= .earg.a)
        scale = eta2theta(eta[,2], .link.scale, earg= .earg.scale)
        parg = 1
        qq = aa

        temp1 = log(y/scale)
        temp2 = (y/scale)^aa
        temp3a = digamma(parg)
        temp3b = digamma(qq)

        dl.da = 1/aa + parg * temp1 - (parg+qq) * temp1 / (1+1/temp2)
        dl.dscale = (aa/scale) * (-parg + (parg+qq) / (1+1/temp2))
        da.deta = dtheta.deta(aa, .link.a, earg= .earg.a)
        dscale.deta = dtheta.deta(scale, .link.scale, earg= .earg.scale)
        w * cbind( dl.da * da.deta, dl.dscale * dscale.deta)
    }), list( .link.a=link.a, .link.scale=link.scale,
              .earg.a=earg.a, .earg.scale=earg.scale 
              ))),
    weight=eval(substitute(expression({
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
        wz = w * wz
        wz
    }), list( .link.a=link.a, .link.scale=link.scale,
              .earg.a=earg.a, .earg.scale=earg.scale
              ))))
}


 invparalogistic = function(link.a="loge",
                            link.scale="loge",
                    earg.a=list(), earg.scale=list(), 
                            init.a=1.0, 
                            init.scale=NULL,
                            zero=NULL)
{

    if(mode(link.a) != "character" && mode(link.a) != "name")
        link.a = as.character(substitute(link.a))
    if(mode(link.scale) != "character" && mode(link.scale) != "name")
        link.scale = as.character(substitute(link.scale))
    if(length(zero) && !is.Numeric(zero, integer=TRUE, posit=TRUE))
        stop("bad input for argument \"zero\"")
    if(!is.list(earg.a)) earg.a = list()
    if(!is.list(earg.scale)) earg.scale = list()

    new("vglmff",
    blurb=c("Inverse paralogistic distribution\n\n",
            "Links:    ",
            namesof("a", link.a, earg=earg.a), ", ", 
            namesof("scale", link.scale, earg=earg.scale), "\n", 
               "Mean:     scale*gamma(a + 1/a)*gamma(1 - 1/a)/gamma(a)"),
    constraints=eval(substitute(expression({
        constraints = cm.zero.vgam(constraints, x, .zero, M)
    }), list( .zero=zero ))),
    initialize=eval(substitute(expression({
        if(ncol(cbind(y)) != 1)
            stop("response must be a vector or a one-column matrix")
        predictors.names =
        c(namesof("a", .link.a, earg=.earg.a, tag=FALSE),
          namesof("scale", .link.scale, earg=.earg.scale, tag=FALSE))

        if(!length(.init.a) || !length(.init.scale)) {
            qvec = c(.25, .5, .75)   # Arbitrary; could be made an argument
            init.p = if(length(.init.a)) .init.a else 1
            xvec = log( qvec^(-1/ init.p ) - 1 )
            fit0 = lsfit(x=xvec, y=log(quantile(y, qvec )))
        }

        qq = 1
        if(!length(etastart)) {
            aa = rep(if(length(.init.a)) .init.a else -1/fit0$coef[2], length=n)
            scale = rep(if(length(.init.scale)) .init.scale else
                        exp(fit0$coef[1]), length=n)
            etastart = cbind(theta2eta(aa, .link.a, earg= .earg.a),
                             theta2eta(scale, .link.scale, earg= .earg.scale))
        }
    }), list( .link.a=link.a, .link.scale=link.scale,
              .earg.a=earg.a, .earg.scale=earg.scale,
              .init.a=init.a, .init.scale=init.scale ))),
    inverse=eval(substitute(function(eta, extra=NULL) {
        aa = eta2theta(eta[,1], .link.a, earg= .earg.a)
        scale = eta2theta(eta[,2], .link.scale, earg= .earg.scale)
        parg = aa
        qq = 1
        scale*gamma(parg + 1/aa)*gamma(qq-1/aa)/(gamma(parg)*gamma(qq))
    }, list( .link.a=link.a, .link.scale=link.scale,
             .earg.a=earg.a, .earg.scale=earg.scale ))),
    last=eval(substitute(expression({
        misc$link = c(a= .link.a, scale= .link.scale )
        misc$earg = list(a= .earg.a, scale= .earg.scale )
    }), list( .link.a=link.a, .link.scale=link.scale,
              .earg.a=earg.a, .earg.scale=earg.scale ))),
    loglikelihood=eval(substitute(
            function(mu, y, w, residuals= FALSE,eta, extra=NULL) {
        aa = eta2theta(eta[,1], .link.a, earg= .earg.a)
        scale = eta2theta(eta[,2], .link.scale, earg= .earg.scale)
        parg = aa
        qq = 1
        if(residuals) stop("loglikelihood residuals not implemented yet") else
        sum(w*(log(aa) + (aa*parg-1)*log(y) - aa*parg*log(scale) +
   (if(is.R()) -lbeta(parg, qq) else lgamma(parg+qq)-lgamma(parg)-lgamma(qq))-
            (parg+qq)*log1p((y/scale)^aa)))
    }, list( .link.a=link.a, .link.scale=link.scale,
             .earg.a=earg.a, .earg.scale=earg.scale ))),
    vfamily=c("invparalogistic"),
    deriv=eval(substitute(expression({
        aa = eta2theta(eta[,1], .link.a, earg= .earg.a)
        scale = eta2theta(eta[,2], .link.scale, earg= .earg.scale)
        parg = aa 
        qq = 1

        temp1 = log(y/scale)
        temp2 = (y/scale)^aa
        temp3a = digamma(parg)
        temp3b = digamma(qq)

        dl.da = 1/aa + parg * temp1 - (parg+qq) * temp1 / (1+1/temp2)
        dl.dscale = (aa/scale) * (-parg + (parg+qq) / (1+1/temp2))
        da.deta = dtheta.deta(aa, .link.a, earg= .earg.a)
        dscale.deta = dtheta.deta(scale, .link.scale, earg= .earg.scale)
        w * cbind( dl.da * da.deta, dl.dscale * dscale.deta )
    }), list( .link.a=link.a, .link.scale=link.scale,
              .earg.a=earg.a, .earg.scale=earg.scale ))),
    weight=eval(substitute(expression({
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
        wz = w * wz
        wz
    }), list( .link.a=link.a, .link.scale=link.scale,
              .earg.a=earg.a, .earg.scale=earg.scale ))))
}



if(FALSE)
genlognormal = function(link.sigma="loge", link.r="loge",
                        esigma=list(), er=list(),
                        init.sigma=1, init.r=1, zero=NULL)
{
warning(paste("2/4/04; doesn't work, possibly because first derivs are",
        "not continuous (sign() is used). Certainly, the derivs wrt",
        "mymu are problematic (run with maxit=4:9 and look at weight",
        "matrices). Possibly fundamentally cannot be estimated by IRLS.",
        "Pooling doesn't seem to help."))

    if(mode(link.sigma) != "character" && mode(link.sigma) != "name")
        link.sigma = as.character(substitute(link.sigma))
    if(mode(link.r) != "character" && mode(link.r) != "name")
        link.r = as.character(substitute(link.r))
    if(length(zero) && !is.Numeric(zero, integer=TRUE, posit=TRUE))
        stop("bad input for argument \"zero\"")
    if(!is.list(esigma)) esigma = list()
    if(!is.list(er)) er = list()

    new("vglmff",
    blurb=c("Three-parameter generalized lognormal distribution\n\n",
            "Links:    ",
            "loc; ", namesof("sigma", link.sigma, earg=esigma, tag= TRUE),
            ", ", namesof("r", link.r, earg=er, tag= TRUE)),
    constraints=eval(substitute(expression({
        constraints = cm.zero.vgam(constraints, x, .zero, M)
    }), list( .zero=zero ))),
    initialize=eval(substitute(expression({
        if(ncol(cbind(y)) != 1)
            stop("response must be a vector or a one-column matrix")
        predictors.names = 
        c(namesof("loc", "identity", earg= list(), tag=FALSE),
          namesof("sigma", .link.sigma, earg=.esigma, tag=FALSE),
          namesof("r", .link.r, earg=.er, tag=FALSE))

        if(!length(.init.sigma) || !length(.init.r)) {
            init.r = if(length(.init.r)) .init.r else 1
            sigma.init = (0.5 * sum(abs(log(y) - mean(log(y )))^init.r))^(1/init.r)
        }
        if(any(y <= 0)) stop("y must be positive")

        if(!length(etastart)) {
            sigma.init = rep(if(length( .init.sigma)) .init.sigma else sigma.init, len=n)
            r.init = if(length( .init.r)) .init.r else init.r
            etastart = cbind(mu=rep(log(median(y)), len=n),
                             sigma=sigma.init,
                             r = r.init)
        }
    }), list( .link.sigma=link.sigma, .link.r=link.r,
             .init.sigma=init.sigma, .init.r=init.r ))),
    inverse=eval(substitute(function(eta, extra=NULL) {
        mymu = eta2theta(eta[,1], "identity", earg=list())
        sigma = eta2theta(eta[,2], .link.sigma, earg= .esigma)
        r = eta2theta(eta[,3], .link.r, earg= .er)
        r
    }, list( .link.sigma=link.sigma, .link.r=link.r ))),
    last=eval(substitute(expression({
        misc$link = c(loc="identity", "sigma"= .link.sigma, r= .link.r )
        misc$expected = TRUE
    }), list( .link.sigma=link.sigma, .link.r=link.r ))),
    loglikelihood=eval(substitute(
        function(mu,y,w,residuals= FALSE,eta, extra=NULL) {
        mymu = eta2theta(eta[,1], "identity", earg=list())
        sigma = eta2theta(eta[,2], .link.sigma, earg= .esigma)
        r = eta2theta(eta[,3], .link.r, earg= .er)
        temp89 = (abs(log(y)-mymu)/sigma)^r
        if(residuals) stop("loglikelihood residuals not implemented yet") else
        sum(w * (-log(r^(1/r) * sigma) - lgamma(1+1/r) - temp89/r))
    }, list( .link.sigma=link.sigma, .link.r=link.r ))),
    vfamily=c("genlognormal3"),
    deriv=eval(substitute(expression({
        mymu = eta2theta(eta[,1], "identity", earg=list())
        sigma = eta2theta(eta[,2], .link.sigma, earg= .esigma)
        r = eta2theta(eta[,3], .link.r, earg= .er)
        ss = 1 + 1/r
        temp33 = (abs(log(y)-mymu)/sigma)
        temp33r1 = temp33^(r-1)
        dl.dmymu = temp33r1 * sign(log(y)-mymu) / sigma
        dl.dsigma = (temp33*temp33r1 - 1) / sigma
        dl.dr = (log(r) - 1 + digamma(ss) + temp33*temp33r1)/r^2 -
                temp33r1 * log(temp33r1) / r

        dmymu.deta = dtheta.deta(mymu, "identity", earg=list())
        dsigma.deta = dtheta.deta(sigma, .link.sigma, earg= .esigma)
        dr.deta = dtheta.deta(r, .link.r, earg= .er)
        w * cbind(dl.dmymu * dmymu.deta, 
                  dl.dsigma * dsigma.deta, 
                  dl.dr * dr.deta)
    }), list( .link.sigma=link.sigma, .link.r=link.r ))),
    weight=expression({
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
        wz = w * wz
        wz
    }))
}


betaprime = function(link="loge", earg=list(), i1=2, i2=NULL, zero=NULL)
{
    if(mode(link) != "character" && mode(link) != "name")
        link = as.character(substitute(link))
    if(!is.list(earg)) earg = list()

    new("vglmff",
    blurb=c("Beta-prime distribution\n",
            "y^(shape1-1) * (1+y)^(-shape1-shape2) / Beta(shape1,shape2),",
            " y>0, shape1>0, shape2>0\n\n",
            "Links:    ",
            namesof("shape1", link, earg=earg),  ", ",
            namesof("shape2", link, earg=earg), "\n",
            "Mean:     shape1/(shape2-1) provided shape2>1"),
    constraints=eval(substitute(expression({
        constraints = cm.zero.vgam(constraints, x, .zero, M)
    }), list( .zero=zero ))),
    initialize=eval(substitute(expression({
        if(ncol(y <- as.matrix(y)) > 1)
            stop("betaprime cannot handle matrix responses yet")
        if(min(y) <= 0)
            stop("response must be positive")
        predictors.names = c(namesof("shape1", .link, earg=.earg, short=TRUE),
                             namesof("shape2", .link, earg=.earg, short=TRUE))
        if(is.numeric(.i1) && is.numeric(.i2)) {
            vec = c(.i1, .i2)
            vec = c(theta2eta(vec[1], .link, earg= .earg),
                    theta2eta(vec[2], .link, earg= .earg))
            etastart = matrix(vec, n, 2, byrow= TRUE)
        }
        if(!length(etastart)) {
            init1 = if(length( .i1)) rep( .i1, len=n) else rep(1, len=n)
            init2 = if(length( .i2)) rep( .i2, len=n) else 1 + init1 / (y + 0.1)
            etastart = matrix(theta2eta(c(init1, init2), .link, earg= .earg),
                              n,2,byrow=TRUE)
        }
    }), list( .link=link, .earg=earg, .i1=i1, .i2=i2 ))), 
    inverse=eval(substitute(function(eta, extra=NULL) {
        shapes = eta2theta(eta, .link, earg= .earg)
        ifelse(shapes[,2] > 1, shapes[,1]/(shapes[,2]-1), NA)
    }, list( .link=link, .earg=earg ))),
    last=eval(substitute(expression({
        misc$link = c(shape1 = .link, shape2 = .link)
        misc$earg = list(shape1 = .earg, shape2 = .earg)
    }), list( .link=link, .earg=earg ))),
    loglikelihood=eval(substitute(
         function(mu, y, w, residuals= FALSE, eta, extra=NULL){
        shapes = eta2theta(eta, .link, earg= .earg)
        temp = if(is.R()) lbeta(shapes[,1], shapes[,2]) else
               lgamma(shapes[,1]) + lgamma(shapes[,2]) -
               lgamma(shapes[,1]+shapes[,2])
        if(residuals) stop("loglikelihood residuals not implemented yet") else
        sum(w *((shapes[,1]-1)*log(y)-(shapes[,2]+shapes[,1])*log1p(y)-temp))
    }, list( .link=link, .earg=earg ))),
    vfamily="betaprime",
    deriv=eval(substitute(expression({
        shapes = eta2theta(eta, .link, earg= .earg)
        dshapes.deta = dtheta.deta(shapes, .link, earg= .earg)
        dl.dshapes = cbind(log(y) - log1p(y) - digamma(shapes[,1]) + 
                           digamma(shapes[,1]+shapes[,2]),
                           - log1p(y) - digamma(shapes[,2]) + 
                           digamma(shapes[,1]+shapes[,2]))
        w * dl.dshapes * dshapes.deta
    }), list( .link=link, .earg=earg ))),
    weight=expression({
        temp2 = trigamma(shapes[,1]+shapes[,2])
        d2l.dshape12 = temp2 - trigamma(shapes[,1])
        d2l.dshape22 = temp2 - trigamma(shapes[,2])
        d2l.dshape1shape2 = temp2

        wz = matrix(as.numeric(NA), n, dimm(M))   #3=dimm(M)
        wz[,iam(1,1,M)] = d2l.dshape12 * dshapes.deta[,1]^2
        wz[,iam(2,2,M)] = d2l.dshape22 * dshapes.deta[,2]^2
        wz[,iam(1,2,M)] = d2l.dshape1shape2 * dshapes.deta[,1] * dshapes.deta[,2]

        -w * wz
    }))
}



maxwell = function(link="loge", earg=list()) {
    if(mode(link) != "character" && mode(link) != "name")
        link = as.character(substitute(link))
    if(!is.list(earg)) earg = list()

    new("vglmff",
    blurb=c("Maxwell distribution f(y) = sqrt(2/pi) * a^(3/2) * y^2 *",
            " exp(-0.5*a*y^2), y>0, a>0\n",
            "Link:    ", namesof("a", link, earg=earg), "\n", "\n",
            "Mean:    sqrt(8 / (a * pi))"),
    initialize=eval(substitute(expression({
        if(ncol(cbind(y)) != 1)
            stop("response must be a vector or a one-column matrix")
        predictors.names = namesof("a", .link, earg=.earg, tag=FALSE) 
        if(!length(etastart)) {
            a.init = rep(8 / (pi*(y+0.1)^2), length=length(y))
            etastart = theta2eta(a.init, .link, earg= .earg)
        }
    }), list( .link=link, .earg=earg ))),
    inverse=eval(substitute(function(eta, extra=NULL) {
        a = eta2theta(eta, .link, earg= .earg)
        sqrt(8 / (a * pi))
    }, list( .link=link, .earg=earg ))),
    last=eval(substitute(expression({
        misc$link = c(a= .link)
        misc$earg = list(a = .earg)
    }), list( .link=link, .earg=earg ))),
    loglikelihood=eval(substitute(
        function(mu,y,w,residuals= FALSE,eta, extra=NULL) {
        a = eta2theta(eta, .link, earg= .earg)
        if(residuals) stop("loglikelihood residuals not implemented yet") else
        sum(w * (1.5 * log(a) + 2 * log(y) - 0.5 * a * y^2 + 0.5*log(2/pi)))
    }, list( .link=link, .earg=earg ))),
    vfamily=c("maxwell"),
    deriv=eval(substitute(expression({
        a = eta2theta(eta, .link, earg= .earg)
        dl.da = 1.5 / a - 0.5 * y^2
        da.deta = dtheta.deta(a, .link, earg= .earg)
        w * dl.da * da.deta
    }), list( .link=link, .earg=earg ))),
    weight=eval(substitute(expression({
        ed2l.da2 = 1.5 / a^2
        wz = w * da.deta^2 * ed2l.da2
        wz
    }), list( .link=link, .earg=earg ))))
}


dmaxwell = function(x, a) {
    if(any(a <= 0)) stop("argument \"a\" must be positive")
    L = max(length(x), length(a)) 
    x = rep(x, len=L); a = rep(a, len=L); 
    ifelse(x>0, sqrt(2/pi) * a^(1.5) * x^2 * exp(-0.5*a*x^2), 0)
}

pmaxwell = function(q, a) {
    if(any(a <= 0)) stop("argument \"a\" must be positive")
    L = max(length(q), length(a)) 
    q = rep(q, len=L); a = rep(a, len=L); 
    ifelse(q > 0, erf(q*sqrt(a/2)) - q*exp(-0.5*a*q^2) * sqrt(2*a/pi), 0)
}

rmaxwell = function(n, a) {
    if(!is.Numeric(n, posit=TRUE, allow=1)) 
        stop("bad input for argument \"n\"")
    if(any(a <= 0)) stop("argument \"a\" must be positive")
    sqrt(2 * rgamma(n=n, 1.5) / a)
}

qmaxwell = function(p, a) {
    if(!is.Numeric(p, posit=TRUE) || any(p>=1)) 
        stop("bad input for argument \"p\"")
    if(any(a <= 0)) stop("argument \"a\" must be positive")
    N = max(length(p), length(a)); p = rep(p, len=N); a = rep(a, len=N)
    sqrt(2 * qgamma(p=p, 1.5) / a)
}




dnaka = function(x, shape, scale=1) {
    L = max(length(x), length(shape), length(scale))
    x = rep(x, len=L); shape = rep(shape, len=L); scale = rep(scale, len=L);
    ifelse(x <= 0, 0, dgamma(x=x^2, shape=shape, scale=scale/shape) * 2 * x)
}


pnaka = function(q, shape, scale=1) {
    if(!is.Numeric(q))
        stop("bad input for argument \"q\"")
    if(!is.Numeric(shape, posit=TRUE))
        stop("bad input for argument \"shape\"")
    if(!is.Numeric(scale, posit=TRUE))
        stop("bad input for argument \"scale\"")
    L = max(length(q), length(shape), length(scale))
    q = rep(q, len=L); shape = rep(shape, len=L); scale = rep(scale, len=L);
    ifelse(q <= 0, 0, pgamma(shape * q^2 / scale, shape))
}


qnaka = function(p, shape, scale=1, ...) {
    if(!is.Numeric(p, posit=TRUE) || max(p) >= 1)
        stop("bad input for argument \"p\"")
    if(!is.Numeric(shape, posit=TRUE))
        stop("bad input for argument \"shape\"")
    if(!is.Numeric(scale, posit=TRUE))
        stop("bad input for argument \"scale\"")
    L = max(length(p), length(shape), length(scale))
    p = rep(p, len=L); shape = rep(shape, len=L); scale = rep(scale, len=L);
    ans = rep(0.0, len=L)
    myfun = function(x, shape, scale=1, p)
        pnaka(q=x, shape=shape, scale=scale) - p
    for(i in 1:L) {
        EY = sqrt(scale[i]/shape[i]) * gamma(shape[i]+0.5) / gamma(shape[i])
        Upper = 5 * EY
        while(pnaka(q=Upper, shape=shape[i], scale=scale[i]) < p[i])
            Upper = Upper + scale[i]
        ans[i] = uniroot(f=myfun, lower=0, upper=Upper,
                         shape=shape[i], scale=scale[i], p=p[i], ...)$root
    }
    ans
}


rnaka = function(n, shape, scale=1, Smallno=1.0e-6) {
    if(!is.Numeric(n, posit=TRUE, integ=TRUE))
        stop("bad input for argument \"n\"")
    if(!is.Numeric(scale, posit=TRUE, allow=1))
        stop("bad input for argument \"scale\"")
    if(!is.Numeric(shape, posit=TRUE, allow=1))
        stop("bad input for argument \"shape\"")
    if(!is.Numeric(Smallno, posit=TRUE, allow=1) || Smallno > 0.01 ||
       Smallno < 2 * .Machine$double.eps)
        stop("bad input for argument \"Smallno\"")
    ans = rep(0.0, len=n)

    ptr1 = 1; ptr2 = 0
    ymax = dnaka(x=sqrt(scale*(1 - 0.5/shape)), shape=shape, scale=scale)
    while(ptr2 < n) {
        EY = sqrt(scale/shape) * gamma(shape+0.5) / gamma(shape)
        Upper = EY + 5 * scale
        while(pnaka(q=Upper, shape=shape, scale=scale) < 1-Smallno)
            Upper = Upper + scale
        x = runif(2*n, min=0, max=Upper)
        index = runif(2*n, max=ymax) < dnaka(x, shape=shape, scale=scale)
        sindex = sum(index)
        if(sindex) {
            ptr2 = min(n, ptr1 + sindex - 1)
            ans[ptr1:ptr2] = (x[index])[1:(1+ptr2-ptr1)]
            ptr1 = ptr2 + 1
        }
    }
    ans
}








nakagami = function(lshape="loge", lscale="loge",
                    eshape=list(), escale=list(), ishape=NULL, iscale=1) {
    if(mode(lshape) != "character" && mode(lshape) != "name")
        lshape = as.character(substitute(lshape))
    if(mode(lscale) != "character" && mode(lscale) != "name")
        lscale = as.character(substitute(lscale))
    if(!is.null(iscale) && !is.Numeric(iscale, positi=TRUE))
        stop("argument \"iscale\" must be a positive number or NULL")
    if(!is.list(eshape)) eshape = list()
    if(!is.list(escale)) escale = list()

    new("vglmff",
    blurb=c("Nakagami distribution f(y) = 2 * (shape/scale)^shape *\n",
            "                             ",
            "y^(2*shape-1) * exp(-shape*y^2/scale) / gamma(shape),\n",
            "                             ",
            "y>0, shape>0, scale>0\n",
            "Links:    ",
            namesof("shape", lshape, earg=eshape), ", ",
            namesof("scale", lscale, earg=escale),
            "\n",
            "\n",
            "Mean:    sqrt(scale/shape) * gamma(shape+0.5) / gamma(shape)"),
    initialize=eval(substitute(expression({
        if(ncol(cbind(y)) != 1)
            stop("response must be a vector or a one-column matrix")
        predictors.names = c(namesof("shape", .lshape, earg=.eshape, tag=FALSE),
                             namesof("scale", .lscale, earg=.escale, tag=FALSE))
        if(!length(etastart)) {
            init2 = if(is.Numeric( .iscale, posit=TRUE))
                        rep( .iscale, len=n) else rep(1, len=n)
            init1 = if(is.Numeric( .ishape, posit=TRUE))
                        rep( .ishape, len=n) else
                    rep(init2 / (y+1/8)^2, len=n)
            etastart = cbind(theta2eta(init1, .lshape, earg= .eshape),
                             theta2eta(init2, .lscale, earg= .escale))
        }
    }), list( .lscale=lscale, .lshape=lshape,
              .escale=escale, .eshape=eshape,
              .ishape=ishape, .iscale=iscale ))),
    inverse=eval(substitute(function(eta, extra=NULL) {
        shape = eta2theta(eta[,1], .lshape, earg= .eshape)
        scale = eta2theta(eta[,2], .lscale, earg= .escale)
        sqrt(scale/shape) * gamma(shape+0.5) / gamma(shape)
    }, list( .lscale=lscale, .lshape=lshape,
             .escale=escale, .eshape=eshape ))),
    last=eval(substitute(expression({
        misc$link = c(shape= .lshape, scale= .lscale)
        misc$earg = list(shape = .eshape, scale = .escale)
        misc$expected = TRUE
    }), list( .lscale=lscale, .lshape=lshape,
              .escale=escale, .eshape=eshape ))),
    loglikelihood=eval(substitute(
        function(mu,y,w,residuals= FALSE,eta, extra=NULL) {
        shape = eta2theta(eta[,1], .lshape, earg= .eshape)
        scale = eta2theta(eta[,2], .lscale, earg= .escale)
        if(residuals) stop("loglikelihood residuals not implemented yet") else
        sum(w * (log(2) + shape * log(shape/scale) - lgamma(shape) +
                 (2*shape-1) * log(y) - shape * y^2 / scale))
    }, list( .lscale=lscale, .lshape=lshape,
             .escale=escale, .eshape=eshape ))),
    vfamily=c("nakagami"),
    deriv=eval(substitute(expression({
        shape = eta2theta(eta[,1], .lshape, earg= .eshape)
        Scale = eta2theta(eta[,2], .lscale, earg= .escale)
        dl.dshape = 1 + log(shape/Scale) - digamma(shape) +
                    2 * log(y) - y^2 / Scale
        dl.dscale = -shape/Scale + shape * (y/Scale)^2
        dshape.deta = dtheta.deta(shape, .lshape, earg= .eshape)
        dscale.deta = dtheta.deta(Scale, .lscale, earg= .escale)
        w * cbind(dl.dshape * dshape.deta, dl.dscale * dscale.deta)
    }), list( .lscale=lscale, .lshape=lshape,
              .escale=escale, .eshape=eshape ))),
    weight=eval(substitute(expression({
        d2l.dshape2 = trigamma(shape) - 1/shape
        d2l.dscale2 = shape / Scale^2
        wz = matrix(as.numeric(NA), n, M)  # diagonal
        wz[,iam(1,1,M)] = d2l.dshape2 * dshape.deta^2
        wz[,iam(2,2,M)] = d2l.dscale2 * dscale.deta^2
        w * wz
    }), list( .lscale=lscale, .lshape=lshape,
              .escale=escale, .eshape=eshape ))))
}




rayleigh = function(link="loge", earg=list(), nrfs=1/3+0.01) {
    if(mode(link) != "character" && mode(link) != "name")
        link = as.character(substitute(link))
    if(!is.list(earg)) earg = list()
    if(!is.Numeric(nrfs, allow=1) || nrfs<0 || nrfs > 1)
        stop("bad input for 'nrfs'")

    new("vglmff",
    blurb=c("Rayleigh distribution f(y) = y*exp(-0.5*(y/a)^2)/a^2, y>0, a>0\n",
            "Link:    ",
            namesof("a", link, earg=earg), "\n\n",
            "Mean:    a * sqrt(pi / 2)"),
    initialize=eval(substitute(expression({
        if(ncol(cbind(y)) != 1)
            stop("response must be a vector or a one-column matrix")
        predictors.names = namesof("a", .link, earg=.earg, tag=FALSE) 
        if(!length(etastart)) {
            a.init = (y+1/8) / sqrt(pi/2)
            etastart = theta2eta(a.init, .link, earg= .earg)
        }
    }), list( .link=link, .earg=earg ))),
    inverse=eval(substitute(function(eta, extra=NULL) {
        a = eta2theta(eta, .link, earg= .earg)
        a * sqrt(pi/2)
    }, list( .link=link, .earg=earg ))),
    last=eval(substitute(expression({
        misc$link = c(a= .link)
        misc$earg = list(a = .earg)
        misc$nrfs = .nrfs
    }), list( .link=link, .earg=earg, .nrfs=nrfs  ))),
    loglikelihood=eval(substitute(
        function(mu,y,w,residuals= FALSE,eta, extra=NULL) {
        a = eta2theta(eta, .link, earg= .earg)
        if(residuals) stop("loglikelihood residuals not implemented yet") else
        sum(w * (log(y) - 2 * log(a) - 0.5 * (y/a)^2))
    }, list( .link=link, .earg=earg ))),
    vfamily=c("rayleigh"),
    deriv=eval(substitute(expression({
        a = eta2theta(eta, .link, earg= .earg)
        dl.da = ((y/a)^2 - 2) / a
        da.deta = dtheta.deta(a, .link, earg= .earg)
        w * dl.da * da.deta
    }), list( .link=link, .earg=earg ))),
    weight=eval(substitute(expression({
        d2l.da2 = (3 * (y/a)^2 - 2) / a^2
        ed2l.da2 = 4 / a^2
        wz = w * da.deta^2 * ((1- .nrfs) * d2l.da2 + .nrfs * ed2l.da2)
        wz
    }), list( .link=link, .earg=earg, .nrfs=nrfs ))))
}


drayleigh = function(x, a) {
    if(any(a <= 0)) stop("argument \"a\" must be positive")
    L = max(length(x), length(a)) 
    x = rep(x, len=L); a = rep(a, len=L); 
    ifelse(x>0, x*exp(-0.5*(x/a)^2)/a^2, 0)
}

prayleigh = function(q, a) {
    if(any(a <= 0)) stop("argument \"a\" must be positive")
    L = max(length(q), length(a)) 
    q = rep(q, len=L); a = rep(a, len=L); 
    ifelse(q > 0, 1 - exp(-0.5*(q/a)^2), 0)
}

qrayleigh = function(p, a) {
   if(any(a <= 0)) stop("argument \"a\" must be positive")
   if(any(p <= 0) || any(p >= 1)) stop("argument \"p\" must be between 0 and 1")
   a * sqrt(-2 * log1p(-p))
}

rrayleigh = function(n, a) {
    if(!is.Numeric(n, posit=TRUE, allow=1, integ=TRUE))
        stop("bad input for argument \"n\"")
    if(any(a <= 0)) stop("a must be positive")
    a * sqrt(-2 * log(runif(n )))
}





dparetoIV = function(x, location=0, scale=1, inequality=1, shape=1) {
    if(!is.Numeric(x)) stop("bad input for argument \"x\"")
    if(!is.Numeric(scale, posit=TRUE)) stop("bad input for argument \"scale\"")
    if(!is.Numeric(inequality, posi=TRUE)) 
        stop("bad input for argument \"inequality\"")
    if(!is.Numeric(shape, posit=TRUE)) stop("bad input for argument \"shape\"")
    N = max(length(x), length(location), length(scale), length(inequality),
            length(shape))
    x = rep(x, len=N); location = rep(location, len=N)
    scale = rep(scale, len=N); inequality = rep(inequality, len=N)
    shape = rep(shape, len=N)
    answer = x * 0
    ii = x > location
    zedd = (x[ii] - location[ii]) / scale[ii]
    answer[ii] = (shape[ii]/(scale[ii]*inequality[ii])) *
        zedd^(1/inequality[ii]-1) / (1+zedd^(1/inequality[ii]))^(shape[ii]+1)
    answer
}

pparetoIV = function(q, location=0, scale=1, inequality=1, shape=1) {
    if(!is.Numeric(q)) stop("bad input for argument \"q\"")
    if(!is.Numeric(scale, posit=TRUE)) 
        stop("bad input for argument \"scale\"")
    if(!is.Numeric(inequality, posi=TRUE)) 
        stop("bad input for argument \"inequality\"")
    if(!is.Numeric(shape, posit=TRUE)) 
        stop("bad input for argument \"shape\"")
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

qparetoIV = function(p, location=0, scale=1, inequality=1, shape=1) {
    if(!is.Numeric(p, posit=TRUE) || any(p >= 1)) 
        stop("bad input for argument \"p\"")
    if(!is.Numeric(scale, posit=TRUE)) 
        stop("bad input for argument \"scale\"")
    if(!is.Numeric(inequality, posi=TRUE)) 
        stop("bad input for argument \"inequality\"")
    if(!is.Numeric(shape, posit=TRUE)) 
        stop("bad input for argument \"shape\"")
    location + scale * (-1 + (1-p)^(-1/shape))^inequality
}

rparetoIV = function(n, location=0, scale=1, inequality=1, shape=1) {
    if(!is.Numeric(n, posit=TRUE, integ=TRUE, allow=1)) 
        stop("bad input for argument n")
    if(!is.Numeric(scale, posit=TRUE)) stop("bad input for argument \"scale\"")
    if(!is.Numeric(inequality, posi=TRUE)) 
        stop("bad input for argument \"inequality\"")
    if(!is.Numeric(shape, posit=TRUE)) stop("bad input for argument \"shape\"")
    location + scale * (-1 + runif(n)^(-1/shape))^inequality
}


dparetoIII = function(x, location=0, scale=1, inequality=1)
    dparetoIV(x=x, location=location, scale=scale, inequality=inequality, shape=1)

pparetoIII = function(q, location=0, scale=1, inequality=1)
    pparetoIV(q=q, location=location, scale=scale, inequality=inequality, shape=1)

qparetoIII = function(p, location=0, scale=1, inequality=1)
    qparetoIV(p=p, location=location, scale=scale, inequality=inequality, shape=1)

rparetoIII = function(n, location=0, scale=1, inequality=1)
    rparetoIV(n=n, location=location, scale=scale, inequality=inequality, shape=1)



dparetoII = function(x, location=0, scale=1, shape=1)
    dparetoIV(x=x, location=location, scale=scale, inequality=1, shape=shape)

pparetoII = function(q, location=0, scale=1, shape=1)
    pparetoIV(q=q, location=location, scale=scale, inequality=1, shape=shape)

qparetoII = function(p, location=0, scale=1, shape=1)
    qparetoIV(p=p, location=location, scale=scale, inequality=1, shape=shape)

rparetoII = function(n, location=0, scale=1, shape=1)
    rparetoIV(n=n, location=location, scale=scale, inequality=1, shape=shape)


dparetoI = function(x, scale=1, shape=1)
    dparetoIV(x=x, location=scale, scale=scale, inequality=1, shape=shape)

pparetoI = function(q, scale=1, shape=1)
    pparetoIV(q=q, location=scale, scale=scale, inequality=1, shape=shape)

qparetoI = function(p, scale=1, shape=1)
    qparetoIV(p=p, location=scale, scale=scale, inequality=1, shape=shape)

rparetoI = function(n, scale=1, shape=1)
    rparetoIV(n=n, location=scale, scale=scale, inequality=1, shape=shape)



paretoIV = function(location=0,
                    lscale="loge",
                    linequality="loge",
                    lshape="loge",
                    escale=list(), einequality=list(), eshape=list(),
                    iscale=1, iinequality=1, ishape=NULL,
                    method.init=1) {
    if(mode(lscale) != "character" && mode(lscale) != "name")
        lscale = as.character(substitute(lscale))
    if(mode(linequality) != "character" && mode(linequality) != "name")
        linequality = as.character(substitute(linequality))
    if(mode(lshape) != "character" && mode(lshape) != "name")
        lshape = as.character(substitute(lshape))
    if(!is.Numeric(location))
        stop("argument \"location\" must be numeric")
    if(is.Numeric(iscale) && any(iscale <= 0))
        stop("argument \"iscale\" must be positive")
    if(is.Numeric(iinequality) && any(iinequality <= 0))
        stop("argument \"iinequality\" must be positive")
    if(is.Numeric(ishape) && any(ishape <= 0))
        stop("argument \"ishape\" must be positive")
    if(!is.Numeric(method.init, allow=1, integ=TRUE) || method.init>2)
        stop("bad input for argument \"method.init\"")
    if(linequality == "nloge" && location != 0)
        warning("The Burr distribution has location=0 and linequality=nloge")
    if(!is.list(escale)) escale = list()
    if(!is.list(einequality)) einequality = list()
    if(!is.list(eshape)) eshape = list()

    new("vglmff",
    blurb=c("Pareto(IV) distribution F(y)=1-[1+((y - ", location,
            ")/scale)^(1/inequality)]^(-shape),",
            "\n", "         y > ",
            location, ", scale > 0, inequality > 0, shape > 0,\n",
            "Links:    ", namesof("scale", lscale, earg=escale ), ", ",
                          namesof("inequality", linequality, earg=einequality ), ", ",
                          namesof("shape", lshape, earg=eshape ), "\n",
            "Mean:    location + scale * NA"),  # zz
    initialize=eval(substitute(expression({
        if(ncol(cbind(y)) != 1)
            stop("response must be a vector or a one-column matrix")
        predictors.names = 
        c(namesof("scale", .lscale, earg=.escale, tag=FALSE),
          namesof("inequality", .linequality, earg=.einequality, tag=FALSE),
          namesof("shape", .lshape, earg=.eshape, tag=FALSE))
        extra$location = location = .location
        if(any(y <= location))
        stop("the response must have values > than the \"location\" argument")
        if(!length(etastart)) {
            inequality.init = if(length(.iinequality)) .iinequality else  1
            scale.init = if(length( .iscale)) .iscale else 1
            shape.init = if(length( .ishape)) .ishape else NULL
            if(!length(shape.init)) {
                zedd = (y - location) / scale.init
                if( .method.init == 1) {
                    A1 = weighted.mean(1/(1 + zedd^(1/inequality.init)), w=w)
                    A2 = weighted.mean(1/(1 + zedd^(1/inequality.init))^2, w=w)
                } else {
                    A1 = median(1/(1 + zedd^(1/inequality.init )))
                    A2 = median(1/(1 + zedd^(1/inequality.init))^2)
                }
                shape.init = max(0.01, (2*A2-A1)/(A1-A2))
            }
            etastart=cbind(
              theta2eta(rep(scale.init, len=n), .lscale, earg= .escale),
              theta2eta(rep(inequality.init, len=n), .linequality, earg= .einequality),
              theta2eta(rep(shape.init, len=n), .lshape, earg= .eshape))
        }
    }), list( .location=location, .lscale=lscale,
             .linequality=linequality, .lshape=lshape, .method.init=method.init,
             .escale=escale, .einequality=einequality, .eshape=eshape,
             .iscale=iscale, .iinequality=iinequality, .ishape=ishape ))),
    inverse=eval(substitute(function(eta, extra=NULL) {
        location = extra$location
        Scale = eta2theta(eta[,1], .lscale, earg= .escale)
        inequality = eta2theta(eta[,2], .linequality, earg= .einequality)
        shape = eta2theta(eta[,3], .lshape, earg= .eshape)
        location + Scale * NA
    }, list( .lscale=lscale, .linequality=linequality, .lshape=lshape,
             .escale=escale, .einequality=einequality, .eshape=eshape ))),
    last=eval(substitute(expression({
        misc$link=c("scale"= .lscale, "inequality"= .linequality,
                    "shape"= .lshape)
        misc$earg = list(scale = .escale, inequality= .einequality,
                         shape = .eshape)
        misc$location = extra$location # Use this for prediction
    }), list( .lscale=lscale,  .linequality=linequality, .lshape=lshape,
              .escale=escale, .einequality=einequality, .eshape=eshape ))),
    loglikelihood=eval(substitute(
        function(mu,y,w,residuals= FALSE,eta, extra=NULL) {
        location = extra$location
        Scale = eta2theta(eta[,1], .lscale, earg= .escale)
        inequality = eta2theta(eta[,2], .linequality, earg= .einequality)
        shape = eta2theta(eta[,3], .lshape, earg= .eshape)
        zedd = (y - location) / Scale
        if(residuals) stop("loglikelihood residuals not implemented yet") else
        sum(w * (log(shape) - log(inequality) - log(Scale) + (1/inequality -1) *
                 log(zedd) - (shape+1) * log1p(zedd^(1/inequality))))
    }, list( .lscale=lscale,  .linequality=linequality, .lshape=lshape,
             .escale=escale, .einequality=einequality, .eshape=eshape ))),
    vfamily=c("paretoIV"),
    deriv=eval(substitute(expression({
        location = extra$location
        Scale = eta2theta(eta[,1], .lscale, earg= .escale)
        inequality = eta2theta(eta[,2], .linequality, earg= .einequality)
        shape = eta2theta(eta[,3], .lshape, earg= .eshape)
        zedd = (y - location) / Scale
        temp100 = 1 + zedd^(1/inequality)
        dl.dscale = (shape  - (1+shape) / temp100) / (inequality * Scale)
        dl.dinequality = ((log(zedd) * (shape - (1+shape)/temp100)) /
                         inequality - 1) / inequality
        dl.dshape = -log(temp100) + 1/shape
        dscale.deta = dtheta.deta(Scale, .lscale, earg= .escale)
        dinequality.deta = dtheta.deta(inequality, .linequality, earg= .einequality)
        dshape.deta = dtheta.deta(shape, .lshape, earg= .eshape)
        w * cbind(dl.dscale * dscale.deta,
                  dl.dinequality * dinequality.deta, 
                  dl.dshape * dshape.deta)
    }), list( .lscale=lscale,  .linequality=linequality, .lshape=lshape,
              .escale=escale, .einequality=einequality, .eshape=eshape ))),
    weight=eval(substitute(expression({
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
        w * wz
    }), list( .lscale=lscale,  .linequality=linequality, .lshape=lshape,
              .escale=escale, .einequality=einequality, .eshape=eshape ))))
}

paretoIII = function(location=0,
                     lscale="loge",
                     linequality="loge",
                     escale=list(), einequality=list(),
                     iscale=NULL, iinequality=NULL) {
    if(mode(lscale) != "character" && mode(lscale) != "name")
        lscale = as.character(substitute(lscale))
    if(mode(linequality) != "character" && mode(linequality) != "name")
        linequality = as.character(substitute(linequality))
    if(!is.Numeric(location))
        stop("argument \"location\" must be numeric")
    if(is.Numeric(iscale) && any(iscale <= 0))
        stop("argument \"iscale\" must be positive")
    if(is.Numeric(iinequality) && any(iinequality <= 0))
        stop("argument \"iinequality\" must be positive")
    if(!is.list(escale)) escale = list()
    if(!is.list(einequality)) einequality = list()

    new("vglmff",
    blurb=c("Pareto(III) distribution F(y)=1-[1+((y - ", location,
            ")/scale)^(1/inequality)]^(-1),",
            "\n", "         y > ",
            location, ", scale > 0, inequality > 0, \n",
            "Links:    ", namesof("scale", lscale, earg=escale ), ", ",
                          namesof("inequality", linequality, earg=einequality ), "\n",
            "Mean:    location + scale * NA"),  # zz
    initialize=eval(substitute(expression({
        if(ncol(cbind(y)) != 1)
            stop("the response must be a vector or a one-column matrix")
        predictors.names = 
        c(namesof("scale", .lscale, earg=.escale, tag=FALSE),
          namesof("inequality", .linequality, earg=.einequality, tag=FALSE))
        extra$location = location = .location
        if(any(y <= location))
        stop("the response must have values > than the \"location\" argument")
        if(!length(etastart)) {
            inequality.init = if(length(.iinequality)) .iinequality else  NULL
            scale.init = if(length( .iscale)) .iscale else NULL
            if(!length(inequality.init) || !length(scale.init)) {
                probs = (1:4)/5
                ytemp = quantile(x=log(y-location), probs=probs)
                fittemp = lsfit(x=logit(probs), y=ytemp, int=TRUE)
                if(!length(inequality.init))
                    inequality.init = max(fittemp$coef["X"], 0.01)
                if(!length(scale.init))
                    scale.init = exp(fittemp$coef["Intercept"])
            }
            etastart=cbind(
            theta2eta(rep(scale.init, len=n), .lscale, earg= .escale),
            theta2eta(rep(inequality.init, len=n), .linequality, earg= .einequality))
        }
    }), list( .location=location, .lscale=lscale, .linequality=linequality,
              .escale=escale, .einequality=einequality,
              .iscale=iscale, .iinequality=iinequality ))),
    inverse=eval(substitute(function(eta, extra=NULL) {
        location = extra$location
        Scale = eta2theta(eta[,1], .lscale, earg= .escale)
        inequality = eta2theta(eta[,2], .linequality, earg= .einequality)
        location + Scale * NA
    }, list( .lscale=lscale, .linequality=linequality,
             .escale=escale, .einequality=einequality ))),
    last=eval(substitute(expression({
        misc$link=c("scale"= .lscale, "inequality"= .linequality)
        misc$earg = list(scale = .escale, inequality= .einequality)
        misc$location = extra$location # Use this for prediction
    }), list( .lscale=lscale, .linequality=linequality,
              .escale=escale, .einequality=einequality ))),
    loglikelihood=eval(substitute(
        function(mu,y,w,residuals= FALSE,eta, extra=NULL) {
        location = extra$location
        Scale = eta2theta(eta[,1], .lscale, earg= .escale)
        inequality = eta2theta(eta[,2], .linequality, earg= .einequality)
        zedd = (y - location) / Scale
        if(residuals) stop("loglikelihood residuals not implemented yet") else
        sum(w * (-log(inequality) - log(Scale) + (1/inequality -1) *
                 log(zedd) - (1+1) * log1p(zedd^(1/inequality))))
    }, list( .lscale=lscale, .linequality=linequality,
             .escale=escale, .einequality=einequality ))),
    vfamily=c("paretoIII"),
    deriv=eval(substitute(expression({
        location = extra$location
        Scale = eta2theta(eta[,1], .lscale, earg= .escale)
        inequality = eta2theta(eta[,2], .linequality, earg= .einequality)
        shape = 1
        zedd = (y - location) / Scale
        temp100 = 1 + zedd^(1/inequality)
        dl.dscale = (shape  - (1+shape) / temp100) / (inequality * Scale)
        dl.dinequality = ((log(zedd) * (shape - (1+shape)/temp100)) /
                         inequality - 1) / inequality
        dscale.deta = dtheta.deta(Scale, .lscale, earg= .escale)
        dinequality.deta = dtheta.deta(inequality, .linequality, earg= .einequality)
        w * cbind(dl.dscale * dscale.deta,
                  dl.dinequality * dinequality.deta)
    }), list( .lscale=lscale, .linequality=linequality,
              .escale=escale, .einequality=einequality ))),
    weight=eval(substitute(expression({
        d2scale.deta2 = 1 / ((inequality*Scale)^2 * 3)
        d2inequality.deta2 = (1 + 2* trigamma(1)) / (inequality^2 * 3)
        wz = matrix(0, n, M) # It is diagonal
        wz[,iam(1,1,M)] = dscale.deta^2 * d2scale.deta2
        wz[,iam(2,2,M)] = dinequality.deta^2 * d2inequality.deta2
        w * wz
    }), list( .lscale=lscale,  .linequality=linequality,
              .escale=escale, .einequality=einequality ))))
}


paretoII = function(location=0,
                    lscale="loge",
                    lshape="loge",
                    escale=list(), eshape=list(),
                    iscale=NULL, ishape=NULL) {
    if(mode(lscale) != "character" && mode(lscale) != "name")
        lscale = as.character(substitute(lscale))
    if(mode(lshape) != "character" && mode(lshape) != "name")
        lshape = as.character(substitute(lshape))
    if(!is.Numeric(location))
        stop("argument \"location\" must be numeric")
    if(is.Numeric(iscale) && any(iscale <= 0))
        stop("argument \"iscale\" must be positive")
    if(is.Numeric(ishape) && any(ishape <= 0))
        stop("argument \"ishape\" must be positive")
    if(!is.list(escale)) escale = list()
    if(!is.list(eshape)) eshape = list()

    new("vglmff",
    blurb=c("Pareto(II) distribution F(y)=1-[1+(y - ", location,
            ")/scale]^(-shape),",
            "\n", "         y > ",
            location, ", scale > 0,  shape > 0,\n",
            "Links:    ", namesof("scale", lscale, earg=escale ), ", ",
                          namesof("shape", lshape, earg=eshape ), "\n",
            "Mean:    location + scale * NA"),  # zz
    initialize=eval(substitute(expression({
        if(ncol(cbind(y)) != 1)
            stop("the response must be a vector or a one-column matrix")
        predictors.names = 
        c(namesof("scale", .lscale, earg=.escale, tag=FALSE),
          namesof("shape", .lshape, earg=.eshape, tag=FALSE))
        extra$location = location = .location
        if(any(y <= location))
        stop("the response must have values > than the \"location\" argument")
        if(!length(etastart)) {
            scale.init = if(length( .iscale)) .iscale else NULL
            shape.init = if(length( .ishape)) .ishape else  NULL
            if(!length(shape.init) || !length(scale.init)) {
                probs = (1:4)/5
                scale.init.0 = 1  # zz; have to put some value here...
                ytemp = quantile(x=log(y-location+scale.init.0), probs=probs)
                fittemp = lsfit(x=log1p(-probs), y=ytemp, int=TRUE)
                if(!length(shape.init))
                    shape.init = max(-1/fittemp$coef["X"], 0.01)
                if(!length(scale.init))
                    scale.init = exp(fittemp$coef["Intercept"])
            }
            etastart=cbind(
            theta2eta(rep(scale.init, len=n), .lscale, earg= .escale),
            theta2eta(rep(shape.init, len=n), .lshape, earg= .eshape))
        }
    }), list( .location=location, .lscale=lscale,
              .escale=escale, .eshape=eshape, 
              .lshape=lshape, .iscale=iscale, .ishape=ishape ))),
    inverse=eval(substitute(function(eta, extra=NULL) {
        location = extra$location
        Scale = eta2theta(eta[,1], .lscale, earg= .escale)
        shape = eta2theta(eta[,2], .lshape, earg= .eshape)
        location + Scale * NA
    }, list( .lscale=lscale, .lshape=lshape,
             .escale=escale, .eshape=eshape ))),
    last=eval(substitute(expression({
        misc$link=c("scale"= .lscale, "shape"= .lshape)
        misc$earg = list(scale = .escale, shape= .eshape)
        misc$location = extra$location # Use this for prediction
    }), list( .lscale=lscale, .lshape=lshape,
              .escale=escale, .eshape=eshape ))),
    loglikelihood=eval(substitute(
        function(mu,y,w,residuals= FALSE,eta, extra=NULL) {
        location = extra$location
        Scale = eta2theta(eta[,1], .lscale, earg= .escale)
        shape = eta2theta(eta[,2], .lshape, earg= .eshape)
        zedd = (y - location) / Scale
        if(residuals) stop("loglikelihood residuals not implemented yet") else
        sum(w * (log(shape) - log(Scale) - (shape+1) * log1p(zedd)))
    }, list( .lscale=lscale, .lshape=lshape,
             .escale=escale, .eshape=eshape ))),
    vfamily=c("paretoII"),
    deriv=eval(substitute(expression({
        location = extra$location
        Scale = eta2theta(eta[,1], .lscale, earg= .escale)
        shape = eta2theta(eta[,2], .lshape, earg= .eshape)
        zedd = (y - location) / Scale
        temp100 = 1 + zedd
        dl.dscale = (shape  - (1+shape) / temp100) / (1 * Scale)
        dl.dshape = -log(temp100) + 1/shape
        dscale.deta = dtheta.deta(Scale, .lscale, earg= .escale)
        dshape.deta = dtheta.deta(shape, .lshape, earg= .eshape)
        w * cbind(dl.dscale * dscale.deta,
                  dl.dshape * dshape.deta)
    }), list( .lscale=lscale, .lshape=lshape,
              .escale=escale, .eshape=eshape ))),
    weight=eval(substitute(expression({
        d2scale.deta2 = shape / (Scale^2 * (shape+2))
        d2shape.deta2 = 1 / shape^2
        d2ss.deta2 = -1 / (Scale * (shape+1))
        wz = matrix(0, n, dimm(M))
        wz[,iam(1,1,M)] = dscale.deta^2 * d2scale.deta2
        wz[,iam(2,2,M)] = dshape.deta^2 * d2shape.deta2
        wz[,iam(1,2,M)] = dscale.deta * dshape.deta * d2ss.deta2
        w * wz
    }), list( .lscale=lscale,  .lshape=lshape,
              .escale=escale, .eshape=eshape ))))
}



pareto1 = function(lshape="loge", earg=list(), location=NULL) {
    if(mode(lshape) != "character" && mode(lshape) != "name")
        lshape = as.character(substitute(lshape))
    if(is.Numeric(location) && location <= 0)
        stop("argument \"location\" must be positive")
    if(!is.list(earg)) earg = list()

    new("vglmff",
    blurb=c("Pareto distribution f(y) = shape * location^shape / y^(shape+1),",
            " y>location>0, shape>0\n",
            "Link:    ", namesof("shape", lshape, earg=earg), "\n", "\n",
            "Mean:    location*shape/(shape-1) for shape>1"),
    initialize=eval(substitute(expression({
        if(ncol(cbind(y)) != 1)
            stop("response must be a vector or a one-column matrix")
        predictors.names = namesof("shape", .lshape, earg=.earg, tag=FALSE) 
        locationhat = if(!length( .location)) {
            locationEstimated = TRUE
            min(y)
        } else {
            locationEstimated = FALSE
            .location
        }
        if(any(y < locationhat))
            stop("the value of location is too high (requires 0 < location < min(y))")
        extra$location = locationhat
        extra$locationEstimated = locationEstimated
        if(!length(etastart)) {
            k.init = (y + 1/8) / (y - locationhat + 1/8)
            etastart = theta2eta(k.init, .lshape, earg= .earg)
        }
    }), list( .lshape=lshape, .earg=earg, .location=location ))),
    inverse=eval(substitute(function(eta, extra=NULL) {
        k = eta2theta(eta, .lshape, earg= .earg)
        location = extra$location
        ifelse(k>1, k * location / (k-1), NA)
    }, list( .lshape=lshape, .earg=earg ))),
    last=eval(substitute(expression({
        misc$link = c(k= .lshape)
        misc$earg = list(k = .earg)
        misc$location = extra$location # Use this for prediction
    }), list( .lshape=lshape, .earg=earg ))),
    loglikelihood=eval(substitute(
        function(mu,y,w,residuals= FALSE,eta, extra=NULL) {
        k = eta2theta(eta, .lshape, earg= .earg)
        location = extra$location
        if(residuals) stop("loglikelihood residuals not implemented yet") else
        sum(w * (log(k) + k * log(location) - (k+1) * log(y )))
    }, list( .lshape=lshape, .earg=earg ))),
    vfamily=c("pareto1"),
    deriv=eval(substitute(expression({
        location = extra$location
        k = eta2theta(eta, .lshape, earg= .earg)
        dl.dk = 1/k + log(location/y)
        dk.deta = dtheta.deta(k, .lshape, earg= .earg)
        w * dl.dk * dk.deta
    }), list( .lshape=lshape, .earg=earg ))),
    weight=eval(substitute(expression({
        ed2l.dk2 = 1 / k^2
        wz = w * dk.deta^2 * ed2l.dk2
        wz
    }), list( .lshape=lshape, .earg=earg ))))
}


dpareto = function(x, location, shape) {
    if(any(location <= 0)) stop("argument \"location\" must be positive")
    if(any(shape <= 0)) stop("argument \"shape\" must be positive")
    L = max(length(x), length(location), length(shape)) 
    x = rep(x, len=L); location = rep(location, len=L); shape= rep(shape, len=L)
    ifelse(x>location, shape * location^shape / x^(shape+1), 0)
}

ppareto = function(q, location, shape) {
    if(any(location <= 0)) stop("argument \"location\" must be positive")
    if(any(shape <= 0)) stop("argument \"shape\" must be positive")
    L = max(length(q), length(location), length(shape))
    q = rep(q, len=L); location = rep(location, len=L); shape= rep(shape, len=L)
    ifelse(q > location, 1 - (location/q)^shape, 0)
}

qpareto = function(p, location, shape) {
    if(any(location <= 0)) stop("argument \"location\" must be positive")
    if(any(shape <= 0)) stop("argument \"shape\" must be positive")
   if(any(p <= 0) || any(p >= 1)) stop("argument \"p\" must be between 0 and 1")
    location / (1 - p)^(1/shape)
}

rpareto = function(n, location, shape) {
    if(!is.Numeric(n, posit=TRUE, integ=TRUE, allow=1)) 
        stop("bad input for argument \"n\"")
    if(any(location <= 0)) stop("argument \"location\" must be positive")
    if(any(shape <= 0)) stop("argument \"shape\" must be positive")
    location / runif(n)^(1/shape)
}


tpareto1 = function(lower, upper, lshape="loge", earg=list(), ishape=NULL,
                    method.init=1) {
    if(mode(lshape) != "character" && mode(lshape) != "name")
        lshape = as.character(substitute(lshape))
    if(!is.Numeric(lower, posit=TRUE, allow=1))
        stop("bad input for argument \"lower\"")
    if(!is.Numeric(upper, posit=TRUE, allow=1))
        stop("bad input for argument \"upper\"")
    if(lower >= upper)
        stop("lower < upper is required")
    if(length(ishape) && !is.Numeric(ishape, posit=TRUE))
        stop("bad input for argument \"ishape\"")
    if(!is.list(earg)) earg = list()
    if(!is.Numeric(method.init, allow=1, integ=TRUE, posit=TRUE) ||
       method.init > 2)
        stop("'method.init' must be 1 or 2")

    new("vglmff",
    blurb=c("Truncated Pareto distribution f(y) = shape * lower^shape /",
            "(y^(shape+1) * (1-(lower/upper)^shape)),",
            " 0 < lower < y < upper < Inf, shape>0\n",
            "Link:    ", namesof("shape", lshape, earg=earg), "\n", "\n",
            "Mean:    shape*lower^shape*(upper^(1-shape)-lower^(1-shape)) /",
                      " ((1-shape) * (1-(lower/upper)^shape))"),
    initialize=eval(substitute(expression({
        if(ncol(cbind(y)) != 1)
            stop("response must be a vector or a one-column matrix")
        predictors.names = namesof("shape", .lshape, earg=.earg, tag=FALSE) 
        if(any(y <= .lower))
            stop(paste("the value of argument \"lower\" is too high",
                 "(requires 0 < lower < min(y))"))
        extra$lower = .lower
        if(any(y >= .upper))
            stop(paste("the value of argument \"upper\" is too low",
                 "(requires max(y) < upper)"))
        extra$upper = .upper
        if(!length(etastart)) {
            shape.init = if(is.Numeric( .ishape)) 0 * y + .ishape else
            if( .method.init == 2) {
                0 * y + median(rep((y + 1/8) / (y - .lower + 1/8), times=w))
            } else {
                tpareto1.Loglikfun = function(shape, y, x, w, extraargs) {
                     myratio = .lower / .upper
                     sum(w * (log(shape) + shape * log( .lower) -
                              (shape+1) * log(y) - log1p(-myratio^shape)))
                 }
                 shape.grid = 2^((-4):4)
                 try.this = getMaxMin(shape.grid, objfun=tpareto1.Loglikfun,
                                      y=y,  x=x, w=w)
                try.this = rep(try.this, len=n)
                try.this
            }
            etastart = theta2eta(shape.init, .lshape, earg= .earg)
        }
    }), list( .ishape=ishape, .earg=earg, .lshape=lshape,
              .method.init=method.init,
              .lower=lower, .upper=upper ))),
    inverse=eval(substitute(function(eta, extra=NULL) {
        shape = eta2theta(eta, .lshape, earg= .earg)
        myratio = .lower / .upper
        constprop = shape * .lower^shape / (1 - myratio^shape)
        constprop * ( .upper^(1-shape) - .lower^(1-shape)) / (1-shape)
    }, list( .lshape=lshape, .earg=earg, .lower=lower, .upper=upper ))),
    last=eval(substitute(expression({
        misc$link = c(shape= .lshape)
        misc$earg = list(shape = .earg)
        misc$lower = extra$lower
        misc$upper = extra$upper
        misc$expected = TRUE
    }), list( .lshape=lshape, .earg=earg, .lower=lower, .upper=upper ))),
    loglikelihood=eval(substitute(
        function(mu,y,w,residuals= FALSE,eta, extra=NULL) {
        shape = eta2theta(eta, .lshape, earg= .earg)
        myratio = .lower / .upper
        if(residuals) stop("loglikelihood residuals not implemented yet") else
        sum(w * (log(shape) + shape * log( .lower) - (shape+1) * log(y) -
                 log1p(-myratio^shape)))
    }, list( .lshape=lshape, .earg=earg, .lower=lower, .upper=upper ))),
    vfamily=c("tpareto1"),
    deriv=eval(substitute(expression({
        shape = eta2theta(eta, .lshape, earg= .earg)
        myratio = .lower / .upper
        myratio2 =  myratio^shape 
        tmp330 = myratio2 * log(myratio) / (1 - myratio2)
        dl.dshape = 1/shape + log( .lower) - log(y) + tmp330 
        dshape.deta = dtheta.deta(shape, .lshape, earg= .earg)
        w * dl.dshape * dshape.deta
    }), list( .lshape=lshape, .earg=earg, .lower=lower, .upper=upper ))),
    weight=eval(substitute(expression({
        ed2l.dshape2 = 1 / shape^2 - tmp330^2 / myratio2
        wz = w * dshape.deta^2 * ed2l.dshape2
        wz
    }), list( .lshape=lshape, .earg=earg, .lower=lower, .upper=upper ))))
}

dtpareto = function(x, lower, upper, shape) {
    if(!is.Numeric(x)) stop("bad input for argument \"x\"")
    if(!is.Numeric(lower, pos=TRUE)) stop("argument \"lower\" must be positive")
    if(!is.Numeric(upper, pos=TRUE)) stop("argument \"upper\" must be positive")
    if(!is.Numeric(shape, pos=TRUE)) stop("argument \"shape\" must be positive")
    L = max(length(x), length(lower), length(upper), length(shape)) 
    x = rep(x, len=L); lower = rep(lower, len=L); upper = rep(upper, len=L);
    shape= rep(shape, len=L)
    ifelse(x > lower & x < upper, 
    shape * lower^shape / (x^(shape+1) * (1-(lower/upper)^shape)), 0)
}

ptpareto = function(q, lower, upper, shape) {
    if(!is.Numeric(q)) stop("bad input for argument \"q\"")
    if(!is.Numeric(lower, pos=TRUE)) stop("argument \"lower\" must be positive")
    if(!is.Numeric(upper, pos=TRUE)) stop("argument \"upper\" must be positive")
    if(!is.Numeric(shape, pos=TRUE)) stop("argument \"shape\" must be positive")
    L = max(length(q), length(lower), length(upper), length(shape)) 
    q = rep(q, len=L); lower = rep(lower, len=L); upper = rep(upper, len=L);
    shape= rep(shape, len=L)
    ans = q * 0
    ans[q > lower & q < upper] = (1-(lower/q)^shape) / (1-(lower/upper)^shape)
    ans[q >= upper] = 1
    ans
}

qtpareto = function(p, lower, upper, shape) {
    if(!is.Numeric(p, posit=TRUE)) stop("bad input for argument \"p\"")
    if(!is.Numeric(lower, pos=TRUE)) stop("argument \"lower\" must be positive")
    if(!is.Numeric(upper, pos=TRUE)) stop("argument \"upper\" must be positive")
    if(!is.Numeric(shape, pos=TRUE)) stop("argument \"shape\" must be positive")
    if(max(p) >= 1) stop("argument \"p\" must be between 0 and 1")
    lower / (1 - p*(1-(lower/upper)^shape))^(1/shape)
}

rtpareto = function(n, lower, upper, shape) {
    if(!is.Numeric(n, posit=TRUE, integ=TRUE, allow=1)) 
        stop("bad input for argument \"n\"")
    if(!is.Numeric(lower, pos=TRUE)) stop("argument \"lower\" must be positive")
    if(!is.Numeric(upper, pos=TRUE)) stop("argument \"upper\" must be positive")
    if(!is.Numeric(shape, pos=TRUE)) stop("argument \"shape\" must be positive")
    lower / (1 - runif(n)*(1-(lower/upper)^shape))^(1/shape)
}


erf = function(x)
    2 * pnorm(x * sqrt(2)) - 1

erfc = function(x)
    2 * pnorm(x * sqrt(2), lower=FALSE)



wald <- function(link.lambda="loge", earg=list(), init.lambda=NULL)
{
    if(mode(link.lambda) != "character" && mode(link.lambda) != "name")
        link.lambda = as.character(substitute(link.lambda))
    if(!is.list(earg)) earg = list()

    new("vglmff",
    blurb=c("Standard Wald distribution\n\n",
           "f(y) = sqrt(lambda/(2*pi*y^3)) * exp(-lambda*(y-1)^2/(2*y)), y&lambda>0",
           "\n", 
           "Link:     ", 
                         namesof("lambda", link.lambda, earg=earg), "\n",
           "Mean:     ", "1\n",
           "Variance: 1 / lambda"),
    initialize=eval(substitute(expression({
        if(ncol(cbind(y)) != 1)
            stop("response must be a vector or a one-column matrix")
        if(any(y <= 0)) stop("Require the response to have positive values")
        predictors.names = 
        namesof("lambda", .link.lambda, earg=.earg, short=TRUE)
        if(!length(etastart)) {
            initlambda = if(length( .init.lambda)) .init.lambda else
                         1 / (0.01 + (y-1)^2)
            initlambda = rep(initlambda, len=n)
            etastart = cbind(theta2eta(initlambda, link=.link.lambda, earg= .earg))
        }
    }), list( .link.lambda=link.lambda, .earg=earg,
             .init.lambda=init.lambda ))),
    inverse=function(eta, extra=NULL) {
        0*eta + 1
    },
    last=eval(substitute(expression({
        misc$link = c(lambda = .link.lambda )
        misc$earg = list(lambda = .earg )
    }), list( .link.lambda=link.lambda, .earg=earg ))),
    loglikelihood=eval(substitute(
             function(mu, y, w, residuals = FALSE, eta, extra=NULL) {
        lambda = eta2theta(eta, link=.link.lambda, earg= .earg)
        if(residuals) stop("loglikelihood residuals not implemented yet") else
        sum(w*(0.5 * log(lambda / (2 * pi * y^3)) -
                lambda *(y-1)^2 / (2* y )))
    }, list( .link.lambda=link.lambda, .earg=earg ))),
    vfamily="wald",
    deriv=eval(substitute(expression({
        lambda = eta2theta(eta, link=.link.lambda, earg= .earg)
        dl.dlambda = 0.5 / lambda + 1 - 0.5 * (y + 1/y)
        dlambda.deta = dtheta.deta(theta=lambda, link=.link.lambda, earg= .earg)
        w * cbind(dl.dlambda * dlambda.deta)
    }), list( .link.lambda=link.lambda, .earg=earg ))),
    weight=eval(substitute(expression({
        d2l.dlambda2 = 0.5 / (lambda^2)
        w * cbind(dlambda.deta^2 * d2l.dlambda2)
    }), list( .link.lambda=link.lambda, .earg=earg ))))
}


expexp = function(lshape="loge", lscale="loge",
                  eshape=list(), escale=list(),
                  ishape=1.1, iscale=NULL,  # ishape cannot be 1
                  tolerance = 1.0e-6,
                  zero=NULL) {

    if(mode(lshape) != "character" && mode(lshape) != "name")
        lshape = as.character(substitute(lshape))
    if(mode(lscale) != "character" && mode(lscale) != "name")
        lscale = as.character(substitute(lscale))
    if(length(zero) && !is.Numeric(zero, integer=TRUE, posit=TRUE))
        stop("bad input for argument \"zero\"")
    if(!is.Numeric(tolerance, posit=TRUE, allow=1) || tolerance>1.0e-2)
        stop("bad input for argument \"tolerance\"")
    if(!is.Numeric(ishape, posit=TRUE))
        stop("bad input for argument \"ishape\"")
    ishape[ishape==1] = 1.1   # Fails in @deriv
    if(!is.list(escale)) escale = list()
    if(!is.list(eshape)) eshape = list()

    new("vglmff",
    blurb=c("Exponentiated Exponential Distribution\n",
           "Links:    ",
           namesof("shape", lshape, earg=eshape), ", ",
           namesof("scale", lscale, earg=escale),"\n",
           "Mean:     (digamma(shape+1)-digamma(1))/scale"),
    constraints=eval(substitute(expression({
        constraints = cm.zero.vgam(constraints, x, .zero, M)
    }), list( .zero=zero ))),
    initialize=eval(substitute(expression({
        if(ncol(cbind(y)) != 1)
            stop("response must be a vector or a one-column matrix")
        predictors.names = 
        c(namesof("shape", .lshape, earg=.eshape, short=TRUE), 
          namesof("scale", .lscale, earg=.escale, short=TRUE))
        if(!length(etastart)) {
            shape.init = if(!is.Numeric( .ishape, posit=TRUE))
                   stop("argument \"ishape\" must be positive") else
                   rep(.ishape, len=n)
            scale.init = if(length( .iscale)) rep(.iscale, len=n) else
                        (digamma(shape.init+1) - digamma(1)) / (y+1/8)
            scale.init = rep(weighted.mean(scale.init, w=w), len=n)
            etastart = cbind(theta2eta(shape.init, .lshape, earg= .eshape),
                             theta2eta(scale.init, .lscale, earg= .escale))
        }
    }), list( .lshape=lshape, .lscale=lscale, .iscale=iscale, .ishape=ishape,
              .eshape=eshape, .escale=escale ))),
    inverse=eval(substitute(function(eta, extra=NULL) {
        shape = eta2theta(eta[,1], .lshape, earg= .eshape)
        scale = eta2theta(eta[,2], .lscale, earg= .escale)
        (digamma(shape+1)-digamma(1)) / scale
    }, list( .lshape=lshape, .lscale=lscale,
             .eshape=eshape, .escale=escale ))),
    last=eval(substitute(expression({
        misc$link = c("shape"= .lshape, "scale"= .lscale)
        misc$earg = list(shape= .eshape, scale= .escale)
        misc$expected = TRUE
    }), list( .lshape=lshape, .lscale=lscale,
              .eshape=eshape, .escale=escale ))),
    loglikelihood= eval(substitute(
        function(mu, y, w, residuals = FALSE, eta, extra=NULL) {
        shape = eta2theta(eta[,1], .lshape, earg= .eshape)
        scale = eta2theta(eta[,2], .lscale, earg= .escale)
        if(residuals) stop("loglikelihood residuals not implemented yet") else
        sum(w * (log(shape) + log(scale) + 
                 (shape-1)*log1p(-exp(-scale*y)) - scale*y))
    }, list( .lscale=lscale, .lshape=lshape,
             .eshape=eshape, .escale=escale ))),
    vfamily=c("expexp"),
    deriv=eval(substitute(expression({
        shape = eta2theta(eta[,1], .lshape, earg= .eshape)
        scale = eta2theta(eta[,2], .lscale, earg= .escale)
        dl.dscale = 1/scale + (shape-1)*y*exp(-scale*y) / (1-exp(-scale*y)) - y
        dl.dshape = 1/shape + log1p(-exp(-scale*y))
        dscale.deta = dtheta.deta(scale, .lscale, earg= .escale)
        dshape.deta = dtheta.deta(shape, .lshape, earg= .eshape)
        w * cbind(dl.dshape * dshape.deta, dl.dscale * dscale.deta)
    }), list( .lshape=lshape, .lscale=lscale,
              .eshape=eshape, .escale=escale ))),
    weight=eval(substitute(expression({
        d11 = 1 / shape^2  # True for all shape
        d22 = d12 = rep(as.numeric(NA), len=n)
        index2 = abs(shape - 2) > .tolerance  # index2 = shape != 1
        largeno = 10000
        if(any(index2)) {
            Shape = shape[index2]
            Shape[abs(Shape-1) < .tolerance] = 1.001 # digamma(0) is undefined
            Scale = scale[index2]
            tmp200 = trigamma(1)-trigamma(Shape-1) +
                  (digamma(Shape-1)-digamma(1))^2    # Fails when Shape==1
            tmp300 = trigamma(1)-digamma(Shape)+(digamma(Shape)-digamma(1))^2
            d22[index2] = (1 + Shape*(Shape-1)*tmp200/(Shape-2)) / Scale^2 +
                          Shape*tmp300 / Scale^2
        }
        if(any(!index2)) {
            Scale = scale[!index2]
            d22[!index2] = (1 + 4 * sum(1/(2 + (0:largeno))^3)) / Scale^2
        }

        index1 = abs(shape - 1) > .tolerance  # index1 = shape != 1
        if(any(index1)) {
            Shape = shape[index1]
            Scale = scale[index1]
            d12[index1] = -(Shape*(digamma(Shape)-digamma(1))/(Shape-1) -
                          digamma(Shape+1) + digamma(1)) / Scale
        }
        if(any(!index1)) {
            Scale = scale[!index1]
            d12[!index1] = -sum(1/(2 + (0:largeno))^2) / Scale
        }
        wz = matrix(0, n, dimm(M))
        wz[,iam(1,1,M)] = dshape.deta^2 * d11
        wz[,iam(2,2,M)] = dscale.deta^2 * d22
        wz[,iam(1,2,M)] = dscale.deta * dshape.deta * d12
        w * wz
    }), list( .tolerance=tolerance ))))
}





expexp1 = function(lscale="loge",
                   escale=list(),
                   iscale=NULL,
                   ishape=1) {
    if(mode(lscale) != "character" && mode(lscale) != "name")
        lscale = as.character(substitute(lscale))
    if(!is.list(escale)) escale = list()

    new("vglmff",
    blurb=c("Exponentiated Exponential Distribution",
            " (profile likelihood estimation)\n",
           "Links:    ",
           namesof("scale", lscale, earg=escale), "\n",
           "Mean:     (digamma(shape+1)-digamma(1))/scale"),
    initialize=eval(substitute(expression({
        if(ncol(cbind(y)) != 1)
            stop("response must be a vector or a one-column matrix")
        predictors.names = namesof("scale", .lscale, earg=.escale, short=TRUE)
        if(length(w) != n || !is.Numeric(w, integer=TRUE, posit=TRUE))
            stop("weights must be a vector of positive integers")
        if(!intercept.only)
  stop("this family function only works for an intercept-only, i.e., y ~ 1")
        extra$yvector = y
        extra$sumw = sum(w)
        extra$w = w
        if(!length(etastart)) {
            shape.init = if(!is.Numeric( .ishape, posit=TRUE))
                   stop("argument \"ishape\" must be positive") else
                   rep(.ishape, len=n)
            scaleinit = if(length( .iscale)) rep(.iscale, len=n) else
                        (digamma(shape.init+1) - digamma(1)) / (y+1/8)  
            etastart = cbind(theta2eta(scaleinit, .lscale, earg= .escale))
        }
    }), list( .lscale=lscale, .iscale=iscale, .ishape=ishape,
              .escale=escale ))),
    inverse=eval(substitute(function(eta, extra=NULL) {
        scale = eta2theta(eta, .lscale, earg= .escale)
        temp7 = 1 - exp(-scale*extra$yvector)
        shape = -extra$sumw / sum(extra$w*log(temp7)) # \gamma(\theta)
        (digamma(shape+1)-digamma(1)) / scale
    }, list( .lscale=lscale,
             .escale=escale ))),
    last=eval(substitute(expression({
        misc$link = c("scale"= .lscale)
        misc$earg = list(scale= .escale)
        temp7 = 1 - exp(-scale*y)
        shape = -extra$sumw / sum(w*log(temp7)) # \gamma(\theta)
        misc$shape = shape   # Store the ML estimate here
        misc$pooled.weight = pooled.weight
    }), list( .lscale=lscale, .escale=escale ))),
    loglikelihood= eval(substitute(
        function(mu, y, w, residuals = FALSE, eta, extra=NULL) {
        scale = eta2theta(eta, .lscale, earg= .escale)
        temp7 = 1 - exp(-scale*y)
        shape = -extra$sumw / sum(w*log(temp7)) # \gamma(\theta)
        if(residuals) stop("loglikelihood residuals not implemented yet") else
        sum(w * (log(shape) + log(scale) + 
                 (shape-1)*log1p(-exp(-scale*y)) - scale*y))
    }, list( .lscale=lscale, .escale=escale ))),
    vfamily=c("expexp1"),
    deriv=eval(substitute(expression({
        scale = eta2theta(eta, .lscale, earg= .escale)
        temp6 = exp(-scale*y)
        temp7 = 1-temp6
        shape = -extra$sumw / sum(w*log(temp7)) # \gamma(\theta)
        d1 = 1/scale + (shape-1)*y*temp6/temp7 - y
        w * cbind(d1 * dtheta.deta(scale, .lscale, earg= .escale))
    }), list( .lscale=lscale, .escale=escale ))),
    weight=eval(substitute(expression({
        d11 = 1/scale^2  + y*(temp6/temp7^2) * ((shape-1) *
              (y*temp7+temp6) - y*temp6 / (log(temp7))^2)
        wz = matrix(0, n, dimm(M))
        wz[,iam(1,1,M)] = dtheta.deta(scale, .lscale, earg= .escale)^2 * d11 -
                          d2theta.deta2(scale, .lscale, earg= .escale) * d1

        if(FALSE && intercept.only) {
            sumw = sum(w)
            for(i in 1:ncol(wz))
                wz[,i] = sum(wz[,i]) / sumw
            pooled.weight = TRUE
            wz = w * wz   # Put back the weights
        } else
            pooled.weight = FALSE
        w * wz
    }), list( .lscale=lscale, .escale=escale ))))
}



betaffqn.control <- function(save.weight=TRUE, ...)
{
    list(save.weight=save.weight)
}



betaffqn = function(link="loge", earg=list(),
                    i1=NULL, i2=NULL, trim=0.05, A=0, B=1)
{
    if(mode(link) != "character" && mode(link) != "name")
        link = as.character(substitute(link))

    if(!is.Numeric(A, allow=1) || !is.Numeric(B, allow=1) || A >= B)
        stop("A must be < B, and both must be of length one")
    stdbeta = (A==0 && B==1)  # stdbeta==T iff standard beta distribution
    if(!is.list(earg)) earg = list()

    new("vglmff",
    blurb=c("Two-parameter Beta distribution\n",
            if(stdbeta)
            "y^(shape1-1) * (1-y)^(shape2-1), 0<=y<=1, shape1>0, shape2>0\n\n"
            else
            paste("(y-",A,")^(shape1-1) * (",B,
            "-y)^(shape2-1), ",A,"<=y<=",B," shape1>0, shape2>0\n\n", sep=""),
            "Links:    ",
            namesof("shape1", link, earg=earg),  ", ",
            namesof("shape2", link, earg=earg)),
    initialize=eval(substitute(expression({
        if(ncol(cbind(y)) != 1)
            stop("response must be a vector or a one-column matrix")
        if(min(y) <= .A || max(y) >= .B)
            stop("data not within (A, B)")
        predictors.names = c(namesof("shape1", .link, earg=.earg, short=TRUE),
                             namesof("shape2", .link, earg=.earg, short=TRUE))
        if(is.numeric(.i1) && is.numeric(.i2)) {
            vec = c(.i1, .i2)
            vec = c(theta2eta(vec[1], .link, earg= .earg),
                    theta2eta(vec[2], .link, earg= .earg))
            etastart = matrix(vec, n, 2, byrow= TRUE)
        }

        # For QN update below
        if(length(w) != n || !is.Numeric(w, posit=TRUE))
            stop("weights must be a vector of positive weights")

        if(!length(etastart)) {
            mu1d = mean(y, trim=.trim)
            uu = (mu1d-.A) / (.B - .A) 
            DD = (.B - .A)^2 
            pinit = uu^2 * (1-uu)*DD/var(y) - uu   # But var(y) is not robust
            qinit = pinit * (1-uu) / uu
            etastart = matrix(theta2eta(c(pinit,qinit), .link, earg= .earg),
                              n,2,byrow=TRUE)
        }
    }), list( .link=link, .earg=earg, .i1=i1, .i2=i2, .trim=trim, .A=A, .B=B ))), 
    inverse=eval(substitute(function(eta, extra=NULL) {
        shapes = eta2theta(eta, .link, earg= .earg)
        .A + (.B-.A) * shapes[,1] / (shapes[,1] + shapes[,2])
    }, list( .link=link, .earg=earg, .A=A, .B=B ))),
    last=eval(substitute(expression({
        misc$link = c(shape1 = .link, shape2 = .link)
        misc$earg = list(shape1 = .earg, shape2 = .earg)
        misc$limits = c(.A, .B)
        misc$expected = FALSE
        misc$BFGS = TRUE
    }), list( .link=link, .earg=earg, .A=A, .B=B ))),
    loglikelihood=eval(substitute(
         function(mu, y, w, residuals= FALSE, eta, extra=NULL){
        shapes = eta2theta(eta, .link, earg= .earg)
        temp = if(is.R()) lbeta(shapes[,1], shapes[,2]) else
               lgamma(shapes[,1]) + lgamma(shapes[,2]) -
               lgamma(shapes[,1]+shapes[,2])
        if(residuals) stop("loglikelihood residuals not implemented yet") else
        sum(w * ((shapes[,1]-1)*log(y-.A) + (shapes[,2]-1)*log(.B-y) - temp -
            (shapes[,1]+shapes[,2]-1)*log(.B-.A )))
    }, list( .link=link, .earg=earg, .A=A, .B=B ))),
    vfamily="betaffqn",
    deriv=eval(substitute(expression({
        shapes = eta2theta(eta, .link, earg= .earg)
        dshapes.deta = dtheta.deta(shapes, .link, earg= .earg)
        dl.dshapes = cbind(log(y-.A), log(.B-y)) - digamma(shapes) +
                     digamma(shapes[,1] + shapes[,2]) - log(.B - .A)
        if(iter == 1) {
            etanew = eta
        } else {
            derivold = derivnew
            etaold = etanew
            etanew = eta
        }
        derivnew = w * dl.dshapes * dshapes.deta
        derivnew
    }), list( .link=link, .earg=earg, .A=A, .B=B ))),
    weight=expression({
        if(iter == 1) {
            wznew = cbind(matrix(w, n, M), matrix(0, n, dimm(M)-M))
        } else {
            wzold = wznew
            wznew = qnupdate(w=w, wzold=wzold, dderiv=(derivold - derivnew),
                             deta=etanew-etaold, M=M,
                             trace=trace)  # weights incorporated in args
        }
        wznew
    }))
}



logistic2 = function(llocation="identity",
                     lscale="loge",
                     elocation=list(),
                     escale=list(),
                     ilocation=NULL, iscale=NULL,
                     method.init=1, zero=NULL) {
    if(mode(llocation) != "character" && mode(llocation) != "name")
        llocation = as.character(substitute(llocation))
    if(mode(lscale) != "character" && mode(lscale) != "name")
        lscale = as.character(substitute(lscale))
    if(!is.Numeric(method.init, allow=1, integ=TRUE, posit=TRUE) ||
       method.init > 2) stop("argument \"method.init\" must be 1 or 2")
    if(length(zero) && !is.Numeric(zero, integer=TRUE, posit=TRUE))
        stop("bad input for argument \"zero\"")
    if(!is.list(elocation)) elocation = list()
    if(!is.list(escale)) escale = list()

    new("vglmff",
    blurb=c("Two-parameter logistic distribution\n\n",
            "Links:    ",
            namesof("location", llocation, earg=elocation), ", ",
            namesof("scale", lscale, earg=escale),
            "\n", "\n",
            "Mean:     location", "\n",
            "Variance: (pi*scale)^2 / 3"),
    constraints=eval(substitute(expression({
        constraints = cm.zero.vgam(constraints, x, .zero, M)
    }), list( .zero=zero ))),
    initialize=eval(substitute(expression({
        if(ncol(cbind(y)) != 1)
            stop("response must be a vector or a one-column matrix")
        predictors.names = 
        c(namesof("location", .llocation, earg=.elocation, tag=FALSE),
          namesof("scale", .lscale, earg=.escale, tag=FALSE))
        if(!length(etastart)) {
            if( .method.init == 1) {
                location.init = y
                scale.init = sqrt(3) * sd(y) / pi
            } else {
                location.init = median(rep(y, w))
                scale.init = sqrt(3) * sum(w*(y-location.init)^2) / (sum(w)*pi)
            }
            location.init = if(length(.ilocation)) rep(.ilocation, len=n) else
                             rep(location.init, len=n)
            if(.llocation == "loge") location.init = abs(location.init) + 0.001
            scale.init = if(length(.iscale)) rep(.iscale, len=n) else
                             rep(1, len=n)
            etastart = cbind(
            theta2eta(location.init, .llocation, earg= .elocation),
            theta2eta(scale.init, .lscale, earg= .escale))
        }
    }), list( .method.init=method.init, .ilocation=ilocation,
              .elocation=elocation, .escale=escale,
              .llocation=llocation, .iscale=iscale, .lscale=lscale ))),
    inverse=eval(substitute(function(eta, extra=NULL) {
        eta2theta(eta[,1], .llocation, earg= .elocation)
    }, list( .llocation=llocation,
             .elocation=elocation, .escale=escale ))),
    last=eval(substitute(expression({
        misc$link = c(location=.llocation, scale= .lscale)
        misc$earg = list(location= .elocation, scale= .escale)
    }), list( .llocation=llocation, .lscale=lscale,
              .elocation=elocation, .escale=escale ))),
    loglikelihood=eval(substitute(
        function(mu,y,w,residuals= FALSE,eta, extra=NULL) {
        location = eta2theta(eta[,1], .llocation, earg= .elocation)
        Scale = eta2theta(eta[,2], .lscale, earg= .escale)
        zedd = (y-location) / Scale
        if(residuals) stop("loglikelihood residuals not implemented yet") else
        sum(w * (-zedd - 2 * log1p(exp(-zedd)) - log(Scale )))
    }, list( .llocation=llocation, .lscale=lscale,
             .elocation=elocation, .escale=escale ))),
    vfamily=c("logistic2"),
    deriv=eval(substitute(expression({
        location = eta2theta(eta[,1], .llocation, earg= .elocation)
        Scale = eta2theta(eta[,2], .lscale, earg= .escale)
        zedd = (y-location) / Scale
        ezedd = exp(-zedd)
        dl.dlocation = (1-ezedd) / ((1 + ezedd) * Scale)
        dlocation.deta = dtheta.deta(location, .llocation, earg= .elocation)
        dl.dscale =  zedd * (1-ezedd) / ((1 + ezedd) * Scale) - 1/Scale
        dscale.deta = dtheta.deta(Scale, .lscale, earg= .escale)
        w * cbind(dl.dlocation * dlocation.deta,
                  dl.dscale * dscale.deta)
    }), list( .llocation=llocation, .lscale=lscale,
              .elocation=elocation, .escale=escale ))),
    weight=eval(substitute(expression({
        d2l.location2 = 1 / (3*Scale^2)
        d2l.dscale2 = (3 + pi^2) / (9*Scale^2)
        wz = matrix(as.numeric(NA), nrow=n, ncol=M) # diagonal
        wz[,iam(1,1,M)] = d2l.location2 * dlocation.deta^2
        wz[,iam(2,2,M)] = d2l.dscale2 * dscale.deta^2
        w * wz
    }), list( .llocation=llocation, .lscale=lscale,
              .elocation=elocation, .escale=escale ))))
}






alaplace = function(llocation="identity", lscale="loge",
                    lkappa="loge",
                    elocation=list(), escale=list(),
                    ekappa=list(),
                    ilocation=NULL, iscale=NULL, ikappa=1.0,
                    method.init=1, zero=NULL) {
    if(mode(llocation) != "character" && mode(llocation) != "name")
        llocation = as.character(substitute(llocation))
    if(mode(lscale) != "character" && mode(lscale) != "name")
        lscale = as.character(substitute(lscale))
    if(mode(lkappa) != "character" && mode(lkappa) != "name")
        lkappa = as.character(substitute(lkappa))
    if(!is.Numeric(method.init, allow=1, integ=TRUE, posit=TRUE) ||
       method.init > 2) stop("argument \"method.init\" must be 1 or 2")
    if(length(zero) && !is.Numeric(zero, integer=TRUE, posit=TRUE))
        stop("bad input for argument \"zero\"")
    if(!is.list(elocation)) elocation = list()
    if(!is.list(escale)) escale = list()
    if(!is.list(ekappa)) ekappa = list()

    new("vglmff",
    blurb=c("Three-parameter asymmetric Laplace distribution\n\n",
            "Links:    ",
            namesof("location", llocation, earg=elocation), ", ",
            namesof("scale", lscale, earg=escale), ", ",
            namesof("kappa", lkappa, earg=ekappa),
            "\n", "\n",
            "Mean:     location + scale * (1/kappa - kappa) / sqrt(2)",
            "\n",
            "Variance: mean^2 + scale^2"),
    constraints=eval(substitute(expression({
        constraints = cm.zero.vgam(constraints, x, .zero, M)
    }), list( .zero=zero ))),
    initialize=eval(substitute(expression({
        if(ncol(cbind(y)) != 1)
            stop("response must be a vector or a one-column matrix")
        predictors.names = 
        c(namesof("location", .llocation, earg=.elocation, tag=FALSE),
          namesof("scale",    .lscale,    earg=.escale,    tag=FALSE),
          namesof("kappa",    .lkappa,    earg=.ekappa,    tag=FALSE))
        if(!length(etastart)) {
            kappa.init = if(length( .ikappa)) rep( .ikappa, len=n) else
                         rep( 1.0, len=n)
            if( .method.init == 1) {
                location.init = median(y)
                scale.init = sqrt(var(y) / 2)
            } else {
                location.init = y
                scale.init = sqrt(sum(w*abs(y-median(y ))) / (sum(w) *2))
            }
            location.init = if(length(.ilocation)) rep(.ilocation, len=n) else
                             rep(location.init, len=n)
            scale.init = if(length(.iscale)) rep(.iscale, len=n) else
                             rep(scale.init, len=n)
            etastart =
                cbind(theta2eta(location.init, .llocation, earg= .elocation),
                      theta2eta(scale.init, .lscale, earg= .escale),
                      theta2eta(kappa.init, .lkappa, earg= .ekappa))
        }
    }), list( .method.init=method.init,
              .elocation=elocation, .escale=escale, .ekappa=ekappa,
              .llocation=llocation, .lscale=lscale, .lkappa=lkappa,
              .ilocation=ilocation, .iscale=iscale, .ikappa=ikappa ))),
    inverse=eval(substitute(function(eta, extra=NULL) {
        location = eta2theta(eta[,1], .llocation, earg= .elocation)
        Scale = eta2theta(eta[,2], .lscale, earg= .escale)
        kappa = eta2theta(eta[,3], .lkappa, earg= .ekappa)
        location + Scale * (1/kappa - kappa) / sqrt(2)
    }, list( .elocation=elocation, .llocation=llocation,
             .escale=ekappa, .lscale=lkappa,
             .ekappa=ekappa, .lkappa=lkappa ))),
    last=eval(substitute(expression({
        misc$link = c(location= .llocation, scale= .lscale, kappa= .lkappa)
        misc$earg = list(location= .elocation, scale= .escale, kappa= .ekappa)
        misc$expected = TRUE
    }), list( .elocation=elocation, .llocation=llocation,
              .escale=escale, .lscale=lscale,
              .ekappa=ekappa, .lkappa=lkappa ))),
    loglikelihood=eval(substitute(
        function(mu,y,w,residuals= FALSE,eta, extra=NULL) {
        location = eta2theta(eta[,1], .llocation, earg= .elocation)
        Scale = eta2theta(eta[,2], .lscale, earg= .escale)
        kappa = eta2theta(eta[,3], .lkappa, earg= .ekappa)
        zedd = ifelse(y >= location, kappa, 1/kappa) * sqrt(2) *
               abs(y-location) / Scale
        if(residuals) stop("loglikelihood residuals not implemented yet") else
        sum(w * (-zedd - log(Scale) - log(2)/2 + log(kappa) - log1p(kappa^2)))
    }, list( .elocation=elocation, .llocation=llocation,
             .escale=escale, .lscale=lscale,
             .ekappa=ekappa, .lkappa=lkappa ))),
    vfamily=c("alaplace"),
    deriv=eval(substitute(expression({
        location = eta2theta(eta[,1], .llocation, earg= .elocation)
        Scale = eta2theta(eta[,2], .lscale, earg= .escale)
        kappa = eta2theta(eta[,3], .lkappa, earg= .ekappa)
        zedd = abs(y-location) / Scale
        dl.dlocation = sqrt(2) * ifelse(y >= location, kappa, 1/kappa) *
                       sign(y-location) / Scale
        dl.dscale =  sqrt(2) * ifelse(y >= location, kappa, 1/kappa) *
                     zedd / Scale - 1 / Scale
        dl.dkappa =  1 / kappa - 2 * kappa / (1+kappa^2) -
                     (sqrt(2) / Scale) *
                     ifelse(y > location, 1, -1/kappa^2) * abs(y-location)  
        dlocation.deta = dtheta.deta(location, .llocation, earg= .elocation)
        dscale.deta = dtheta.deta(Scale, .lscale, earg= .escale)
        dkappa.deta = dtheta.deta(kappa, .lkappa, earg= .ekappa)
        w * cbind(dl.dlocation * dlocation.deta,
                  dl.dscale * dscale.deta,
                  dl.dkappa * dkappa.deta)
    }), list( .escale=escale, .lscale=lscale,
              .elocation=elocation, .llocation=llocation,
              .ekappa=ekappa, .lkappa=lkappa ))),
    weight=eval(substitute(expression({
        d2l.dlocation2 = 2 / Scale^2
        d2l.dscale2 = 1 / Scale^2
        d2l.dkappa2 = 1 / kappa^2 + 4 / (1+kappa^2)^2
        d2l.dkappadloc = -sqrt(8) / ((1+kappa^2) * Scale)
        d2l.dkappadscale = -(1-kappa^2) / ((1+kappa^2) * kappa * Scale)
        wz = matrix(0, nrow=n, dimm(M))
        wz[,iam(1,1,M)] = d2l.dlocation2 * dlocation.deta^2
        wz[,iam(2,2,M)] = d2l.dscale2 * dscale.deta^2
        wz[,iam(3,3,M)] = d2l.dkappa2 * dkappa.deta^2
        wz[,iam(1,3,M)] = d2l.dkappadloc * dkappa.deta * dlocation.deta
        wz[,iam(2,3,M)] = d2l.dkappadscale  * dkappa.deta * dscale.deta
        w * wz
    }), list( .escale=escale, .lscale=lscale,
              .elocation=elocation, .llocation=llocation ))))
}


laplace = function(llocation="identity", lscale="loge",
                   elocation=list(), escale=list(),
                   ilocation=NULL, iscale=NULL,
                   method.init=1, zero=NULL) {
    if(mode(llocation) != "character" && mode(llocation) != "name")
        llocation = as.character(substitute(llocation))
    if(mode(lscale) != "character" && mode(lscale) != "name")
        lscale = as.character(substitute(lscale))
    if(!is.Numeric(method.init, allow=1, integ=TRUE, posit=TRUE) ||
       method.init > 2) stop("argument \"method.init\" must be 1 or 2")
    if(length(zero) && !is.Numeric(zero, integer=TRUE, posit=TRUE))
        stop("bad input for argument \"zero\"")
    if(!is.list(elocation)) elocation = list()
    if(!is.list(escale)) escale = list()

    new("vglmff",
    blurb=c("Two-parameter Laplace distribution\n\n",
            "Links:    ",
            namesof("location", llocation, earg=elocation), ", ",
            namesof("scale", lscale, earg=escale),
            "\n", "\n",
            "Mean:     location", "\n",
            "Variance: 2*scale^2"),
    constraints=eval(substitute(expression({
        constraints = cm.zero.vgam(constraints, x, .zero, M)
    }), list( .zero=zero ))),
    initialize=eval(substitute(expression({
        if(ncol(cbind(y)) != 1)
            stop("response must be a vector or a one-column matrix")
        predictors.names = 
        c(namesof("location", .llocation, earg=.elocation, tag=FALSE),
          namesof("scale",    .lscale,    earg=.escale,    tag=FALSE))
        if(!length(etastart)) {
            if( .method.init == 1) {
                location.init = median(y)
                scale.init = sqrt(var(y) / 2)
            } else {
                location.init = y
                scale.init = sqrt(sum(w*abs(y-median(y ))) / (sum(w) *2))
            }
            location.init = if(length(.ilocation)) rep(.ilocation, len=n) else
                             rep(location.init, len=n)
            scale.init = if(length(.iscale)) rep(.iscale, len=n) else
                             rep(scale.init, len=n)
            etastart =
                cbind(theta2eta(location.init, .llocation, earg= .elocation),
                      theta2eta(scale.init, .lscale, earg= .escale))
        }
    }), list( .method.init=method.init,
             .elocation=elocation, .escale=escale,
             .llocation=llocation, .lscale=lscale,
             .ilocation=ilocation, .iscale=iscale ))),
    inverse=eval(substitute(function(eta, extra=NULL) {
        eta2theta(eta[,1], .llocation, earg= .elocation)
    }, list( .elocation=elocation, .llocation=llocation ))),
    last=eval(substitute(expression({
        misc$link = c(location= .llocation, scale= .lscale)
        misc$earg = list(location= .elocation, scale= .escale)
        misc$expected = TRUE
        misc$RegCondOK = FALSE # Save this for later
    }), list( .escale=escale, .lscale=lscale,
              .elocation=elocation, .llocation=llocation ))),
    loglikelihood=eval(substitute(
        function(mu,y,w,residuals= FALSE,eta, extra=NULL) {
        location = eta2theta(eta[,1], .llocation, earg= .elocation)
        Scale = eta2theta(eta[,2], .lscale, earg= .escale)
        zedd = abs(y-location) / Scale
        if(residuals) stop("loglikelihood residuals not implemented yet") else
        sum(w * (-zedd - log(Scale) - log(2)))
    }, list( .escale=escale, .lscale=lscale,
             .elocation=elocation, .llocation=llocation ))),
    vfamily=c("laplace"),
    deriv=eval(substitute(expression({
        location = eta2theta(eta[,1], .llocation, earg= .elocation)
        Scale = eta2theta(eta[,2], .lscale, earg= .escale)
        zedd = abs(y-location) / Scale
        dl.dlocation = sign(y-location) / Scale
        dl.dscale =  zedd / Scale - 1/Scale
        dlocation.deta = dtheta.deta(location, .llocation, earg= .elocation)
        dscale.deta = dtheta.deta(Scale, .lscale, earg= .escale)
        w * cbind(dl.dlocation * dlocation.deta, dl.dscale * dscale.deta)
    }), list( .escale=escale, .lscale=lscale,
              .elocation=elocation, .llocation=llocation ))),
    weight=eval(substitute(expression({
        d2l.dlocation2 = 1 / Scale^2
        d2l.dscale2 = 1 / Scale^2
        wz = matrix(0, nrow=n, ncol=2) # diagonal
        wz[,iam(1,1,M)] = d2l.dlocation2 * dlocation.deta^2
        wz[,iam(2,2,M)] = d2l.dscale2 * dscale.deta^2
        w * wz
    }), list( .escale=escale, .lscale=lscale,
              .elocation=elocation, .llocation=llocation ))))
}

dlaplace = function(x, location=0, scale=1) {
    if(!is.Numeric(scale, posit=TRUE)) 
        stop("argument \"scale\" must be positive")
    exp(-abs(x-location)/scale) / (2*scale)
}

plaplace = function(q, location=0, scale=1) {
    if(!is.Numeric(scale, posit=TRUE)) 
        stop("argument \"scale\" must be positive")
    zedd = (q-location) / scale
    L = max(length(q), length(location), length(scale))
    q = rep(q, len=L); location = rep(location, len=L); scale= rep(scale, len=L)
    ifelse(q < location, 0.5*exp(zedd), 1-0.5*exp(-zedd))
}

qlaplace = function(p, location=0, scale=1) {
    if(!is.Numeric(scale, posit=TRUE)) 
        stop("argument \"scale\" must be positive")
    L = max(length(p), length(location), length(scale))
    p = rep(p, len=L); location = rep(location, len=L); scale= rep(scale, len=L)
    location - sign(p-0.5) * scale * log(2*ifelse(p < 0.5, p, 1-p))
}

rlaplace = function(n, location=0, scale=1) {
    if(!is.Numeric(n, posit=TRUE, integ=TRUE, allow=1)) 
        stop("bad input for argument \"n\"")
    if(!is.Numeric(scale, posit=TRUE)) stop("\"scale\" must be positive")
    location = rep(location, len=n); scale= rep(scale, len=n)
    r = runif(n)
    location - sign(r-0.5) * scale * log(2*ifelse(r < 0.5, r, 1-r))
}




fff.control <- function(save.weight=TRUE, ...)
{
    list(save.weight=save.weight)
}

fff = function(link="loge", earg=list(),
               idf1=NULL, idf2=NULL,
               method.init=1, zero=NULL) {
    if(mode(link) != "character" && mode(link) != "name")
        link = as.character(substitute(link))
    if(!is.Numeric(method.init, allow=1, integ=TRUE, posit=TRUE) ||
       method.init > 2) stop("argument \"method.init\" must be 1 or 2")
    if(length(zero) && !is.Numeric(zero, integer=TRUE, posit=TRUE))
        stop("bad input for argument \"zero\"")
    if(!is.list(earg)) earg = list()

    new("vglmff",
    blurb=c("F-distribution\n\n",
            "Links:    ",
            namesof("df1", link, earg=earg), ", ",
            namesof("df2", link, earg=earg),
            "\n", "\n",
            "Mean:     df2/(df2-2) provided df2>2", "\n",
      "Variance: 2*df2^2*(df1+df2-2)/(df1*(df2-2)^2*(df2-4)) provided df2>4"),
    constraints=eval(substitute(expression({
        constraints = cm.zero.vgam(constraints, x, .zero, M)
    }), list( .zero=zero ))),
    initialize=eval(substitute(expression({
        if(ncol(cbind(y)) != 1)
            stop("response must be a vector or a one-column matrix")
        predictors.names = c(namesof("df1", .link, earg=.earg, tag=FALSE),
                             namesof("df2", .link, earg=.earg, tag=FALSE))
        if(!length(etastart)) {
            if( .method.init == 1) {
                df2.init = b = 2*mean(y) / (mean(y)-1)
                df1.init = 2*b^2*(b-2)/(var(y)*(b-2)^2 * (b-4) - 2*b^2)
                if(df2.init > 4) df2.init = 5
                if(df1.init > 2) df1.init = 3
            } else {
                df2.init = b = 2*median(y) / (median(y)-1)
                summy = summary(y)
                var.est = summy[5] - summy[2]
                df1.init = 2*b^2*(b-2)/(var.est*(b-2)^2 * (b-4) - 2*b^2)
            }
            df1.init = if(length(.idf1)) rep(.idf1, len=n) else
                           rep(df1.init, len=n)
            df2.init = if(length(.idf2)) rep(.idf2, len=n) else rep(1, len=n)
            etastart = cbind(theta2eta(df1.init, .link, earg= .earg),
                             theta2eta(df2.init, .link, earg= .earg))
        }
    }), list( .method.init=method.init, .idf1=idf1, .earg=earg,
             .idf2=idf2, .link=link ))),
    inverse=eval(substitute(function(eta, extra=NULL) {
        df2 = eta2theta(eta[,2], .link, earg= .earg)
        ans = df2 * NA
        ans[df2>2] = df2[df2>2] / (df2[df2>2]-2)
        ans
    }, list( .link=link, .earg=earg ))),
    last=eval(substitute(expression({
        misc$link = c(df1= .link, df2= .link)
        misc$earg = list(df1= .earg, df2= .earg)
    }), list( .link=link, .earg=earg ))),
    loglikelihood=eval(substitute(
        function(mu,y,w,residuals= FALSE,eta, extra=NULL) {
        df1 = eta2theta(eta[,1], .link, earg= .earg)
        df2 = eta2theta(eta[,2], .link, earg= .earg)
        if(residuals) stop("loglikelihood residuals not implemented yet") else
        sum(w * (lgamma(0.5*(df1+df2)) + 0.5*df1*log(df1/df2) +
                 0.5*(df1-2) * log(y) - lgamma(df1/2) - lgamma(df2/2) -
                 0.5*(df1+df2)*log1p(df1*y/df2)))
    }, list( .link=link, .earg=earg ))),
    vfamily=c("fff"),
    deriv=eval(substitute(expression({
        df1 = eta2theta(eta[,1], .link, earg= .earg)
        df2 = eta2theta(eta[,2], .link, earg= .earg)
        dl.ddf1 = 0.5*digamma(0.5*(df1+df2)) + 0.5 + 0.5*log(df1/df2) +
                  0.5*log(y) - 0.5*digamma(0.5*df1) -
                  0.5*(df1+df2)*(y/df2) / (1 + df1*y/df2) -
                  0.5*log1p(df1*y/df2)
        ddf1.deta = dtheta.deta(df1, .link, earg= .earg)
        dl.ddf2 = 0.5*digamma(0.5*(df1+df2)) - 0.5*df1/df2 - 
                  0.5*digamma(0.5*df2) -
                  0.5*(df1+df2) * (-df1*y/df2^2) / (1 + df1*y/df2) -
                  0.5*log1p(df1*y/df2)
        ddf2.deta = dtheta.deta(df2, .link, earg= .earg)
        if(iter == 1) {
            etanew = eta
        } else {
            derivold = derivnew
            etaold = etanew
            etanew = eta
        }
        derivnew = w * cbind(dl.ddf1 * ddf1.deta,
                             dl.ddf2 * ddf2.deta)
        derivnew
    }), list( .link=link, .earg=earg ))),
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
    }), list( .link=link, .earg=earg ))))
}



vonmises = function(llocation="elogit",
                    lscale="loge",
      elocation=if(llocation=="elogit") list(min=0, max=2*pi) else list(),
      escale=list(),
                    ilocation=NULL, iscale=NULL,
                    method.init=1, zero=NULL) {
    if(mode(llocation) != "character" && mode(llocation) != "name")
        llocation = as.character(substitute(llocation))
    if(mode(lscale) != "character" && mode(lscale) != "name")
        lscale = as.character(substitute(lscale))
    if(!is.Numeric(method.init, allow=1, integ=TRUE, posit=TRUE) ||
       method.init > 2) stop("argument \"method.init\" must be 1 or 2")
    if(length(zero) && !is.Numeric(zero, integer=TRUE, posit=TRUE))
        stop("bad input for argument \"zero\"")
    if(!is.list(escale)) escale = list()

    new("vglmff",
    blurb=c("Von Mises distribution\n\n",
            "Links:    ",
            namesof("location", llocation, earg= elocation), ", ",
            namesof("scale", lscale, earg=escale),
            "\n", "\n",
            "Mean:     location"),
    constraints=eval(substitute(expression({
        constraints = cm.zero.vgam(constraints, x, .zero, M)
    }), list( .zero=zero ))),
    initialize=eval(substitute(expression({
        if(ncol(cbind(y)) != 1)
            stop("response must be a vector or a one-column matrix")
        predictors.names = 
        c(namesof("location", .llocation, earg= .elocation, tag=FALSE),
          namesof("scale", .lscale, earg=.escale, tag=FALSE))
        if(!length(etastart)) {
            if( .method.init == 1) {
                location.init = mean(y)
                rat10 = sqrt((sum(w*cos(y )))^2 + sum(w*sin(y))^2) / sum(w)
                scale.init = sqrt(1 - rat10)
            } else {
                location.init = median(y)
                scale.init = sqrt(sum(w*abs(y - location.init)) / sum(w))
            }
            location.init = if(length(.ilocation)) rep(.ilocation, len=n) else
                           rep(location.init, len=n)
            scale.init= if(length(.iscale)) rep(.iscale,len=n) else rep(1,len=n)
            etastart = cbind(
                theta2eta(location.init, .llocation, earg= .elocation),
                theta2eta(scale.init, .lscale, earg= .escale))
        }
        y = y %% (2*pi) # Coerce after initial values have been computed
    }), list( .method.init=method.init, .ilocation=ilocation,
              .escale=escale, .iscale=iscale,
              .lscale=lscale, .llocation=llocation, .elocation=elocation ))),
    inverse=eval(substitute(function(eta, extra=NULL) {
        eta2theta(eta[,1], .llocation, earg= .elocation) %% (2*pi)
    }, list( .escale=escale, .lscale=lscale,
             .llocation=llocation, .elocation=elocation ))),
    last=eval(substitute(expression({
        misc$link = c(location= .llocation, scale= .lscale)
        misc$earg = list(location= .elocation, scale= .escale )
    }), list( .escale=escale, .lscale=lscale,
              .llocation=llocation, .elocation=elocation ))),
    loglikelihood=eval(substitute(
        function(mu,y,w,residuals= FALSE,eta, extra=NULL) {
        location = eta2theta(eta[,1], .llocation, earg= .elocation)
        Scale = eta2theta(eta[,2], .lscale, earg= .escale)
        if(residuals) stop("loglikelihood residuals not implemented yet") else
        sum(w * (Scale * cos(y - location) -
                 log(mbesselI0(x=Scale ))))
    }, list( .escale=escale, .lscale=lscale,
             .llocation=llocation, .elocation=elocation ))),
    vfamily=c("vonmises"),
    deriv=eval(substitute(expression({
        location = eta2theta(eta[,1], .llocation, earg= .elocation)
        Scale = eta2theta(eta[,2], .lscale, earg= .escale)
        tmp6 = mbesselI0(x=Scale, deriv=2)
        dl.dlocation = Scale * sin(y - location)
        dlocation.deta = dtheta.deta(location, .llocation, earg= .elocation)
        dl.dscale = cos(y - location) - tmp6[,2] / tmp6[,1]
        dscale.deta = dtheta.deta(Scale, .lscale, earg= .escale)
        w * cbind(dl.dlocation * dlocation.deta,
                  dl.dscale * dscale.deta)
    }), list( .escale=escale, .lscale=lscale,
              .llocation=llocation, .elocation=elocation ))),
    weight=eval(substitute(expression({
        d2l.location2 = Scale * tmp6[,2] / tmp6[,1]
        d2l.dscale2 = tmp6[,3] / tmp6[,1] - (tmp6[,2] / tmp6[,1])^2
        wz = matrix(as.numeric(NA), nrow=n, ncol=2) # diagonal
        wz[,iam(1,1,M)] = d2l.location2 * dlocation.deta^2
        wz[,iam(2,2,M)] = d2l.dscale2 * dscale.deta^2
        w * wz
    }), list( .escale=escale, .lscale=lscale,
              .llocation=llocation, .elocation=elocation ))))
}



hyperg = function(N=NULL, D=NULL,
                 lprob="logit", earg=list(),
                 iprob=NULL) {
    if(mode(lprob) != "character" && mode(lprob) != "name")
        lprob = as.character(substitute(lprob))
    inputN = is.Numeric(N, positive=TRUE)
    inputD = is.Numeric(D, positive=TRUE)
    if(inputD && inputN)
        stop("only one of \"N\" and \"D\" is to be inputted")
    if(!inputD && !inputN)
        stop("one of \"N\" and \"D\" needs to be inputted")
    if(!is.list(earg)) earg = list()

    new("vglmff",
    blurb=c("Hypergeometric distribution\n\n",
            "Link:     ",
            namesof("prob", lprob, earg=earg), "\n",
            "Mean:     D/N\n"),
    initialize=eval(substitute(expression({
            NCOL = function (x)
                if(is.array(x) && length(dim(x)) > 1 ||
                is.data.frame(x)) ncol(x) else as.integer(1)
            if(NCOL(y) == 1) {
                if(is.factor(y)) y = y != levels(y)[1]
                nn = rep(1, len=n)
                if(!all(y >= 0 & y <= 1))
                    stop("response values must be in [0, 1]")
                mustart = (0.5 + w * y) / (1 + w)
                no.successes = w * y
                if(any(abs(no.successes - round(no.successes)) > 0.001))
                    stop("Number of successes must be integer-valued")
            } else if(NCOL(y) == 2) {
                if(any(abs(y - round(y)) > 0.001))
                    stop("Count data must be integer-valued")
                nn = y[,1] + y[,2]
                y = ifelse(nn > 0, y[,1]/nn, 0)
                w = w * nn
                mustart = (0.5 + nn * y) / (1 + nn)
            } else
                 stop("Response not of the right form")

        predictors.names = namesof("prob", .lprob, earg=.earg, tag=FALSE)
        extra$Nvector = .N
        extra$Dvector = .D
        extra$Nunknown = length(extra$Nvector) == 0
        if(!length(etastart)) {
            init.prob = if(length( .iprob)) rep( .iprob, len=n) else mustart
            etastart = matrix(init.prob, n, ncol(cbind(y )))
        }
    }), list( .lprob=lprob, .earg=earg, .N=N, .D=D, .iprob=iprob ))), 
    inverse=eval(substitute(function(eta, extra=NULL) {
        eta2theta(eta, .lprob, earg= .earg)
    }, list( .lprob=lprob, .earg=earg ))),
    last=eval(substitute(expression({
        misc$link = c("prob"= .lprob) 
        misc$earg = list("prob"= .earg) 
        misc$Dvector = .D
        misc$Nvector = .N
    }), list( .N=N, .D=D, .lprob=lprob, .earg=earg ))),
    link=eval(substitute(function(mu, extra=NULL) {
        theta2eta(mu, .lprob, earg= .earg)
    }, list( .lprob=lprob, .earg=earg ))),
    loglikelihood=eval(substitute(
        function(mu, y, w, residuals = FALSE, eta, extra=NULL) {
        N = extra$Nvector
        Dvec = extra$Dvector
        prob = mu
        yvec = w * y
        if(residuals) stop("loglikelihood residuals not implemented yet") else {
        if(extra$Nunknown) {
            tmp12 = Dvec * (1-prob) / prob
            sum(lgamma(1+tmp12) + lgamma(1+Dvec/prob-w) -
                   lgamma(1+tmp12-w+yvec) - lgamma(1+Dvec/prob))
        } else
            sum(lgamma(1+N*prob) + lgamma(1+N*(1-prob)) -
                   lgamma(1+N*prob-yvec) - lgamma(1+N*(1-prob) -w + yvec))
        }
    }, list( .lprob=lprob, .earg=earg ))), 
    vfamily=c("hyperg"),
    deriv=eval(substitute(expression({
        prob = mu   # equivalently, eta2theta(eta, .lprob, earg= .earg)
        dprob.deta = dtheta.deta(prob, .lprob, earg= .earg)
        Dvec = extra$Dvector
        Nvec = extra$Nvector
        yvec = w * y
        if(extra$Nunknown) {
            tmp72 = -Dvec / prob^2
            tmp12 =  Dvec * (1-prob) / prob
            dl.dprob = tmp72 * (digamma(1 + tmp12) + digamma(1 + Dvec/prob -w) -
                       digamma(1 + tmp12-w+yvec) - digamma(1 + Dvec/prob))
        } else {
            dl.dprob = Nvec * (digamma(1+Nvec*prob) - digamma(1+Nvec*(1-prob)) -
                digamma(1+Nvec*prob-yvec) + digamma(1+Nvec*(1-prob)-w+yvec))
        }
        w * dl.dprob * dprob.deta
    }), list( .lprob=lprob, .earg=earg ))),
    weight=eval(substitute(expression({
        if(extra$Nunknown) {
            tmp722 = tmp72^2
            tmp13  = 2*Dvec / prob^3
            d2l.dprob2 = tmp722 * (trigamma(1 + tmp12) + 
                         trigamma(1 + Dvec/prob - w) -
                         trigamma(1 + tmp12 - w + yvec) -
                         trigamma(1 + Dvec/prob)) +
                         tmp13 * (digamma(1 + tmp12) +
                         digamma(1 + Dvec/prob - w) -
                         digamma(1 + tmp12 - w + yvec) -
                         digamma(1 + Dvec/prob))
        } else {
            d2l.dprob2 = Nvec^2 * (trigamma(1+Nvec*prob) +
                         trigamma(1+Nvec*(1-prob)) -
                         trigamma(1+Nvec*prob-yvec) -
                         trigamma(1+Nvec*(1-prob)-w+yvec))
        }
        d2prob.deta2 = d2theta.deta2(prob, .lprob, earg= .earg)
        wz = -(dprob.deta^2) * d2l.dprob2 - d2prob.deta2 * dl.dprob
        wz = w * wz
        wz[wz < .Machine$double.eps] = .Machine$double.eps
        wz
    }), list( .lprob=lprob, .earg=earg ))))
}



dbenini = function(x, shape, y0) {
    if(!is.Numeric(x)) stop("bad input for argument \"x\"")
    if(!is.Numeric(shape, posit=TRUE)) stop("bad input for argument \"shape\"")
    if(!is.Numeric(y0, posit=TRUE)) stop("bad input for argument \"y0\"")
    N = max(length(x), length(shape), length(y0))
    x = rep(x, len=N); shape = rep(shape, len=N); y0 = rep(y0, len=N); 
    ok = x > y0
    temp = log(x[ok]/y0[ok])
    ans = y0 * 0
    ans[ok] = 2 * shape[ok] * exp(-shape[ok] * temp^2) * temp / x[ok]
    ans
}

pbenini = function(q, shape, y0) {
    if(!is.Numeric(q)) stop("bad input for argument \"q\"")
    if(!is.Numeric(shape, posit=TRUE)) stop("bad input for argument \"shape\"")
    if(!is.Numeric(y0, posit=TRUE)) stop("bad input for argument \"y0\"")
    N = max(length(q), length(shape), length(y0))
    q = rep(q, len=N); shape = rep(shape, len=N); y0 = rep(y0, len=N); 
    ans = y0 * 0
    ok = q > y0
    ans[ok] = 1 - exp(-shape[ok] * (log(q[ok]/y0[ok]))^2)
    ans
}

qbenini = function(p, shape, y0) {
    if(!is.Numeric(p, posit=TRUE) || any(p >= 1)) 
        stop("bad input for argument \"p\"")
    if(!is.Numeric(shape, posit=TRUE)) stop("bad input for argument \"shape\"")
    if(!is.Numeric(y0, posit=TRUE)) stop("bad input for argument \"y0\"")
    y0 * exp(sqrt(-log1p(-p) / shape))
}

rbenini = function(n, shape, y0) {
    if(!is.Numeric(n, posit=TRUE, integ=TRUE, allow=1)) 
        stop("bad input for argument \"n\"")
    if(!is.Numeric(shape, posit=TRUE)) stop("bad input for argument \"shape\"")
    if(!is.Numeric(y0, posit=TRUE)) stop("bad input for argument \"y0\"")
    y0 * exp(sqrt(-log(runif(n)) / shape))
}

benini = function(y0=stop("argument \"y0\" must be specified"),
                  lshape="loge", earg=list(),
                  ishape=NULL, method.init=1) {
    if(mode(lshape) != "character" && mode(lshape) != "name")
        lshape = as.character(substitute(lshape))
    if(!is.Numeric(method.init, allow=1, integ=TRUE, posit=TRUE) ||
       method.init > 2) stop("argument \"method.init\" must be 1 or 2")
    if(!is.Numeric(y0, allow=1, posit=TRUE))
       stop("bad input for argument \"y0\"")
    if(!is.list(earg)) earg = list()

    new("vglmff",
    blurb=c("1-parameter Benini distribution\n\n",
            "Link:    ",
            namesof("shape", lshape, earg=earg),
            "\n", "\n"),
    initialize=eval(substitute(expression({
        if(ncol(cbind(y)) != 1)
            stop("response must be a vector or a one-column matrix")
        predictors.names = c(namesof("shape", .lshape, earg=.earg, tag=FALSE))
        extra$y0 = .y0
        if(min(y) <= extra$y0) stop("argument \"y0\" is too large")
        if(!length(etastart)) {
            probs = (1:3) / 4
            qofy= quantile(rep(y, times=w), probs=probs) # fails if w != integer
            if( .method.init == 1) {
                shape.init = mean(-log1p(-probs) / (log(qofy))^2)
            } else {
                shape.init = median(-log1p(-probs) / (log(qofy))^2)
            }
            shape.init = if(length(.ishape)) rep(.ishape, len=n) else
                         rep(shape.init, len=n)
            etastart = cbind(theta2eta(shape.init, .lshape, earg= .earg))
        }
    }), list( .method.init=method.init, .ishape=ishape, .lshape=lshape, .earg=earg,
             .y0=y0 ))),
    inverse=eval(substitute(function(eta, extra=NULL) {
        shape = eta2theta(eta, .lshape, earg= .earg)
        temp = 1/(4*shape)
        extra$y0 * exp(temp) *
        ((sqrt(pi) * (1 - pgamma(temp, 0.5 ))) / (2*sqrt(shape)) +
                     1 - pgamma(temp, 1))
    }, list( .lshape=lshape, .earg=earg ))),
    last=eval(substitute(expression({
        misc$link = c(shape= .lshape)
        misc$earg = list(shape= .earg )
    }), list( .lshape=lshape, .earg=earg ))),
    loglikelihood=eval(substitute(
        function(mu,y,w,residuals= FALSE,eta, extra=NULL) {
        shape = eta2theta(eta, .lshape, earg= .earg)
        y0 = extra$y0
        if(residuals) stop("loglikelihood residuals not implemented yet") else
        sum(w * (log(shape) - log(y) - shape*(log(y/y0))^2 + log(log(y/y0 ))))
    }, list( .lshape=lshape, .earg=earg ))),
    vfamily=c("benini"),
    deriv=eval(substitute(expression({
        shape = eta2theta(eta, .lshape, earg= .earg)
        y0 = extra$y0
        dl.dshape = 1/shape - (log(y/y0))^2
        dshape.deta = dtheta.deta(shape, .lshape, earg= .earg)
        w * dl.dshape * dshape.deta
    }), list( .lshape=lshape, .earg=earg ))),
    weight=eval(substitute(expression({
        d2l.dshape2 = 1 / shape^2
        wz = d2l.dshape2 * dshape.deta^2
        w * wz
    }), list( .lshape=lshape, .earg=earg ))))
}




dpolono = function(x, meanlog=0, sdlog=1, ...) {
    if(!is.Numeric(x)) stop("bad input for argument \"x\"")
    if(!is.Numeric(meanlog)) stop("bad input for argument \"meanlog\"")
    if(!is.Numeric(sdlog, posit=TRUE)) stop("bad input for argument \"sdlog\"")
    N = max(length(x), length(meanlog), length(sdlog))
    x = rep(x, len=N); meanlog = rep(meanlog, len=N); sdlog = rep(sdlog, len=N);
    ans = x * 0
    integrand = function(t, x, meanlog, sdlog)
        exp(t*x - exp(t) - 0.5*((t-meanlog)/sdlog)^2)
    for(i in 1:N) {
        if(x[i] == round(x[i]) && x[i] >= 0) {
            temp = integrate(f=integrand, lower=-Inf, upper=Inf,
                             x=x[i], meanlog=meanlog[i], sdlog=sdlog[i], ...)
            if(temp$message == "OK") ans[i] = temp$value else {
            warning(paste("could not integrate (numerically) observation", i))
                ans[i] = NA
            }
        }
    }
    ans = ans / (sqrt(2*pi) * sdlog * gamma(x+1))
    ifelse(x == round(x) & x >= 0, ans, 0)
}


rpolono = function(n, meanlog=0, sdlog=1) {
    if(!is.Numeric(n, integ=TRUE,allow=1)) stop("bad input for argument \"n\"")
    if(!is.Numeric(meanlog)) stop("bad input for argument \"meanlog\"")
    if(!is.Numeric(sdlog)) stop("bad input for argument \"sdlog\"")
    meanlog = rep(meanlog, len=n); sdlog = rep(sdlog, len=n);
    lambda = if(is.R()) rlnorm(n=n, meanlog=meanlog, sdlog=sdlog) else
             stop("suppressing a warning message")
    rpois(n=n, lambda=lambda)
}









dtriangle = function(x, theta, lower=0, upper=1) {
    if(!is.Numeric(x)) stop("bad input for argument \"x\"")
    if(!is.Numeric(theta)) stop("bad input for argument \"theta\"")
    if(!is.Numeric(lower)) stop("bad input for argument \"lower\"")
    if(!is.Numeric(upper)) stop("bad input for argument \"upper\"")
    if(!all(lower < theta & theta < upper))
        stop("lower < theta < upper values are required")
    N = max(length(x), length(theta), length(lower), length(upper))
    x = rep(x, len=N); lower = rep(lower, len=N); upper = rep(upper, len=N);
    theta = rep(theta, len=N)
    ans = x * 0
    neg = (lower <= x) & (x <= theta)
    pos = (theta <= x) & (x <= upper)
    denom1 = ((upper-lower)*(theta-lower))
    denom2 = ((upper-lower)*(upper-theta))
    ans[neg] = pmax(2 * (x-lower) / denom1, 0)[neg]
    ans[pos] = pmax(2 * (upper-x) / denom2, 0)[pos]
    ans
}


rtriangle = function(n, theta, lower=0, upper=1) {
    if(!is.Numeric(n, integ=TRUE,allow=1)) stop("bad input for argument \"n\"")
    if(!is.Numeric(theta)) stop("bad input for argument \"theta\"")
    if(!is.Numeric(lower)) stop("bad input for argument \"lower\"")
    if(!is.Numeric(upper)) stop("bad input for argument \"upper\"")
    if(!all(lower < theta & theta < upper))
        stop("lower < theta < upper values are required")
    N = n
    lower = rep(lower, len=N); upper = rep(upper, len=N);
    theta = rep(theta, len=N)
    t1 = sqrt(runif(n))
    t2 = sqrt(runif(n))
    ifelse(runif(n) < (theta-lower)/(upper-lower),
           lower + (theta-lower)*t1,
           upper - (upper-theta)*t2)
}


qtriangle = function(p, theta, lower=0, upper=1) {
    if(!is.Numeric(p, posit=TRUE)) stop("bad input for argument \"p\"")
    if(!is.Numeric(theta)) stop("bad input for argument \"theta\"")
    if(!is.Numeric(lower)) stop("bad input for argument \"lower\"")
    if(!is.Numeric(upper)) stop("bad input for argument \"upper\"")
    if(!all(lower < theta & theta < upper))
        stop("lower < theta < upper values are required")

    N = max(length(p), length(theta), length(lower), length(upper))
    p = rep(p, len=N); lower = rep(lower, len=N); upper = rep(upper, len=N);
    theta = rep(theta, len=N)

    bad = (p < 0) | (p > 1)
    if(any(bad))
        stop("bad input for 'p'")

    Neg = (p <= (theta - lower)/(upper - lower))
    ans = as.numeric(NA) * p
    temp1 = p * (upper-lower) * (theta-lower)
    ans[ Neg] = lower[ Neg] + sqrt(temp1[ Neg])

    Pos = (p >= (theta - lower)/(upper - lower))
    if(any(Pos)) {
        pstar = (p - (theta-lower)/(upper-lower)) / (1 -
                (theta-lower)/(upper-lower))
        qstar = cbind(1 - sqrt(1-pstar), 1 + sqrt(1-pstar))
        qstar = qstar[Pos,]
        qstar = ifelse(qstar[,1] >= 0 & qstar[,1] <= 1, qstar[,1], qstar[,2])
        ans[Pos] = theta[Pos] + qstar * (upper-theta)[Pos]
    }
    ans
}


ptriangle = function(q, theta, lower=0, upper=1) {
    if(!is.Numeric(q)) stop("bad input for argument \"q\"")
    if(!is.Numeric(theta)) stop("bad input for argument \"theta\"")
    if(!is.Numeric(lower)) stop("bad input for argument \"lower\"")
    if(!is.Numeric(upper)) stop("bad input for argument \"upper\"")
    if(!all(lower < theta & theta < upper))
        stop("lower < theta < upper values are required")

    N = max(length(q), length(theta), length(lower), length(upper))
    q = rep(q, len=N); lower = rep(lower, len=N); upper = rep(upper, len=N);
    theta = rep(theta, len=N)
    ans = q * 0

    qstar = (q - lower)^2 / ((upper-lower) * (theta-lower))
    Neg = (lower <= q & q <= theta)
    ans[Neg] = (qstar)[Neg]

    Pos = (theta <= q & q <= upper)
    qstar = (q - theta) / (upper-theta)
    ans[Pos] = ((theta-lower)/(upper-lower))[Pos] +
               (qstar * (2-qstar) * (upper-theta) / (upper - lower))[Pos]
    ans[q >= upper] = 1
    ans
}



triangle = function(lower=0, upper=1,
                    link="elogit", earg=if(link=="elogit") 
                    list(min = lower, max = upper) else list(), itheta=NULL)
{
    if(!is.Numeric(lower)) stop("bad input for argument \"lower\"")
    if(!is.Numeric(upper)) stop("bad input for argument \"upper\"")
    if(!all(lower < upper))
        stop("lower < upper values are required")
    if(length(itheta) && !is.Numeric(itheta))
        stop("bad input for 'itheta'")

    if(mode(link) != "character" && mode(link) != "name")
        link = as.character(substitute(link))
    if(!is.list(earg)) earg = list()

    new("vglmff",
    blurb=c(
    "Triangle distribution\n\n",
            "Link:    ",
            namesof("theta", link, earg=earg)),
    initialize=eval(substitute(expression({
        y = as.numeric(y)
        if(ncol(cbind(y)) != 1)
            stop("response must be a vector or a one-column matrix")
        extra$lower = rep( .lower, len=n)
        extra$upper = rep( .upper, len=n)

        if(any(y <= extra$lower | y >= extra$upper))
            stop("some y values in [lower,upper] detected")
        predictors.names = namesof("theta", .link, earg= .earg, tag=FALSE)
        if(!length(etastart)) {
            Theta.init = if(length( .itheta)) .itheta else {
                weighted.mean(y, w)
            }
            Theta.init = rep(Theta.init, length=n)
            etastart = theta2eta(Theta.init, .link, earg= .earg )
        }
    }), list( .link=link, .earg=earg, .itheta=itheta,
              .upper=upper, .lower=lower ))),
    inverse=eval(substitute(function(eta, extra=NULL) {
        Theta = eta2theta(eta, .link, earg= .earg )
        lower = extra$lower
        upper = extra$upper
        mu =  ((Theta^3 / 3 - lower * Theta^2 / 2 +
              lower^3 / 6) / (Theta - lower) + 
               ((Theta^3 / 3 - upper * Theta^2 / 2 +
              upper^3 / 6) / (upper - Theta))) * 2  / (upper-lower)
        mu
    }, list( .link=link, .earg=earg ))),
    last=eval(substitute(expression({
        misc$link = c(theta = .link)
        misc$earg = list(theta = .earg)
        misc$expected = TRUE
    }), list( .link=link, .earg=earg ))),
    loglikelihood=eval(substitute(
        function(mu,y,w,residuals= FALSE,eta, extra=NULL) {
        Theta = eta2theta(eta, .link, earg= .earg )
        lower = extra$lower
        upper = extra$upper
        if(residuals) stop("loglikelihood residuals not implemented yet") else {
            pos = y >= Theta
            neg = y <  Theta
            sum(w * (log(2) - log(upper-lower))) +
            sum(w[neg]*(log(y[neg]-lower[neg]) - log(Theta[neg]-lower[neg]))) +
            sum(w[pos]*(log(upper[pos]-y[pos]) - log(upper[pos]-Theta[pos])))
        }
    }, list( .link=link, .earg=earg ))),
    vfamily=c("triangle"),
    deriv=eval(substitute(expression({
        Theta = eta2theta(eta, .link, earg= .earg ) 
        dTheta.deta = dtheta.deta(Theta, .link, earg= .earg )
        pos = y > Theta
        neg = y < Theta
        lower = extra$lower
        upper = extra$upper
        dl.dTheta =  0 * y
        dl.dTheta[neg] =  -1 / (Theta[neg]-lower[neg])
        dl.dTheta[pos] =   1 / (upper[pos]-Theta[pos])
        dl.dTheta * dTheta.deta
    }), list( .link=link, .earg=earg ))),
    weight=eval(substitute(expression({
        d2l.dTheta2 =  1 / ((Theta-lower)*(upper-Theta))
        wz = dTheta.deta^2 * d2l.dTheta2
        w * wz
    }), list( .link=link, .earg=earg ))))
}




# These functions are
# Copyright (C) 1998-2009 T.W. Yee, University of Auckland. All rights reserved.











lms.bcn.control <-
lms.bcg.control <-
lms.yjn.control <- function(trace=TRUE, ...)
   list(trace=trace) 





lms.bcn <- function(percentiles=c(25,50,75),
                    zero=c(1,3),
                    llambda="identity",
                    lmu="identity",
                    lsigma="loge",
                    elambda=list(), emu=list(), esigma=list(),
                    dfmu.init=4,
                    dfsigma.init=2,
                    ilambda=1,
                    isigma=NULL)
{
    if(mode(llambda) != "character" && mode(llambda) != "name")
        llambda = as.character(substitute(llambda))
    if(mode(lmu) != "character" && mode(lmu) != "name")
        lmu = as.character(substitute(lmu))
    if(mode(lsigma) != "character" && mode(lsigma) != "name")
        lsigma = as.character(substitute(lsigma))
    if(!is.list(elambda)) elambda = list()
    if(!is.list(emu)) emu = list()
    if(!is.list(esigma)) esigma = list()
    if(!is.Numeric(ilambda))
        stop("bad input for argument 'ilambda'")
    if(length(isigma) && !is.Numeric(isigma, posit=TRUE))
        stop("bad input for argument 'isigma'")

    new("vglmff",
        blurb=c("LMS Quantile Regression ",
                "(Box-Cox transformation to normality)\n",
            "Links:    ",
            namesof("lambda", link=llambda, earg= elambda), ", ",
            namesof("mu", link=lmu, earg= emu), ", ",
            namesof("sigma", link=lsigma, earg= esigma)),
    constraints=eval(substitute(expression({
        constraints = cm.zero.vgam(constraints, x, .zero, M)
    }), list(.zero=zero))),
    initialize=eval(substitute(expression({
        if(ncol(cbind(y)) != 1)
            stop("response must be a vector or a one-column matrix")
        if(any(y<0, na.rm = TRUE))
            stop("negative responses not allowed")

        predictors.names =
            c(namesof("lambda", .llambda, earg= .elambda,  short= TRUE),
              namesof("mu",  .lmu, earg= .emu,  short= TRUE),
              namesof("sigma",  .lsigma, earg= .esigma,  short= TRUE))
 
        if(!length(etastart)) {

            fit500=vsmooth.spline(x=x[,min(ncol(x),2)],y=y,w=w, df= .dfmu.init)
            fv.init = c(predict(fit500, x=x[,min(ncol(x),2)])$y)

            lambda.init = if(is.Numeric( .ilambda)) .ilambda else 1.0
            sigma.init = if(is.null(.isigma)) {
                myratio = ((y/fv.init)^lambda.init - 1) / lambda.init
                if(is.Numeric( .dfsigma.init)) {
                    fit600 = vsmooth.spline(x=x[,min(ncol(x),2)], y=myratio^2,
                                            w=w, df= .dfsigma.init)
                    sqrt(c(abs(predict(fit600, x=x[,min(ncol(x),2)])$y)))
                } else 
                    sqrt(var(myratio))
            } else .isigma
 
            etastart = cbind(theta2eta(lambda.init, .llambda, earg= .elambda),
                             theta2eta(fv.init,     .lmu, earg= .emu),
                             theta2eta(sigma.init,  .lsigma, earg= .esigma))
        }
    }), list( .llambda=llambda, .lmu=lmu, .lsigma=lsigma,
              .elambda=elambda, .emu=emu, .esigma=esigma, 
              .dfmu.init=dfmu.init,
              .dfsigma.init=dfsigma.init,
              .ilambda=ilambda, .isigma=isigma ))),
    inverse=eval(substitute(function(eta, extra=NULL) {
        eta[,1] = eta2theta(eta[,1], .llambda, earg= .elambda)
        eta[,2] = eta2theta(eta[,2], .lmu, earg= .emu)
        eta[,3] = eta2theta(eta[,3], .lsigma, earg= .esigma)
        qtplot.lms.bcn(percentiles= .percentiles, eta=eta)
    }, list( .llambda=llambda, .lmu=lmu, .lsigma=lsigma,
             .elambda=elambda, .emu=emu, .esigma=esigma, 
             .percentiles=percentiles ))),
    last=eval(substitute(expression({
        misc$percentiles = .percentiles
        misc$links = c(lambda = .llambda, mu = .lmu, sigma = .lsigma)
        misc$earg = list(lambda = .elambda, mu = .emu, sigma = .esigma)
        misc$true.mu = FALSE    # $fitted is not a true mu
        if(control$cdf) {
            post$cdf = cdf.lms.bcn(y, eta0=matrix(c(lambda,mymu,sigma), 
                ncol=3, dimnames=list(dimnames(x)[[1]], NULL)))
        }
    }), list( .llambda=llambda, .lmu=lmu, .lsigma=lsigma,
              .elambda=elambda, .emu=emu, .esigma=esigma, 
              .percentiles=percentiles ))),
    loglikelihood=eval(substitute(
        function(mu,y,w, residuals= FALSE, eta, extra=NULL) {
            lambda = eta2theta(eta[,1], .llambda, earg= .elambda)
            mu = eta2theta(eta[,2], .lmu, earg= .emu)
            sigma = eta2theta(eta[,3], .lsigma, earg= .esigma)
            z = ((y/mu)^lambda - 1) / (lambda * sigma)
         if(residuals) stop("loglikelihood residuals not implemented yet") else
            sum(w * (lambda * log(y/mu) - log(sigma) - 0.5*z^2))
        }, list( .llambda=llambda, .lmu=lmu, .lsigma=lsigma,
                 .elambda=elambda, .emu=emu, .esigma=esigma ))),
    vfamily=c("lms.bcn", "lmscreg"),
    deriv=eval(substitute(expression({
        lambda = eta2theta(eta[,1], .llambda, earg= .elambda)
        mymu   = eta2theta(eta[,2], .lmu, earg= .emu)
        sigma  = eta2theta(eta[,3], .lsigma, earg= .esigma)
        z = ((y/mymu)^lambda - 1) / (lambda * sigma)
        z2m1 = z * z - 1
        dl.dlambda = z*(z - log(y/mymu) / sigma) / lambda - z2m1 * log(y/mymu)
        dl.dmu = z / (mymu * sigma) + z2m1 * lambda / mymu
        dl.dsigma = z2m1 / sigma
        dlambda.deta  = dtheta.deta(lambda, .llambda, earg= .elambda)
        dmu.deta  = dtheta.deta(mymu, .lmu, earg= .emu)
        dsigma.deta = dtheta.deta(sigma, .lsigma, earg= .esigma)
        w * cbind(dl.dlambda * dlambda.deta,
                  dl.dmu * dmu.deta,
                  dl.dsigma * dsigma.deta)
    }), list( .llambda=llambda, .lmu=lmu, .lsigma=lsigma,
              .elambda=elambda, .emu=emu, .esigma=esigma ))),
    weight=eval(substitute(expression({
        wz = matrix(as.numeric(NA), n, 6)
        wz[,iam(1,1,M)] = (7 * sigma^2 / 4) * dlambda.deta^2
        wz[,iam(2,2,M)] = (1 + 2*(lambda*sigma)^2)/(mymu*sigma)^2 * dmu.deta^2
        wz[,iam(3,3,M)] = (2 / sigma^2) * dsigma.deta^2
        wz[,iam(1,2,M)] = (-1 / (2 * mymu)) * dlambda.deta * dmu.deta
        wz[,iam(1,3,M)] = (lambda * sigma) * dlambda.deta * dsigma.deta
        wz[,iam(2,3,M)] = (2*lambda/(mymu * sigma)) * dmu.deta * dsigma.deta
        wz * w
    }), list( .llambda=llambda, .lmu=lmu, .lsigma=lsigma,
              .elambda=elambda, .emu=emu, .esigma=esigma ))))
}



lms.bcg = function(percentiles=c(25,50,75),
                   zero=c(1,3),
                   llambda="identity",
                   lmu="identity",
                   lsigma="loge",
                   elambda=list(), emu=list(), esigma=list(),
                   dfmu.init=4,
                   dfsigma.init=2,
                   ilambda=1,
                   isigma=NULL)
{
    if(mode(llambda) != "character" && mode(llambda) != "name")
        llambda = as.character(substitute(llambda))
    if(mode(lmu) != "character" && mode(lmu) != "name")
        lmu = as.character(substitute(lmu))
    if(mode(lsigma) != "character" && mode(lsigma) != "name")
        lsigma = as.character(substitute(lsigma))
    if(!is.list(elambda)) elambda = list()
    if(!is.list(emu)) emu = list()
    if(!is.list(esigma)) esigma = list()
    if(!is.Numeric(ilambda))
        stop("bad input for argument 'ilambda'")
    if(length(isigma) && !is.Numeric(isigma, posit=TRUE))
        stop("bad input for argument 'isigma'")

    new("vglmff",
    blurb=c("LMS Quantile Regression ",
            "(Box-Cox transformation to a Gamma distribution)\n",
            "Links:    ",
            namesof("lambda", link=llambda, earg= elambda), ", ",
            namesof("mu", link=lmu, earg= emu), ", ",
            namesof("sigma", link=lsigma, earg= esigma)),
    constraints=eval(substitute(expression({
        constraints = cm.zero.vgam(constraints, x, .zero, M)
    }), list(.zero=zero))),
    initialize=eval(substitute(expression({
      if(ncol(cbind(y)) != 1)
          stop("response must be a vector or a one-column matrix")
      if(any(y<0, na.rm = TRUE))
            stop("negative responses not allowed")

        predictors.names = c(
            namesof("lambda", .llambda, earg= .elambda,  short=TRUE),
            namesof("mu",     .lmu, earg= .emu,  short=TRUE),
            namesof("sigma",  .lsigma, earg= .esigma, short=TRUE))

        if(!length(etastart)) {

            fit500=vsmooth.spline(x=x[,min(ncol(x),2)],y=y,w=w, df= .dfmu.init)
            fv.init = c(predict(fit500, x=x[,min(ncol(x),2)])$y)

            lambda.init = if(is.Numeric( .ilambda)) .ilambda else 1.0

            sigma.init = if(is.null(.isigma)) {
               myratio=((y/fv.init)^lambda.init-1)/lambda.init #~(0,var=sigma^2)
                if(is.numeric( .dfsigma.init) && is.finite( .dfsigma.init)) {
                    fit600 = vsmooth.spline(x=x[,min(ncol(x),2)],
                                            y=(myratio)^2,
                                            w=w, df= .dfsigma.init)
                    sqrt(c(abs(predict(fit600, x=x[,min(ncol(x),2)])$y)))
                } else 
                    sqrt(var(myratio))
            } else .isigma

            etastart = cbind(theta2eta(lambda.init,  .llambda, earg= .elambda),
                             theta2eta(fv.init,      .lmu, earg= .emu),
                             theta2eta(sigma.init,   .lsigma, earg= .esigma))
        }
    }), list( .llambda=llambda, .lmu=lmu, .lsigma=lsigma,
              .elambda=elambda, .emu=emu, .esigma=esigma, 
              .dfmu.init=dfmu.init,
              .dfsigma.init=dfsigma.init,
              .ilambda=ilambda, .isigma=isigma ))),
    inverse=eval(substitute(function(eta, extra=NULL) {
        eta[,1] = eta2theta(eta[,1], .llambda, earg= .elambda)
        eta[,2] = eta2theta(eta[,2], .lmu, earg= .emu)
        eta[,3] = eta2theta(eta[,3], .lsigma, earg= .esigma)
        qtplot.lms.bcg(percentiles= .percentiles, eta=eta)
    }, list( .llambda=llambda, .lmu=lmu, .lsigma=lsigma,
             .elambda=elambda, .emu=emu, .esigma=esigma, 
             .percentiles=percentiles ))),
    last=eval(substitute(expression({
        misc$percentiles = .percentiles
        misc$link = c(lambda = .llambda, mu = .lmu, sigma = .lsigma)
        misc$earg = list(lambda = .elambda, mu = .emu, sigma = .esigma)
        misc$true.mu = FALSE    # $fitted is not a true mu
        if(control$cdf) {
            post$cdf = cdf.lms.bcg(y, eta0=matrix(c(lambda,mymu,sigma), 
                ncol=3, dimnames=list(dimnames(x)[[1]], NULL)))
        }
    }), list( .llambda=llambda, .lmu=lmu, .lsigma=lsigma,
              .elambda=elambda, .emu=emu, .esigma=esigma, 
              .percentiles=percentiles ))),
    loglikelihood=eval(substitute(
        function(mu,y,w, residuals= FALSE, eta, extra=NULL) {
            lambda = eta2theta(eta[,1], .llambda, earg= .elambda)
            mu = eta2theta(eta[,2], .lmu, earg= .emu)
            sigma = eta2theta(eta[,3], .lsigma, earg= .esigma)
            g = (y/mu)^lambda
            theta = 1 / (sigma * lambda)^2
         if(residuals) stop("loglikelihood residuals not implemented yet") else
            sum(w * (log(abs(lambda)) + theta*(log(theta)+log(g)-g) - 
                     lgamma(theta) - log(y)))
        }, list( .llambda=llambda, .lmu=lmu, .lsigma=lsigma,
                 .elambda=elambda, .emu=emu, .esigma=esigma ))),
    vfamily=c("lms.bcg", "lmscreg"),
    deriv=eval(substitute(expression({
        lambda = eta2theta(eta[,1], .llambda, earg= .elambda)
        mymu = eta2theta(eta[,2], .lmu, earg= .emu)
        sigma = eta2theta(eta[,3], .lsigma, earg= .esigma)

        g = (y/mymu)^lambda
        theta = 1 / (sigma * lambda)^2
        dd = digamma(theta)

        dl.dlambda = (1 + 2*theta*(dd+g-1-log(theta)-0.5*(g+1)*log(g)))/lambda
        dl.dmu = lambda * theta * (g-1) / mymu
        dl.dsigma = 2*theta*(dd+g-log(theta * g)-1) / sigma
        dlambda.deta = dtheta.deta(lambda, link=.llambda, earg= .elambda)
        dmu.deta = dtheta.deta(mymu, link=.lmu, earg= .emu)
        dsigma.deta = dtheta.deta(sigma, link=.lsigma, earg= .esigma)

        cbind(dl.dlambda * dlambda.deta,
              dl.dmu     * dmu.deta,
              dl.dsigma  * dsigma.deta) * w
    }), list( .llambda=llambda, .lmu=lmu, .lsigma=lsigma,
              .elambda=elambda, .emu=emu, .esigma=esigma ))),
    weight=eval(substitute(expression({
        tritheta = trigamma(theta)
        wz = matrix(0, n, 6)

        if(TRUE) {
            part2 = dd + 2/theta - 2*log(theta)
            wz[,iam(1,1,M)] = ((1 + theta*(tritheta*(1+4*theta) -
                               4*(1+1/theta) - log(theta)*(2/theta -
                               log(theta)) + dd*part2)) / lambda^2) *
                              dlambda.deta^2
        } else {
            temp = mean( g*(log(g))^2 )
            wz[,iam(1,1,M)] = ((4*theta*(theta*tritheta-1) -1 +
                               theta*temp)/lambda^2) *
                              dlambda.deta^2
        }

        wz[,iam(2,2,M)] = dmu.deta^2 / (mymu*sigma)^2
        wz[,iam(3,3,M)] = (4*theta*(theta*tritheta-1)/sigma^2) * dsigma.deta^2
        wz[,iam(1,2,M)] = (-theta * (dd + 1/theta - log(theta)) / mymu) *
                          dlambda.deta * dmu.deta
        wz[,iam(1,3,M)] = 2 * theta^1.5 * (2 * theta * tritheta - 2 -
                           1/theta) * dlambda.deta * dsigma.deta
        wz * w
    }), list( .llambda=llambda, .lmu=lmu, .lsigma=lsigma,
              .elambda=elambda, .emu=emu, .esigma=esigma ))))
}


dy.dpsi.yeojohnson = function(psi, lambda) {

    L = max(length(psi), length(lambda))
    psi = rep(psi, len=L); lambda = rep(lambda, len=L);
    ifelse(psi>0, (1 + psi * lambda)^(1/lambda - 1),
                  (1 - (2-lambda) * psi)^((lambda - 1)/(2-lambda)))
}

dyj.dy.yeojohnson = function(y, lambda) {
    L = max(length(y), length(lambda))
    y = rep(y, len=L); lambda = rep(lambda, len=L);
    ifelse(y>0, (1 + y)^(lambda - 1), (1 - y)^(1 - lambda))
}

yeo.johnson = function(y, lambda, derivative=0,
                        epsilon=sqrt(.Machine$double.eps), inverse= FALSE)
{

    if(!is.Numeric(derivative, allow=1, integ=TRUE) || derivative<0)
        stop("'derivative' must be a non-negative integer")
    ans = y
    if(!is.Numeric(epsilon, allow=1, posit=TRUE))
        stop("'epsilon' must be a single positive number")
    L = max(length(lambda), length(y))
    if(length(y) != L) y = rep(y, len=L)
    if(length(lambda) != L) lambda = rep(lambda, len=L)  # lambda may be of length 1

    if(inverse) {
        if(derivative!=0)
            stop("derivative must 0 when inverse=TRUE")
        if(any(index <- y >= 0 & abs(lambda) > epsilon))
            ans[index] = (y[index]*lambda[index] + 1)^(1/lambda[index]) - 1
        if(any(index <- y >= 0 & abs(lambda) <= epsilon))
            ans[index] = expm1(y[index])
        if(any(index <- y <  0 & abs(lambda-2) > epsilon))
            ans[index] = 1-(-(2-lambda[index])*y[index]+1)^(1/(2-lambda[index]))
        if(any(index <- y <  0 & abs(lambda-2) <= epsilon))
            ans[index] = -expm1(-y[index])
        return(ans)
    }
    if(derivative==0) {
        if(any(index <- y >= 0 & abs(lambda) > epsilon))
            ans[index] = ((y[index]+1)^(lambda[index]) - 1) / lambda[index]
        if(any(index <- y >= 0 & abs(lambda) <= epsilon))
            ans[index] = log1p(y[index])
        if(any(index <- y <  0 & abs(lambda-2) > epsilon))
            ans[index] = -((-y[index]+1)^(2-lambda[index])-1)/(2-lambda[index])
        if(any(index <- y <  0 & abs(lambda-2) <= epsilon))
            ans[index] = -log1p(-y[index])
    } else {
        psi <- Recall(y=y, lambda=lambda, derivative=derivative-1,
                      epsilon=epsilon, inverse=inverse)
        if(any(index <- y >= 0 & abs(lambda) > epsilon))
            ans[index] = ( (y[index]+1)^(lambda[index]) *
                          (log1p(y[index]))^(derivative) - derivative *
                          psi[index] ) / lambda[index]
        if(any(index <- y >= 0 & abs(lambda) <= epsilon))
            ans[index] = (log1p(y[index]))^(derivative + 1) / (derivative + 1)
        if(any(index <- y <  0 & abs(lambda-2) > epsilon))
            ans[index] = -( (-y[index]+1)^(2-lambda[index]) *
                          (-log1p(-y[index]))^(derivative) - derivative *
                          psi[index] ) / (2-lambda[index])
        if(any(index <- y <  0 & abs(lambda-2) <= epsilon))
            ans[index] = (-log1p(-y[index]))^(derivative + 1) / (derivative + 1)
    }
    ans
}


dpsi.dlambda.yjn = function(psi, lambda, mymu, sigma,
                            derivative=0, smallno=1.0e-8) {

    if(!is.Numeric(derivative, allow=1, integ=TRUE) || derivative<0)
        stop("'derivative' must be a non-negative integer")
    if(!is.Numeric(smallno, allow=1, posit=TRUE))
        stop("'smallno' must be a single positive number")

    L = max(length(psi), length(lambda), length(mymu), length(sigma))
    if(length(psi) != L) psi = rep(psi, len=L)
    if(length(lambda) != L) lambda = rep(lambda, len=L)
    if(length(mymu) != L) mymu = rep(mymu, len=L)
    if(length(sigma) != L) sigma = rep(sigma, len=L)

    answer = matrix(as.numeric(NA), L, derivative+1)
    CC = psi >= 0
    BB = ifelse(CC, lambda, -2+lambda)
    AA = psi * BB 
    temp8 = if(derivative > 0) {
        answer[,1:derivative] =
            Recall(psi=psi, lambda=lambda, mymu=mymu, sigma=sigma,
                   derivative=derivative-1, smallno=smallno) 
        answer[,derivative] * derivative
    } else { 
        0
    }
    answer[,1+derivative] = ((AA+1) * (log1p(AA)/BB)^derivative - temp8) / BB

    pos = (CC & abs(lambda) <= smallno) | (!CC & abs(lambda-2) <= smallno)
    if(any(pos)) 
        answer[pos,1+derivative] = (answer[pos,1]^(1+derivative))/(derivative+1)
    answer
}

gh.weight.yjn.11 = function(z, lambda, mymu, sigma, derivmat=NULL) {


    if(length(derivmat)) {
        ((derivmat[,2]/sigma)^2 + sqrt(2) * z * derivmat[,3] / sigma) / sqrt(pi)
    } else {
        # Long-winded way 
        psi = mymu + sqrt(2) * sigma * z
        (1 / sqrt(pi)) *
        (dpsi.dlambda.yjn(psi, lambda, mymu, sigma, derivative=1)[,2]^2 +
        (psi - mymu) * 
        dpsi.dlambda.yjn(psi, lambda, mymu, sigma, derivative=2)[,3]) / sigma^2
    }
}

gh.weight.yjn.12 = function(z, lambda, mymu, sigma, derivmat=NULL) {
    if(length(derivmat)) {
        (-derivmat[,2]) / (sqrt(pi) * sigma^2)
    } else {
        psi = mymu + sqrt(2) * sigma * z
        (1 / sqrt(pi)) *
        (- dpsi.dlambda.yjn(psi, lambda, mymu, sigma, derivative=1)[,2]) / sigma^2
    }
}

gh.weight.yjn.13 = function(z, lambda, mymu, sigma, derivmat=NULL) {
    if(length(derivmat)) {
        sqrt(8 / pi) * (-derivmat[,2]) * z / sigma^2
    } else {
        psi = mymu + sqrt(2) * sigma * z
        (1 / sqrt(pi)) *
        (-2 * dpsi.dlambda.yjn(psi, lambda, mymu, sigma, derivative=1)[,2]) *
        (psi - mymu) / sigma^3
    }
}


glag.weight.yjn.11 = function(z, lambda, mymu, sigma, derivmat=NULL) {


    if(length(derivmat)) {
        derivmat[,4] * (derivmat[,2]^2 + sqrt(2) * sigma * z * derivmat[,3])
    } else {
        psi = mymu + sqrt(2) * sigma * z
        discontinuity = -mymu / (sqrt(2) * sigma)
        (1 / (2 * sqrt((z-discontinuity^2)^2 + discontinuity^2))) *
        (1 / sqrt(pi)) *
        (dpsi.dlambda.yjn(psi, lambda, mymu, sigma, derivative=1)[,2]^2 +
        (psi - mymu) * 
        dpsi.dlambda.yjn(psi, lambda, mymu, sigma, derivative=2)[,3]) / sigma^2
    }
}

glag.weight.yjn.12 = function(z, lambda, mymu, sigma, derivmat=NULL) {
    discontinuity = -mymu / (sqrt(2) * sigma)
    if(length(derivmat)) {
        derivmat[,4] * (-derivmat[,2])
    } else {
        psi = mymu + sqrt(2) * sigma * z
        (1 / (2 * sqrt((z-discontinuity^2)^2 + discontinuity^2))) *
        (1 / sqrt(pi)) *
        (- dpsi.dlambda.yjn(psi, lambda, mymu, sigma, derivative=1)[,2]) / sigma^2
    }
}

glag.weight.yjn.13 = function(z, lambda, mymu, sigma, derivmat=NULL) {
    if(length(derivmat)) {
        derivmat[,4] * (-derivmat[,2]) * sqrt(8) * z
    } else {
        psi = mymu + sqrt(2) * sigma * z
        discontinuity = -mymu / (sqrt(2) * sigma)
        (1 / (2 * sqrt((z-discontinuity^2)^2 + discontinuity^2))) *
        (1 / sqrt(pi)) *
        (-2 * dpsi.dlambda.yjn(psi, lambda, mymu, sigma, derivative=1)[,2]) *
        (psi - mymu) / sigma^3
    }
}


gleg.weight.yjn.11 = function(z, lambda, mymu, sigma, derivmat=NULL) {




    if(length(derivmat)) {
        derivmat[,4] * (derivmat[,2]^2 + sqrt(2) * sigma * z * derivmat[,3])
    } else {
        psi = mymu + sqrt(2) * sigma * z
        (exp(-z^2) / sqrt(pi)) *
        (dpsi.dlambda.yjn(psi, lambda, mymu, sigma, derivative=1)[,2]^2 +
        (psi - mymu) * 
        dpsi.dlambda.yjn(psi, lambda, mymu, sigma, derivative=2)[,3]) / sigma^2
    }
}

gleg.weight.yjn.12 = function(z, lambda, mymu, sigma, derivmat=NULL) {
    if(length(derivmat)) {
        derivmat[,4] * (- derivmat[,2])
    } else {
        psi = mymu + sqrt(2) * sigma * z
        (exp(-z^2) / sqrt(pi)) *
        (- dpsi.dlambda.yjn(psi, lambda, mymu, sigma, derivative=1)[,2]) / sigma^2
    }
}

gleg.weight.yjn.13 = function(z, lambda, mymu, sigma, derivmat=NULL) {
    if(length(derivmat)) {
        derivmat[,4] * (-derivmat[,2]) * sqrt(8) * z
    } else {
        psi = mymu + sqrt(2) * sigma * z
        (exp(-z^2) / sqrt(pi)) *
        (-2 * dpsi.dlambda.yjn(psi, lambda, mymu, sigma, derivative=1)[,2]) *
        (psi - mymu) / sigma^3
    }
}



lms.yjn2.control <- function(save.weight=TRUE, ...)
{
    list(save.weight=save.weight)
}

lms.yjn2 = function(percentiles=c(25,50,75),
                    zero=c(1,3),
                    llambda="identity",
                    lmu="identity",
                    lsigma="loge",
                    elambda=list(), emu=list(), esigma=list(),
                    dfmu.init=4,
                    dfsigma.init=2,
                    ilambda=1.0,
                    isigma=NULL,
                    yoffset=NULL,
                    nsimEIM=250)
{

    if(mode(llambda) != "character" && mode(llambda) != "name")
        llambda = as.character(substitute(llambda))
    if(mode(lmu) != "character" && mode(lmu) != "name")
        lmu = as.character(substitute(lmu))
    if(mode(lsigma) != "character" && mode(lsigma) != "name")
        lsigma = as.character(substitute(lsigma))
    if(!is.list(elambda)) elambda = list()
    if(!is.list(emu)) emu = list()
    if(!is.list(esigma)) esigma = list()
    if(!is.Numeric(ilambda))
        stop("bad input for argument 'ilambda'")
    if(length(isigma) && !is.Numeric(isigma, posit=TRUE))
        stop("bad input for argument 'isigma'")

    new("vglmff",
    blurb=c("LMS Quantile Regression (Yeo-Johnson transformation",
            " to normality)\n",
            "Links:    ",
            namesof("lambda", link=llambda, earg= elambda),
            ", ",
            namesof("mu", link=lmu, earg= emu),
            ", ",
            namesof("sigma", link=lsigma, earg= esigma)),
    constraints=eval(substitute(expression({
        constraints = cm.zero.vgam(constraints, x, .zero, M)
    }), list(.zero=zero))),
    initialize=eval(substitute(expression({
      if(ncol(cbind(y)) != 1)
          stop("response must be a vector or a one-column matrix")
        predictors.names =
          c(namesof("lambda", .llambda, earg= .elambda, short= TRUE),
            namesof("mu",     .lmu,     earg= .emu,     short= TRUE),
            namesof("sigma",  .lsigma, earg= .esigma,  short= TRUE))

        y.save = y
        yoff = if(is.Numeric( .yoffset)) .yoffset else -median(y) 
        extra$yoffset = yoff
        y = y + yoff

        if(!length(etastart)) {
            lambda.init = if(is.Numeric( .ilambda)) .ilambda else 1.

            y.tx = yeo.johnson(y, lambda.init)
            fv.init = 
            if(smoothok <- (length(unique(sort(x[,min(ncol(x),2)]))) > 7)) {
                fit700=vsmooth.spline(x=x[,min(ncol(x),2)],
                                      y=y.tx, w=w, df= .dfmu.init)
                c(predict(fit700, x=x[,min(ncol(x),2)])$y)
            } else {
                rep(weighted.mean(y, w), len=n)
            }

            sigma.init = if(!is.Numeric(.isigma)) {
                              if(is.Numeric( .dfsigma.init) && smoothok) {
                                   fit710 = vsmooth.spline(x=x[,min(ncol(x),2)],
                                            y=(y.tx - fv.init)^2,
                                            w=w, df= .dfsigma.init)
                                   sqrt(c(abs(predict(fit710,
                                        x=x[,min(ncol(x),2)])$y)))
                              } else {
                                   sqrt( sum( w * (y.tx - fv.init)^2 ) / sum(w) )
                              }
                          } else
                              .isigma

            etastart = matrix(0, n, 3)
            etastart[,1] = theta2eta(lambda.init, .llambda, earg=.elambda)
            etastart[,2] = theta2eta(fv.init, .lmu, earg=.emu)
            etastart[,3] = theta2eta(sigma.init, .lsigma, earg=.esigma)

        }
    }), list(.llambda=llambda, .lmu=lmu, .lsigma=lsigma,
             .elambda=elambda, .emu=emu, .esigma=esigma, 
             .dfmu.init=dfmu.init,
             .dfsigma.init=dfsigma.init,
             .ilambda=ilambda,
             .yoffset=yoffset,
             .isigma=isigma))),
    inverse=eval(substitute(function(eta, extra=NULL) {
        eta[,1] = eta2theta(eta[,1], .llambda, earg= .elambda)
        eta[,3] = eta2theta(eta[,3], .lsigma, earg= .esigma)
        qtplot.lms.yjn(percentiles= .percentiles, eta=eta, yoffset= extra$yoff)
    }, list(.percentiles=percentiles,
            .esigma=esigma, .elambda=elambda,
            .llambda=llambda,
            .lsigma=lsigma))),
    last=eval(substitute(expression({
        misc$expected = TRUE
        misc$nsimEIM = .nsimEIM
        misc$percentiles = .percentiles
        misc$link = c(lambda= .llambda, mu= .lmu, sigma= .lsigma)
        misc$earg = list(lambda = .elambda, mu = .emu, sigma = .esigma)
        misc$true.mu = FALSE # $fitted is not a true mu
        misc[["yoffset"]] = extra$yoffset

        y = y.save   # Restore back the value; to be attached to object

        if(control$cdf) {
            post$cdf = cdf.lms.yjn(y + misc$yoffset,
                eta0=matrix(c(lambda,mymu,sigma), 
                ncol=3, dimnames=list(dimnames(x)[[1]], NULL)))
        }
    }), list(.percentiles=percentiles,
             .elambda=elambda, .emu=emu, .esigma=esigma, 
             .nsimEIM=nsimEIM,
             .llambda=llambda, .lmu=lmu, .lsigma=lsigma ))),
    loglikelihood=eval(substitute(
        function(mu,y,w, residuals= FALSE, eta, extra=NULL) {
            lambda = eta2theta(eta[,1], .llambda, earg= .elambda)
            mu = eta2theta(eta[,2], .lmu, earg= .emu)
            sigma = eta2theta(eta[,3], .lsigma, earg= .esigma)
            psi = yeo.johnson(y, lambda)
         if(residuals) stop("loglikelihood residuals not implemented yet") else
            sum(w * (-log(sigma) - 0.5 * ((psi-mu)/sigma)^2 +
                     (lambda-1) * sign(y) * log1p(abs(y))))
        }, list( .elambda=elambda, .emu=emu, .esigma=esigma, 
                 .llambda=llambda, .lmu=lmu,
                 .lsigma=lsigma ))),
    vfamily=c("lms.yjn2", "lmscreg"),
    deriv=eval(substitute(expression({
        lambda = eta2theta(eta[,1], .llambda, earg= .elambda)
        mymu = eta2theta(eta[,2], .lmu, earg= .emu)
        sigma = eta2theta(eta[,3], .lsigma, earg= .esigma)
        dlambda.deta = dtheta.deta(lambda, link=.llambda, earg= .elambda)
        dmu.deta = dtheta.deta(mymu, link=.lmu, earg= .emu)
        dsigma.deta = dtheta.deta(sigma, link=.lsigma, earg= .esigma)

        psi = yeo.johnson(y, lambda)
        d1 = yeo.johnson(y, lambda, deriv=1)
        AA = (psi - mymu) / sigma 
        dl.dlambda = -AA * d1 /sigma + sign(y) * log1p(abs(y))
        dl.dmu = AA / sigma 
        dl.dsigma = (AA^2 -1) / sigma
        dthetas.detas = cbind(dlambda.deta, dmu.deta, dsigma.deta)
        w * cbind(dl.dlambda, dl.dmu, dl.dsigma) * dthetas.detas
    }), list( .elambda=elambda, .emu=emu, .esigma=esigma, 
              .llambda=llambda, .lmu=lmu,
                 .lsigma=lsigma ))),
    weight=eval(substitute(expression({


        run.varcov = 0
        ind1 = iam(NA, NA, M=M, both=TRUE, diag=TRUE)
        for(ii in 1:( .nsimEIM )) {
            psi = rnorm(n, mymu, sigma)
            ysim = yeo.johnson(y=psi, lam=lambda, inv=TRUE)
            d1 = yeo.johnson(ysim, lambda, deriv=1)
            AA = (psi - mymu) / sigma 
            dl.dlambda = -AA * d1 /sigma + sign(ysim) * log1p(abs(ysim))
            dl.dmu = AA / sigma 
            dl.dsigma = (AA^2 -1) / sigma
            rm(ysim)
            temp3 = cbind(dl.dlambda, dl.dmu, dl.dsigma)
            run.varcov = ((ii-1) * run.varcov +
                       temp3[,ind1$row.index]*temp3[,ind1$col.index]) / ii
        }

        if(intercept.only)
            run.varcov = matrix(colMeans(run.varcov),
                                nr=n, nc=ncol(run.varcov), byrow=TRUE)


        wz = run.varcov * dthetas.detas[,ind1$row] * dthetas.detas[,ind1$col]
        dimnames(wz) = list(rownames(wz), NULL)  # Remove the colnames
        wz * w
    }), list(.lsigma=lsigma,
             .esigma=esigma, .elambda=elambda,
             .nsimEIM=nsimEIM,
             .llambda=llambda))))
}


lms.yjn <- function(percentiles=c(25,50,75),
                    zero=c(1,3),
                    llambda="identity",
                    lsigma="loge",
                    elambda=list(), esigma=list(),
                    dfmu.init=4,
                    dfsigma.init=2,
                    ilambda=1.0,
                    isigma=NULL,
                    rule=c(10,5),
                    yoffset=NULL,
                    diagW=FALSE, iters.diagW=6)
{



    if(mode(lsigma) != "character" && mode(lsigma) != "name")
        lsigma = as.character(substitute(lsigma))
    if(mode(llambda) != "character" && mode(llambda) != "name")
        llambda = as.character(substitute(llambda))
    if(!is.list(elambda)) elambda = list()
    if(!is.list(esigma)) esigma = list()

    rule = rule[1] # Number of points (common) for all the quadrature schemes
    if(rule != 5 && rule != 10)
        stop("only rule=5 or 10 is supported")

    new("vglmff",
    blurb=c("LMS Quantile Regression ",
            "(Yeo-Johnson transformation to normality)\n",
            "Links:    ",
            namesof("lambda", link=llambda, earg= elambda),
            ", mu, ",
            namesof("sigma", link=lsigma, earg= esigma)),
    constraints=eval(substitute(expression({
        constraints = cm.zero.vgam(constraints, x, .zero, M)
    }), list(.zero=zero))),
    initialize=eval(substitute(expression({
      if(ncol(cbind(y)) != 1)
          stop("response must be a vector or a one-column matrix")
        predictors.names =
          c(namesof("lambda", .llambda, earg= .elambda, short= TRUE),
                "mu",
            namesof("sigma",  .lsigma, earg= .esigma,  short= TRUE))

        y.save = y
        yoff = if(is.Numeric( .yoffset)) .yoffset else -median(y) 
        extra$yoffset = yoff
        y = y + yoff

        if(!length(etastart)) {

            lambda.init = if(is.Numeric( .ilambda)) .ilambda else 1.0

            y.tx = yeo.johnson(y, lambda.init)
            if(smoothok <- (length(unique(sort(x[,min(ncol(x),2)]))) > 7)) {
                fit700=vsmooth.spline(x=x[,min(ncol(x),2)],
                                      y=y.tx, w=w, df= .dfmu.init)
                fv.init = c(predict(fit700, x=x[,min(ncol(x),2)])$y)
            } else {
                fv.init = rep(weighted.mean(y, w), len=n)
            }

            sigma.init = if(!is.Numeric(.isigma)) {
                              if(is.Numeric( .dfsigma.init) && smoothok) {
                                   fit710 = vsmooth.spline(x=x[,min(ncol(x),2)],
                                            y=(y.tx - fv.init)^2,
                                            w=w, df= .dfsigma.init)
                                   sqrt(c(abs(predict(fit710,
                                        x=x[,min(ncol(x),2)])$y)))
                              } else {
                                   sqrt( sum( w * (y.tx - fv.init)^2 ) / sum(w) )
                              }
                          } else
                              .isigma

            etastart = cbind(theta2eta(lambda.init,.llambda, earg=.elambda),
                             fv.init,
                             theta2eta(sigma.init, .lsigma, earg=.esigma))

        }
    }), list(.lsigma=lsigma,
             .llambda=llambda,
             .esigma=esigma, .elambda=elambda,
             .dfmu.init=dfmu.init,
             .dfsigma.init=dfsigma.init,
             .ilambda=ilambda,
             .yoffset=yoffset,
             .isigma=isigma))),
    inverse=eval(substitute(function(eta, extra=NULL) {
        eta[,1] = eta2theta(eta[,1], .llambda, earg= .elambda)
        eta[,3] = eta2theta(eta[,3], .lsigma, earg= .esigma)
        qtplot.lms.yjn(percentiles= .percentiles, eta=eta, yoffset= extra$yoff)
    }, list(.percentiles=percentiles,
             .esigma=esigma, .elambda=elambda,
            .llambda=llambda,
            .lsigma=lsigma))),
    last=eval(substitute(expression({
        misc$percentiles = .percentiles
        misc$link = c(lambda= .llambda, mu= "identity", sigma= .lsigma)
        misc$earg = list(lambda = .elambda, mu = list(), sigma = .esigma)
        misc$true.mu = FALSE    # $fitted is not a true mu
        misc[["yoffset"]] = extra$yoff

        y = y.save   # Restore back the value; to be attached to object

        if(control$cdf) {
            post$cdf = cdf.lms.yjn(y + misc$yoffset,
                eta0=matrix(c(lambda,mymu,sigma), 
                ncol=3, dimnames=list(dimnames(x)[[1]], NULL)))
        }
    }), list(.percentiles=percentiles,
             .esigma=esigma, .elambda=elambda,
            .llambda=llambda,
            .lsigma=lsigma))),
    loglikelihood=eval(substitute(
        function(mu,y,w, residuals= FALSE, eta, extra=NULL) {
            lambda = eta2theta(eta[,1], .llambda, earg= .elambda)
            mu = eta[,2]
            sigma = eta2theta(eta[,3], .lsigma, earg= .esigma)
            psi = yeo.johnson(y, lambda)
         if(residuals) stop("loglikelihood residuals not implemented yet") else
            sum(w * (-log(sigma) - 0.5 * ((psi-mu)/sigma)^2 +
                     (lambda-1) * sign(y) * log1p(abs(y))))
        }, list( .esigma=esigma, .elambda=elambda,
                 .lsigma=lsigma, .llambda=llambda))),
    vfamily=c("lms.yjn", "lmscreg"),
    deriv=eval(substitute(expression({
        lambda = eta2theta(eta[,1], .llambda, earg= .elambda)
        mymu = eta[,2]
        sigma = eta2theta(eta[,3], .lsigma, earg= .esigma)

        psi = yeo.johnson(y, lambda)
        d1 = yeo.johnson(y, lambda, deriv=1)
        AA = (psi - mymu) / sigma 

        dl.dlambda = -AA * d1 /sigma + sign(y) * log1p(abs(y))
        dl.dmu = AA / sigma 
        dl.dsigma = (AA^2 -1) / sigma
        dlambda.deta = dtheta.deta(lambda, link=.llambda, earg= .elambda)
        dsigma.deta = dtheta.deta(sigma, link=.lsigma, earg= .esigma)

        cbind(dl.dlambda * dlambda.deta,
              dl.dmu,
              dl.dsigma * dsigma.deta) * w
    }), list( .esigma=esigma, .elambda=elambda,
              .lsigma=lsigma, .llambda=llambda ))),
    weight=eval(substitute(expression({
        wz = matrix(0, n, 6)


        wz[,iam(2,2,M)] = 1 / sigma^2
        wz[,iam(3,3,M)] = 2 * wz[,iam(2,2,M)]   # 2 / sigma^2


        if(.rule == 10) {
        glag.abs=c(0.13779347054,0.729454549503,1.80834290174,3.40143369785,
                     5.55249614006,8.33015274676,11.8437858379,16.2792578314,
                     21.996585812, 29.9206970123)
        glag.wts = c(0.308441115765, 0.401119929155, 0.218068287612,
                     0.0620874560987, 0.00950151697517, 0.000753008388588, 
                     2.82592334963e-5,
                     4.24931398502e-7, 1.83956482398e-9, 9.91182721958e-13)
        } else {
        glag.abs = c(0.2635603197180449, 1.4134030591060496, 3.5964257710396850,
                     7.0858100058570503, 12.6408008442729685)
        glag.wts=c(5.217556105826727e-01,3.986668110832433e-01,7.594244968176882e-02,
                     3.611758679927785e-03, 2.336997238583738e-05)
        }

        if(.rule == 10) {
        sgh.abs = c(0.03873852801690856, 0.19823332465268367, 0.46520116404433082,
                    0.81686197962535023, 1.23454146277833154, 1.70679833036403172,
                    2.22994030591819214, 2.80910399394755972, 3.46387269067033854,
                    4.25536209637269280)
        sgh.wts=c(9.855210713854302e-02,2.086780884700499e-01,2.520517066468666e-01,
             1.986843323208932e-01,9.719839905023238e-02,2.702440190640464e-02,
             3.804646170194185e-03, 2.288859354675587e-04, 4.345336765471935e-06,
             1.247734096219375e-08)
        } else {
      sgh.abs = c(0.1002421519682381, 0.4828139660462573, 1.0609498215257607,
                  1.7797294185202606, 2.6697603560875995)
      sgh.wts=c(0.2484061520284881475,0.3923310666523834311,0.2114181930760276606,
                0.0332466603513424663, 0.0008248533445158026)
        }

        if(.rule == 10) {
            gleg.abs = c(-0.973906528517, -0.865063366689, -0.679409568299,
                         -0.433395394129, -0.148874338982)
            gleg.abs = c(gleg.abs, rev(-gleg.abs))
            gleg.wts = c(0.0666713443087, 0.149451349151, 0.219086362516,
                         0.26926671931, 0.295524224715)
            gleg.wts = c(gleg.wts, rev(gleg.wts))
        } else {
            gleg.abs = c(-0.9061798459386643,-0.5384693101056820, 0,
                          0.5384693101056828, 0.9061798459386635)
            gleg.wts=c(0.2369268850561853,0.4786286704993680,0.5688888888888889,
                       0.4786286704993661, 0.2369268850561916)
        }


        discontinuity = -mymu/(sqrt(2)*sigma) # Needs to be near 0, eg within 4


        LL = pmin(discontinuity, 0)
        UU = pmax(discontinuity, 0)
        if(FALSE) {
            AA = (UU-LL)/2
            for(kk in 1:length(gleg.wts)) {
                temp1 = AA * gleg.wts[kk] 
                abscissae = (UU+LL)/2 + AA * gleg.abs[kk]
                psi = mymu + sqrt(2) * sigma * abscissae
                temp9 = dpsi.dlambda.yjn(psi, lambda, mymu, sigma, derivative=2)
                temp9 = cbind(temp9, exp(-abscissae^2) / (sqrt(pi) * sigma^2))
    
                wz[,iam(1,1,M)] = wz[,iam(1,1,M)] + temp1 *
                    gleg.weight.yjn.11(abscissae, lambda, mymu, sigma, temp9)
                wz[,iam(1,2,M)] = wz[,iam(1,2,M)] + temp1 *
                    gleg.weight.yjn.12(abscissae, lambda, mymu, sigma, temp9)
                wz[,iam(1,3,M)] = wz[,iam(1,3,M)] + temp1 *
                    gleg.weight.yjn.13(abscissae, lambda, mymu, sigma, temp9)
            }
        } else {
            temp9 = dotFortran(name="yjngintf", as.double(LL), as.double(UU),
                     as.double(gleg.abs), as.double(gleg.wts), as.integer(n),
                     as.integer(length(gleg.abs)), as.double(lambda),
                     as.double(mymu), as.double(sigma), answer=double(3*n),
                     eps=as.double(1.0e-5))$ans
            dim(temp9) = c(3,n)
            wz[,iam(1,1,M)] = temp9[1,]
            wz[,iam(1,2,M)] = temp9[2,]
            wz[,iam(1,3,M)] = temp9[3,]
        }



        for(kk in 1:length(sgh.wts)) {

            abscissae = sign(-discontinuity) * sgh.abs[kk]
            psi = mymu + sqrt(2) * sigma * abscissae   # abscissae = z
            temp9 = dpsi.dlambda.yjn(psi, lambda, mymu, sigma, derivative=2)
            wz[,iam(1,1,M)] = wz[,iam(1,1,M)] + sgh.wts[kk] * 
                gh.weight.yjn.11(abscissae, lambda, mymu, sigma, temp9)
            wz[,iam(1,2,M)] = wz[,iam(1,2,M)] + sgh.wts[kk] * 
                gh.weight.yjn.12(abscissae, lambda, mymu, sigma, temp9)
            wz[,iam(1,3,M)] = wz[,iam(1,3,M)] + sgh.wts[kk] * 
                gh.weight.yjn.13(abscissae, lambda, mymu, sigma, temp9)
        }

        temp1 = exp(-discontinuity^2)
        for(kk in 1:length(glag.wts)) {
            abscissae = sign(discontinuity) * sqrt(glag.abs[kk]) + discontinuity^2
            psi = mymu + sqrt(2) * sigma * abscissae
            temp9 = dpsi.dlambda.yjn(psi, lambda, mymu, sigma, derivative=2)
            temp9 = cbind(temp9, 
                   1 / (2 * sqrt((abscissae-discontinuity^2)^2 + discontinuity^2) *
                        sqrt(pi) * sigma^2))
            temp7 = temp1 * glag.wts[kk]
            wz[,iam(1,1,M)] = wz[,iam(1,1,M)] + temp7 * 
                glag.weight.yjn.11(abscissae, lambda, mymu, sigma, temp9)
            wz[,iam(1,2,M)] = wz[,iam(1,2,M)] + temp7 * 
                glag.weight.yjn.12(abscissae, lambda, mymu, sigma, temp9)
            wz[,iam(1,3,M)] = wz[,iam(1,3,M)] + temp7 * 
                glag.weight.yjn.13(abscissae, lambda, mymu, sigma, temp9)
        }

        wz[,iam(1,1,M)] = wz[,iam(1,1,M)] * dlambda.deta^2
        wz[,iam(1,2,M)] = wz[,iam(1,2,M)] * dlambda.deta
        wz[,iam(1,3,M)] = wz[,iam(1,3,M)] * dsigma.deta * dlambda.deta
        if( .diagW && iter <= .iters.diagW) {
            wz[,iam(1,2,M)] = wz[,iam(1,3,M)] = 0
        }
        wz[,iam(2,3,M)] = wz[,iam(2,3,M)] * dsigma.deta
        wz[,iam(3,3,M)] = wz[,iam(3,3,M)] * dsigma.deta^2

        wz = wz * w
        wz
    }), list(.lsigma=lsigma,
             .esigma=esigma, .elambda=elambda,
             .rule=rule,
             .diagW=diagW,
             .iters.diagW=iters.diagW,
             .llambda=llambda))))
}



lmscreg.control <- function(cdf= TRUE, at.arg=NULL, x0=NULL, ...)
{

    if(!is.logical(cdf)) {
        warning("'cdf' is not logical; using TRUE instead")
        cdf = TRUE
    }
    list(cdf=cdf, at.arg=at.arg, x0=x0)
}





Wr1 <- function(r, w) ifelse(r <= 0, 1, w)


Wr2 <- function(r, w) (r <= 0) * 1 + (r > 0) * w


amlnormal.deviance = function(mu, y, w, residuals = FALSE, eta, extra=NULL) {

    M <- length(extra$w.als)

    if(M > 1) y = matrix(y,extra$n,extra$M)

    devi =  cbind((y - mu)^2)
    if(residuals) {
        stop("not sure here")
        wz = VGAM.weights.function(w = w, M = extra$M, n = extra$n)
        return((y - mu) * sqrt(wz) * matrix(extra$w.als,extra$n,extra$M))
    } else {
        all.deviances = numeric(M)
        myresid = matrix(y,extra$n,extra$M) - cbind(mu)
        for(ii in 1:M) all.deviances[ii] = sum(w * devi[,ii] *
                               Wr1(myresid[,ii], w=extra$w.als[ii]))
    }
    if(is.logical(extra$individual) && extra$individual)
        all.deviances else sum(all.deviances)
}



 amlnormal <- function(w.als=1, parallel=FALSE,
                       lexpectile = "identity", eexpectile = list(),
                       iexpectile = NULL,
                       method.init=1, digw=4)
{


    if(!is.Numeric(w.als, posit=TRUE))
        stop("'w.als' must be a vector of positive values")
    if(!is.Numeric(method.init, allow=1, integ=TRUE, posit=TRUE) ||
       method.init > 3) stop("argument 'method.init' must be 1, 2 or 3")
    if(mode(lexpectile) != "character" && mode(lexpectile) != "name")
        lexpectile = as.character(substitute(lexpectile))
    if(!is.list(eexpectile)) eexpectile = list()
    if(length(iexpectile) && !is.Numeric(iexpectile))
        stop("bad input for argument 'iexpectile'")

    new("vglmff",
        blurb=c("Asymmetric least squares quantile regression\n\n",
                "Links:    ",
                namesof("expectile", link=lexpectile, earg= eexpectile)),
    constraints=eval(substitute(expression({
        constraints = cm.vgam(matrix(1,M,1), x, .parallel, constraints)
    }), list( .parallel=parallel ))),
    deviance= function(mu, y, w, residuals = FALSE, eta, extra=NULL) {
        amlnormal.deviance(mu=mu, y=y, w=w, residuals=residuals,
                           eta=eta, extra=extra)
    },
    initialize=eval(substitute(expression({
        extra$w.als = .w.als
        if(ncol(y <- cbind(y)) != 1)
            stop("response must be a vector or a one-column matrix")
        extra$M = M = length(extra$w.als)  # Recycle if necessary
        extra$n = n
        extra$y.names = y.names = paste("w.als=", round(extra$w.als, dig=.digw),
                                        sep="")
        predictors.names = c(namesof(
            paste("expectile(",y.names,")", sep=""), .lexpectile,
                   earg=.eexpectile, tag=FALSE))

        if(!length(etastart)) {
            mean.init = if( .method.init == 1)
                    rep(median(y), length=n) else
                if( .method.init == 2)
                    rep(weighted.mean(y, w), length=n) else {
                        junk = if(is.R()) lm.wfit(x=x, y=y, w=w) else
                               lm.wfit(x=x, y=y, w=w, method="qr")
                        junk$fitted
                    }
            if(length( .iexpectile))
                mean.init = matrix( .iexpectile, n, M, byrow=TRUE)
            etastart = matrix(theta2eta(mean.init, .lexpectile,
                                        earg= .eexpectile), n, M)
        }
    }), list( .lexpectile=lexpectile, .eexpectile=eexpectile,
              .iexpectile=iexpectile,
              .method.init=method.init, .digw = digw, .w.als=w.als ))),
    inverse=eval(substitute(function(eta, extra=NULL) {
        ans = eta = as.matrix(eta)
        for(ii in 1:ncol(eta))
            ans[,ii] = eta2theta(eta[,ii], .lexpectile, earg= .eexpectile)
        dimnames(ans) = list(dimnames(eta)[[1]], extra$y.names)
        ans
    }, list( .lexpectile=lexpectile, .eexpectile=eexpectile ))),
    last=eval(substitute(expression({
        misc$link = rep(.lexpectile, length=M)
        names(misc$link) = extra$y.names
        misc$earg = vector("list", M)
        names(misc$earg) = names(misc$link)
        misc$parallel = .parallel
        misc$expected = TRUE
        extra$percentile = numeric(M)
        for(ii in 1:M)
            extra$percentile[ii] = 100 * weighted.mean(myresid[,ii] <= 0, w)
        names(extra$percentile) = names(misc$link)

        extra$individual = TRUE
        extra$deviance = amlnormal.deviance(mu=mu, y=y, w=w, residuals=FALSE,
                                            eta=eta, extra=extra)
        names(extra$deviance) = extra$y.names
    }), list( .lexpectile=lexpectile, .eexpectile=eexpectile,
              .parallel=parallel ))),
    vfamily=c("amlnormal"),
    deriv=eval(substitute(expression({
        mymu = eta2theta(eta, .lexpectile, earg= .eexpectile)
        dexpectile.deta = dtheta.deta(mymu, .lexpectile, earg= .eexpectile)
        myresid = matrix(y,extra$n,extra$M) - cbind(mu)
        wor1 = Wr2(myresid, w= matrix(extra$w.als, extra$n, extra$M,
                                       byrow=TRUE))
        w * myresid * wor1 * dexpectile.deta
    }), list( .lexpectile=lexpectile, .eexpectile=eexpectile ))),
    weight=eval(substitute(expression({
        wz = w * wor1 * dexpectile.deta^2
        wz
    }), list( .lexpectile=lexpectile, .eexpectile=eexpectile ))))
}










amlpoisson.deviance = function(mu, y, w, residuals = FALSE, eta, extra=NULL) {

    M <- length(extra$w.aml)

    if(M > 1) y = matrix(y,extra$n,extra$M)

    nz <- y > 0
    devi =  cbind(-(y - mu))
    devi[nz] = devi[nz] + y[nz] * log(y[nz]/mu[nz])
    if(residuals) {
        stop("not sure here")
        return(sign(y - mu) * sqrt(2 * abs(devi) * w) *
               matrix(extra$w,extra$n,extra$M))
    } else {
        all.deviances = numeric(M)
        myresid = matrix(y,extra$n,extra$M) - cbind(mu)
        for(ii in 1:M) all.deviances[ii] = 2 * sum(w * devi[,ii] *
                               Wr1(myresid[,ii], w=extra$w.aml[ii]))
    }
    if(is.logical(extra$individual) && extra$individual)
        all.deviances else sum(all.deviances)
}


amlpoisson <- function(w.aml=1, parallel=FALSE, method.init=1, digw=4,
                       link="loge", earg=list())
{
    if(!is.Numeric(w.aml, posit=TRUE))
        stop("'w.aml' must be a vector of positive values")
    if(mode(link)!= "character" && mode(link)!= "name")
        link = as.character(substitute(link))
    if(!is.list(earg)) earg = list()

    new("vglmff",
        blurb=c("Poisson expectile regression by",
                " asymmetric maximum likelihood estimation\n\n",
           "Link:     ", namesof("expectile", link, earg= earg)),
    constraints=eval(substitute(expression({
        constraints = cm.vgam(matrix(1,M,1), x, .parallel, constraints)
    }), list( .parallel=parallel ))),
    deviance= function(mu, y, w, residuals = FALSE, eta, extra=NULL) {
        amlpoisson.deviance(mu=mu, y=y, w=w, residuals=residuals,
                            eta=eta, extra=extra)
    },
    initialize=eval(substitute(expression({
        extra$w.aml = .w.aml
        if(ncol(y <- cbind(y)) != 1)
            stop("response must be a vector or a one-column matrix")
        extra$M = M = length(extra$w.aml)  # Recycle if necessary
        extra$n = n
        extra$y.names = y.names = paste("w.aml=", round(extra$w.aml, dig=.digw),
                                        sep="")
        extra$individual = FALSE
        predictors.names = c(namesof(paste("expectile(",y.names,")", sep=""),
                                     .link, earg=.earg, tag=FALSE))

        if(!length(etastart)) {
            mean.init = if( .method.init == 2)
                    rep(median(y), length=n) else
                if( .method.init == 1)
                    rep(weighted.mean(y, w), length=n) else {
                        junk = if(is.R()) lm.wfit(x=x, y=y, w=w) else
                               lm.wfit(x=x, y=y, w=w, method="qr")
                        abs(junk$fitted)
                    }
            etastart = matrix(theta2eta(mean.init, .link, earg= .earg), n, M)
        }
    }), list( .link=link, .earg=earg, .method.init=method.init,
              .digw = digw, .w.aml=w.aml ))),
    inverse=eval(substitute(function(eta, extra=NULL) {
        mu.ans = eta = as.matrix(eta)
        for(ii in 1:ncol(eta))
            mu.ans[,ii] = eta2theta(eta[,ii], .link, earg= .earg)
        dimnames(mu.ans) = list(dimnames(eta)[[1]], extra$y.names)
        mu.ans
    }, list( .link=link, .earg=earg ))),
    last=eval(substitute(expression({
        misc$link = rep(.link, length=M)
        names(misc$link) = extra$y.names
        misc$earg = vector("list", M)
        names(misc$earg) = names(misc$link)
        misc$parallel = .parallel
        misc$expected = TRUE
        extra$percentile = numeric(M)
        for(ii in 1:M)
            extra$percentile[ii] = 100 * weighted.mean(myresid[,ii] <= 0, w)
        names(extra$percentile) = names(misc$link)

        extra$individual = TRUE
        extra$deviance = amlpoisson.deviance(mu=mu, y=y, w=w,
                         residuals=FALSE, eta=eta, extra=extra)
        names(extra$deviance) = extra$y.names
    }), list( .link=link, .earg=earg, .parallel=parallel ))),
    link=eval(substitute(function(mu, extra=NULL) {
        theta2eta(mu, link= .link, earg= .earg)
    }, list( .link=link, .earg=earg ))),
    vfamily=c("amlpoisson"),
    deriv=eval(substitute(expression({
        mymu = eta2theta(eta, .link, earg= .earg)
        dexpectile.deta = dtheta.deta(mymu, .link, earg=.earg)
        myresid = matrix(y,extra$n,extra$M) - cbind(mu)
        wor1 = Wr2(myresid, w= matrix(extra$w.aml, extra$n, extra$M,
                                       byrow=TRUE))
        w * myresid * wor1 * (dexpectile.deta / mymu)
    }), list( .link=link, .earg=earg ))),
    weight=eval(substitute(expression({
        use.mu = mymu
        use.mu[use.mu < .Machine$double.eps^(3/4)] = .Machine$double.eps^(3/4)
        wz = w * wor1 * use.mu * (dexpectile.deta / mymu)^2
        wz
    }), list( .link=link, .earg=earg ))))
}





amlbinomial.deviance = function(mu, y, w, residuals = FALSE, eta, extra=NULL) {

    M <- length(extra$w.aml)


    if(M > 1) y = matrix(y,extra$n,extra$M)


    devy <- y
    nz <- y != 0
    devy[nz] <- y[nz] * log(y[nz])
    nz <- (1 - y) != 0
    devy[nz] <- devy[nz] + (1 - y[nz]) * log1p(-y[nz])
    devmu <- y * log(mu) + (1 - y) * log1p(-mu)
    if(any(small <- mu * (1 - mu) < .Machine$double.eps)) {
        warning("fitted values close to 0 or 1")
        smu <- mu[small]
        sy <- y[small]
        smu <- ifelse(smu < .Machine$double.eps, .Machine$double.eps, smu)
        onemsmu <- ifelse((1 - smu) < .Machine$double.eps,
                          .Machine$double.eps, 1 - smu)
        devmu[small] <- sy * log(smu) + (1 - sy) * log(onemsmu)
    }
    devi <- 2 * (devy - devmu)
    if(residuals) {
        stop("not sure here")
        return(sign(y - mu) * sqrt(abs(devi) * w))
    } else {
        all.deviances = numeric(M)
        myresid = matrix(y,extra$n,extra$M) - matrix(mu,extra$n,extra$M)
        for(ii in 1:M) all.deviances[ii] = sum(w * devi[,ii] *
                               Wr1(myresid[,ii], w=extra$w.aml[ii]))
    }
    if(is.logical(extra$individual) && extra$individual)
        all.deviances else sum(all.deviances)
}


amlbinomial <- function(w.aml=1, parallel=FALSE, digw=4,
                       link="logit", earg=list())
{
    if(!is.Numeric(w.aml, posit=TRUE))
        stop("'w.aml' must be a vector of positive values")
    if(mode(link)!= "character" && mode(link)!= "name")
        link = as.character(substitute(link))
    if(!is.list(earg)) earg = list()

    new("vglmff",
        blurb=c("Logistic expectile regression by ",
                "asymmetric maximum likelihood estimation\n\n",
         "Link:     ", namesof("expectile", link, earg= earg)),
    constraints=eval(substitute(expression({
        constraints = cm.vgam(matrix(1,M,1), x, .parallel, constraints)
    }), list( .parallel=parallel ))),
    deviance= function(mu, y, w, residuals = FALSE, eta, extra=NULL) {
        amlbinomial.deviance(mu=mu, y=y, w=w, residuals=residuals,
                            eta=eta, extra=extra)
    },
    initialize=eval(substitute(expression({
        {
            NCOL = function (x)
                if(is.array(x) && length(dim(x)) > 1 ||
                is.data.frame(x)) ncol(x) else as.integer(1)

            if(NCOL(y) == 1) {
                if(is.factor(y)) y = y != levels(y)[1]
                nn = rep(1, n)
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
        }

        mustart = matrix(mustart, n, length( .w.aml ))

        extra$w.aml = .w.aml
        if(ncol(y <- cbind(y)) != 1)
            stop("response must be a vector or a one-column matrix")
        extra$M = M = length(extra$w.aml)  # Recycle if necessary
        extra$n = n
        extra$y.names = y.names = paste("w.aml=", round(extra$w.aml, dig=.digw),
                                        sep="")
        extra$individual = FALSE
        predictors.names = c(namesof(paste("expectile(",y.names,")", sep=""),
                                     .link, earg=.earg, tag=FALSE))

    }), list( .link=link, .earg=earg,
              .digw = digw, .w.aml=w.aml ))),
    inverse=eval(substitute(function(eta, extra=NULL) {
        mu.ans = eta = as.matrix(eta)
        for(ii in 1:ncol(eta))
            mu.ans[,ii] = eta2theta(eta[,ii], .link, earg= .earg)
        dimnames(mu.ans) = list(dimnames(eta)[[1]], extra$y.names)
        mu.ans
    }, list( .link=link, .earg=earg ))),
    last=eval(substitute(expression({
        misc$link = rep(.link, length=M)
        names(misc$link) = extra$y.names
        misc$earg = vector("list", M)
        names(misc$earg) = names(misc$link)
        misc$parallel = .parallel
        misc$expected = TRUE
        extra$percentile = numeric(M)
        for(ii in 1:M)
            extra$percentile[ii] = 100 * weighted.mean(myresid[,ii] <= 0, w)
        names(extra$percentile) = names(misc$link)

        extra$individual = TRUE
        extra$deviance = amlbinomial.deviance(mu=mu, y=y, w=w,
                         residuals=FALSE, eta=eta, extra=extra)
        names(extra$deviance) = extra$y.names
    }), list( .link=link, .earg=earg, .parallel=parallel ))),
    link=eval(substitute(function(mu, extra=NULL) {
        theta2eta(mu, link= .link, earg= .earg)
    }, list( .link=link, .earg=earg ))),
    vfamily=c("amlbinomial"),
    deriv=eval(substitute(expression({
        mymu = eta2theta(eta, .link, earg= .earg)
        use.mu = mymu
        use.mu[use.mu < .Machine$double.eps^(3/4)] = .Machine$double.eps^(3/4)
        dexpectile.deta = dtheta.deta(use.mu, .link, earg=.earg)
        myresid = matrix(y,extra$n,extra$M) - cbind(mu)
        wor1 = Wr2(myresid, w= matrix(extra$w.aml, extra$n, extra$M,
                                       byrow=TRUE))
        w * myresid * wor1 * (dexpectile.deta / (use.mu * (1-use.mu)))
    }), list( .link=link, .earg=earg ))),
    weight=eval(substitute(expression({
        wz = w * wor1 * (dexpectile.deta^2 / (use.mu * (1-use.mu)))
        wz
    }), list( .link=link, .earg=earg ))))
}










amlexponential.deviance = function(mu, y, w, residuals = FALSE, eta, extra=NULL) {

    M <- length(extra$w.aml)

    if(M > 1) y = matrix(y,extra$n,extra$M)

    devy =  cbind(-log(y) - 1)
    devi =  cbind(-log(mu) - y / mu)
    if(residuals) {
        stop("not sure here")
        return(sign(y - mu) * sqrt(2 * abs(devi) * w) *
               matrix(extra$w,extra$n,extra$M))
    } else {
        all.deviances = numeric(M)
        myresid = matrix(y,extra$n,extra$M) - cbind(mu)
        for(ii in 1:M) all.deviances[ii] = 2 * sum(w *
                               (devy[,ii] - devi[,ii]) *
                               Wr1(myresid[,ii], w=extra$w.aml[ii]))
    }
    if(is.logical(extra$individual) && extra$individual)
        all.deviances else sum(all.deviances)
}


amlexponential <- function(w.aml=1, parallel=FALSE, method.init=1, digw=4,
                           link="loge", earg=list())
{
    if(!is.Numeric(w.aml, posit=TRUE))
        stop("'w.aml' must be a vector of positive values")
    if(mode(link)!= "character" && mode(link)!= "name")
        link = as.character(substitute(link))
    if(!is.list(earg)) earg = list()
    if(!is.Numeric(method.init, allow=1, integ=TRUE, posit=TRUE) ||
       method.init > 3) stop("argument 'method.init' must be 1, 2 or 3")

    y.names = paste("w.aml=", round(w.aml, dig=digw), sep="")
    predictors.names = c(namesof(
        paste("expectile(", y.names,")", sep=""), link, earg=earg))
    predictors.names = paste(predictors.names, collapse=", ")

    new("vglmff",
        blurb=c("Exponential expectile regression by",
                " asymmetric maximum likelihood estimation\n\n",
           "Link:     ", predictors.names),
    constraints=eval(substitute(expression({
        constraints = cm.vgam(matrix(1,M,1), x, .parallel, constraints)
    }), list( .parallel=parallel ))),
    deviance= function(mu, y, w, residuals = FALSE, eta, extra=NULL) {
        amlexponential.deviance(mu=mu, y=y, w=w, residuals=residuals,
                            eta=eta, extra=extra)
    },
    initialize=eval(substitute(expression({
        extra$w.aml = .w.aml
        if(ncol(y <- cbind(y)) != 1)
            stop("response must be a vector or a one-column matrix")
        if(any(y <= 0.0))
            stop("all responses must be positive")
        extra$M = M = length(extra$w.aml)  # Recycle if necessary
        extra$n = n
        extra$y.names = y.names = paste("w.aml=", round(extra$w.aml, dig=.digw),
                                        sep="")
        extra$individual = FALSE
        predictors.names = c(namesof(
            paste("expectile(",y.names,")", sep=""), .link, earg=.earg, tag=FALSE))

        if(!length(etastart)) {
            mean.init = if( .method.init == 1)
                    rep(median(y), length=n) else
                if( .method.init == 2)
                    rep(weighted.mean(y, w), length=n) else {
                        1 / (y + 1)
                    }
            etastart = matrix(theta2eta(mean.init, .link, earg= .earg), n, M)
        }
    }), list( .link=link, .earg=earg, .method.init=method.init,
              .digw = digw, .w.aml=w.aml ))),
    inverse=eval(substitute(function(eta, extra=NULL) {
        mu.ans = eta = as.matrix(eta)
        for(ii in 1:ncol(eta))
            mu.ans[,ii] = eta2theta(eta[,ii], .link, earg= .earg)
        dimnames(mu.ans) = list(dimnames(eta)[[1]], extra$y.names)
        mu.ans
    }, list( .link=link, .earg=earg ))),
    last=eval(substitute(expression({
        misc$link = rep(.link, length=M)
        names(misc$link) = extra$y.names
        misc$earg = vector("list", M)
        names(misc$earg) = names(misc$link)
        misc$parallel = .parallel
        misc$expected = TRUE
        extra$percentile = numeric(M)
        for(ii in 1:M)
            extra$percentile[ii] = 100 * weighted.mean(myresid[,ii] <= 0, w)
        names(extra$percentile) = names(misc$link)

        extra$individual = TRUE
        extra$deviance = amlexponential.deviance(mu=mu, y=y, w=w,
                         residuals=FALSE, eta=eta, extra=extra)
        names(extra$deviance) = extra$y.names
    }), list( .link=link, .earg=earg, .parallel=parallel ))),
    link=eval(substitute(function(mu, extra=NULL) {
        theta2eta(mu, link= .link, earg= .earg)
    }, list( .link=link, .earg=earg ))),
    vfamily=c("amlexponential"),
    deriv=eval(substitute(expression({
        mymu = eta2theta(eta, .link, earg= .earg)
        bigy = matrix(y,extra$n,extra$M)
        dl.dmu = (bigy - mymu) / mymu^2
        dmu.deta = dtheta.deta(mymu, .link, earg=.earg)
        myresid = bigy - cbind(mymu)
        wor1 = Wr2(myresid, w= matrix(extra$w.aml, extra$n, extra$M,
                                       byrow=TRUE))
        w * wor1 * dl.dmu * dmu.deta
    }), list( .link=link, .earg=earg ))),
    weight=eval(substitute(expression({
        ned2l.dmu2 = 1 / mymu^2
        wz = w * wor1 * ned2l.dmu2 * dmu.deta^2
        wz
    }), list( .link=link, .earg=earg ))))
}






rho1check = function(u, tau=0.5)
    u * (tau - (u <= 0))

dalap = function(x, location=0, scale=1, tau=0.5,
                 kappa=sqrt(tau/(1-tau)), log=FALSE) {
    if(!is.logical(log.arg <- log)) stop("bad input for argument 'log'")
    rm(log)

    NN = max(length(x), length(location), length(scale), length(kappa))
    location = rep(location, len=NN); scale= rep(scale, len=NN)
    kappa = rep(kappa, len=NN); x = rep(x, len=NN)
    tau = rep(tau, len=NN)

    logconst = 0.5 * log(2) - log(scale) + log(kappa) - log1p(kappa^2)
    exponent = -(sqrt(2) / scale) * abs(x - location) *
               ifelse(x >= location, kappa, 1/kappa)

    indexTF = (scale > 0) & (tau > 0) & (tau < 1) & (kappa > 0) # &
    logconst[!indexTF] = NaN

    if(log.arg) logconst + exponent else exp(logconst + exponent)
}


ralap = function(n, location=0, scale=1, tau=0.5, kappa=sqrt(tau/(1-tau))) {
    use.n = if((length.n <- length(n)) > 1) length.n else
            if(!is.Numeric(n, integ=TRUE, allow=1, posit=TRUE))
                stop("bad input for argument 'n'") else n

    location = rep(location, len=use.n); scale= rep(scale, len=use.n)
    tau = rep(tau, len=use.n); kappa = rep(kappa, len=use.n);
    ans = location + scale *
          log(runif(use.n)^kappa / runif(use.n)^(1/kappa)) / sqrt(2)
    indexTF = (scale > 0) & (tau > 0) & (tau < 1) & (kappa > 0) # &
    ans[!indexTF] = NaN
    ans
}


palap = function(q, location=0, scale=1, tau=0.5, kappa=sqrt(tau/(1-tau))) {
    NN = max(length(q), length(location), length(scale), length(kappa))
    location = rep(location, len=NN); scale= rep(scale, len=NN)
    kappa = rep(kappa, len=NN); q= rep(q, len=NN)
    tau = rep(tau, len=NN);

    exponent = -(sqrt(2) / scale) * abs(q - location) *
               ifelse(q >= location, kappa, 1/kappa)
    temp5 = exp(exponent) / (1 + kappa^2)
    ans = 1 - temp5
    index1 = (q < location)
    ans[index1] = (kappa[index1])^2 * temp5[index1]

    indexTF = (scale > 0) & (tau > 0) & (tau < 1) & (kappa > 0) # &
    ans[!indexTF] = NaN
    ans
}


qalap = function(p, location=0, scale=1, tau=0.5, kappa=sqrt(tau/(1-tau))) {
    NN = max(length(p), length(location), length(scale), length(kappa))
    location = rep(location, len=NN); scale= rep(scale, len=NN)
    kappa = rep(kappa, len=NN); p = rep(p, len=NN)
    tau = rep(tau, len=NN)
    ans = p
    temp5 = kappa^2 / (1 + kappa^2)
    index1 = (p <= temp5)
    exponent = p[index1] / temp5[index1]
    ans[index1] = location[index1] + (scale[index1] * kappa[index1]) *
                  log(exponent) / sqrt(2)
    ans[!index1] = location[!index1] - (scale[!index1] / kappa[!index1]) *
                   (log1p((kappa[!index1])^2) + log1p(-p[!index1])) / sqrt(2)

    indexTF = (scale > 0) & (tau > 0) & (tau < 1) & (kappa > 0) &
              (p >= 0) & (p <= 1)
    ans[!indexTF] = NaN
    ans[p == 0 & indexTF] = -Inf
    ans[p == 1 & indexTF] = Inf
    ans
}




if(FALSE)
dqregal = function(x, tau=0.5, location=0, scale=1) {
    if(!is.Numeric(scale, posit=TRUE)) stop("'scale' must be positive")
    if(!is.Numeric(tau, posit=TRUE) || max(tau) >= 1)
        stop("'tau' must have values in (0,1)")
    const = tau * (1-tau) / scale
    const * exp(-rho1check((x-location)/scale, tau=tau))
}



if(FALSE)
rqregal = function(n, tau=0.5, location=0, scale=1) {
    if(!is.Numeric(n, posit=TRUE, integ=TRUE, allow=1))
        stop("bad input for argument 'n'")
    if(!is.Numeric(scale, posit=TRUE)) stop("'scale' must be positive")
    if(!is.Numeric(tau, posit=TRUE) || max(tau) >= 1)
        stop("'tau' must have values in (0,1)")
    location = rep(location, len=n); scale= rep(scale, len=n)
    r = runif(n)
    location - sign(r-tau) * scale * log(2*ifelse(r < tau, r, 1-r))
}



if(FALSE)
pqregal = function(q, tau=0.5, location=0, scale=1) {
    if(!all(scale == 1))
        stop("currently can only handle scale == 1")
    if(!is.Numeric(q))
        stop("bad input for argument 'q'")
    if(!is.Numeric(location))
        stop("bad input for argument 'location'")
    if(!is.Numeric(scale, posit=TRUE)) stop("'scale' must be positive")
    if(!is.Numeric(tau, posit=TRUE) || max(tau) >= 1)
        stop("'tau' must have values in (0,1)")
    N = max(length(q), length(tau), length(location), length(scale))
    location = rep(location, len=N); scale= rep(scale, len=N)
    tau = rep(tau, len=N); q= rep(q, len=N)
    ans = tau * exp(-(location - q) * (1 - tau))
    index1 = (q > location)
    ans[index1] = (1 - (1-tau) * exp(-tau * (q - location)))[index1]
    ans
}

if(FALSE)
qregal = function(tau=c(0.25, 0.5, 0.75),
                  llocation="identity",
                  elocation=list(),
                  lscale="loge", escale=list(),
                  ilocation=NULL,
                  parallel=FALSE, method.init=1, digt=4) {
    if(mode(llocation) != "character" && mode(llocation) != "name")
        llocation = as.character(substitute(llocation))
    if(!is.Numeric(method.init, allow=1, integ=TRUE, posit=TRUE) ||
       method.init > 2) stop("argument 'method.init' must be 1 or 2")
    if(!is.Numeric(tau, posit=TRUE) || max(tau) >= 1)
        stop("bad input for argument 'tau'")
    if(!is.list(elocation)) elocation = list()
    if(mode(lscale) != "character" && mode(lscale) != "name")
        lscale = as.character(substitute(lscale))
    if(!is.list(escale)) escale = list()

    new("vglmff",
    blurb=c("Quantile REGression via an Asymmetric Laplace distribution\n\n",
            "Links:    ",
            namesof("scale", lscale, earg=escale), ", ",
            namesof("location", llocation, earg=elocation)),
    constraints=eval(substitute(expression({
        constraints = cm.vgam(matrix(1,M,1), x, .parallel, constraints)
    }), list( .parallel=parallel ))),
    initialize=eval(substitute(expression({
        extra$tau = .tau
        if(ncol(y <- cbind(y)) != 1)
            stop("response must be a vector or a one-column matrix")
        extra$M = M = 1 + length(extra$tau)
        extra$n = n
        extra$y.names = y.names = paste("tau=", round(extra$tau, dig=.digt),
                                        sep="")
        extra$individual = FALSE
        predictors.names = c(
                  namesof("scale",    .lscale,    earg=.escale,    tag=FALSE),
                  namesof(paste("quantile(",y.names,")", sep=""),
                  link = .llocation, earg=.elocation, tag=FALSE))

        if(!length(etastart)) {
            if( .method.init == 1) {
                location.init = median(y)
            } else {
                location.init = y
            }
            location.init = if(length(.ilocation)) {
                matrix( .ilocation, n, M-1, byrow=TRUE)
            } else {
                rep(location.init, len=n)
            }
            scale.init = rep(1.0, len=n)
            etastart = cbind(
                theta2eta(scale.init,    .lscale, earg = .escale),
                matrix(
                theta2eta(location.init, .llocation, earg= .elocation), n, M-1))
        }
    }), list( .method.init=method.init, .tau=tau, .digt=digt,
              .elocation=elocation, .escale=escale,
              .llocation=llocation, .lscale=lscale,
              .ilocation=ilocation ))),
    inverse=eval(substitute(function(eta, extra=NULL) {
        eta = as.matrix(eta)
        xi.ans = matrix(0, nrow(eta), ncol(eta)-1)
        for(ii in 1:(ncol(eta)-1))
            xi.ans[,ii] = eta2theta(eta[,ii+1], .llocation, earg= .elocation)
        dimnames(xi.ans) = list(dimnames(eta)[[1]], extra$y.names)
        xi.ans
    }, list( .elocation=elocation, .llocation=llocation, .tau=tau,
             .escale=escale, .lscale=lscale ))),
    last=eval(substitute(expression({
        misc$link = rep( .llocation, length=M)
        names(misc$link) = extra$y.names
        misc$earg = vector("list", M)
        names(misc$earg) = names(misc$link)

        extra$percentile = numeric(M)
        for(ii in 1:M)
            extra$percentile[ii] = 100 *
                weighted.mean(ymat[,ii] - mu[,ii] <= 0, w)
        names(extra$percentile) = names(misc$link)

        misc$expected = TRUE
        misc$RegCondOK = FALSE # Save this for later
        misc$tau = .tau
    }), list( .elocation=elocation, .llocation=llocation, .tau=tau,
             .escale=escale, .lscale=lscale ))),
    loglikelihood=eval(substitute(
        function(mu,y,w,residuals= FALSE,eta, extra=NULL) {
        locmat = eta2theta(eta[,-1,drop=FALSE], .llocation, earg= .elocation)
        scalemat = matrix(eta2theta(eta[,1,drop=FALSE], .lscale,
                          earg= .escale), nrow=extra$n, ncol=extra$M - 1)
        taumat = matrix(extra$tau, nrow=extra$n, ncol=extra$M - 1, byrow=TRUE)
        ymat = matrix(y, nrow=extra$n, ncol=extra$M - 1)
        if(residuals) stop("loglikelihood residuals not implemented yet") else
        sum(w * (-log(scalemat) + log(taumat) + log1p(-taumat) -
                 rho1check((ymat-locmat)/scalemat, tau=taumat)))
    }, list( .elocation=elocation, .llocation=llocation,
             .escale=escale, .lscale=lscale, .tau=tau ))),
    vfamily=c("qregal"),
    deriv=eval(substitute(expression({
        ymat = matrix(y, nrow=extra$n, ncol=extra$M - 1)
        taumat = matrix(extra$tau, nrow=extra$n, ncol=extra$M - 1, byrow=TRUE)
        scalemat = matrix(eta2theta(eta[,1,drop=FALSE], .lscale,
                          earg= .escale), nrow=extra$n, ncol=extra$M - 1)
        locmat = eta2theta(eta[,-1,drop=FALSE], .llocation, earg= .elocation)
        dl.dlocation = taumat /scalemat
        index1 = (ymat < locmat)
        dl.dlocation[index1] = ((taumat - 1)/scalemat)[index1]
        dlocation.deta = dtheta.deta(locmat, .llocation, earg= .elocation)
        dscale.deta = dtheta.deta(scalemat, .lscale, earg= .escale)
        w * cbind(dl.dlocation * dlocation.deta)
    }), list( .tau=tau, .elocation=elocation, .llocation=llocation,
             .escale=escale, .lscale=lscale ))),
    weight=eval(substitute(expression({
        wz = matrix(0, nrow=n, M)  # Diagonal
        ed2l.dlocation2 = taumat * (1 - taumat) / scalemat^2
        ed2l.dscale2 = 2 * (3*taumat^2 - 3*taumat+1) / (scalemat^2 *
                       taumat * (1-taumat))
        wz[,iam(1,1,M)] = ed2l.dscale2 * dscale.deta^2
        wz[,-1] = ed2l.dlocation2 * dlocation.deta^2
        w * wz
    }), list( .tau=tau, .elocation=elocation, .llocation=llocation,
             .escale=escale, .lscale=lscale ))))
}





rloglap = function(n, location.ald=0, scale.ald=1, tau=0.5,
                       kappa=sqrt(tau/(1-tau))) {
    use.n = if((length.n <- length(n)) > 1) length.n else
            if(!is.Numeric(n, integ=TRUE, allow=1, posit=TRUE))
                stop("bad input for argument 'n'") else n
    location.ald = rep(location.ald, len=use.n);
    scale.ald= rep(scale.ald, len=use.n)
    tau = rep(tau, len=use.n);
    kappa = rep(kappa, len=use.n);
    ans = exp(location.ald) *
          (runif(use.n)^kappa / runif(use.n)^(1/kappa))^(scale.ald / sqrt(2))
    indexTF = (scale.ald > 0) & (tau > 0) & (tau < 1) & (kappa > 0) # &
    ans[!indexTF] = NaN
    ans
}


dloglap = function(x, location.ald=0, scale.ald=1, tau=0.5,
                       kappa=sqrt(tau/(1-tau)), log=FALSE) {
    if(!is.logical(log.arg <- log)) stop("bad input for argument 'log'")
    rm(log)

    NN = max(length(x), length(location.ald), length(scale.ald), length(kappa))
    location = rep(location.ald, len=NN); scale= rep(scale.ald, len=NN)
    kappa = rep(kappa, len=NN); x = rep(x, len=NN)
    tau = rep(tau, len=NN)

    Alpha = sqrt(2) * kappa / scale.ald
    Beta  = sqrt(2) / (scale.ald * kappa)
    Delta = exp(location.ald)
    exponent = ifelse(x >= Delta, -(Alpha+1), (Beta-1)) * (log(x) - location.ald)
    logdensity = -location.ald + log(Alpha) + log(Beta) -
                 log(Alpha + Beta) + exponent
    indexTF = (scale.ald > 0) & (tau > 0) & (tau < 1) & (kappa > 0) # &
    logdensity[!indexTF] = NaN
    logdensity[x <  0 & indexTF] = -Inf
    if(log.arg) logdensity else exp(logdensity)
}


qloglap = function(p, location.ald=0, scale.ald=1,
                       tau=0.5, kappa=sqrt(tau/(1-tau))) {
    NN = max(length(p), length(location.ald), length(scale.ald), length(kappa))
    location = rep(location.ald, len=NN); scale= rep(scale.ald, len=NN)
    kappa = rep(kappa, len=NN); p = rep(p, len=NN)
    tau = rep(tau, len=NN)

    Alpha = sqrt(2) * kappa / scale.ald
    Beta  = sqrt(2) / (scale.ald * kappa)
    Delta = exp(location.ald)

    temp9 = Alpha + Beta
    ans = Delta * (p * temp9 / Alpha)^(1/Beta)
    index1 = (p > Alpha / temp9)
    ans[index1] = (Delta * ((1-p) * temp9 / Beta)^(-1/Alpha))[index1]
    ans[p == 0] = 0
    ans[p == 1] = Inf

    indexTF = (scale.ald > 0) & (tau > 0) & (tau < 1) & (kappa > 0)
              (p >= 0) & (p <= 1) # &
    ans[!indexTF] = NaN
    ans
}



ploglap = function(q, location.ald=0, scale.ald=1,
                       tau=0.5, kappa=sqrt(tau/(1-tau))) {
    NN = max(length(q), length(location.ald), length(scale.ald), length(kappa))
    location = rep(location.ald, len=NN); scale = rep(scale.ald, len=NN)
    kappa = rep(kappa, len=NN); q = rep(q, len=NN)
    tau = rep(tau, len=NN)

    Alpha = sqrt(2) * kappa / scale.ald
    Beta  = sqrt(2) / (scale.ald * kappa)
    Delta = exp(location.ald)

    temp9 = Alpha + Beta
    ans = (Alpha / temp9) * (q / Delta)^(Beta)
    ans[q <= 0] = 0
    index1 = (q >= Delta)
    ans[index1] = (1 - (Beta/temp9) * (Delta/q)^(Alpha))[index1]

    indexTF = (scale.ald > 0) & (tau > 0) & (tau < 1) & (kappa > 0) # &
    ans[!indexTF] = NaN
    ans
}




rlogitlap = function(n, location.ald=0, scale.ald=1, tau=0.5,
                         kappa=sqrt(tau/(1-tau)), earg=list()) {
    logit(ralap(n=n, location=location.ald, scale=scale.ald,
                tau=tau, kappa=kappa), inverse=TRUE, earg=earg)
}


dlogitlap = function(x, location.ald=0, scale.ald=1, tau=0.5,
                         kappa=sqrt(tau/(1-tau)), log=FALSE, earg=list()) {
    if(!is.logical(log.arg <- log)) stop("bad input for argument 'log'")
    rm(log)

    NN = max(length(x), length(location.ald), length(scale.ald), length(kappa))
    location = rep(location.ald, len=NN); scale= rep(scale.ald, len=NN)
    kappa = rep(kappa, len=NN); x = rep(x, len=NN)
    tau = rep(tau, len=NN)

    Alpha = sqrt(2) * kappa / scale.ald
    Beta  = sqrt(2) / (scale.ald * kappa)
    Delta = logit(location.ald, inverse=TRUE, earg=earg)

    exponent = ifelse(x >= Delta, -Alpha, Beta) *
               (logit(x, earg=earg) - location.ald)
    logdensity = log(Alpha) + log(Beta) - log(Alpha + Beta) -
                 log(x) - log1p(-x) + exponent
    indexTF = (scale.ald > 0) & (tau > 0) & (tau < 1) & (kappa > 0) # &
    logdensity[!indexTF] = NaN
    logdensity[x <  0 & indexTF] = -Inf
    logdensity[x >  1 & indexTF] = -Inf
    if(log.arg) logdensity else exp(logdensity)
}


qlogitlap = function(p, location.ald=0, scale.ald=1,
                         tau=0.5, kappa=sqrt(tau/(1-tau)), earg=list()) {
    qqq = qalap(p=p, location=location.ald, scale=scale.ald,
                tau=tau, kappa=kappa)
    ans = logit(qqq, inverse=TRUE, earg=earg)
    ans[(p < 0) | (p > 1)] = NaN
    ans[p == 0] = 0
    ans[p == 1] = 1
    ans
}



plogitlap = function(q, location.ald=0, scale.ald=1,
                         tau=0.5, kappa=sqrt(tau/(1-tau)), earg=list()) {
    NN = max(length(q), length(location.ald), length(scale.ald),
             length(kappa))
    location.ald = rep(location.ald, len=NN); scale.ald= rep(scale.ald, len=NN)
    kappa = rep(kappa, len=NN); q= rep(q, len=NN)
    tau = rep(tau, len=NN);

    indexTF = (q > 0) & (q < 1)
    qqq = logit(q[indexTF], earg=earg)
    ans = q
    ans[indexTF] = palap(q=qqq, location=location.ald[indexTF],
                         scale=scale.ald[indexTF],
                         tau=tau[indexTF], kappa=kappa[indexTF])
    ans[q >= 1] = 1
    ans[q <= 0] = 0
    ans
}



rprobitlap = function(n, location.ald=0, scale.ald=1, tau=0.5,
                          kappa=sqrt(tau/(1-tau)), earg=list()) {
    probit(ralap(n=n, location=location.ald, scale=scale.ald,
                 tau=tau, kappa=kappa), inverse=TRUE, earg=earg)
}


dprobitlap = function(x, location.ald=0, scale.ald=1, tau=0.5,
                          kappa=sqrt(tau/(1-tau)), log=FALSE,
                          earg=list(), meth2=TRUE) {
    if(!is.logical(log.arg <- log)) stop("bad input for argument 'log'")
    rm(log)

    NN = max(length(x), length(location.ald), length(scale.ald), length(kappa))
    location.ald = rep(location.ald, len=NN); scale.ald= rep(scale.ald, len=NN)
    kappa = rep(kappa, len=NN); x = rep(x, len=NN)
    tau = rep(tau, len=NN)

    logdensity = x * NaN
    index1 = (x > 0) & (x < 1)
    indexTF = (scale.ald > 0) & (tau > 0) & (tau < 1) & (kappa > 0) # &
    if(meth2) {
        dx.dy = x
        use.x = probit(x[index1], earg=earg)
        logdensity[index1] = dalap(x=use.x, location=location.ald[index1],
                                   scale=scale.ald[index1], tau=tau[index1],
                                   kappa=kappa[index1], log=TRUE)
    } else {
        Alpha = sqrt(2) * kappa / scale.ald
        Beta  = sqrt(2) / (scale.ald * kappa)
        Delta = pnorm(location.ald)
        use.x  = qnorm(x) # qnorm(x[index1])
        log.dy.dw = dnorm(use.x, log=TRUE)

        exponent = ifelse(x >= Delta, -Alpha, Beta) * (use.x - location.ald) -
                   log.dy.dw

        logdensity[index1] = (log(Alpha) + log(Beta) -
                             log(Alpha + Beta) + exponent)[index1]
    }
    logdensity[!indexTF] = NaN
    logdensity[x <  0 & indexTF] = -Inf
    logdensity[x >  1 & indexTF] = -Inf

    if(meth2) {
        dx.dy[index1] = probit(x[index1], earg=earg, inverse=FALSE, deriv=1)
        dx.dy[!index1] = 0 # zz 0 seems to work. -Inf
        dx.dy[!indexTF] = NaN
        if(log.arg) logdensity - log(abs(dx.dy)) else exp(logdensity) / abs(dx.dy)
    } else {
        if(log.arg) logdensity else exp(logdensity)
    }
}


qprobitlap = function(p, location.ald=0, scale.ald=1,
                          tau=0.5, kappa=sqrt(tau/(1-tau)), earg=list()) {
    qqq = qalap(p=p, location=location.ald, scale=scale.ald,
                tau=tau, kappa=kappa)
    ans = probit(qqq, inverse=TRUE, earg=earg)
    ans[(p < 0) | (p > 1)] = NaN
    ans[p == 0] = 0
    ans[p == 1] = 1
    ans
}



pprobitlap = function(q, location.ald=0, scale.ald=1,
                          tau=0.5, kappa=sqrt(tau/(1-tau)), earg=list()) {
    NN = max(length(q), length(location.ald), length(scale.ald),
             length(kappa))
    location.ald = rep(location.ald, len=NN); scale.ald= rep(scale.ald, len=NN)
    kappa = rep(kappa, len=NN); q= rep(q, len=NN)
    tau = rep(tau, len=NN);

    indexTF = (q > 0) & (q < 1)
    qqq = probit(q[indexTF], earg=earg)
    ans = q
    ans[indexTF] = palap(q=qqq, location=location.ald[indexTF],
                         scale=scale.ald[indexTF],
                         tau=tau[indexTF], kappa=kappa[indexTF])
    ans[q >= 1] = 1
    ans[q <= 0] = 0
    ans
}


rclogloglap = function(n, location.ald=0, scale.ald=1, tau=0.5,
                          kappa=sqrt(tau/(1-tau)), earg=list()) {
    cloglog(ralap(n=n, location=location.ald, scale=scale.ald,
                  tau=tau, kappa=kappa), inverse=TRUE, earg=earg)
}


dclogloglap = function(x, location.ald=0, scale.ald=1, tau=0.5,
                           kappa=sqrt(tau/(1-tau)), log=FALSE,
                           earg=list(), meth2=TRUE) {
    if(!is.logical(log.arg <- log)) stop("bad input for argument 'log'")
    rm(log)

    NN = max(length(x), length(location.ald), length(scale.ald), length(kappa))
    location.ald = rep(location.ald, len=NN); scale.ald= rep(scale.ald, len=NN)
    kappa = rep(kappa, len=NN); x = rep(x, len=NN)
    tau = rep(tau, len=NN)

    logdensity = x * NaN
    index1 = (x > 0) & (x < 1)
    indexTF = (scale.ald > 0) & (tau > 0) & (tau < 1) & (kappa > 0) # &
    if(meth2) {
        dx.dy = x
        use.w = cloglog(x[index1], earg=earg)
        logdensity[index1] = dalap(x=use.w, location=location.ald[index1],
                                   scale=scale.ald[index1], tau=tau[index1],
                                   kappa=kappa[index1], log=TRUE)

    } else {
        Alpha = sqrt(2) * kappa / scale.ald
        Beta  = sqrt(2) / (scale.ald * kappa)
        Delta = cloglog(location.ald, inverse=TRUE)

        exponent = ifelse(x >= Delta, -(Alpha+1), Beta-1) * log(-log1p(-x)) +
                   ifelse(x >= Delta, Alpha, -Beta) * location.ald
        logdensity[index1] = (log(Alpha) + log(Beta) -
                             log(Alpha + Beta) - log1p(-x) + exponent)[index1]
    }
    logdensity[!indexTF] = NaN
    logdensity[x <  0 & indexTF] = -Inf
    logdensity[x >  1 & indexTF] = -Inf

    if(meth2) {
        dx.dy[index1] = cloglog(x[index1], earg=earg, inverse=FALSE, deriv=1)
        dx.dy[!index1] = 0 # zz 0 seems to work. -Inf
        dx.dy[!indexTF] = NaN
        if(log.arg) logdensity - log(abs(dx.dy)) else exp(logdensity) / abs(dx.dy)
    } else {
        if(log.arg) logdensity else exp(logdensity)
    }
}


qclogloglap = function(p, location.ald=0, scale.ald=1,
                          tau=0.5, kappa=sqrt(tau/(1-tau)), earg=list()) {
    qqq = qalap(p=p, location=location.ald, scale=scale.ald,
                tau=tau, kappa=kappa)
    ans = cloglog(qqq, inverse=TRUE, earg=earg)
    ans[(p < 0) | (p > 1)] = NaN
    ans[p == 0] = 0
    ans[p == 1] = 1
    ans
}



pclogloglap = function(q, location.ald=0, scale.ald=1,
                           tau=0.5, kappa=sqrt(tau/(1-tau)), earg=list()) {
    NN = max(length(q), length(location.ald), length(scale.ald),
             length(kappa))
    location.ald = rep(location.ald, len=NN); scale.ald= rep(scale.ald, len=NN)
    kappa = rep(kappa, len=NN); q= rep(q, len=NN)
    tau = rep(tau, len=NN);

    indexTF = (q > 0) & (q < 1)
    qqq = cloglog(q[indexTF], earg=earg)
    ans = q
    ans[indexTF] = palap(q=qqq, location=location.ald[indexTF],
                         scale=scale.ald[indexTF],
                         tau=tau[indexTF], kappa=kappa[indexTF])
    ans[q >= 1] = 1
    ans[q <= 0] = 0
    ans
}


















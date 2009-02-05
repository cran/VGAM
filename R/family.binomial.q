# These functions are
# Copyright (C) 1998-2009 T.W. Yee, University of Auckland. All rights reserved.



process.binomial2.data.vgam <- expression({


    if(!is.matrix(y)) {
        yf <- as.factor(y)
        lev <- levels(yf)
        llev <- length(lev)
        if(llev != 4)
            stop("response must have 4 levels")
        nn <- length(yf)
        y <- matrix(0, nn, llev)
        y[cbind(1:nn,as.vector(unclass(yf)))] <- 1
        colnamesy <- paste(lev, ":", c("00","01","10","11"), sep="")
        dimnames(y) <- list(names(yf), colnamesy)
        input.type <- 1
    } else if(ncol(y)==2) {
        if(!all(y==0 | y==1))
            stop("response must contains 0's and 1's only")
        col.index <- y[,2] + 2*y[,1] + 1    # 1:4
        nn <- nrow(y)
        y <- matrix(0, nn, 4)
        y[cbind(1:nn,col.index)] <- 1
        dimnames(y) <- list(dimnames(y)[[1]], c("00","01","10","11"))
        input.type <- 2
    } else if(ncol(y)==4) {
        input.type <- 3
    } else
        stop("response unrecognized")


    nvec <- drop(y %*% rep(1,ncol(y)))

    w <- w * nvec
    y <- y / nvec             # Convert to proportions

    mu <- y + (1/ncol(y) - y)/nvec
    dimnames(mu) <- dimnames(y)

})










betabinomial <- function(lmu="logit", lrho="logit",
                         emu=list(), erho=list(),
                         irho=NULL, method.init=1, zero=2)
{
    if(mode(lmu) != "character" && mode(lmu) != "name")
        lmu = as.character(substitute(lmu))
    if(mode(lrho) != "character" && mode(lrho) != "name")
        lrho = as.character(substitute(lrho))
    if(!is.list(emu )) emu  = list()
    if(!is.list(erho)) erho = list()
    if(!is.Numeric(method.init, allow=1, integ=TRUE, posit=TRUE) ||
       method.init > 2) stop("argument \"method.init\" must be 1 or 2")

    new("vglmff",
    blurb=c("Beta-binomial model\n",
           "Links:      ",
           namesof("mu", lmu, earg= emu), ", ",
           namesof("rho", lrho, earg= erho), "\n",
           "Variance:   mu*(1-mu)*(1+(w-1)*rho)/w"),
    constraints=eval(substitute(expression({
        constraints <- cm.zero.vgam(constraints, x, .zero, M)
    }), list( .zero=zero ))),
    initialize=eval(substitute(expression({
        eval(binomialff()@initialize)   # Note: n,w,y,mustart is changed 
        ycounts = y * w   # Convert proportions to counts
        if(max(abs(ycounts-round(ycounts))) > 1.0e-6)
           stop("the response (as counts) does not appear to be integer-valued")
        predictors.names = c(namesof("mu",  .lmu,  earg= .emu,  tag=FALSE),
                             namesof("rho", .lrho, earg= .erho, tag=FALSE))
        if(!length(etastart)) {
            if(is.Numeric( .irho )) {
                init.rho = rep( .irho, length=n)
            } else {
                betabinomial.Loglikfun = function(rhoval, y, x, w, extraargs) {
                    shape1 = extraargs$mustart*(1-rhoval)/rhoval
                    shape2 = (1-extraargs$mustart)*(1-rhoval)/rhoval
                    ycounts = extraargs$ycounts
                    nvec = extraargs$nvec
                    if(is.R()) sum(lbeta(shape1+ycounts, shape2+nvec-ycounts) -
                                   lbeta(shape1, shape2)) else
                    sum(lgamma(shape1+ycounts) + lgamma(shape2+nvec-ycounts) -
                        lgamma(shape1+shape2+nvec) -
                        (lgamma(shape1) + lgamma(shape2) -
                         lgamma(shape1+shape2)))
                }
                rho.grid = rvar = seq(0.05, 0.95, len=21)  # 
                mustart.use = if( .method.init == 2) {
                    mustart
                } else {
                    y.matrix = cbind(y)
                    mat.temp = matrix(apply(y.matrix, 2, mean), nrow(y.matrix),
                                      ncol(y.matrix), byrow=TRUE)
                    0.5 * mustart + 0.5 * mat.temp
                }
                try.this = getMaxMin(rho.grid, objfun=betabinomial.Loglikfun,
                                     y=y,  x=x, w=w, extraargs=list(
                                     ycounts=ycounts, nvec=w,
                                     mustart=mustart.use))
                init.rho = rep(try.this, len=n)
            }
            etastart = cbind(theta2eta(mustart.use,  .lmu, earg= .emu),
                             theta2eta(init.rho,     .lrho, earg= .erho))
          }
    }), list( .lmu=lmu, .lrho=lrho,
              .emu=emu, .erho=erho,
              .method.init=method.init,
              .irho=irho ))),
    inverse=eval(substitute(function(eta, extra=NULL)
        eta2theta(eta[,1], .lmu, earg= .emu), 
    list( .lmu=lmu, .emu=emu ))),
    last=eval(substitute(expression({
        misc$link <- c(mu = .lmu, rho = .lrho)
        misc$earg <- list(mu = .emu, rho = .erho)
        misc$zero <- .zero
        misc$expected <- TRUE
    }), list( .lmu=lmu, .lrho=lrho,
              .emu=emu, .erho=erho,
              .zero=zero ))),
    loglikelihood=eval(substitute(
        function(mu,y,w,residuals=FALSE, eta, extra=NULL) {
        ycounts = y * w   # Convert proportions to counts
        mymu = eta2theta(eta[,1], .lmu, earg= .emu)
        rho  = eta2theta(eta[,2], .lrho, earg= .erho)
        smallno = 100 * .Machine$double.eps
        rho  = pmax(rho, smallno)
        rho  = pmin(rho, 1-smallno)
        shape1 = mymu * (1 - rho) / rho
        shape2 = (1-mymu) * (1 - rho) / rho
        nvec = w
        if(residuals) stop("loglikelihood residuals not implemented yet") else {
            if(is.R()) sum(lbeta(shape1+ycounts, shape2+nvec-ycounts) -
                           lbeta(shape1, shape2)) else
            sum(lgamma(shape1+ycounts) + lgamma(shape2+nvec-ycounts) -
                lgamma(shape1+shape2+nvec) -
                (lgamma(shape1) + lgamma(shape2) - lgamma(shape1+shape2)))
        }
    }, list( .lmu=lmu,
             .emu=emu, .erho=erho,
             .lrho=lrho ))),
    vfamily=c("betabinomial"),
    deriv=eval(substitute(expression({
        nvec = w  # extra$nvec # for summary()
        ycounts = y * w   # Convert proportions to counts
        mymu = eta2theta(eta[,1], .lmu, earg= .emu)
        rho  = eta2theta(eta[,2], .lrho, earg= .erho)
        smallno = 100 * .Machine$double.eps
        rho  = pmax(rho, smallno)
        rho  = pmin(rho, 1-smallno)
        shape1 = mymu * (1 - rho) / rho
        shape2 = (1-mymu) * (1 - rho) / rho
        dshape1.dmu =  (1 - rho) / rho
        dshape2.dmu = -(1 - rho) / rho
        dshape1.drho = -mymu / rho^2
        dshape2.drho =  -(1 - mymu) / rho^2
        dmu.deta  = dtheta.deta(mymu, .lmu, earg= .emu)
        drho.deta = dtheta.deta(rho,  .lrho, earg= .erho)
        dl.dmu = dshape1.dmu * (digamma(shape1+ycounts) -
                     digamma(shape2+nvec-ycounts) -
                     digamma(shape1) + digamma(shape2))
        dl.drho = (-1/rho^2) * (mymu * digamma(shape1+ycounts) +
                     (1-mymu) * digamma(shape2+nvec-ycounts) -
                     digamma(shape1+shape2+nvec) - 
                     mymu * digamma(shape1) -
                     (1-mymu)*digamma(shape2) + digamma(shape1+shape2))
        temp5 = cbind(dl.dmu * dmu.deta, dl.drho * drho.deta)
        temp5
    }), list( .lmu=lmu,
              .emu=emu, .erho=erho,
              .lrho=lrho ))),
    weight=eval(substitute(expression({
        wz = matrix(as.numeric(NA), n, dimm(M))  #3=dimm(2)
        wz11 = -(expected.betabin.ab(nvec, shape1, shape2, TRUE) -
                            trigamma(shape1+shape2+nvec) -
                            trigamma(shape1) + trigamma(shape1+shape2))
        wz22 = -(expected.betabin.ab(nvec, shape1, shape2, FALSE) -
                            trigamma(shape1+shape2+nvec) -
                            trigamma(shape2) + trigamma(shape1+shape2))
        wz21 = -(trigamma(shape1+shape2) - trigamma(shape1+shape2+nvec))
        wz[,iam(1,1,M)] = dmu.deta^2 * (wz11 * dshape1.dmu^2 +
                                        wz22 * dshape2.dmu^2 +
                           2 * wz21 * dshape1.dmu * dshape2.dmu)
        wz[,iam(2,2,M)] = drho.deta^2 * (wz11 * dshape1.drho^2 +
                                         wz22 * dshape2.drho^2 +
                           2 * wz21 * dshape1.drho * dshape2.drho)
        wz[,iam(2,1,M)] = dmu.deta * drho.deta *
                          (dshape1.dmu*(wz11*dshape1.drho + wz21*dshape2.drho) +
                           dshape2.dmu*(wz21*dshape1.drho + wz22*dshape2.drho))
        wz
    }), list( .lmu=lmu,
              .emu=emu, .erho=erho,
              .lrho=lrho ))))
}






dbinom2.or = function(mu1,
                      mu2=if(exchangeable) mu1 else stop("'mu2' not specified"),
                      oratio=1,
                      exchangeable=FALSE,
                      tol=0.001,
                      colnames=c("00", "01", "10", "11"),
                      ErrorCheck=TRUE)
{
    if(ErrorCheck) {
        if(!is.Numeric(mu1, positive=TRUE) || max(mu1) >= 1)
            stop("bad input for argument 'mu1'") 
        if(!is.Numeric(mu2, positive=TRUE) || max(mu2) >= 1)
            stop("bad input for argument 'mu2'") 
        if(!is.Numeric(oratio, positive=TRUE))
            stop("bad input for argument 'oratio'") 
        if(!is.Numeric(tol, positive=TRUE, allow=1) || tol > 0.1)
            stop("bad input for argument \"tol\"") 
        if(exchangeable && max(abs(mu1 - mu2)) > 0.00001)
            stop("argument 'exchangeable' is TRUE but 'mu1' and 'mu2' differ") 
    }

    n = max(length(mu1), length(mu2), length(oratio))
    oratio = rep(oratio, len=n)
    mu1 = rep(mu1, len=n)
    mu2 = rep(mu2, len=n)

    a.temp = 1 + (mu1+mu2)*(oratio-1)
    b.temp = -4 * oratio * (oratio-1) * mu1 * mu2
    temp = sqrt(a.temp^2 + b.temp)
    p11 = ifelse(abs(oratio-1) < tol, mu1*mu2, (a.temp-temp)/(2*(oratio-1)))
    p01 = mu2 - p11
    p10 = mu1 - p11
    p00 = 1 - p11 - p01 - p10
    matrix(c(p00,p01,p10,p11), n, 4, dimnames=list(NULL,colnames))
}




rbinom2.or = function(n, mu1,
                      mu2=if(exchangeable) mu1 else stop("'mu2' not specified"),
                      oratio=1,
                      exchangeable=FALSE,
                      tol=0.001,
                      twoCols=TRUE,
            colnames=if(twoCols) c("y1","y2") else c("00", "01", "10", "11"),
                      ErrorCheck=TRUE)
{
    if(ErrorCheck) {
        if(!is.Numeric(n, integer=TRUE, posit=TRUE, allow=1))
            stop("bad input for argument 'n'")
        if(!is.Numeric(mu1, positive=TRUE) || max(mu1) >= 1)
            stop("bad input for argument 'mu1'") 
        if(!is.Numeric(mu2, positive=TRUE) || max(mu2) >= 1)
            stop("bad input for argument 'mu2'") 
        if(!is.Numeric(oratio, positive=TRUE))
            stop("bad input for argument 'oratio'") 
        if(!is.Numeric(tol, positive=TRUE, allow=1) || tol > 0.1)
            stop("bad input for argument \"tol\"") 
        if(exchangeable && max(abs(mu1 - mu2)) > 0.00001)
            stop("argument 'exchangeable' is TRUE but 'mu1' and 'mu2' differ") 
    }

    dmat = dbinom2.or(mu1=mu1, mu2=mu2, oratio=oratio, exchang=exchangeable,
                      tol=tol, ErrorCheck=ErrorCheck)

    answer = matrix(0, n, 2, dimnames=list(NULL, if(twoCols) colnames else NULL))
    yy = runif(n)
    cs1 = dmat[,"00"] + dmat[,"01"]
    cs2 = cs1 + dmat[,"10"]
    index = (dmat[,"00"] < yy) & (yy <= cs1)
    answer[index,2] = 1
    index = (cs1 < yy) & (yy <= cs2)
    answer[index,1] = 1
    index = (yy > cs2)
    answer[index,] = 1
    if(twoCols) answer else {
        answer4 = matrix(0, n, 4, dimnames=list(NULL, colnames))
        answer4[cbind(1:n, 1 + 2*answer[,1] + answer[,2])] = 1
        answer4
    }
}




binom2.or = function(lmu="logit", lmu1=lmu, lmu2=lmu, loratio="loge",
                     emu=list(), emu1=emu, emu2=emu, eoratio=list(),
                     imu1=NULL, imu2=NULL, ioratio = NULL,
                     zero=3, exchangeable=FALSE, tol=0.001, morerobust=FALSE)
{
    if(mode(lmu) != "character" && mode(lmu) != "name")
        lmu = as.character(substitute(lmu))
    if(mode(lmu1) != "character" && mode(lmu1) != "name")
        lmu1 = as.character(substitute(lmu1))
    if(mode(lmu2) != "character" && mode(lmu2) != "name")
        lmu2 = as.character(substitute(lmu2))
    if(mode(loratio) != "character" && mode(loratio) != "name")
        loratio = as.character(substitute(loratio))
    if(is.logical(exchangeable) && exchangeable && ((lmu1 != lmu2) ||
       !all.equal(emu1, emu2)))
        stop("exchangeable=TRUE but marginal links are not equal") 
    if(!is.Numeric(tol, positive=TRUE, allow=1) || tol > 0.1)
        stop("bad input for argument \"tol\"") 
    if(!is.list(emu1)) emu1  = list()
    if(!is.list(emu2)) emu2  = list()
    if(!is.list(eoratio)) eoratio = list()

    new("vglmff",
    blurb=c("Bivariate binomial regression with an odds ratio\n",
            "Links:    ",
            namesof("mu1", lmu1, earg=emu1), ", ",
            namesof("mu2", lmu2, earg=emu2), "; ",
            namesof("oratio", loratio, earg=eoratio)),
    constraints=eval(substitute(expression({
        constraints = cm.vgam(matrix(c(1,1,0,0,0,1),3,2), x, 
                              .exchangeable, constraints,
                              intercept.apply=TRUE)
        constraints = cm.zero.vgam(constraints, x, .zero, M)
    }), list( .exchangeable=exchangeable, .zero=zero ))),
    deviance=Deviance.categorical.data.vgam,
    initialize=eval(substitute(expression({
        eval(process.binomial2.data.vgam)
        predictors.names = c(namesof("mu1", .lmu1, earg= .emu1, short=TRUE), 
                 namesof("mu2", .lmu2, earg= .emu2, short=TRUE), 
                 namesof("oratio",  .loratio, earg= .eoratio, short=TRUE))

        if(!length(etastart)) {
            pmargin = cbind(mu[,3]+mu[,4], mu[,2]+mu[,4])
            ioratio = if(length( .ioratio)) rep( .ioratio, len=n) else
                      mu[,4]*mu[,1]/(mu[,2]*mu[,3])
            if(length( .imu1)) pmargin[,1] = .imu1
            if(length( .imu2)) pmargin[,2] = .imu2
            etastart = cbind(theta2eta(pmargin[,1], .lmu1, earg= .emu1),
                             theta2eta(pmargin[,2], .lmu2, earg= .emu2), 
                             theta2eta(ioratio, .loratio, earg= .eoratio))
        }
    }), list( .lmu1=lmu1, .lmu2=lmu2, .loratio=loratio,
              .emu1=emu1, .emu2=emu2, .eoratio=eoratio,
              .imu1=imu1, .imu2=imu2, .ioratio=ioratio ))),
    inverse=eval(substitute(function(eta, extra=NULL) {
        pmargin = cbind(eta2theta(eta[,1], .lmu1, earg= .emu1),
                        eta2theta(eta[,2], .lmu2, earg= .emu2))
        oratio = eta2theta(eta[,3], .loratio, earg= .eoratio)
        a.temp = 1 + (pmargin[,1]+pmargin[,2])*(oratio-1)
        b.temp = -4 * oratio * (oratio-1) * pmargin[,1] * pmargin[,2]
        temp = sqrt(a.temp^2 + b.temp)
        pj4 = ifelse(abs(oratio-1) < .tol, pmargin[,1]*pmargin[,2],
                     (a.temp-temp)/(2*(oratio-1)))
        pj2 = pmargin[,2] - pj4
        pj3 = pmargin[,1] - pj4
        cbind("00" = 1-pj4-pj2-pj3, "01" = pj2, "10" = pj3, "11" = pj4)
    }, list( .tol=tol, .lmu1=lmu1, .lmu2=lmu2,
             .emu1=emu1, .emu2=emu2, .eoratio=eoratio,
             .loratio=loratio ))),
    last=eval(substitute(expression({
        misc$link = c("mu1"= .lmu1, "mu2"= .lmu2, "oratio"= .loratio)
        misc$earg = list(mu1 = .emu1, mu2 = .emu2, oratio = .eoratio)
        misc$tol = .tol
        misc$expected = TRUE
    }), list( .tol=tol, .lmu1=lmu1, .lmu2=lmu2,
              .emu1=emu1, .emu2=emu2, .eoratio=eoratio,
              .loratio=loratio ))),
    link=eval(substitute(function(mu, extra=NULL) {
        pmargin = cbind(mu[,3]+mu[,4], mu[,2]+mu[,4])
        oratio = mu[,4]*mu[,1] / (mu[,2]*mu[,3])
        cbind(theta2eta(pmargin[,1], .lmu1, earg= .emu1),
              theta2eta(pmargin[,2], .lmu2, earg= .emu2), 
              theta2eta(oratio, .loratio, earg= .eoratio))
    }, list( .lmu1=lmu1, .lmu2=lmu2,
             .emu1=emu1, .emu2=emu2, .eoratio=eoratio,
             .loratio=loratio ))),
    loglikelihood=eval(substitute(
        function(mu, y, w, residuals = FALSE, eta, extra=NULL) {
            if(residuals) stop("loglikelihood residuals not implemented yet") else {
                if( .morerobust) {
                    vsmallno =  1.0e4 * .Machine$double.xmin
                    mu.use = mu
                    mu.use[mu.use < vsmallno] = vsmallno
                    sum(w * y * log(mu.use))
                } else
                    sum(w * y * log(mu))
        }
        }, list( .morerobust=morerobust ))),
    vfamily=c("binom2.or", "binom2"),
    deriv=eval(substitute(expression({
        smallno = 1.0e4 * .Machine$double.eps
iii = c(46,55,63)
iii = c(39)
iii = 1:n
        mu.use = mu
        mu.use[mu.use < smallno] = smallno
        mu.use[mu.use > 1-smallno] = 1-smallno
        pmargin = cbind(mu.use[,3]+mu.use[,4], mu.use[,2]+mu.use[,4])
        pmargin[,1] = pmax(smallno,   pmargin[,1])
        pmargin[,1] = pmin(1-smallno, pmargin[,1])
        pmargin[,2] = pmax(smallno,   pmargin[,2])
        pmargin[,2] = pmin(1-smallno, pmargin[,2])

        oratio = mu.use[,4]*mu.use[,1] / (mu.use[,2]*mu.use[,3])
        use.oratio = pmax(smallno, oratio)
        a.temp = 1 + (pmargin[,1]+pmargin[,2])*(oratio-1)
        b.temp = -4 * oratio * (oratio-1) * pmargin[,1] * pmargin[,2]
        temp9 = sqrt(a.temp^2 + b.temp)

        coeff12 = -0.5 + (2*oratio*pmargin - a.temp) / (2*temp9)
        dl.dmu1 = coeff12[,2] * (y[,1]/mu.use[,1]-y[,3]/mu.use[,3]) -
           (1+coeff12[,2]) * (y[,2]/mu.use[,2]-y[,4]/mu.use[,4])
    
        dl.dmu2 = coeff12[,1] * (y[,1]/mu.use[,1]-y[,2]/mu.use[,2]) -
           (1+coeff12[,1]) * (y[,3]/mu.use[,3]-y[,4]/mu.use[,4])
    
        coeff3 = (y[,1]/mu.use[,1] - y[,2]/mu.use[,2] - y[,3]/mu.use[,3] + y[,4]/mu.use[,4])
        Vab = pmax(smallno, 1 / (1/mu.use[,1] + 1/mu.use[,2] + 1/mu.use[,3] + 1/mu.use[,4]))
        dp11.doratio = Vab / use.oratio
        dl.doratio = coeff3 * dp11.doratio

        w * cbind(dl.dmu1 * dtheta.deta(pmargin[,1], .lmu1, earg= .emu1),
                  dl.dmu2 * dtheta.deta(pmargin[,2], .lmu2, earg= .emu2),
                  dl.doratio * dtheta.deta(oratio, .loratio, earg= .eoratio))
    }), list( .lmu1=lmu1, .lmu2=lmu2,
              .emu1=emu1, .emu2=emu2, .eoratio=eoratio,
              .loratio=loratio ))),
    weight=eval(substitute(expression({
        Deltapi = mu.use[,3]*mu.use[,2] - mu.use[,4]*mu.use[,1]
        myDelta  = pmax(smallno, mu.use[,1] * mu.use[,2] * mu.use[,3] * mu.use[,4])
        pqmargin = pmargin * (1-pmargin)
        pqmargin[pqmargin < smallno] = smallno

        wz = matrix(0, n, 4)
        wz[,iam(1,1,M)] = (pqmargin[,2] * Vab / myDelta) *
                          dtheta.deta(pmargin[,1], .lmu1, earg= .emu1)^2
        wz[,iam(2,2,M)] = (pqmargin[,1] * Vab / myDelta) *
                          dtheta.deta(pmargin[,2], .lmu2, earg= .emu2)^2
        wz[,iam(3,3,M)] = (Vab / use.oratio^2) *
                          dtheta.deta(use.oratio, .loratio, earg= .eoratio)^2
        wz[,iam(1,2,M)] = (Vab * Deltapi / myDelta) *
                          dtheta.deta(pmargin[,1], .lmu1, earg= .emu1) *
                          dtheta.deta(pmargin[,2], .lmu2, earg= .emu2)
        w * wz
    }), list( .lmu1=lmu1, .lmu2=lmu2,
              .emu1=emu1, .emu2=emu2, .eoratio=eoratio,
              .loratio=loratio ))))
}


dbinom2.rho = function(mu1,
                      mu2=if(exchangeable) mu1 else stop("'mu2' not specified"),
                       rho=0,
                       exchangeable=FALSE,
                       colnames=c("00", "01", "10", "11"),
                       ErrorCheck=TRUE)
{
    if(ErrorCheck) {
        if(!is.Numeric(mu1, positive=TRUE) || max(mu1) >= 1)
            stop("bad input for argument 'mu1'") 
        if(!is.Numeric(mu2, positive=TRUE) || max(mu2) >= 1)
            stop("bad input for argument 'mu2'") 
        if(!is.Numeric(rho) || min(rho) <= -1 || max(rho) >= 1)
            stop("bad input for argument 'rho'") 
        if(exchangeable && max(abs(mu1 - mu2)) > 0.00001)
            stop("argument 'exchangeable' is TRUE but 'mu1' and 'mu2' differ") 
    }

    n = max(length(mu1), length(mu2), length(rho))
    rho = rep(rho, len=n)
    mu1 = rep(mu1, len=n)
    mu2 = rep(mu2, len=n)
    eta1 = qnorm(mu1)
    eta2 = qnorm(mu2)
    p11 = pnorm2(eta1, eta2, rho)
    p01 = mu2 - p11
    p10 = mu1 - p11
    p00 = 1 - p01 - p10 - p11
    matrix(c(p00,p01,p10,p11), n, 4, dimnames=list(NULL,colnames))
}



rbinom2.rho = function(n, mu1,
                      mu2=if(exchangeable) mu1 else stop("'mu2' not specified"),
                       rho=0,
                       exchangeable=FALSE,
                       twoCols=TRUE,
             colnames=if(twoCols) c("y1","y2") else c("00", "01", "10", "11"),
                       ErrorCheck=TRUE)
{
    if(ErrorCheck) {
        if(!is.Numeric(n, integer=TRUE, posit=TRUE, allow=1))
            stop("bad input for argument 'n'")
        if(!is.Numeric(mu1, positive=TRUE) || max(mu1) >= 1)
            stop("bad input for argument 'mu1'") 
        if(!is.Numeric(mu2, positive=TRUE) || max(mu2) >= 1)
            stop("bad input for argument 'mu2'") 
        if(!is.Numeric(rho) || min(rho) <= -1 || max(rho) >= 1)
            stop("bad input for argument 'rho'") 
        if(exchangeable && max(abs(mu1 - mu2)) > 0.00001)
            stop("argument 'exchangeable' is TRUE but 'mu1' and 'mu2' differ") 
    }

    dmat = dbinom2.rho(mu1=mu1, mu2=mu2, rho=rho, exchang=exchangeable,
                       ErrorCheck=ErrorCheck)

    answer = matrix(0, n, 2, dimnames=list(NULL, if(twoCols) colnames else NULL))
    yy = runif(n)
    cs1 = dmat[,"00"] + dmat[,"01"]
    cs2 = cs1 + dmat[,"10"]
    index = (dmat[,"00"] < yy) & (yy <= cs1)
    answer[index,2] = 1
    index = (cs1 < yy) & (yy <= cs2)
    answer[index,1] = 1
    index = (yy > cs2)
    answer[index,] = 1
    if(twoCols) answer else {
        answer4 = matrix(0, n, 4, dimnames=list(NULL, colnames))
        answer4[cbind(1:n, 1 + 2*answer[,1] + answer[,2])] = 1
        answer4
    }
}



binom2.rho.control <- function(save.weight=TRUE, ...)
{
    list(save.weight=save.weight)
}



binom2.rho = function(lrho="rhobit", erho=list(),
                      imu1=NULL, imu2=NULL, 
                      init.rho=NULL,
                      zero=3, exchangeable=FALSE, nsimEIM=NULL)
{

    if(mode(lrho) != "character" && mode(lrho) != "name")
        lrho = as.character(substitute(lrho))
    if(!is.list(erho)) erho = list()
    lmu12 = "probit"
    emu12 = list()
    if(is.Numeric(nsimEIM)) {
        if(!is.Numeric(nsimEIM, allow=1, integ=TRUE))
            stop("bad input for argument 'nsimEIM'")
        if(nsimEIM <= 100)
            warning("'nsimEIM' should be an integer greater than 100")
    }

    new("vglmff",
    blurb=c("Bivariate probit model\n",
           "Links:    ",
            namesof("mu1", lmu12, earg= emu12), ", ",
            namesof("mu2", lmu12, earg= emu12), ", ",
            namesof("rho", lrho, earg= erho)),
    constraints=eval(substitute(expression({
        constraints = cm.vgam(matrix(c(1,1,0,0,0,1),3,2), x, 
                              .exchangeable, constraints, intercept.apply=TRUE)
        constraints = cm.zero.vgam(constraints, x, .zero, M)
    }), list( .exchangeable=exchangeable, .zero=zero ))),
    deviance=Deviance.categorical.data.vgam,
    initialize=eval(substitute(expression({
        eval(process.binomial2.data.vgam)
        predictors.names = c(
                      namesof("mu1", .lmu12, earg= .emu12, short=TRUE),
                      namesof("mu2", .lmu12, earg= .emu12, short=TRUE),
                      namesof("rho", .lrho,  earg= .erho,  short=TRUE))

        if(is.null( .nsimEIM)) {
             save.weight <- control$save.weight <- FALSE
        }
        if(is.null(etastart)) {
            mu1.init= if(is.Numeric(.imu1)) rep(.imu1, len=n) else mu[,3]+mu[,4]
            mu2.init= if(is.Numeric(.imu2)) rep(.imu2, len=n) else mu[,2]+mu[,4]
            rho.init = if(is.Numeric(.init.rho)) rep( .init.rho, len=n) else {
                temp4 = oratio = mu[,1] * mu[,4] / (mu[,2] * mu[,3])
                temp4[oratio <= 0.1] = -0.6
                temp4[oratio >  0.1] = -0.4
                temp4[oratio >  0.5] = -0.2
                temp4[oratio >  0.9] = -0.1
                temp4[oratio >  1.1] =  0.1
                temp4[oratio >  2.0] =  0.3
                temp4[oratio >  6.0] =  0.6
                temp4[oratio > 15.0] =  0.8
                temp4
            }
            etastart = cbind(theta2eta(mu1.init, .lmu12, earg= .emu12),
                             theta2eta(mu2.init, .lmu12, earg= .emu12),
                             theta2eta(rho.init, .lrho, earg= .erho))
        }
    }), list( .lmu12=lmu12, .emu12=emu12, .nsimEIM=nsimEIM,
              .lrho=lrho, .erho=erho, 
              .imu1=imu1, .imu2=imu2, .init.rho=init.rho ))),
    inverse=eval(substitute(function(eta, extra=NULL) {
        pmargin = cbind(eta2theta(eta[,1], .lmu12, earg= .emu12),
                        eta2theta(eta[,2], .lmu12, earg= .emu12))
        rho = eta2theta(eta[,3], .lrho, earg= .erho)
        p11 = pnorm2(eta[,1], eta[,2], rho)
        p01 = pmargin[,2]-p11
        p10 = pmargin[,1]-p11
        p00 = 1-p01-p10-p11
        cbind("00"=p00, "01"=p01, "10"=p10, "11"=p11)
    }, list( .lmu12=lmu12, .emu12=emu12, .lrho=lrho, .erho=erho ))),
    last=eval(substitute(expression({
        misc$link = c(mu1 = .lmu12, mu2 = .lmu12, rho = .lrho)
        misc$earg = list(mu1 = .emu12, mu2 = .emu12, rho = .erho)
        misc$nsimEIM = .nsimEIM
        misc$expected = TRUE
    }), list( .lmu12=lmu12, .emu12=emu12, .lrho=lrho, .erho=erho,
              .nsimEIM=nsimEIM ))),
    vfamily=c("binom2.rho", "binom2"),
    deriv=eval(substitute(expression({
        pmargin = cbind(eta2theta(eta[,1], .lmu12, earg= .emu12),
                        eta2theta(eta[,2], .lmu12, earg= .emu12))
        rhovec = eta2theta(eta[,3], .lrho, earg= .erho)
        p11 = pnorm2(eta[,1], eta[,2], rhovec)
        p01 = pmargin[,2]-p11
        p10 = pmargin[,1]-p11
        p00 = 1-p01-p10-p11

        ABmat = (eta[,1:2] - rhovec*eta[,2:1]) / sqrt(1-rhovec^2)
        PhiA = pnorm(ABmat[,1])
        PhiB = pnorm(ABmat[,2])
        onemPhiA = pnorm(ABmat[,1], lower.tail=FALSE)
        onemPhiB = pnorm(ABmat[,2], lower.tail=FALSE)

        smallno = 1000 * .Machine$double.eps
        p00[p00 < smallno] = smallno
        p01[p01 < smallno] = smallno
        p10[p10 < smallno] = smallno
        p11[p11 < smallno] = smallno

        dprob00 = dnorm2(eta[,1], eta[,2], rhovec)
        dl.dprob1 = PhiB*(y[,4]/p11-y[,2]/p01) + onemPhiB*(y[,3]/p10-y[,1]/p00)
        dl.dprob2 = PhiA*(y[,4]/p11-y[,3]/p10) + onemPhiA*(y[,2]/p01-y[,1]/p00)
        dl.drho = (y[,4]/p11-y[,3]/p10-y[,2]/p01+y[,1]/p00) * dprob00
        dprob1.deta = dtheta.deta(pmargin[,1], .lmu12, earg= .emu12)
        dprob2.deta = dtheta.deta(pmargin[,2], .lmu12, earg= .emu12)
        drho.deta = dtheta.deta(rhovec, .lrho, earg= .erho)
        dthetas.detas = cbind(dprob1.deta, dprob2.deta, drho.deta)

        w * cbind(dl.dprob1, dl.dprob2, dl.drho) * dthetas.detas
    }), list( .lmu12=lmu12, .emu12=emu12, .lrho=lrho, .erho=erho ))),
    weight=eval(substitute(expression({
        if(is.null( .nsimEIM)) {
            d2l.dprob1prob1 = PhiB^2 *(1/p11+1/p01) + onemPhiB^2 *(1/p10+1/p00)
            d2l.dprob2prob2 = PhiA^2 *(1/p11+1/p10) + onemPhiA^2 *(1/p01+1/p00)
            d2l.dprob1prob2 = PhiA * (PhiB/p11 - onemPhiB/p10) +
                              onemPhiA * (onemPhiB/p00 - PhiB/p01)
            d2l.dprob1rho = (PhiB*(1/p11+1/p01) -
                             onemPhiB*(1/p10+1/p00)) * dprob00
            d2l.dprob2rho = (PhiA*(1/p11+1/p10) -
                             onemPhiA*(1/p01+1/p00)) * dprob00
            d2l.drho2 = (1/p11+1/p01+1/p10+1/p00) * dprob00^2
            wz = matrix(0, n, dimm(M))  # 6=dimm(M)
            wz[,iam(1,1,M)] = d2l.dprob1prob1 * dprob1.deta^2
            wz[,iam(2,2,M)] = d2l.dprob2prob2 * dprob2.deta^2
            wz[,iam(1,2,M)] = d2l.dprob1prob2 * dprob1.deta * dprob2.deta
            wz[,iam(1,3,M)] = d2l.dprob1rho * dprob1.deta * drho.deta
            wz[,iam(2,3,M)] = d2l.dprob2rho * dprob2.deta * drho.deta
            wz[,iam(3,3,M)] = d2l.drho2 * drho.deta^2
        } else {
            run.varcov = 0
            ind1 = iam(NA, NA, M=M, both=TRUE, diag=TRUE)
            for(ii in 1:( .nsimEIM )) {
                ysim = rbinom2.rho(n=n, mu1=pmargin[,1], mu2=pmargin[,2],
                                   twoCols=FALSE, rho=rhovec)
                dl.dprob1 = PhiB * (ysim[,4]/p11-ysim[,2]/p01) +
                            onemPhiB * (ysim[,3]/p10-ysim[,1]/p00)
                dl.dprob2 = PhiA * (ysim[,4]/p11-ysim[,3]/p10) +
                            onemPhiA * (ysim[,2]/p01-ysim[,1]/p00)
                dl.drho = (ysim[,4]/p11-ysim[,3]/p10 -
                           ysim[,2]/p01+ysim[,1]/p00) * dprob00

                rm(ysim)
                temp3 = cbind(dl.dprob1, dl.dprob2, dl.drho)
                run.varcov = ((ii-1) * run.varcov +
                           temp3[,ind1$row.index] * temp3[,ind1$col.index]) / ii
            }
            wz = if(intercept.only)
                matrix(apply(run.varcov, 2, mean),
                       n, ncol(run.varcov), byrow=TRUE) else run.varcov
    
            wz = wz * dthetas.detas[,ind1$row] * dthetas.detas[,ind1$col]
        }
        w * wz
    }), list( .nsimEIM=nsimEIM, .lmu12=lmu12, .emu12=emu12, .lrho=lrho,
              .erho=erho ))))
}



dnorm2 <- function(x, y, r) 
    exp(-0.5*(x^2+y^2-2*x*y*r)/(1-r^2)) / (2*pi*sqrt(1-r^2))


pnorm2 <- function(ah, ak, r) { 

    ans <- ah 
    size <- length(ah)
    singler <- ifelse(length(r)==1,1,0)
    dotC(name="pnorm2", ah=as.double(-ah), ak=as.double(-ak), r=as.double(r),
       size=as.integer(size), singler=as.integer(singler), 
       ans=as.double(ans))$ans
}





my.dbinom <- function(x,
                      size = stop("no size arg"),
                      prob = stop("no prob arg"))
{

    exp( lgamma(size+1) - lgamma(size-x+1) - lgamma(x+1) +
              x * log(prob/(1-prob)) + size * log1p(-prob) )
}



size.binomial <- function(prob=0.5, link="loge", earg=list())
{
    if(any(prob<=0 || prob>=1))
        stop("some values of prob out of range")
    if(!missing(link)) link <- as.character(substitute(link))
    if(!is.list(earg)) earg = list()

    new("vglmff",
    blurb=c("Binomial with n unknown, prob known (prob=",prob,")\n",
           "Links:    ",
           namesof("size", link, tag=TRUE),
           " (treated as real-valued)\n",
           "Variance:  Var(Y) = size * prob * (1-prob);",
           " Var(size) is intractable"),
    initialize=eval(substitute(expression({
        predictors.names <- "size"
        extra$temp2 <- rep( .prob , length=n)
        if(is.null(etastart)) {
            nvec <- (y+0.1)/extra$temp2
            etastart <- theta2eta(nvec, .link)
        }
    }), list( .prob =prob, .link=link ))),
    inverse=eval(substitute(function(eta, extra=NULL) {
        nvec <- eta2theta(eta, .link)
        nvec*extra$temp2
    }, list( .link=link ))),
    last=eval(substitute(expression({
        misc$link <- c(size = .link)
        misc$prob <- extra$temp2
    }), list( .link=link ))),
    link=eval(substitute(function(mu, extra=NULL) {
        nvec <- mu/extra$temp2
        theta2eta(nvec, .link)
    }, list( .link=link ))),
    loglikelihood=eval(substitute(
        function(mu, y, w, res=FALSE,eta, extra=NULL) {
        nvec <- mu/extra$temp2
        sum(w * (lgamma(nvec+1) - lgamma(y+1) - lgamma(nvec-y+1) +
            y * log(.prob / (1- .prob)) + nvec * log1p(- .prob)))
    }, list( .prob=prob ))),
    vfamily=c("size.binomial"),
    deriv=eval(substitute(expression({
        nvec <- mu/extra$temp2
        dldnvec = digamma(nvec+1) - digamma(nvec-y+1) + log1p(-extra$temp2)
        dnvecdeta <- dtheta.deta(nvec, .link)
        w * cbind(dldnvec * dnvecdeta)
    }), list( .link=link ))),
    weight=eval(substitute(expression({
        d2ldnvec2 <- trigamma(nvec+1) - trigamma(nvec-y+1)
        # Note: if y==0 then d2ldnvec2 is 0. Below is a quick fix.
        d2ldnvec2[y==0] = -sqrt(.Machine$double.eps)
        wz = -w * dnvecdeta^2 * d2ldnvec2
        wz
    }), list( .link=link ))))
}




dbetabin.ab = function(x, size, shape1, shape2, log = FALSE) {
    log.arg = log
    rm(log)
    if(!is.Numeric(x)) stop("bad input for argument \"x\"")
    if(!is.Numeric(size, posit=TRUE, integer=TRUE))
        stop("bad input for argument \"size\"")
    if(!is.Numeric(shape1, pos=TRUE)) stop("bad input for argument \"shape1\"")
    if(!is.Numeric(shape2, pos=TRUE)) stop("bad input for argument \"shape2\"")
    L = max(length(x), length(size), length(shape1), length(shape2))
    x = rep(x, len=L); size = rep(size, len=L);
    shape1 = rep(shape1, len=L); shape2 = rep(shape2, len=L);
    answer = 0 * x
    ok <- round(x) == x & x >= 0 & x <= size
    answer[ok] = lchoose(size[ok], x[ok]) +
                 lbeta(shape1[ok]+x[ok], shape2[ok]+size[ok]-x[ok]) -
                 lbeta(shape1[ok], shape2[ok])
    if(log.arg) {
        answer[!ok] = log(0.0)
    } else {
        answer[ok] = exp(answer[ok])
    }
    answer
}


pbetabin.ab = function(q, size, shape1, shape2, log.p=FALSE) {
    if(!is.Numeric(q)) stop("bad input for argument \"q\"")
    if(!is.Numeric(size, posit=TRUE, integer=TRUE))
        stop("bad input for argument \"size\"")
    if(!is.Numeric(shape1, pos=TRUE)) stop("bad input for argument \"shape1\"")
    if(!is.Numeric(shape2, pos=TRUE)) stop("bad input for argument \"shape2\"")
    N = max(length(q), length(size), length(shape1), length(shape2))
    q = rep(q, len=N); shape1 = rep(shape1, len=N); shape2 = rep(shape2, len=N)
    size = rep(size, len=N);
    ans = q * 0  # Retains names(q)
    if(max(abs(size-size[1])) < 1.0e-08 &&
       max(abs(shape1-shape1[1])) < 1.0e-08 &&
       max(abs(shape2-shape2[1])) < 1.0e-08) {
        qstar = floor(q)
        temp = if(max(qstar) >= 0) dbetabin.ab(0:max(qstar), 
               size=size[1], shape1=shape1[1], shape2=shape2[1]) else 0*qstar
        unq = unique(qstar)
        for(i in unq) {
            index = qstar == i
            ans[index] = if(i >= 0) sum(temp[1:(1+i)]) else 0
        }
    } else
    for(i in 1:N) {
        qstar = floor(q[i])
        ans[i] = if(qstar >= 0) sum(dbetabin.ab(x=0:qstar, size=size[i],
                 shape1=shape1[i], shape2=shape2[i])) else 0
    }
    if(log.p) log(ans) else ans
}

rbetabin.ab = function(n, size, shape1, shape2) {
    if(!is.Numeric(n,integ=TRUE, allow=1)) stop("bad input for argument \"n\"")
    if(!is.Numeric(size, posit=TRUE, integer=TRUE))
        stop("bad input for argument \"size\"")
    if(!is.Numeric(shape1, pos=TRUE)) stop("bad input for argument \"shape1\"")
    if(!is.Numeric(shape2, pos=TRUE)) stop("bad input for argument \"shape2\"")
    size = rep(size, len=n);
    shape1 = rep(shape1, len=n); shape2 = rep(shape2, len=n);
    rbinom(n=n, size=size, prob=rbeta(n=n, shape1=shape1, shape2=shape2))
}


dbetabin = function(x, size, prob, rho, log = FALSE) {
    dbetabin.ab(x=x, size=size, shape1=prob*(1-rho)/rho,
                shape2=(1-prob)*(1-rho)/rho, log=log)
}

pbetabin = function(q, size, prob, rho, log.p=FALSE) {
    pbetabin.ab(q=q, size=size, shape1=prob*(1-rho)/rho,
                shape2=(1-prob)*(1-rho)/rho, log.p=log.p)
}

rbetabin = function(n, size, prob, rho) {
    rbetabin.ab(n=n, size=size, shape1=prob*(1-rho)/rho,
                shape2=(1-prob)*(1-rho)/rho)
}


expected.betabin.ab = function(nvec, shape1, shape2, first) {

    n = length(nvec)
    ans = rep(0.0, len=n)
    if(!is.R()) {
        lchoose = function(a,b) log(choose(a,b))
        lbeta = function(a,b) lgamma(a) + lgamma(b) - lgamma(a+b)
    }
    if(first) {
        for(i in 1:n) {
            temp639 = lbeta(shape1[i], shape2[i])
            for(y in 0:nvec[i])
                ans[i] = ans[i] + trigamma(shape1[i]+y) *
                         exp(lchoose(nvec[i], y) +
                         lbeta(shape1[i]+y, shape2[i]+nvec[i]-y) - temp639)
        }
    } else {
        for(i in 1:n) {
            temp639 = lbeta(shape1[i], shape2[i])
            for(y in 0:nvec[i])
                ans[i] = ans[i] + trigamma(nvec[i]+shape2[i]-y) *
                         exp(lchoose(nvec[i], y) +
                             lbeta(shape1[i]+y, shape2[i]+nvec[i]-y) - temp639)
        }
    }
    ans
}




betabin.ab = function(link.shape12="loge", earg = list(),
                      i1=1, i2=NULL, zero=NULL)
{
    if(mode(link.shape12) != "character" && mode(link.shape12) != "name")
        link.shape12 = as.character(substitute(link.shape12))
    if(!is.Numeric(i1, positive=TRUE)) stop("bad input for argument \"i1\"")
    if(length(i2) && !is.Numeric(i2, pos=TRUE))
        stop("bad input for argument \"i2\"")
    if(!is.list(earg)) earg = list()

    new("vglmff",
    blurb=c("Beta-binomial model\n",
           "Links:    ",
           namesof("shape1", link.shape12, earg= earg), ", ",
           namesof("shape2", link.shape12, earg= earg), "\n",
           "Variance: mu*(1-mu)[1+(w-1)*rho]/w where mu=alpha/(alpha+beta)"),
    constraints=eval(substitute(expression({
        constraints <- cm.zero.vgam(constraints, x, .zero, M)
    }), list( .zero=zero ))),
    initialize=eval(substitute(expression({
        # Compute initial values for mustart -------
        eval(binomialff()@initialize)   # Note: n,w,y,mustart is changed 
        predictors.names = c(namesof("shape1", .link.shape12, earg= .earg, tag=FALSE),
                             namesof("shape2", .link.shape12, earg= .earg, short=FALSE))

        if(!length(etastart)) {
            shape1 = rep( .i1, len=n)
            shape2 = if(length( .i2)) rep( .i2,len=n) else shape1*(1/mustart-1)
            ycounts = y * w   # Convert proportions to counts
            if(max(abs(ycounts-round(ycounts))) > 1.0e-6)
                stop("ycounts not integer")
            ycounts = round(ycounts) # Make sure it is an integer
            etastart = cbind(theta2eta(shape1, .link.shape12, earg= .earg),
                             theta2eta(shape2, .link.shape12, earg= .earg))
        }
    }), list( .link.shape12=link.shape12, .earg=earg, .i1=i1 , .i2=i2 ))),
    inverse=eval(substitute(function(eta, extra=NULL) {
        shape1 = eta2theta(eta[,1], .link.shape12, earg= .earg)
        shape2 = eta2theta(eta[,2], .link.shape12, earg= .earg)
        shape1 / (shape1 + shape2)
    }, list( .link.shape12=link.shape12, .earg=earg ))),
    last=eval(substitute(expression({
        misc$link = c("shape1" = .link.shape12, "shape2" = .link.shape12)
        misc$earg <- list(shape1 = .earg, shape2 = .earg)
        shape1 = eta2theta(eta[,1], .link.shape12, earg= .earg)
        shape2 = eta2theta(eta[,2], .link.shape12, earg= .earg)
        misc$rho = 1 / (shape1 + shape2 + 1)
        misc$expected = TRUE
    }), list( .link.shape12=link.shape12, .earg=earg ))),
    loglikelihood=eval(substitute(
        function(mu,y,w,residuals=FALSE,eta, extra=NULL) {
        ycounts = y * w   # Convert proportions to counts
        shape1 = eta2theta(eta[,1], .link.shape12, earg= .earg)
        shape2 = eta2theta(eta[,2], .link.shape12, earg= .earg)
        nvec = w
        if(residuals) stop("loglikelihood residuals not implemented yet") else {
            if(is.R()) sum(lbeta(shape1+ycounts, shape2+nvec-ycounts) -
                           lbeta(shape1, shape2)) else
            sum(lgamma(shape1+ycounts) + lgamma(shape2+nvec-ycounts) -
                lgamma(shape1+shape2+nvec) -
                (lgamma(shape1) + lgamma(shape2) - lgamma(shape1+shape2)))
        }
    }, list( .link.shape12=link.shape12, .earg=earg ))),
    vfamily=c("betabin.ab"),
    deriv=eval(substitute(expression({
        nvec = w  # extra$nvec # for summary()
        ycounts = y * w   # Convert proportions to counts
        shape1 = eta2theta(eta[,1], .link.shape12, earg= .earg)
        shape2 = eta2theta(eta[,2], .link.shape12, earg= .earg)
        dshape1.deta = dtheta.deta(shape1, .link.shape12, earg= .earg)
        dshape2.deta = dtheta.deta(shape2, .link.shape12, earg= .earg)
        dl.dshape1 = digamma(shape1+ycounts) - digamma(shape1+shape2+nvec) -
                     digamma(shape1) + digamma(shape1+shape2)
        dl.dshape2 = digamma(nvec+shape2-ycounts) -
                     digamma(shape1+shape2+nvec) -
                     digamma(shape2) + digamma(shape1+shape2)
        cbind(dl.dshape1 * dshape1.deta, dl.dshape2 * dshape2.deta)
    }), list( .link.shape12=link.shape12, .earg=earg ))),
    weight=eval(substitute(expression({
        wz = matrix(as.numeric(NA), n, dimm(M))  #3=dimm(2)
        wz[,iam(1,1,M)] = -(expected.betabin.ab(nvec, shape1, shape2, TRUE) -
                            trigamma(shape1+shape2+nvec) -
                            trigamma(shape1) + trigamma(shape1+shape2)) *
                          dshape1.deta^2
        wz[,iam(2,2,M)] = -(expected.betabin.ab(nvec, shape1, shape2, FALSE) -
                            trigamma(shape1+shape2+nvec) -
                            trigamma(shape2) + trigamma(shape1+shape2)) *
                          dshape2.deta^2
        wz[,iam(2,1,M)] = -(trigamma(shape1+shape2) -
                           trigamma(shape1+shape2+nvec)) *
                          dshape1.deta * dshape2.deta
        wz
    }), list( .link.shape12=link.shape12, .earg=earg ))))
}



betageometric = function(lprob="logit", lshape="loge",
                         eprob=list(), eshape=list(),
                         iprob = NULL, ishape = 0.1,
                         moreSummation=c(2,100), tolerance=1.0e-10, zero=NULL)
{
    if(mode(lprob) != "character" && mode(lprob) != "name")
        lprob = as.character(substitute(lprob))
    if(mode(lshape) != "character" && mode(lshape) != "name")
        lshape = as.character(substitute(lshape))
    if(!is.Numeric(ishape, positive=TRUE))
        stop("bad input for argument \"ishape\"")
    if(!is.Numeric(moreSummation, positive=TRUE, allow=2, integ=TRUE))
        stop("bad input for argument \"moreSummation\"")
    if(!is.Numeric(tolerance, positive=TRUE, allow=1) || 1.0-tolerance >= 1.0)
        stop("bad input for argument \"tolerance\"")
    if(!is.list(eprob)) eprob = list()
    if(!is.list(eshape)) eshape = list()

    new("vglmff",
    blurb=c("Beta-geometric distribution\n",
           "Links:    ", namesof("prob", lprob, earg= eprob), ", ",
                         namesof("shape", lshape, earg= eshape)),
    constraints=eval(substitute(expression({
        constraints <- cm.zero.vgam(constraints, x, .zero, M)
    }), list( .zero=zero ))),
    initialize=eval(substitute(expression({
        eval(geometric()@initialize)
        predictors.names = c(namesof("prob",  .lprob, earg= .eprob,  tag=FALSE),
                             namesof("shape", .lshape, earg= .eshape, short=FALSE))
        if(length( .iprob))
            prob.init = rep( .iprob, len=n)
        if(!length(etastart) || ncol(cbind(etastart)) != 2) {
            shape.init = rep( .ishape, len=n)
            etastart = cbind(theta2eta(prob.init,  .lprob, earg= .eprob),
                             theta2eta(shape.init, .lshape, earg= .eshape))
        }
    }), list( .iprob=iprob, .ishape=ishape, .lprob=lprob,
              .eprob=eprob, .eshape=eshape,
              .lshape=lshape ))),
    inverse=eval(substitute(function(eta, extra=NULL) {
        prob  = eta2theta(eta[,1], .lprob, earg= .eprob)
        shape = eta2theta(eta[,2], .lshape, earg= .eshape)
        mymu = (1-prob) / (prob - shape)
        ifelse(mymu >= 0, mymu, NA)
    }, list( .lprob=lprob, .lshape=lshape,
             .eprob=eprob, .eshape=eshape ))),
    last=eval(substitute(expression({
        misc$link = c("prob" = .lprob, "shape" = .lshape)
        misc$earg <- list(prob = .eprob, shape = .eshape)
        if(intercept.only) {
            misc$shape1 = shape1[1]  # These quantities computed in @deriv
            misc$shape2 = shape2[1]
        }
        misc$expected = TRUE
        misc$tolerance = .tolerance
        misc$zero = .zero
        misc$moreSummation = .moreSummation
    }), list( .lprob=lprob, .lshape=lshape, .tolerance=tolerance,
              .eprob=eprob, .eshape=eshape,
              .moreSummation=moreSummation, .zero=zero ))),
    loglikelihood=eval(substitute(
        function(mu,y,w,residuals=FALSE,eta, extra=NULL) {
        prob  = eta2theta(eta[,1], .lprob, earg= .eprob)
        shape = eta2theta(eta[,2], .lshape, earg= .eshape)
        ans = log(prob)
        maxy = max(y)
        if(residuals) stop("loglikelihood residuals not implemented yet") else {
            for(ii in 1:maxy) {
                index = ii <= y
                ans[index]=ans[index] + log1p(-prob[index]+(ii-1)*shape[index])-
                           log1p((ii-1)*shape[index])
            }
            ans = ans - log1p((y+1-1)*shape)
            sum(w * ans)
        }
    }, list( .lprob=lprob, .lshape=lshape,
             .eprob=eprob, .eshape=eshape ))),
    vfamily=c("betageometric"),
    deriv=eval(substitute(expression({
        prob  = eta2theta(eta[,1], .lprob, earg= .eprob)
        shape = eta2theta(eta[,2], .lshape, earg= .eshape)
        shape1 = prob / shape; shape2 = (1-prob) / shape;
        dprob.deta  = dtheta.deta(prob,  .lprob, earg= .eprob)
        dshape.deta = dtheta.deta(shape, .lshape, earg= .eshape)
        dl.dprob = 1 / prob
        dl.dshape = 0 * y
        maxy = max(y)
        for(ii in 1:maxy) {
            index = ii <= y
            dl.dprob[index] = dl.dprob[index] -
                              1/(1-prob[index]+(ii-1)*shape[index])
            dl.dshape[index] = dl.dshape[index] +
                              (ii-1)/(1-prob[index]+(ii-1)*shape[index]) -
                              (ii-1)/(1+(ii-1)*shape[index])
        }
        dl.dshape = dl.dshape - (y+1 -1)/(1+(y+1 -1)*shape)
        w * cbind(dl.dprob * dprob.deta, dl.dshape * dshape.deta)
    }), list( .lprob=lprob, .lshape=lshape,
              .eprob=eprob, .eshape=eshape ))),
    weight=eval(substitute(expression({
        wz = matrix(0, n, dimm(M))  #3=dimm(2)
        wz[,iam(1,1,M)] = 1 / prob^2
        moresum = .moreSummation
        maxsummation = round(maxy * moresum[1] + moresum[2])
        for(ii in 3:maxsummation) {
            temp7 = 1 - pbetageom(q=ii-1-1, shape1=shape1, shape2=shape2)
            denom1 = (1-prob+(ii-2)*shape)^2
            denom2 = (1+(ii-2)*shape)^2
            wz[,iam(1,1,M)] = wz[,iam(1,1,M)] + temp7 / denom1
            wz[,iam(1,2,M)] = wz[,iam(1,2,M)] - (ii-2) * temp7 / denom1
            wz[,iam(2,2,M)] = wz[,iam(2,2,M)] + (ii-2)^2 * temp7 / denom1 -
                              (ii-1)^2 * temp7 / denom2
            if(max(temp7) < .tolerance ) break;
        }
        ii = 2
        temp7 = 1 - pbetageom(q=ii-1-1, shape1=shape1, shape2=shape2)
        denom1 = (1-prob+(ii-2)*shape)^2
        denom2 = (1+(ii-2)*shape)^2
        wz[,iam(1,1,M)] = wz[,iam(1,1,M)] + temp7 / denom1
        wz[,iam(2,2,M)] = wz[,iam(2,2,M)] - (ii-1)^2 * temp7 / denom2
        wz[,iam(1,1,M)] = wz[,iam(1,1,M)] * dprob.deta^2
        wz[,iam(2,2,M)] = wz[,iam(2,2,M)] * dshape.deta^2
        wz[,iam(2,1,M)] = wz[,iam(2,1,M)] * dprob.deta * dshape.deta
        w * wz
    }), list( .lprob=lprob, .lshape=lshape, .moreSummation=moreSummation,
              .eprob=eprob, .eshape=eshape,
              .tolerance=tolerance ))))
}




seq2binomial = function(lprob1="logit", lprob2="logit",
                        eprob1=list(), eprob2=list(),
                        iprob1 = NULL, iprob2 = NULL,
                        zero=NULL)
{
    if(mode(lprob1) != "character" && mode(lprob1) != "name")
        lprob1 = as.character(substitute(lprob1))
    if(mode(lprob2) != "character" && mode(lprob2) != "name")
        lprob2 = as.character(substitute(lprob2))
    if(length(iprob1) &&
       (!is.Numeric(iprob1, positive=TRUE) || max(iprob1) >= 1))
        stop("bad input for argument \"iprob1\"")
    if(length(iprob2) &&
       (!is.Numeric(iprob2, positive=TRUE) || max(iprob2) >= 1))
        stop("bad input for argument \"iprob2\"")
    if(!is.list(eprob1)) eprob1 = list()
    if(!is.list(eprob2)) eprob2 = list()

    new("vglmff",
    blurb=c("Sequential binomial distribution (Crowder and Sweeting, 1989)\n",
           "Links:    ", namesof("prob1", lprob1, earg= eprob1), ", ",
                         namesof("prob2", lprob2, earg= eprob2)),
    constraints=eval(substitute(expression({
        constraints <- cm.zero.vgam(constraints, x, .zero, M)
    }), list( .zero=zero ))),
    initialize=eval(substitute(expression({
        if(!is.vector(w))
            stop("the 'weights' argument must be a vector")
        if(any(w != round(w)))
            warning("the 'weights' argument should be integer-valued")
        if(ncol(y <- cbind(y)) != 2)
            stop("the response must be a 2-column matrix")
        if(any(y < 0 | y > 1))
            stop("the response must have values between 0 and 1")
        rvector = w * y[,1]
        if(any(abs(rvector - round(rvector)) > 1.0e-8))
        warning("number of successes in column one should be integer-valued")
        svector = rvector * y[,2]
        if(any(abs(svector - round(svector)) > 1.0e-8))
        warning("number of successes in column two should be integer-valued")
        predictors.names = c(namesof("prob1", .lprob1,earg= .eprob1, tag=FALSE),
                             namesof("prob2", .lprob2,earg= .eprob2, tag=FALSE))
        prob1.init = if(is.Numeric( .iprob1)) rep( .iprob1, len=n) else
                     rep(weighted.mean(y[,1], w=w), len=n)
        prob2.init = if(is.Numeric( .iprob2)) rep( .iprob2, len=n) else
                     rep(weighted.mean(y[,2], w=w*y[,1]), len=n)
        if(!length(etastart)) {
            etastart = cbind(theta2eta(prob1.init, .lprob1, earg= .eprob1),
                             theta2eta(prob2.init, .lprob2, earg= .eprob2))
        }
    }), list( .iprob1=iprob1, .iprob2=iprob2, .lprob1=lprob1,
              .eprob1=eprob1, .eprob2=eprob2,
              .lprob2=lprob2 ))),
    inverse=eval(substitute(function(eta, extra=NULL) {
        prob1 = eta2theta(eta[,1], .lprob1, earg= .eprob1)
        prob2 = eta2theta(eta[,2], .lprob2, earg= .eprob2)
        cbind(prob1, prob2)
    }, list( .lprob1=lprob1, .lprob2=lprob2,
             .eprob1=eprob1, .eprob2=eprob2 ))),
    last=eval(substitute(expression({
        misc$link = c("prob1" = .lprob1, "prob2" = .lprob2)
        misc$earg <- list(prob1 = .eprob1, prob2 = .eprob2)
        misc$expected = TRUE
        misc$zero = .zero
    }), list( .lprob1=lprob1, .lprob2=lprob2,
              .eprob1=eprob1, .eprob2=eprob2,
              .zero=zero ))),
    loglikelihood=eval(substitute(
        function(mu,y,w,residuals=FALSE,eta, extra=NULL) {
        prob1 = eta2theta(eta[,1], .lprob1, earg= .eprob1)
        prob2 = eta2theta(eta[,2], .lprob2, earg= .eprob2)
        smallno = 100 * .Machine$double.eps
        prob1 = pmax(prob1, smallno)
        prob1 = pmin(prob1, 1-smallno)
        prob2 = pmax(prob2, smallno)
        prob2 = pmin(prob2, 1-smallno)
        rvector = w * y[,1]
        svector = rvector * y[,2]
        if(residuals) stop("loglikelihood residuals not implemented yet") else {
            sum(rvector * log(prob1) + (mvector-rvector)*log1p(-prob1) +
                svector * log(prob2) + (rvector-svector)*log1p(-prob2))
        }
    }, list( .lprob1=lprob1, .lprob2=lprob2,
             .eprob1=eprob1, .eprob2=eprob2 ))),
    vfamily=c("seq2binomial"),
    deriv=eval(substitute(expression({
        prob1 = eta2theta(eta[,1], .lprob1, earg= .eprob1)
        prob2 = eta2theta(eta[,2], .lprob2, earg= .eprob2)
        smallno = 100 * .Machine$double.eps
        prob1 = pmax(prob1, smallno)
        prob1 = pmin(prob1, 1-smallno)
        prob2 = pmax(prob2, smallno)
        prob2 = pmin(prob2, 1-smallno)
        dprob1.deta = dtheta.deta(prob1, .lprob1, earg= .eprob1)
        dprob2.deta = dtheta.deta(prob2, .lprob2, earg= .eprob2)
        rvector = w * y[,1]
        svector = rvector * y[,2]
        dl.dprob1 = rvector / prob1 - (mvector-rvector) / (1-prob1)
        dl.dprob2 = svector / prob2 - (rvector-svector) / (1-prob2)
        cbind(dl.dprob1 * dprob1.deta, dl.dprob2 * dprob2.deta)
    }), list( .lprob1=lprob1, .lprob2=lprob2,
              .eprob1=eprob1, .eprob2=eprob2 ))),
    weight=eval(substitute(expression({
        wz = matrix(0, n, M)
        wz[,iam(1,1,M)] = (dprob1.deta^2) / (prob1 * (1-prob1))
        wz[,iam(2,2,M)] = (dprob2.deta^2) * prob1 / (prob2 * (1-prob2))
        w * wz
    }), list( .lprob1=lprob1, .lprob2=lprob2,
              .eprob1=eprob1, .eprob2=eprob2 ))))
}



zipebcom   = function(lmu12="cloglog", lphi12="logit", loratio="loge",
                      emu12=list(), ephi12=list(), eoratio=list(), 
                      imu12=NULL, iphi12=NULL, ioratio = NULL, 
                      zero=2:3, tol=0.001,
                      addRidge=0.001)
{

    if(mode(lphi12) != "character" && mode(lphi12) != "name")
        lphi12 = as.character(substitute(lphi12))
    if(mode(loratio) != "character" && mode(loratio) != "name")
        loratio = as.character(substitute(loratio))
    if(!is.Numeric(tol, positive=TRUE, allow=1) || tol > 0.1)
        stop("bad input for argument \"tol\"") 
    if(!is.Numeric(addRidge, allow=1, posit=TRUE) || addRidge > 0.5)
        stop("bad input for argument \"addRidge\"") 
    if(!is.list(emu12)) emu12  = list()
    if(!is.list(ephi12)) ephi12  = list()
    if(!is.list(eoratio)) eoratio = list()
    if(lmu12 != "cloglog")
        warning("argument 'lmu12' should be 'cloglog'")

    new("vglmff",
    blurb=c("Exchangeable bivariate ", lmu12, " odds-ratio model based on\n",
            "a zero-inflated Poisson distribution\n\n",
            "Links:    ",
            namesof("mu12", lmu12, earg=emu12), ", ",
            namesof("phi12", lphi12, earg=ephi12), ", ",
            namesof("oratio", loratio, earg=eoratio)),
    constraints=eval(substitute(expression({
        constraints = cm.zero.vgam(constraints, x, .zero, M)
    }), list( .zero=zero ))),
    initialize=eval(substitute(expression({
        eval(process.binomial2.data.vgam)
        predictors.names = c(
                 namesof("mu12",   .lmu12, earg= .emu12, short=TRUE), 
                 namesof("phi12",  .lphi12, earg= .ephi12, short=TRUE),
                 namesof("oratio", .loratio, earg= .eoratio, short=TRUE))

        propY1.eq.0 = weighted.mean(y[,'00'], w) + weighted.mean(y[,'01'], w)
        propY2.eq.0 = weighted.mean(y[,'00'], w) + weighted.mean(y[,'10'], w)
        if(length( .iphi12) && any(.iphi12 > propY1.eq.0))
            warning("iphi12 must be less than the sample proportion of Y1==0")
        if(length( .iphi12) && any(.iphi12 > propY2.eq.0))
            warning("iphi12 must be less than the sample proportion of Y2==0")

        if(!length(etastart)) {
            pstar.init = ((mu[,3]+mu[,4]) + (mu[,2]+mu[,4])) / 2
            phi.init = if(length(.iphi12)) rep(.iphi12, len=n) else
                min(propY1.eq.0 * 0.95, propY2.eq.0 * 0.95, pstar.init/1.5)
            oratio.init = if(length( .ioratio)) rep( .ioratio, len=n) else
                      mu[,4]*mu[,1]/(mu[,2]*mu[,3])
            mu12.init = if(length(.imu12)) rep(.imu12, len=n) else
                pstar.init / (1-phi.init)
            etastart = cbind(
                theta2eta(mu12.init, .lmu12, earg= .emu12),
                theta2eta(phi.init, .lphi12, earg= .ephi12),
                theta2eta(oratio.init, .loratio, earg= .eoratio))
        }
    }), list( .lmu12=lmu12, .lphi12=lphi12, .loratio=loratio,
              .emu12=emu12, .ephi12=ephi12, .eoratio=eoratio,
              .imu12=imu12, .iphi12=iphi12, .ioratio=ioratio ))),
    inverse=eval(substitute(function(eta, extra=NULL) {
        A1vec  = eta2theta(eta[,1], .lmu12,  earg= .emu12)
        phivec = eta2theta(eta[,2], .lphi12, earg= .ephi12)
        pmargin = matrix((1 - phivec) * A1vec, nrow(eta), 2)
        oratio = eta2theta(eta[,3], .loratio, earg= .eoratio)
        a.temp = 1 + (pmargin[,1]+pmargin[,2])*(oratio-1)
        b.temp = -4 * oratio * (oratio-1) * pmargin[,1] * pmargin[,2]
        temp = sqrt(a.temp^2 + b.temp)
        pj4 = ifelse(abs(oratio-1) < .tol, pmargin[,1]*pmargin[,2],
                     (a.temp-temp)/(2*(oratio-1)))
        pj2 = pmargin[,2] - pj4
        pj3 = pmargin[,1] - pj4
        cbind("00" = 1-pj4-pj2-pj3, "01" = pj2, "10" = pj3, "11" = pj4)
    }, list( .tol=tol,
             .lmu12=lmu12, .lphi12=lphi12, .loratio=loratio,
             .emu12=emu12, .ephi12=ephi12, .eoratio=eoratio ))),
    last=eval(substitute(expression({
        misc$link = c("mu12"= .lmu12, "phi12" = .lphi12, "oratio"= .loratio)
        misc$earg = list("mu12"= .emu12, "phi12"= .ephi12, "oratio"= .eoratio)
        misc$tol = .tol
        misc$expected = TRUE
        misc$addRidge = .addRidge
    }), list( .tol=tol, .addRidge = addRidge,
              .lmu12=lmu12, .lphi12=lphi12, .loratio=loratio,
              .emu12=emu12, .ephi12=ephi12, .eoratio=eoratio ))),
    loglikelihood= function(mu, y, w, residuals = FALSE, eta, extra=NULL)
        if(residuals) stop("loglikelihood residuals not implemented yet") else
        sum(w * y * log(mu)),
    vfamily=c("zipebcom"),
    deriv=eval(substitute(expression({
        A1vec  = eta2theta(eta[,1], .lmu12,  earg= .emu12)
        smallno = .Machine$double.eps^(2/4)
        A1vec[A1vec > 1.0 -smallno] = 1.0 - smallno

        phivec = eta2theta(eta[,2], .lphi12, earg= .ephi12)
        pmargin = matrix((1 - phivec) * A1vec, n, 2)
        oratio = eta2theta(eta[,3], .loratio, earg= .eoratio)

        Vab = 1 / (1/mu[,1] + 1/mu[,2] + 1/mu[,3] + 1/mu[,4])
        Vabc = 1/mu[,1] + 1/mu[,2]
        denom3 = 2 * oratio * mu[,2] + mu[,1] + mu[,4]
        temp1 = oratio * mu[,2] + mu[,4]
        dp11star.dp1unstar = 2*(1-phivec)*Vab * Vabc
        dp11star.dphi1 = -2 * A1vec * Vab * Vabc
        dp11star.doratio = Vab / oratio 
        yandmu = (y[,1]/mu[,1] - y[,2]/mu[,2] - y[,3]/mu[,3] + y[,4]/mu[,4])
        dp11.doratio = Vab / oratio
        check.dl.doratio = yandmu * dp11.doratio

        cyandmu = (y[,2]+y[,3])/mu[,2] - 2 * y[,1]/mu[,1]
        dl.dmu1 = dp11star.dp1unstar * yandmu + (1-phivec) * cyandmu
        dl.dphi1 = dp11star.dphi1 * yandmu - A1vec * cyandmu
        dl.doratio = check.dl.doratio
        dthetas.detas = cbind(dtheta.deta(A1vec,  .lmu12,   earg= .emu12),
                              dtheta.deta(phivec, .lphi12,  earg= .ephi12),
                              dtheta.deta(oratio, .loratio, earg= .eoratio))
        w * cbind(dl.dmu1, dl.dphi1, dl.doratio) * dthetas.detas
    }), list( .lmu12=lmu12, .lphi12=lphi12, .loratio=loratio,
              .emu12=emu12, .ephi12=ephi12, .eoratio=eoratio ))),
    weight=eval(substitute(expression({
        wz = matrix(0, n, 4)
        alternwz11 = 2*(1-phivec)^2 *(2/mu[,1] + 1/mu[,2] - 2*Vab*Vabc^2) *
                     (dthetas.detas[,1])^2
        wz[,iam(1,1,M)] = alternwz11

        alternwz22 = 2* A1vec^2 *(2/mu[,1] + 1/mu[,2] - 2*Vab*Vabc^2) *
                     (dthetas.detas[,2])^2
        wz[,iam(2,2,M)] = alternwz22

        alternwz12 = -2*A1vec*(1-phivec)*(2/mu[,1] + 1/mu[,2] - 2*Vab*Vabc^2) *
                     dthetas.detas[,1] * dthetas.detas[,2]
        wz[,iam(1,2,M)] = alternwz12

        alternwz33 = (Vab / oratio^2) * dthetas.detas[,3]^2
        wz[,iam(3,3,M)] = alternwz33

        wz[,1:2] = wz[,1:2] * (1 + .addRidge)
        w * wz
    }), list( .addRidge = addRidge ))))
}





# These functions are
# Copyright (C) 1998-2006 T.W. Yee, University of Auckland. All rights reserved.






mix2normal1.control <- function(save.weight=TRUE, ...)
{
    list(save.weight=save.weight)
}

mix2normal1 = function(lphi="logit",
                       lmu="identity",
                       lsd="loge",
                       ephi=list(), emu1=list(), emu2=list(), esd1=list(), esd2=list(),
                       iphi=0.5, imu1=NULL, imu2=NULL, isd1=NULL, isd2=NULL,
                       qmu=c(0.2, 0.8),
                       esd=FALSE,
                       zero=1)
{
    if(mode(lphi) != "character" && mode(lphi) != "name")
        lphi = as.character(substitute(lphi))
    if(mode(lmu) != "character" && mode(lmu) != "name")
        lmu = as.character(substitute(lmu))
    if(mode(lsd) != "character" && mode(lsd) != "name")
        lsd = as.character(substitute(lsd))
    if(!is.Numeric(qmu, allow=2, positive=TRUE) || any(qmu >= 1))
        stop("bad input for argument \"qmu\"")
    if(length(iphi) && (!is.Numeric(iphi, allow=1, positive=TRUE) || iphi>= 1))
        stop("bad input for argument \"iphi\"")
    if(length(imu1) && !is.Numeric(imu1))
        stop("bad input for argument \"imu1\"")
    if(length(imu2) && !is.Numeric(imu2))
        stop("bad input for argument \"imu2\"")
    if(length(isd1) && !is.Numeric(isd1, positive=TRUE))
        stop("bad input for argument \"isd1\"")
    if(length(isd2) && !is.Numeric(isd2, positive=TRUE))
        stop("bad input for argument \"isd2\"")
    if(!is.list(ephi)) ephi = list()
    if(!is.list(emu1)) emu1 = list()
    if(!is.list(emu2)) emu2 = list()
    if(!is.list(esd1)) esd1 = list()
    if(!is.list(esd2)) esd2 = list()

    new("vglmff",
    blurb=c("Mixture of two univariate normals\n\n",
           "Links:    ",
           namesof("phi",lphi, earg= ephi), ", ", 
           namesof("mu1", lmu, earg= emu1, tag=FALSE), ", ",
           namesof("sd1", lsd, earg= esd1, tag=FALSE), ", ",
           namesof("mu2", lmu, earg= emu2, tag=FALSE), ", ",
           namesof("sd2", lsd, earg= esd2, tag=FALSE), "\n",
           "Mean:     phi*mu1 + (1-phi)*mu2\n",
           "Variance: phi*sd1^2 + (1-phi)*sd2^2 + phi*(1-phi)*(mu1-mu2)^2"),
    constraints=eval(substitute(expression({
        constraints = cm.vgam(rbind(diag(4), c(0,0,1,0)), x, .esd,
                              constraints, int=TRUE)
        constraints = cm.zero.vgam(constraints, x, .zero, M)
    }), list(.zero=zero, .esd=esd))),
    initialize=eval(substitute(expression({
        if(ncol(y <- cbind(y)) != 1)
            stop("the response must be a vector or one-column matrix")
        predictors.names = c(
            namesof("phi", .lphi, tag=FALSE),
            namesof("mu1", .lmu, earg= .emu1, tag=FALSE),
            namesof("sd1", .lsd, earg= .esd1, tag=FALSE),
            namesof("mu2", .lmu, earg= .emu2, tag=FALSE),
            namesof("sd2", .lsd, earg= .esd2, tag=FALSE))
        if(!length(etastart)) {
            qy = quantile(y, prob= .qmu)
            init.phi = if(length(.iphi)) rep(.iphi, length=n) else {
                0.5
            }
            init.mu1 = if(length(.imu1)) rep(.imu1, length=n) else {
                rep(qy[1], length=n)
            }
            init.mu2 = if(length(.imu2)) rep(.imu2, length=n) else {
                rep(qy[2], length=n)
            }
            ind.1 = if(init.mu1[1] < init.mu2[1]) 1:round(n* init.phi[1]) else
                round(n* init.phi[1]):n
            ind.2 = if(init.mu1[1] < init.mu2[1]) round(n* init.phi[1]):n else
                1:round(n* init.phi[1])
            sorty = sort(y)
            init.sd1 = if(length(.isd1)) rep(.isd1, length=n) else {
                sd(sorty[ind.1])
            }
            init.sd2 = if(length(.isd2)) rep(.isd2, length=n) else {
                sd(sorty[ind.2])
            }
            etastart = cbind(theta2eta(init.phi, .lphi, earg= .ephi),
                             theta2eta(init.mu1, .lmu, earg= .emu1),
                             theta2eta(init.sd1, .lsd, earg= .esd1),
                             theta2eta(init.mu2, .lmu, earg= .emu2),
                             theta2eta(init.sd2, .lsd, earg= .esd2))
        }
    }), list(.lphi=lphi, .lmu=lmu, .iphi=iphi, .imu1=imu1, .imu2=imu2,
             .ephi=ephi, .emu1=emu1, .emu2=emu2, .esd1=esd1, .esd2=esd2,
             .lsd=lsd, .isd1=isd1, .isd2=isd2, .qmu=qmu))),
    inverse=eval(substitute(function(eta, extra=NULL){
        phi = eta2theta(eta[,1], link= .lphi, earg= .ephi)
        mu1 = eta2theta(eta[,2], link= .lmu, earg= .emu1)
        mu2 = eta2theta(eta[,4], link= .lmu, earg= .emu2)
        phi*mu1 + (1-phi)*mu2
    }, list(.lphi=lphi, .lmu=lmu,
             .ephi=ephi, .emu1=emu1, .emu2=emu2, .esd1=esd1, .esd2=esd2 ))),
    last=eval(substitute(expression({
        misc$link = c("phi"= .lphi, "mu1"= .lmu,
                      "sd1"= .lsd, "mu2"= .lmu, "sd2"= .lsd)
        misc$earg = list("phi"= .ephi, "mu1"= .emu1,
                         "sd1"= .esd1, "mu2"= .emu2, "sd2"= .esd2)
        misc$expected = FALSE
        misc$esd = .esd
        misc$BFGS = TRUE
    }), list(.lphi=lphi, .lmu=lmu, .lsd=lsd, .esd=esd,
             .ephi=ephi, .emu1=emu1, .emu2=emu2, .esd1=esd1, .esd2=esd2 ))),
    loglikelihood=eval(substitute(
            function(mu,y,w,residuals=FALSE,eta,extra=NULL) {
        phi = eta2theta(eta[,1], link= .lphi, earg= .ephi)
        mu1 = eta2theta(eta[,2], link= .lmu, earg= .emu1)
        sd1 = eta2theta(eta[,3], link= .lsd, earg= .esd1)
        mu2 = eta2theta(eta[,4], link= .lmu, earg= .emu2)
        sd2 = eta2theta(eta[,5], link= .lsd, earg= .esd2)
        f1 = dnorm(y, mean=mu1, sd=sd1)
        f2 = dnorm(y, mean=mu2, sd=sd2)
        if(residuals) stop("loglikelihood residuals not implemented yet") else
        sum(w * log(phi*f1 + (1-phi)*f2))
    }, list(.lphi=lphi, .lmu=lmu,
            .ephi=ephi, .emu1=emu1, .emu2=emu2, .esd1=esd1, .esd2=esd2,
            .lsd=lsd ))),
    vfamily=c("mix2normal1"),
    deriv=eval(substitute(expression({
        phi = eta2theta(eta[,1], link= .lphi, earg= .ephi)
        mu1 = eta2theta(eta[,2], link= .lmu, earg= .emu1)
        sd1 = eta2theta(eta[,3], link= .lsd, earg= .esd1)
        mu2 = eta2theta(eta[,4], link= .lmu, earg= .emu2)
        sd2 = eta2theta(eta[,5], link= .lsd, earg= .esd2)
        f1 = dnorm(y, mean=mu1, sd=sd1)
        f2 = dnorm(y, mean=mu2, sd=sd2)
        pdf = phi*f1 + (1-phi)*f2
        df1.dmu1 = (y-mu1) * f1 / sd1^2
        df2.dmu2 = (y-mu2) * f2 / sd2^2
        dl.dphi = (f1-f2) / pdf
        dl.dmu1 = phi * df1.dmu1 / pdf
        dl.dmu2 = (1-phi) * df2.dmu2 / pdf
        dl.dsd1 = phi * f1 * (((y-mu1)/sd1)^2 - 1) / (sd1 * pdf)
        dl.dsd2 = (1-phi) * f2 * (((y-mu2)/sd2)^2 - 1) / (sd2 * pdf)
        dphi.deta = dtheta.deta(phi, link= .lphi, earg= .ephi)
        dmu1.deta = dtheta.deta(mu1, link= .lmu, earg= .emu1)
        dmu2.deta = dtheta.deta(mu2, link= .lmu, earg= .emu2)
        dsd1.deta = dtheta.deta(sd1, link= .lsd, earg= .esd1)
        dsd2.deta = dtheta.deta(sd2, link= .lsd, earg= .esd2)
        if(iter == 1) {
            etanew = eta
        } else {
            derivold = derivnew
            etaold = etanew
            etanew = eta
        }
        derivnew = w * cbind(dl.dphi * dphi.deta,
                             dl.dmu1 * dmu1.deta,
                             dl.dsd1 * dsd1.deta,
                             dl.dmu2 * dmu2.deta,
                             dl.dsd2 * dsd2.deta)
        derivnew
    }), list(.lphi=lphi, .lmu=lmu, .lsd=lsd,
             .ephi=ephi, .emu1=emu1, .emu2=emu2, .esd1=esd1, .esd2=esd2 ))),
    weight = eval(substitute(expression({
        if(iter == 1) {
            wznew = cbind(matrix(w, n, M), matrix(0, n, dimm(M)-M))
        } else {
            wzold = wznew
            wznew = qnupdate(w=w, wzold=wzold, dderiv=(derivold - derivnew),
                             deta=etanew-etaold, M=M,
                             trace=trace)  # weights incorporated in args
        }
        wznew
    }), list(.lphi=lphi, .lmu=lmu))))
}




mix2poisson.control <- function(save.weight=TRUE, ...)
{
    list(save.weight=save.weight)
}


mix2poisson = function(lphi="logit", llambda="loge",
                       ephi=list(), el1=list(), el2=list(),
                       iphi=0.5, il1=NULL, il2=NULL,
                       qmu=c(0.2, 0.8), zero=1)
{
    if(mode(lphi) != "character" && mode(lphi) != "name")
        lphi = as.character(substitute(lphi))
    if(mode(llambda) != "character" && mode(llambda) != "name")
        llambda = as.character(substitute(llambda))
    if(!is.Numeric(qmu, allow=2, positive=TRUE) || any(qmu >= 1))
        stop("bad input for argument \"qmu\"")
    if(length(iphi) && (!is.Numeric(iphi, allow=1, positive=TRUE) || iphi>= 1))
        stop("bad input for argument \"iphi\"")
    if(length(il1) && !is.Numeric(il1))
        stop("bad input for argument \"il1\"")
    if(length(il2) && !is.Numeric(il2))
        stop("bad input for argument \"il2\"")
    if(!is.list(ephi)) ephi = list()
    if(!is.list(el1)) el1 = list()
    if(!is.list(el2)) el2 = list()

    new("vglmff",
    blurb=c("Mixture of two univariate normals\n\n",
           "Links:    ",
           namesof("phi",lphi, earg= .ephi), ", ", 
           namesof("lambda1", llambda, earg= el1, tag=FALSE), ", ",
           namesof("lambda2", llambda, earg= el2, tag=FALSE), "\n",
           "Mean:     phi*lambda1 + (1-phi)*lambda2\n",
           "Variance: phi*lambda1^2 + (1-phi)*lambda2^2 + phi*(1-phi)*(lambda1-lambda2)^2"),
    constraints=eval(substitute(expression({
        constraints = cm.zero.vgam(constraints, x, .zero, M)
    }), list(.zero=zero ))),
    initialize=eval(substitute(expression({
        if(ncol(y <- cbind(y)) != 1)
            stop("the response must be a vector or one-column matrix")
        predictors.names = c(namesof("phi", .lphi, earg= .ephi, tag=FALSE),
                             namesof("lambda1", .llambda, earg= .el1, tag=FALSE),
                             namesof("lambda2", .llambda, earg= .el2, tag=FALSE))
        if(!length(etastart)) {
            qy = quantile(y, prob= .qmu)
            init.phi = if(length(.iphi)) rep(.iphi, length=n) else {
                0.5
            }
            init.lambda1 = if(length(.il1)) rep(.il1, length=n) else {
                rep(qy[1], length=n)
            }
            init.lambda2 = if(length(.il2)) rep(.il2, length=n) else {
                rep(qy[2], length=n)
            }
            etastart = cbind(theta2eta(init.phi, .lphi, earg= .ephi),
                             theta2eta(init.lambda1, .llambda, earg= .el1),
                             theta2eta(init.lambda2, .llambda, earg= .el2))
        }
    }), list(.lphi=lphi, .llambda=llambda, .iphi=iphi, .il1=il1, .il2=il2,
             .ephi=ephi, .el1=el1, .el2=el2,
             .qmu=qmu))),
    inverse=eval(substitute(function(eta, extra=NULL){
        phi = eta2theta(eta[,1], link= .lphi, earg= .ephi)
        lambda1 = eta2theta(eta[,2], link= .llambda, earg= .el1)
        lambda2 = eta2theta(eta[,3], link= .llambda, earg= .el2)
        phi*lambda1 + (1-phi)*lambda2
    }, list(.lphi=lphi, .llambda=llambda,
            .ephi=ephi, .el1=el1, .el2=el2 ))),
    last=eval(substitute(expression({
        misc$link = c("phi"= .lphi, "lambda1"= .llambda, "lambda2"= .llambda)
        misc$earg = list("phi"= .ephi, "lambda1"= .el1, "lambda2"= .el2)
        misc$expected = FALSE
        misc$BFGS = TRUE
    }), list(.lphi=lphi, .llambda=llambda,
             .ephi=ephi, .el1=el1, .el2=el2 ))),
    loglikelihood=eval(substitute(
            function(mu,y,w,residuals=FALSE,eta,extra=NULL) {
        phi = eta2theta(eta[,1], link= .lphi, earg= .ephi)
        lambda1 = eta2theta(eta[,2], link= .llambda, earg= .el1)
        lambda2 = eta2theta(eta[,3], link= .llambda, earg= .el2)
        f1 = dpois(y, lam=lambda1)
        f2 = dpois(y, lam=lambda2)
        if(residuals) stop("loglikelihood residuals not implemented yet") else
        sum(w * log(phi*f1 + (1-phi)*f2))
    }, list(.lphi=lphi, .llambda=llambda,
             .ephi=ephi, .el1=el1, .el2=el2 ))),
    vfamily=c("mix2poisson"),
    deriv=eval(substitute(expression({
        phi = eta2theta(eta[,1], link= .lphi, earg= .ephi)
        lambda1 = eta2theta(eta[,2], link= .llambda, earg= .el1)
        lambda2 = eta2theta(eta[,3], link= .llambda, earg= .el2)
        f1 = dpois(x=y, lam=lambda1)
        f2 = dpois(x=y, lam=lambda2)
        pdf = phi*f1 + (1-phi)*f2
        df1.dlambda1 = dpois(y-1, lam=lambda1) - f1
        df2.dlambda2 = dpois(y-1, lam=lambda2) - f2
        dl.dphi = (f1-f2) / pdf
        dl.dlambda1 = phi * df1.dlambda1 / pdf
        dl.dlambda2 = (1-phi) * df2.dlambda2 / pdf
        dphi.deta = dtheta.deta(phi, link= .lphi, earg= .ephi)
        dlambda1.deta = dtheta.deta(lambda1, link= .llambda, earg= .el1)
        dlambda2.deta = dtheta.deta(lambda2, link= .llambda, earg= .el2)
        if(iter == 1) {
            etanew = eta
        } else {
            derivold = derivnew
            etaold = etanew
            etanew = eta
        }
        derivnew = w * cbind(dl.dphi * dphi.deta,
                             dl.dlambda1 * dlambda1.deta,
                             dl.dlambda2 * dlambda2.deta)
        derivnew
    }), list(.lphi=lphi, .llambda=llambda,
             .ephi=ephi, .el1=el1, .el2=el2 ))),
    weight = eval(substitute(expression({
        if(iter == 1) {
            wznew = cbind(matrix(w, n, M), matrix(0, n, dimm(M)-M))
        } else {
            wzold = wznew
            wznew = qnupdate(w=w, wzold=wzold, dderiv=(derivold - derivnew),
                             deta=etanew-etaold, M=M,
                             trace=trace)  # weights incorporated in args
        }
        wznew
    }), list(.lphi=lphi, .llambda=llambda))))
}



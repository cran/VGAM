# These functions are
# Copyright (C) 1998-2008 T.W. Yee, University of Auckland. All rights reserved.





dcard = function(x, mu, rho) {
    if(!is.Numeric(mu) || any(mu < 0) || any(mu > 2*pi))
        stop("'mu' must be between 0 and 2*pi inclusive")
    if(!is.Numeric(rho) || max(abs(rho) > 0.5))
        stop("'rho' must be between -0.5 and 0.5 inclusive")
    L = max(length(x), length(mu), length(rho))
    x = rep(x, len=L); mu = rep(mu, len=L); rho = rep(rho, len=L);
    ans = (1 + 2 * rho * cos(x-mu)) / (2*pi)
    ans[x > (2*pi)] = 0
    ans[x < 0] = 0
    ans
}

pcard = function(q, mu, rho) {
    if(!is.Numeric(mu) || any(mu < 0) || any(mu > 2*pi))
        stop("'mu' must be between 0 and 2*pi inclusive")
    if(!is.Numeric(rho) || max(abs(rho) > 0.5))
        stop("'rho' must be between -0.5 and 0.5 inclusive")
    ans = (q + 2 * rho * (sin(q-mu) + sin(mu))) / (2*pi)
    ans[q >= (2*pi)] = 1
    ans[q <= 0] = 0
    ans
}

qcard = function(p, mu, rho, tolerance=1.0e-7, maxits=500) {
    if(!is.Numeric(mu) || any(mu < 0) || any(mu > 2*pi))
        stop("'mu' must be between 0 and 2*pi inclusive")
    if(!is.Numeric(rho) || max(abs(rho) > 0.5))
        stop("'rho' must be between -0.5 and 0.5 inclusive")
    if(!is.Numeric(p, positive=TRUE) || any(p > 1))
        stop("'p' must be between 0 and 1")
    nn = max(length(p), length(mu), length(rho))
    p = rep(p, len=nn)
    mu = rep(mu, len=nn)
    rho = rep(rho, len=nn)


    oldans = 2 * pi * p

    for(its in 1:maxits) {
        ans = oldans - (oldans + 2 * rho * (sin(oldans-mu)+sin(mu)) -
              2*pi*p) / (1 + 2 * rho * cos(oldans - mu))
        index = (ans <= 0) | (ans > 2*pi)
        if(any(index)) {
            ans[index] = runif(sum(index), 0, 2*pi)
        }
        if(max(abs(ans - oldans)) < tolerance) break;
        if(its == maxits) {warning("did not converge"); break}
        oldans = ans
    }
    ans
}

rcard = function(n, mu, rho, ...) {
    if(!is.Numeric(mu) || any(mu < 0) || any(mu > 2*pi))
        stop("'mu' must be between 0 and 2*pi inclusive")
    if(!is.Numeric(rho) || max(abs(rho) > 0.5))
        stop("'rho' must be between -0.5 and 0.5 inclusive")
    if(!is.Numeric(n, positive=TRUE, integer=TRUE, allow=1))
        stop("'n' must be a single positive integer")
    mu = rep(mu, len=n)
    rho = rep(rho, len=n)
    qcard(runif(n), mu=mu, rho=rho, ...)
}



cardioid.control <- function(save.weight=TRUE, ...)
{
    list(save.weight=save.weight)
}


cardioid = function(lmu="elogit", lrho="elogit",
                    emu=if(lmu=="elogit") list(min=0, max=2*pi) else list(),
                    erho=if(lmu=="elogit") list(min=-0.5, max=0.5) else list(),
                    imu=NULL, irho=0.3,
                    nsimEIM=100, zero=NULL)
{
    if(mode(lmu) != "character" && mode(lmu) != "name")
        lmu = as.character(substitute(lmu))
    if(mode(lrho) != "character" && mode(lrho) != "name")
        lrho = as.character(substitute(lrho))
    if(length(imu) && (!is.Numeric(imu, positive=TRUE) || any(imu > 2*pi)))
        stop("bad input for argument \"imu\"")
    if(!is.Numeric(irho) || max(abs(irho)) > 0.5)
        stop("bad input for argument \"irho\"")
    if(!is.list(emu)) emu = list()
    if(!is.list(erho)) erho = list()
    if(!is.Numeric(nsimEIM, allow=1, integ=TRUE) || nsimEIM <= 50)
        stop("'nsimEIM' should be an integer greater than 50")

    new("vglmff",
    blurb=c("Cardioid distribution\n\n",
           "Links:    ",
           namesof("mu", lmu, earg= emu), ", ", 
           namesof("rho", lrho, earg= erho, tag=FALSE), "\n",
           "Mean:     ",
           "pi + (rho/pi) *",
           "((2*pi-mu)*sin(2*pi-mu)+cos(2*pi-mu)-mu*sin(mu)-cos(mu))"),
    constraints=eval(substitute(expression({
        constraints = cm.zero.vgam(constraints, x, .zero, M)
    }), list( .zero=zero ))),
    initialize=eval(substitute(expression({
        if(ncol(y <- cbind(y)) != 1)
            stop("the response must be a vector or one-column matrix")
        if(any((y <= 0) | (y >=2*pi)))
            stop("the response must be in (0,2*pi)")
        predictors.names = c(
                       namesof("mu", .lmu, earg= .emu, tag=FALSE),
                       namesof("rho", .lrho, earg= .erho, tag=FALSE))
        if(!length(etastart)) {
            rho.init = rep(if(length(.irho)) .irho else 0.3, length=n)

            cardioid.Loglikfun = function(mu, y, x, w, extraargs) {
                rho = extraargs$irho
                sum(w * (-log(2*pi) + log1p(2*rho*cos(y-mu))))
            }
            mu.grid = seq(0.1, 6.0, len=19)
            mu.init = if(length( .imu )) .imu else
                getMaxMin(mu.grid, objfun=cardioid.Loglikfun, y=y,  x=x, w=w,
                          extraargs=list(irho = rho.init))
            mu.init = rep(mu.init, length=length(y))
            etastart = cbind(theta2eta(mu.init, .lmu, earg= .emu),
                             theta2eta(rho.init, .lrho, earg= .erho))
        }
    }), list( .lmu=lmu, .lrho=lrho,
              .imu=imu, .irho=irho,
              .emu=emu, .erho=erho ))),
    inverse=eval(substitute(function(eta, extra=NULL){
        mu = eta2theta(eta[,1], link= .lmu, earg= .emu)
        rho = eta2theta(eta[,2], link= .lrho, earg= .erho)
        pi + (rho/pi) *
        ((2*pi-mu)*sin(2*pi-mu) + cos(2*pi-mu) - mu*sin(mu) - cos(mu))
    }, list( .lmu=lmu, .lrho=lrho,
             .emu=emu, .erho=erho ))),
    last=eval(substitute(expression({
        misc$link = c("mu"= .lmu, "rho"= .lrho)
        misc$earg = list("mu"= .emu, "rho"= .erho)
        misc$expected = TRUE
        misc$nsimEIM = .nsimEIM
    }), list( .lmu=lmu, .lrho=lrho,
              .emu=emu, .erho=erho, .nsimEIM=nsimEIM ))),
    loglikelihood=eval(substitute(
            function(mu,y,w,residuals=FALSE,eta,extra=NULL) {
        mu = eta2theta(eta[,1], link= .lmu, earg= .emu)
        rho = eta2theta(eta[,2], link= .lrho, earg= .erho)
        if(residuals) stop("loglikelihood residuals not implemented yet") else
        sum(w * (-log(2*pi) + log1p(2*rho*cos(y-mu))))
    }, list( .lmu=lmu, .lrho=lrho,
             .emu=emu, .erho=erho ))),
    vfamily=c("cardioid"),
    deriv=eval(substitute(expression({
        mu = eta2theta(eta[,1], link= .lmu, earg= .emu)
        rho = eta2theta(eta[,2], link= .lrho, earg= .erho)
        dmu.deta = dtheta.deta(mu, link= .lmu, earg= .emu)
        drho.deta = dtheta.deta(rho, link= .lrho, earg= .erho)
        dl.dmu =  2 * rho * sin(y-mu) / (1 + 2 * rho * cos(y-mu))
        dl.drho = 2 * cos(y-mu) / (1 + 2 * rho * cos(y-mu))
        w * cbind(dl.dmu * dmu.deta,
                  dl.drho * drho.deta)
    }), list( .lmu=lmu, .lrho=lrho,
              .emu=emu, .erho=erho, .nsimEIM=nsimEIM ))),
    weight = eval(substitute(expression({
        run.varcov = 0
        ind1 = iam(NA, NA, M=M, both=TRUE, diag=TRUE)
        index0 = iam(NA, NA, M=M, both=TRUE, diag=TRUE)
        for(ii in 1:( .nsimEIM )) {
            ysim = rcard(n, mu=mu, rho=rho)
            dl.dmu =  2 * rho * sin(ysim-mu) / (1 + 2 * rho * cos(ysim-mu))
            dl.drho = 2 * cos(ysim-mu) / (1 + 2 * rho * cos(ysim-mu))
            rm(ysim)
            temp3 = cbind(dl.dmu, dl.drho)
            run.varcov = ((ii-1) * run.varcov +
                       temp3[,ind1$row.index]*temp3[,ind1$col.index]) / ii
        }
        wz = if(intercept.only)
            matrix(apply(run.varcov, 2, mean),
                   n, ncol(run.varcov), byrow=TRUE) else run.varcov

        dtheta.detas = cbind(dmu.deta, drho.deta)
        wz = wz * dtheta.detas[,index0$row] * dtheta.detas[,index0$col]
        w * wz
    }), list( .lmu=lmu, .lrho=lrho,
              .emu=emu, .erho=erho, .nsimEIM=nsimEIM ))))
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




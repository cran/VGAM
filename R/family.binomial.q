# These functions are
# Copyright (C) 1998-2006 T.W. Yee, University of Auckland. All rights reserved.



process.binomial2.data.vgam <- expression({


    if(!is.matrix(y))
    {
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
    } else
    if(ncol(y)==2)
    {
        if(!all(y==0 | y==1))
            stop("response must contains 0's and 1's only")
        col.index <- y[,2] + 2*y[,1] + 1    # 1:4
        nn <- nrow(y)
        y <- matrix(0, nn, 4)
        y[cbind(1:nn,col.index)] <- 1
        dimnames(y) <- list(dimnames(y)[[1]], c("00","01","10","11"))
        input.type <- 2
    } else
    if(ncol(y)==4)
    {
        input.type <- 3
    } else
        stop("response unrecognized")


    nvec <- drop(y %*% rep(1,ncol(y)))

    w <- w * nvec
    y <- y / nvec             # Convert to proportions

    mu <- y + (1/ncol(y) - y)/nvec
    dimnames(mu) <- dimnames(y)

})










betabinomial <- function(lmu="logit", lrho="logit", irho=0.5, zero=2)
{
    if(mode(lmu) != "character" && mode(lmu) != "name")
        lmu = as.character(substitute(lmu))
    if(mode(lrho) != "character" && mode(lrho) != "name")
        lrho = as.character(substitute(lrho))
    if(length(irho) && (!is.Numeric(irho, positive=TRUE) || max(irho) >= 1))
        stop("bad input for argument \"irho\"") 

    new("vglmff",
    blurb=c("Beta-binomial model\n",
           "Links:      ",
           namesof("mu", lmu), ", ",
           namesof("rho", lrho), "\n",
           "Variance:   mu*(1-mu)*(1+(w-1)*rho)/w"),
    constraints=eval(substitute(expression({
        constraints <- cm.zero.vgam(constraints, x, .zero, M)
    }), list( .zero=zero ))),
    initialize=eval(substitute(expression({
        eval(binomialff()@initialize)   # Note: n,w,y,mustart is changed 
        ycounts = y * w   # Convert proportions to counts
        if(max(abs(ycounts-round(ycounts))) > 1.0e-6)
           stop("the response (as counts) does not appear to be integer-valued")
        predictors.names = c(namesof("mu",  .lmu, tag=FALSE),
                             namesof("rho", .lrho, tag=FALSE))
        if(!length(etastart)) {
            if(is.Numeric( .irho )) {
                init.rho = rep( .irho, length=n)
            } else {
                init.rho = rep(0, length=n)
                Loglikfun = function(ycounts, nvec, shape1, shape2)
                if(is.R()) sum(lbeta(shape1+ycounts, shape2+nvec-ycounts) -
                               lbeta(shape1, shape2)) else
                sum(lgamma(shape1+ycounts) + lgamma(shape2+nvec-ycounts) -
                    lgamma(shape1+shape2+nvec) -
                    (lgamma(shape1) + lgamma(shape2) - lgamma(shape1+shape2)))
                rho.grid = rvar = seq(0.05, 0.95, len=11)  # 
                for(ii in 1:length(rho.grid))
                    rvar[ii] = Loglikfun(ycounts=y*w,
                        shape1=mustart*(1-rho.grid[ii])/rho.grid[ii],
                        shape2=(1-mustart)*(1-rho.grid[ii])/rho.grid[ii],
                        nvec=w)
                try.this = rho.grid[rvar == max(rvar)]
                init.rho = rep(try.this, len=n)
            }

            etastart = cbind(theta2eta(mustart,  .lmu),
                             theta2eta(init.rho, .lrho))
          }
    }), list( .lmu=lmu, .lrho=lrho, .irho=irho ))),
    inverse=eval(substitute(function(eta, extra=NULL)
        eta2theta(eta[,1], .lmu), 
    list( .lmu=lmu ))),
    last=eval(substitute(expression({
        misc$link <- c(mu = .lmu, rho = .lrho)
        misc$zero <- .zero
        misc$expected <- TRUE
    }), list( .lmu=lmu, .lrho=lrho, .zero=zero ))),
    loglikelihood=eval(substitute(
        function(mu,y,w,residuals=FALSE, eta, extra=NULL) {
        ycounts = y * w   # Convert proportions to counts
        mymu = eta2theta(eta[,1], .lmu)
        rho  = eta2theta(eta[,2], .lrho)
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
    }, list( .lmu=lmu, .lrho=lrho ))),
    vfamily=c("betabinomial"),
    deriv=eval(substitute(expression({
        nvec = w  # extra$nvec # for summary()
        ycounts = y * w   # Convert proportions to counts
        mymu = eta2theta(eta[,1], .lmu)
        rho  = eta2theta(eta[,2], .lrho)
        shape1 = mymu * (1 - rho) / rho
        shape2 = (1-mymu) * (1 - rho) / rho
        dshape1.dmu =  (1 - rho) / rho
        dshape2.dmu = -(1 - rho) / rho
        dshape1.drho = -mymu / rho^2
        dshape2.drho =  -(1 - mymu) / rho^2
        dmu.deta  = dtheta.deta(mymu, .lmu)
        drho.deta = dtheta.deta(rho,  .lrho)
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
    }), list( .lmu=lmu, .lrho=lrho ))),
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
    }), list( .lmu=lmu, .lrho=lrho ))))
}







binom2.or <- function(lp="logit", lp1=lp, lp2=lp, lor="loge",
                      zero=3, exchangeable=FALSE, tol=0.001)
{
    if(mode(lp1) != "character" && mode(lp1) != "name")
        lp1 <- as.character(substitute(lp1))
    if(mode(lp2) != "character" && mode(lp2) != "name")
        lp2 <- as.character(substitute(lp2))
    if(mode(lor) != "character" && mode(lor) != "name")
        lor <- as.character(substitute(lor))
    if(is.logical(exchangeable) && exchangeable && (lp1 != lp2))
        stop("exchangeable=TRUE but marginal links are not equal") 
    if(!is.Numeric(tol, positive=TRUE, allow=1))
        stop("bad input for argument \"tol\"") 

    new("vglmff",
    blurb=c("Palmgren model\n",
           "Links:    ",
           namesof("mu1", lp1), ", ",
           namesof("mu2", lp2), "; ",
           namesof("OR", lor)),
    constraints=eval(substitute(expression({
        constraints <- cm.vgam(matrix(c(1,1,0,0,0,1),3,2), x, 
                               .exchangeable, constraints,
                               intercept.apply=TRUE)
        constraints <- cm.zero.vgam(constraints, x, .zero, M)
    }), list( .exchangeable=exchangeable, .zero=zero ))),
    deviance=Deviance.categorical.data.vgam,
    initialize=eval(substitute(expression({
        eval(process.binomial2.data.vgam)
        predictors.names <- c(namesof("mu1", .lp1, short=TRUE), 
                              namesof("mu2", .lp2, short=TRUE), 
                              namesof("OR",  .lor, short=TRUE))
    }), list( .lp1=lp1, .lp2=lp2, .lor=lor ))),
    inverse=eval(substitute(function(eta, extra=NULL) {
        pm <- cbind(eta2theta(eta[,1], .lp1), eta2theta(eta[,2], .lp2))
        or <- eta2theta(eta[,3], .lor)
        a <- 1 + (pm[,1]+pm[,2])*(or-1)
        b <- -4 * or * (or-1) * pm[,1] * pm[,2]
        temp <- sqrt(a^2+b)
        pj4 <- ifelse(abs(or-1) < .tol, pm[,1]*pm[,2], (a-temp)/(2*(or-1)))
        pj2 <- pm[,2] - pj4
        pj3 <- pm[,1] - pj4
        cbind("00" = 1-pj4-pj2-pj3, "01" = pj2, "10" = pj3, "11" = pj4)
    }, list( .tol=tol, .lp1=lp1, .lp2=lp2, .lor=lor ))),
    last=eval(substitute(expression({
        misc$link <- c("mu1"= .lp1, "mu2"= .lp2, "OR"= .lor)
        misc$tol <- .tol
    }), list( .tol=tol, .lp1=lp1, .lp2=lp2, .lor=lor ))),
    link=eval(substitute(function(mu, extra=NULL) {
        pm <- cbind(mu[,3]+mu[,4], mu[,2]+mu[,4])
        or <- mu[,4]*mu[,1]/(mu[,2]*mu[,3])
        cbind(theta2eta(pm[,1], .lp1),
              theta2eta(pm[,2], .lp2), 
              theta2eta(or, .lor))
    }, list( .lp1=lp1, .lp2=lp2, .lor=lor ))),
    loglikelihood= function(mu, y, w, residuals = FALSE, eta, extra=NULL)
        if(residuals) stop("loglikelihood residuals not implemented yet") else
        sum(w * y * log(mu)),
    vfamily=c("binom2.or", "binom2"),
    deriv=eval(substitute(expression({
        pm <- cbind(mu[,3]+mu[,4], mu[,2]+mu[,4])
        or <- mu[,4]*mu[,1]/(mu[,2]*mu[,3])

        a <- 1 + (pm[,1]+pm[,2])*(or-1)
        b <- -4 * or * (or-1) * pm[,1] * pm[,2]
        temp <- sqrt(a^2+b)

        coeff <- -0.5 + (2*or*pm[,2]-a)/(2*temp)
        d1 <- coeff*(y[,1]/mu[,1]-y[,3]/mu[,3])-
           (1+coeff)*(y[,2]/mu[,2]-y[,4]/mu[,4])
    
        coeff <- -0.5 + (2*or*pm[,1]-a)/(2*temp)
        d2 <- coeff*(y[,1]/mu[,1]-y[,2]/mu[,2])-
           (1+coeff)*(y[,3]/mu[,3]-y[,4]/mu[,4])
    
        coeff <- (y[,1]/mu[,1]-y[,2]/mu[,2]-y[,3]/mu[,3]+y[,4]/mu[,4])
        d3 <- ifelse(abs(or-1) < .tol,
                 coeff * pm[,1] * (1-pm[,1]) * pm[,2] * (1-pm[,2]),
                 (1/(or-1)) * coeff * ( (pm[,1]+pm[,2])*(1-a/temp)/2 +
                 (2*or-1)*pm[,1]*pm[,2]/temp  - (a-temp)/(2*(or-1)) ))
        w * cbind(d1 * dtheta.deta(pm[,1], .lp1),
                  d2 * dtheta.deta(pm[,2], .lp2),
                  d3 * dtheta.deta(or, .lor))
    }), list( .tol=tol, .lp1=lp1, .lp2=lp2, .lor=lor ))),
    weight=eval(substitute(expression({
        Vab <- 1/(1/mu[,1] + 1/mu[,2] + 1/mu[,3] + 1/mu[,4])
        deltapi <- mu[,3]*mu[,2] - mu[,4]*mu[,1]
        delta  <- mu[,1]*mu[,2]*mu[,3]*mu[,4]
        pq <- pm[,1:2]*(1-pm[,1:2])

        wz <- matrix(0, n, 4)
        wz[,iam(1,1,M)] <- dtheta.deta(pm[,1], .lp1)^2 * pq[,2] * Vab / delta
        wz[,iam(2,2,M)] <- dtheta.deta(pm[,2], .lp2)^2 * pq[,1] * Vab / delta
        wz[,iam(3,3,M)] <- Vab * (dtheta.deta(or, .lor) / 
                                  dtheta.deta(or, "loge"))^2
        wz[,iam(1,2,M)] <- Vab * deltapi * dtheta.deta(pm[,1], .lp1) *
                           dtheta.deta(pm[,2], .lp2) / delta
       w * wz
    }), list( .lp1=lp1, .lp2=lp2, .lor=lor ))))
}






binom2.rho <- function(lrho="rhobit", init.rho=0.4, zero=3, exchangeable=FALSE)
{

    if(mode(lrho) != "character" && mode(lrho) != "name")
        lrho <- as.character(substitute(lrho))

    new("vglmff",
    blurb=c("Bivariate probit model\n",
           "Links:    ",
           "probit(mu1), probit(mu2); ",
            namesof("rho", lrho)),
    constraints=eval(substitute(expression({
        constraints <- cm.vgam(matrix(c(1,1,0,0,0,1),3,2), x, 
                               .exchangeable, constraints, intercept.apply=TRUE)
        constraints <- cm.zero.vgam(constraints, x, .zero, M)
    }), list( .exchangeable=exchangeable, .zero=zero ))),
    deviance=Deviance.categorical.data.vgam,
    initialize=eval(substitute(expression({
        eval(process.binomial2.data.vgam)
        predictors.names <- c("probit(mu1)", "probit(mu2)", 
                      namesof("rho", .lrho, short=TRUE))
        if(is.null(etastart))
            etastart <- cbind(theta2eta(mu[,3]+mu[,4], "probit"), 
                              theta2eta(mu[,2]+mu[,4], "probit"), 
                              theta2eta(.init.rho, .lrho))
    }), list( .lrho=lrho, .init.rho=init.rho ))),
    inverse=eval(substitute(function(eta, extra=NULL) {
        pm <- cbind(pnorm(eta[,1]),pnorm(eta[,2]))
        rho <- eta2theta(eta[,3], .lrho) 
        p11 <- pnorm2(eta[,1], eta[,2], rho)
        p01 <- pm[,2]-p11
        p10 <- pm[,1]-p11
        p00 <- 1-p01-p10-p11
        cbind("00"=p00, "01"=p01, "10"=p10, "11"=p11)
    }, list( .lrho=lrho ))),
    last=eval(substitute(expression({
        misc$link <- c(mu1 = "probit", mu2 = "probit", rho = .lrho)
    }), list( .lrho=lrho ))),
    link=eval(substitute(function(mu, extra=NULL) {
        if(is.null(extra))
            stop("rho must be passed into $link through \"extra\"")
        pm <- cbind(mu[,3]+mu[,4], mu[,2]+mu[,4])
        cbind("probit(mu1)"=qnorm(pm[,1]), 
              "probit(mu2)"=qnorm(pm[,2]),
              "link(rho)"=theta2eta(extra, .lrho))
    }, list( .lrho=lrho ))),
    vfamily=c("binom2.rho", "binom2"),
    deriv=eval(substitute(expression({
        pm <- cbind(pnorm(eta[,1]),pnorm(eta[,2]))
        rho <- eta2theta(eta[,3], .lrho) 
        p11 <- pnorm2(eta[,1], eta[,2], rho)
        p01 <- pm[,2]-p11
        p10 <- pm[,1]-p11
        p00 <- 1-p01-p10-p11

        B <- (eta[,2]-rho*eta[,1])/sqrt(1-rho^2)
        A <- (eta[,1]-rho*eta[,2])/sqrt(1-rho^2)
        phi1 <- dnorm(eta[,1])
        phi2 <- dnorm(eta[,2])
        PhiA <- pnorm(A)
        PhiB <- pnorm(B)

        ff <- dnorm2(eta[,1], eta[,2], rho)
        d1 = phi1*(PhiB*(y[,4]/p11-y[,2]/p01) + (1-PhiB)*(y[,3]/p10-y[,1]/p00))
        d2 = phi2*(PhiA*(y[,4]/p11-y[,3]/p10) + (1-PhiA)*(y[,2]/p01-y[,1]/p00))
        dl.drho <- (y[,4]/p11-y[,3]/p10-y[,2]/p01+y[,1]/p00)* ff
        drho.deta <- dtheta.deta(rho, .lrho)
        w * cbind(d1, d2, dl.drho * drho.deta)
    }), list( .lrho=lrho ))),
    weight=eval(substitute(expression({
        wz <- matrix(as.numeric(NA), n, dimm(M))  # 6=dimm(M)
        wz[,iam(1,1,M)] = phi1^2*(PhiB^2*(1/p11+1/p01)+(1-PhiB)^2*(1/p10+1/p00))
        wz[,iam(2,2,M)] = phi2^2*(PhiA^2*(1/p11+1/p10)+(1-PhiA)^2*(1/p01+1/p00))
        wz[,iam(1,2,M)] = phi1*phi2*(PhiA*PhiB/p11 + (1-PhiA)*(1-PhiB)/p00 -
                    PhiA*(1-PhiB)/p10 - (1-PhiA)*PhiB/p01)
        d2l.drhoeta1 <- ff*phi1*(PhiB*(1/p11+1/p01) - (1-PhiB)*(1/p10+1/p00))
        wz[,iam(1,3,M)] <- d2l.drhoeta1 * drho.deta
        d2l.drhoeta2 <- ff*phi2*(PhiA*(1/p11+1/p10) - (1-PhiA)*(1/p01+1/p00))
        wz[,iam(2,3,M)] <- d2l.drhoeta2 * drho.deta
        d2l.drho2 <- ff^2 * (1/p11+1/p01+1/p10+1/p00)
        wz[,iam(3,3,M)] <-  d2l.drho2 * drho.deta^2
        wz * w
    }), list( .lrho=lrho ))))
}



dnorm2 <- function(x, y, r) 
    exp(-0.5*(x^2+y^2-2*x*y*r)/(1-r^2)) / (2*pi*sqrt(1-r^2))


pnorm2 <- function(ah, ak, r) 
{ 

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
              x * log(prob/(1-prob)) + size * log(1-prob) )
}



size.binomial <- function(prob=0.5, link="loge")
{
    if(any(prob<=0 || prob>=1))
        stop("some values of prob out of range")
    if(!missing(link)) link <- as.character(substitute(link))

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
            y * log(.prob / (1- .prob)) + nvec * log(1- .prob)))
    }, list( .prob=prob ))),
    vfamily=c("size.binomial"),
    deriv=eval(substitute(expression({
        nvec <- mu/extra$temp2
        dldnvec = digamma(nvec+1) - digamma(nvec-y+1) + log(1-extra$temp2)
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
    if(!is.Numeric(x)) stop("bad input for argument \"x\"")
    if(!is.Numeric(size, posit=TRUE, integer=TRUE))
        stop("bad input for argument \"size\"")
    if(!is.Numeric(shape1, pos=TRUE)) stop("bad input for argument \"shape1\"")
    if(!is.Numeric(shape2, pos=TRUE)) stop("bad input for argument \"shape2\"")
    L = max(length(x), length(size), length(shape1), length(shape2))
    x = rep(x, len=L); size = rep(size, len=L);
    shape1 = rep(shape1, len=L); shape2 = rep(shape2, len=L);
    answer = 0 * x
    if(length(ok <- round(x) == x & x >= 0 & x <= size))
        answer[ok] = if(log) lchoose(size[ok], x[ok]) +
                     lbeta(shape1[ok]+x[ok], shape2[ok]+size[ok]-x[ok]) -
                     lbeta(shape1[ok], shape2[ok]) else 
                     choose(size[ok], x[ok]) *
                     beta(shape1[ok]+x[ok], shape2[ok]+size[ok]-x[ok]) /
                     beta(shape1[ok], shape2[ok])
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
    rbetabin.ab(x=x, size=size, shape1=prob*(1-rho)/rho,
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




betabin.ab = function(link.shape12="loge", i1=1, i2=NULL, zero=NULL)
{
    if(mode(link.shape12) != "character" && mode(link.shape12) != "name")
        link.shape12 = as.character(substitute(link.shape12))
    if(!is.Numeric(i1, positive=TRUE)) stop("bad input for argument \"i1\"")
    if(length(i2) && !is.Numeric(i2, pos=TRUE))
        stop("bad input for argument \"i2\"")

    new("vglmff",
    blurb=c("Beta-binomial model\n",
           "Links:    ",
           namesof("shape1", link.shape12), ", ",
           namesof("shape2", link.shape12), "\n",
           "Variance: mu*(1-mu)[1+(w-1)*rho]/w where mu=alpha/(alpha+beta)"),
    constraints=eval(substitute(expression({
        constraints <- cm.zero.vgam(constraints, x, .zero, M)
    }), list( .zero=zero ))),
    initialize=eval(substitute(expression({
        # Compute initial values for mustart -------
        eval(binomialff()@initialize)   # Note: n,w,y,mustart is changed 
        predictors.names = c(namesof("shape1", .link.shape12, tag=FALSE),
                             namesof("shape2", .link.shape12, short=FALSE))

        if(!length(etastart)) {
            shape1 = rep( .i1, len=n)
            shape2 = if(length( .i2)) rep( .i2,len=n) else shape1*(1/mustart-1)
            ycounts = y * w   # Convert proportions to counts
            if(max(abs(ycounts-round(ycounts))) > 1.0e-6)
                stop("ycounts not integer")
            ycounts = round(ycounts) # Make sure it is an integer
            etastart = cbind(theta2eta(shape1, .link.shape12),
                             theta2eta(shape2, .link.shape12))
        }
    }), list( .link.shape12=link.shape12, .i1=i1 , .i2=i2 ))),
    inverse=eval(substitute(function(eta, extra=NULL) {
        shape1 = eta2theta(eta[,1], .link.shape12)
        shape2 = eta2theta(eta[,2], .link.shape12)
        shape1 / (shape1 + shape2)
    }, list( .link.shape12=link.shape12 ))),
    last=eval(substitute(expression({
        misc$link = c("shape1" = .link.shape12, "shape2" = .link.shape12)
        shape1 = eta2theta(eta[,1], .link.shape12)
        shape2 = eta2theta(eta[,2], .link.shape12)
        misc$rho = 1 / (shape1 + shape2 + 1)
        misc$expected = TRUE
    }), list( .link.shape12=link.shape12 ))),
    loglikelihood=eval(substitute(
        function(mu,y,w,residuals=FALSE,eta, extra=NULL) {
        ycounts = y * w   # Convert proportions to counts
        shape1 = eta2theta(eta[,1], .link.shape12)
        shape2 = eta2theta(eta[,2], .link.shape12)
        nvec = w
        if(residuals) stop("loglikelihood residuals not implemented yet") else {
            if(is.R()) sum(lbeta(shape1+ycounts, shape2+nvec-ycounts) -
                           lbeta(shape1, shape2)) else
            sum(lgamma(shape1+ycounts) + lgamma(shape2+nvec-ycounts) -
                lgamma(shape1+shape2+nvec) -
                (lgamma(shape1) + lgamma(shape2) - lgamma(shape1+shape2)))
        }
    }, list( .link.shape12=link.shape12 ))),
    vfamily=c("betabin.ab"),
    deriv=eval(substitute(expression({
        nvec = w  # extra$nvec # for summary()
        ycounts = y * w   # Convert proportions to counts
        shape1 = eta2theta(eta[,1], .link.shape12)
        shape2 = eta2theta(eta[,2], .link.shape12)
        dshape1.deta = dtheta.deta(shape1, .link.shape12) 
        dshape2.deta = dtheta.deta(shape2, .link.shape12) 
        dl.dshape1 = digamma(shape1+ycounts) - digamma(shape1+shape2+nvec) -
                     digamma(shape1) + digamma(shape1+shape2)
        dl.dshape2 = digamma(nvec+shape2-ycounts) -
                     digamma(shape1+shape2+nvec) -
                     digamma(shape2) + digamma(shape1+shape2)
        cbind(dl.dshape1 * dshape1.deta, dl.dshape2 * dshape2.deta)
    }), list( .link.shape12=link.shape12 ))),
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
    }), list( .link.shape12=link.shape12 ))))
}



betageometric = function(lprob="logit", lshape="loge",
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

    new("vglmff",
    blurb=c("Beta-geometric distribution\n",
           "Links:    ", namesof("prob", lprob), ", ",
                         namesof("shape", lshape)),
    constraints=eval(substitute(expression({
        constraints <- cm.zero.vgam(constraints, x, .zero, M)
    }), list( .zero=zero ))),
    initialize=eval(substitute(expression({
        eval(geometric()@initialize)
        predictors.names = c(namesof("prob",  .lprob,  tag=FALSE),
                             namesof("shape", .lshape, short=FALSE))
        if(length( .iprob))
            prob.init = rep( .iprob, len=n)
        if(!length(etastart) || ncol(cbind(etastart)) != 2) {
            shape.init = rep( .ishape, len=n)
            etastart = cbind(theta2eta(prob.init,  .lprob),
                             theta2eta(shape.init, .lshape))
        }
    }), list( .iprob=iprob, .ishape=ishape, .lprob=lprob, .lshape=lshape ))),
    inverse=eval(substitute(function(eta, extra=NULL) {
        prob  = eta2theta(eta[,1], .lprob)
        shape = eta2theta(eta[,2], .lshape)
        mymu = (1-prob) / (prob - shape)
        ifelse(mymu >= 0, mymu, NA)
    }, list( .lprob=lprob, .lshape=lshape ))),
    last=eval(substitute(expression({
        misc$link = c("prob" = .lprob, "shape" = .lshape)
        if(intercept.only) {
            misc$shape1 = shape1[1]  # These quantities computed in @deriv
            misc$shape2 = shape2[1]
        }
        misc$expected = TRUE
        misc$tolerance = .tolerance
        misc$zero = .zero
        misc$moreSummation = .moreSummation
    }), list( .lprob=lprob, .lshape=lshape, .tolerance=tolerance,
              .moreSummation=moreSummation, .zero=zero ))),
    loglikelihood=eval(substitute(
        function(mu,y,w,residuals=FALSE,eta, extra=NULL) {
        prob  = eta2theta(eta[,1], .lprob)
        shape = eta2theta(eta[,2], .lshape)
        ans = log(prob)
        maxy = max(y)
        if(residuals) stop("loglikelihood residuals not implemented yet") else {
            for(ii in 1:maxy) {
                index = ii <= y
                ans[index]=ans[index] + log(1-prob[index]+(ii-1)*shape[index])-
                           log(1+(ii-1)*shape[index])
            }
            ans = ans - log(1+(y+1-1)*shape)
            sum(w * ans)
        }
    }, list( .lprob=lprob, .lshape=lshape ))),
    vfamily=c("betageometric"),
    deriv=eval(substitute(expression({
        prob  = eta2theta(eta[,1], .lprob)
        shape = eta2theta(eta[,2], .lshape)
        shape1 = prob / shape; shape2 = (1-prob) / shape;
        dprob.deta  = dtheta.deta(prob,  .lprob) 
        dshape.deta = dtheta.deta(shape, .lshape) 
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
    }), list( .lprob=lprob, .lshape=lshape ))),
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
              .tolerance=tolerance ))))
}








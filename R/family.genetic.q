# These functions are
# Copyright (C) 1998-2007 T.W. Yee, University of Auckland. All rights reserved.








G1G2G3 <- function(link="logit", earg = list(), ip1=NULL, ip2=NULL, iF=NULL)
{
    if(mode(link) != "character" && mode(link) != "name")
        link <- as.character(substitute(link))
    if(!is.list(earg)) earg = list()

    new("vglmff",
    blurb=c("G1/G2/G3 zzphenotype\n\n",
            "Links:    ",
            namesof("p1", link, earg= earg), ", ", 
            namesof("p2", link, earg= earg), ", ", 
            namesof("f",  link, earg= earg, tag=FALSE),
           "\n",
           "Variance: multinomial type variance"),
    initialize=eval(substitute(expression({
        delete.zero.colns <- FALSE
        eval(process.categorical.data.vgam)
        predictors.names <- c(namesof("p1", .link, earg= .earg, tag=FALSE),
                      namesof("p2", .link, earg= .earg, tag=FALSE),
                      namesof("f",  .link, earg= .earg, tag=FALSE))
        if(is.null(etastart)) {
            p1 <- if(is.numeric(.ip1)) rep(.ip1, n) else
                       sqrt(mustart[,1])
            f <- if(is.numeric(.iF)) rep(.iF, n) else
                rep(0.01, n) # close to zero 
            p2 <- if(is.numeric(.ip2)) rep(.ip2, n) else
                       mustart[,2] / (sqrt(mustart[,1]) * 2)
            if(any(p1 <= 0) || any(p1 >= 1))
                stop("bad initial value for p1")
            if(any(p2 <= 0) || any(p2 >= 1))
                stop("bad initial value for p2")
            etastart <- cbind(theta2eta(p1, .link, earg= .earg),
                              theta2eta(p2, .link, earg= .earg),
                              theta2eta(f,  .link, earg= .earg))
        }
    }), list(.link=link, .ip1=ip1, .ip2=ip2, .iF=iF,
             .earg=earg ))),
    inverse=eval(substitute(function(eta, extra=NULL){
        p1 <- eta2theta(eta[,1], link=.link, earg= .earg)
        p2 <- eta2theta(eta[,2], link=.link, earg= .earg)
        p3 <- 1-p1-p2
        f  <- eta2theta(eta[,3], link=.link, earg= .earg)
        cbind("G1/G1"=f*p1+(1-f)*p1^2,
              "G1/G2"=2*p1*p2*(1-f),
              "G1/G3"=2*p1*p3*(1-f),
              "G2/G2"=f*p2+(1-f)*p2^2,
              "G2/G3"=2*p2*p3*(1-f),
              "G3/G3"=f*p3+(1-f)*p3^2)
    }, list(.link=link,
             .earg=earg ))),
    last=eval(substitute(expression({
        misc$link <- c(p1= .link, p2= .link, f= .link)
        misc$earg = list(p1= .earg, p2= .earg, f= .earg )
    }), list(.link=link,
             .earg=earg ))),
    loglikelihood=function(mu,y,w,residuals=FALSE,eta,extra=NULL)
        if(residuals) stop("loglikelihood residuals not implemented yet") else
        sum((w*y)*log(mu)),
    vfamily=c("G1G2G3", "vgenetic"),
    deriv=eval(substitute(expression({
        p1 <- eta2theta(eta[,1], link=.link, earg= .earg)
        p2 <- eta2theta(eta[,2], link=.link, earg= .earg)
        p3 <- 1-p1-p2
        f  <- eta2theta(eta[,3], link=.link, earg= .earg)
        dP1 <- cbind(f + 2*p1*(1-f), 2*(1-f)*p2, 2*(1-f)*(1-p2-2*p1), 0,
                     -2*(1-f)*p2, -f - 2*p3*(1-f))
        dP2 <- cbind(0, 2*p1*(1-f), -2*(1-f)*p1, f+2*p2*(1-f),
                     2*(1-f)*(1-p1-2*p2), -f - 2*p3*(1-f))
        dP3 <- cbind(p1*(1-p1), -2*p1*p2, -2*p1*p3, p2*(1-p2), -2*p2*p3, 
                     p3*(1-p3))
        dl1 <- apply(y * dP1 / mu, 1, sum)
        dl2 <- apply(y * dP2 / mu, 1, sum)
        dl3 <- apply(y * dP3 / mu, 1, sum)
        dPP.deta <- dtheta.deta(cbind(p1,p2,f), link=.link, earg= .earg)
        w * cbind(dPP.deta[,1] * dl1, dPP.deta[,2] * dl2, 
                  dPP.deta[,3] * dl3)
    }), list(.link=link,
             .earg=earg ))),
    weight=eval(substitute(expression({
        dPP <- array(c(dP1,dP2,dP3), c(n,6,3))

        wz <- matrix(as.numeric(NA), n, dimm(M))   # dimm(M)==6 because M==3
        for(i1 in 1:M)
            for(i2 in i1:M) {
                index <- iam(i1,i2,M)
                wz[,index] <- apply(dPP[,,i1,drop=TRUE] * dPP[,,i2,drop=TRUE] /
                                mu, 1, sum) * dPP.deta[,i1] * dPP.deta[,i2]
        }
        w * wz
    }), list(.link=link,
             .earg=earg ))))
}



AAaa.nohw <- function(link="logit", earg = list(), ipA=NULL, iF=NULL)
{

    if(mode(link) != "character" && mode(link) != "name")
        link <- as.character(substitute(link))
    if(!is.list(earg)) earg = list()

    new("vglmff",
    blurb=c("AA-Aa-aa phenotype (without Hardy-Weinberg assumption)\n\n",
            "Links:    ",
            namesof("pA", link, earg= earg), ", ", 
            namesof("f",  "identity", tag=FALSE),
           "\n",
           "Variance: multinomial type variance"),
    initialize=eval(substitute(expression({
        delete.zero.colns <- FALSE
        eval(process.categorical.data.vgam)
        predictors.names <- c(namesof("pA", .link, earg= .earg, tag=FALSE),
                      namesof("f",  "identity", tag=FALSE))
        if(is.null(etastart)) {
            pA <- if(is.numeric(.ipA)) rep(.ipA, n) else
                       c(sqrt(mustart[,1] - mustart[,2]/2))
            f <- if(is.numeric(.iF)) rep(.iF, n) else
                rep(0.01, n) # 1- mustart[,2]/(2*pA*(1-pA))
            if(any(pA <= 0) || any(pA >= 1))
                stop("bad initial value for pA")
            etastart <- cbind(theta2eta(pA, .link, earg= .earg),
                              theta2eta(f,  "identity"))
        }
    }), list(.link=link, .ipA=ipA, .iF=iF,
             .earg=earg ))),
    inverse=eval(substitute(function(eta, extra=NULL){
        pA <- eta2theta(eta[,1], link=.link, earg= .earg)
        f  <- eta2theta(eta[,2], link="identity")
        cbind(AA=pA^2+pA*(1-pA)*f, Aa=2*pA*(1-pA)*(1-f),
              aa=(1-pA)^2 + pA*(1-pA)*f)
    }, list(.link=link,
             .earg=earg ))),
    last=eval(substitute(expression({
        misc$link <- c(pA= .link, f= "identity")
        misc$earg = list(pA= .earg, f= list() )
    }), list(.link=link,
             .earg=earg ))),
    link=eval(substitute(function(mu, extra=NULL){
        pA <- sqrt(mu[,1] - mu[,2]/2)
        f <- 1- mu[,2]/(2*pA*(1-pA))
        cbind(theta2eta(pA, .link, earg= .earg),
              theta2eta(f,  "identity"))
    }, list(.link=link,
             .earg=earg ))),
    loglikelihood=function(mu,y,w,residuals=FALSE,eta,extra=NULL)
        if(residuals) stop("loglikelihood residuals not implemented yet") else
        sum((w*y)*log(mu)),
    vfamily=c("AAaa.nohw", "vgenetic"),
    deriv=eval(substitute(expression({
        pA <- eta2theta(eta[,1], link=.link, earg= .earg)
        f  <- eta2theta(eta[,2], link="identity")
        dP1 <- cbind(f + 2*pA*(1-f), 2*(1-f)*(1-2*pA), -2*(1-pA) +f*(1-2*pA))
        dP2 <- cbind(pA*(1-pA), -2*pA*(1-pA), pA*(1-pA))
        dl1 <- apply(y * dP1 / mu, 1, sum)
        dl2 <- apply(y * dP2 / mu, 1, sum)
        dPP.deta <- dtheta.deta(pA, link=.link, earg= .earg)
        w * cbind(dPP.deta * dl1, dl2)
    }), list(.link=link,
             .earg=earg ))),
    weight=eval(substitute(expression({
        dPP <- array(c(dP1,dP2), c(n,3,2))
        dPP.deta <- cbind(dtheta.deta(pA, link=.link, earg= .earg),
                          dtheta.deta(f,  link="identity"))
        wz <- matrix(as.numeric(NA), n, dimm(M))   # dimm(M)==3 because M==2
        for(i1 in 1:M)
            for(i2 in i1:M) {
                index <- iam(i1,i2,M)
                wz[,index] <- apply(dPP[,,i1,drop=T] * dPP[,,i2,drop=T] /
                                mu, 1, sum) * dPP.deta[,i1] * dPP.deta[,i2]
        }
        w * wz
    }), list(.link=link,
             .earg=earg ))))
}


AB.Ab.aB.ab2 <- function(link="logit", earg = list(), init.p=NULL)
{
    if(mode(link) != "character" && mode(link) != "name")
        link <- as.character(substitute(link))
    if(!is.list(earg)) earg = list()

    new("vglmff",
    blurb=c("AB-Ab-aB-ab2 phenotype\n\n",
            "Links:    ",
            namesof("p", link, earg= earg),
            "\n",
            "Variance: multinomial type variance"),
    initialize=eval(substitute(expression({
        delete.zero.colns <- FALSE
        eval(process.categorical.data.vgam)
        predictors.names <- namesof("p", .link, earg= .earg, tag=FALSE)
      
        if(is.null(etastart)) {
            p.init <- if(is.numeric(.init.p)) rep(.init.p, n) else
                      c(1 - 2 * sqrt(mustart[,4]))
            etastart <- theta2eta(p.init, .link, earg= .earg)
        }
    }), list(.link=link, .init.p=init.p,
             .earg=earg ))),
    inverse=eval(substitute(function(eta,extra=NULL){
        p <- eta2theta(eta, link=.link, earg= .earg)
        cbind("A-B-"=(2+(1-p)^2),
              "A-bb"=(1-(1-p)^2),
              "aaB-"=(1-(1-p)^2),
              "aabb"=(1-p)^2) / 4
    }, list(.link=link,
             .earg=earg ) )),
    last=eval(substitute(expression({
        misc$link <- c(p = .link)
        misc$earg = list(p= .earg )
    }), list(.link=link,
             .earg=earg ) )),
    link=eval(substitute(function(mu, extra=NULL){
        p <- 1 - 2 * sqrt(mu[,4])
        theta2eta(p, .link, earg= .earg)
    }, list(.link=link,
             .earg=earg ) )),
    loglikelihood= function(mu,y,w,residuals=FALSE,eta,extra=NULL)
        if(residuals) stop("loglikelihood residuals not implemented yet") else
        sum((w*y)*log(mu)),
    vfamily=c("AB.Ab.aB.ab2", "vgenetic"),
    deriv=eval(substitute(expression({
        pp <- eta2theta(eta, link=.link, earg= .earg)
        dP1 <- cbind(-0.5*(1-pp), 0.5*(1-pp), 0.5*(1-pp), -0.5*(1-pp))
        dl1 <- apply(y * dP1 / mu, 1, sum)
        dPP.deta <- dtheta.deta(pp, link=.link, earg= .earg)
        w * dPP.deta * dl1
        }), list(.link=link,
             .earg=earg ) )),
    weight=eval(substitute(expression({
        wz <- apply(dP1 * dP1 / mu, 1, sum) * dPP.deta^2
        w * wz
    }), list(.link=link,
             .earg=earg ) )))
}



A1A2A3 <- function(link="logit", earg = list(), ip1=NULL, ip2=NULL)
{
    if(mode(link) != "character" && mode(link) != "name")
        link <- as.character(substitute(link))
    if(!is.list(earg)) earg = list()

    new("vglmff",
    blurb=c("A1A2A3 Allele System (A1A1, A1A2, A2A2, A1A3, A2A3, A3A3)\n\n",
            "Links:    ",
            namesof("p1",link, earg= earg), ", ", 
            namesof("p2", link, earg= earg, tag=FALSE),
           "\n",
           "Variance: multinomial type variance"),
    initialize=eval(substitute(expression({
        delete.zero.colns <- FALSE
        eval(process.categorical.data.vgam)
        predictors.names <- c(namesof("pA", .link, earg= .earg,tag=FALSE),
                              namesof("pB", .link, earg= .earg,tag=FALSE))
        if(is.null(etastart)) {
            p1 <- if(is.numeric(.ip1)) rep(.ip1, n) else
                       c(sqrt(mustart[,1]))
            p2 <- if(is.numeric(.ip2)) rep(.ip2, n) else
                       c(sqrt(mustart[,3]))
            etastart <- cbind(theta2eta(p1, .link, earg= .earg),
                              theta2eta(p2, .link, earg= .earg))
        }
    }), list(.link=link, .ip1=ip1, .ip2=ip2,
             .earg=earg ))),
    inverse=eval(substitute(function(eta, extra=NULL){
        p1 <- eta2theta(eta[,1], link=.link, earg= .earg)
        p2 <- eta2theta(eta[,2], link=.link, earg= .earg)
        qq <- 1-p1-p2
        cbind(A1A1=p1*p1, A1A2=2*p1*p2, A2A2=p2*p2, A1A3=2*p1*qq,
              A2A3=2*p2*qq, A3A3=qq*qq)
    }, list(.link=link,
             .earg=earg ))),
    last=eval(substitute(expression({
        misc$link <- c(p1= .link, p2= .link)
        misc$earg = list(p1= .earg, p2= .earg )
    }), list(.link=link,
             .earg=earg ))),
    link=eval(substitute(function(mu, extra=NULL){
        p1 <- sqrt(mu[,1])
        p2 <- sqrt(mu[,3])
        qq <- 1-p1-p2
        cbind(theta2eta(p1, .link, earg= .earg),
              theta2eta(p2, .link, earg= .earg))
    }, list(.link=link,
             .earg=earg ))),
    loglikelihood=function(mu,y,w,residuals=FALSE,eta,extra=NULL)
        if(residuals) stop("loglikelihood residuals not implemented yet") else
        sum((w*y)*log(mu)),
    vfamily=c("A1A2A3", "vgenetic"),
    deriv=eval(substitute(expression({
        p1 <- eta2theta(eta[,1], link=.link, earg= .earg)
        p2 <- eta2theta(eta[,2], link=.link, earg= .earg)
        dl.dp1 <- (2*y[,1]+y[,2]+y[,4])/p1 - (2*y[,6]+y[,4]+y[,5])/(1-p1-p2)
        dl.dp2 <- (2*y[,3]+y[,2]+y[,5])/p2 - (2*y[,6]+y[,4]+y[,5])/(1-p1-p2)
        dp1.deta <- dtheta.deta(p1, link=.link, earg= .earg)
        dp2.deta <- dtheta.deta(p2, link=.link, earg= .earg)
        w * cbind(dl.dp1 * dp1.deta, dl.dp2 * dp2.deta)
    }), list(.link=link,
             .earg=earg ))),
    weight=eval(substitute(expression({
        qq <- 1-p1-p2
        wz <- matrix(as.numeric(NA), n, dimm(M))   # dimm(M)==3 because M==2
        ed2l.dp12  <-  2 * (1/p1 + 1/qq)
        ed2l.dp22  <-  2 * (1/p2 + 1/qq)
        ed2l.dp1dp2 <-  2 / qq
        wz[,iam(1,1,M)] <- dp1.deta^2 * ed2l.dp12
        wz[,iam(2,2,M)] <- dp2.deta^2 * ed2l.dp22
        wz[,iam(1,2,M)] <- ed2l.dp1dp2 * dp1.deta * dp2.deta
        w * wz
    }), list(.link=link,
             .earg=earg ))))
}




MNSs <- function(link="logit", earg = list(), imS=NULL, ims=NULL, inS=NULL)
{
    if(mode(link) != "character" && mode(link) != "name")
        link <- as.character(substitute(link))
    if(!is.list(earg)) earg = list()

    new("vglmff",
    blurb=c("MNSs Blood Group System (MS-Ms-MNS-MNs-NS-Ns phenotype)\n\n",
            "Links:    ",
            namesof("mS",link, earg= earg), ", ", 
            namesof("ms",link, earg= earg), ", ", 
            namesof("nS", link, earg= earg, tag=FALSE),
           "\n",
           "Variance: multinomial type variance"),
    initialize=eval(substitute(expression({
        delete.zero.colns <- FALSE
        eval(process.categorical.data.vgam)
        predictors.names <-
           c(namesof("mS", .link, earg= .earg,tag=FALSE),
             namesof("ms", .link, earg= .earg,tag=FALSE),
             namesof("nS", .link, earg= .earg,tag=FALSE))
        if(is.null(etastart)) {
            ms <- if(is.numeric(.ims)) rep(.ims, n) else
                       c(sqrt(mustart[,2]))
            ns <- c(sqrt(mustart[,6]))
            nS <- if(is.numeric(.inS)) rep(.inS, n) else
                c(-ns + sqrt(ns^2 + mustart[,5]))  # Solve a quadratic eqn
            mS <- if(is.numeric(.imS)) rep(.imS, n) else
                    1-ns-ms-nS
            etastart <- cbind(theta2eta(mS, .link, earg= .earg),
                              theta2eta(ms, .link, earg= .earg),
                              theta2eta(nS, .link, earg= .earg))
        }
    }), list(.link=link, .imS=imS, .ims=ims, .inS=inS,
             .earg=earg ))),
    inverse=eval(substitute(function(eta, extra=NULL){
        mS <- eta2theta(eta[,1], link=.link, earg= .earg)
        ms <- eta2theta(eta[,2], link=.link, earg= .earg)
        nS <- eta2theta(eta[,3], link=.link, earg= .earg)
        ns <- 1-mS-ms-nS
       cbind(MS=mS^2+2*mS*ms, Ms=ms^2, MNS=2*(mS*nS+ms*nS+mS*ns),
             MNs=2*ms*ns, NS=nS^2 + 2*nS*ns, Ns=ns^2)
    }, list(.link=link,
             .earg=earg ))),
    last=eval(substitute(expression({
        misc$link <- c(mS= .link, ms= .link, nS= .link)
        misc$earg = list(mS= .earg, ms= .earg, nS= .earg )
    }), list(.link=link,
             .earg=earg ))),
    link=eval(substitute(function(mu, extra=NULL){
        ms <- sqrt(mu[,2])
        ns <- sqrt(mu[,6])
        nS <- c(-nS + sqrt(nS^2 + mu[,5]))
        mS <- 1-ns-ms-nS
        cbind(theta2eta(mS, .link, earg= .earg),
              theta2eta(ms, .link, earg= .earg),
              theta2eta(nS, .link, earg= .earg))
    }, list(.link=link,
             .earg=earg ))),
    loglikelihood=function(mu,y,w,residuals=FALSE,eta,extra=NULL)
        if(residuals) stop("loglikelihood residuals not implemented yet") else
        sum((w*y)*log(mu)),
    vfamily=c("MNSs", "vgenetic"),
    deriv=eval(substitute(expression({
        mS <- eta2theta(eta[,1], link=.link, earg= .earg)
        ms <- eta2theta(eta[,2], link=.link, earg= .earg)
        nS <- eta2theta(eta[,3], link=.link, earg= .earg)
        ns <- 1-mS-ms-nS
        dP1 <- cbind(2*(mS+ms), 0, 2*(nS+ns-mS), -2*ms, -2*nS, -2*ns)
        dP2 <- cbind(2*mS, 2*ms, 2*(nS-mS), 2*(ns-ms), -2*nS, -2*ns)
        dP3 <- cbind(0, 0, 2*ms, -2*ms,  2*ns, -2*ns) # n x 6
        dl1 <- apply(y * dP1 / mu, 1, sum)
        dl2 <- apply(y * dP2 / mu, 1, sum)
        dl3 <- apply(y * dP3 / mu, 1, sum)
        dPP.deta <- dtheta.deta(cbind(mS,ms,nS), link=.link, earg= .earg)
        w * dPP.deta * cbind(dl1, dl2, dl3)
    }), list(.link=link,
             .earg=earg ))),
    weight=eval(substitute(expression({
        dPP <- array(c(dP1,dP2,dP3), c(n,6,3))
        wz <- matrix(as.numeric(NA), n, dimm(M))   # dimm(M)==6 because M==3
        for(i1 in 1:M)
            for(i2 in i1:M) {
                index <- iam(i1,i2,M)
                wz[,index] <- apply(dPP[,,i1,drop=TRUE] * dPP[,,i2,drop=TRUE] /
                                mu, 1, sum) * dPP.deta[,i1] * dPP.deta[,i2]
        }
        w * wz
    }), list(.link=link,
             .earg=earg ))))
}


ABO <- function(link="logit", earg = list(), ir=NULL, ip=NULL)
{
    if(mode(link) != "character" && mode(link) != "name")
        link <- as.character(substitute(link))
    if(!is.list(earg)) earg = list()

    new("vglmff",
    blurb=c("ABO Blood Group System (A-B-AB-O phenotype)\n\n",
            "Links:    ",
            namesof("p",link, earg= earg), ", ", 
            namesof("q", link, earg= earg, tag=FALSE),
           "\n",
           "Variance: multinomial type variance"),
    initialize=eval(substitute(expression({
        delete.zero.colns <- FALSE
        eval(process.categorical.data.vgam)
        predictors.names <-
          c(namesof("pA", .link, earg= .earg,tag=FALSE),
            namesof("pB", .link, earg= .earg,tag=FALSE))
        if(is.null(etastart)) {
            r <- if(is.numeric(.ir)) rep(.ir, n) else
                       c(sqrt(mustart[,4]))
            p <- if(is.numeric(.ip)) rep(.ip, n) else
                c(1-sqrt(mustart[,2]+mustart[,4]))
            q <- 1-p-r
            etastart <- cbind(theta2eta(p, .link, earg= .earg),
                              theta2eta(q, .link, earg= .earg))
        }
    }), list(.link=link, .ir=ir, .ip=ip,
             .earg=earg ))),
    inverse=eval(substitute(function(eta, extra=NULL){
        p <- eta2theta(eta[,1], link=.link, earg= .earg)
        q <- eta2theta(eta[,2], link=.link, earg= .earg)
        r <- 1-p-q
        cbind(A=p*(p+2*r), B=q*(q+2*r), AB=2*p*q, O=r*r) 
    }, list(.link=link,
             .earg=earg ))),
    last=eval(substitute(expression({
        misc$link <- c(p = .link, q = .link)
        misc$earg = list(p= .earg, q= .earg )
    }), list(.link=link,
             .earg=earg ))),
    link=eval(substitute(function(mu, extra=NULL){
         r <- sqrt(mu[,4])
         p1 <- ( (1-r)+sqrt((1-r)^2 + 2*mu[,3]) )/2
         p2 <- ( (1-r)-sqrt((1-r)^2 + 2*mu[,3]) )/2
         index <- p2 >= 0 & p2 <= 1
         p <- p1
         p[index] <- p2[index]
         q <- 1-p-r
         cbind(theta2eta(p, .link, earg= .earg),
               theta2eta(q, .link, earg= .earg))
    }, list(.link=link,
             .earg=earg ))),
 
    loglikelihood=function(mu,y,w,residuals=FALSE,eta,extra=NULL)
        if(residuals) stop("loglikelihood residuals not implemented yet") else
        sum((w*y)*log(mu)),
    vfamily=c("ABO", "vgenetic"),
    deriv=eval(substitute(expression({
        p <- eta2theta(eta[,1], link=.link, earg= .earg)
        q <- eta2theta(eta[,2], link=.link, earg= .earg)
        r <- 1-p-q
        pbar <- 2*r+p
        qbar <- 2*r+q
        na <- y[,1]
        nb <- y[,2]
        nab <- y[,3]
        no <- y[,4]
        dl.dp <- (na+nab)/p - na/pbar - 2*nb/qbar - 2*no/r
        dl.dq <- (nb+nab)/q - 2*na/pbar - nb/qbar - 2*no/r
        dp.deta <- dtheta.deta(p, link=.link, earg= .earg)
        dq.deta <- dtheta.deta(q, link=.link, earg= .earg)
        w * cbind(dl.dp * dp.deta, dl.dq * dq.deta)
    }), list(.link=link,
             .earg=earg ))),
    weight=eval(substitute(expression({
        wz <- matrix(as.numeric(NA), n, dimm(M))   # dimm(M)==3 because M==2
        ed2l.dp2  <-  w * (1 + 2/p + 4*q/qbar + p/pbar)
        ed2l.dq2  <-  w * (1 + 2/q + 4*p/pbar + q/qbar)
        ed2l.dpdq <-  2 * w * (1 + q/qbar + p/pbar)
        wz[,iam(1,1,M)] <- dp.deta^2 * ed2l.dp2
        wz[,iam(2,2,M)] <- dq.deta^2 * ed2l.dq2
        wz[,iam(1,2,M)] <- ed2l.dpdq * dp.deta * dq.deta
        wz
    }), list(.link=link,
             .earg=earg ))))
}




AB.Ab.aB.ab <- function(link="logit", earg = list(), init.p=NULL)
{
    if(mode(link) != "character" && mode(link) != "name")
        link <- as.character(substitute(link))
    if(!is.list(earg)) earg = list()

    new("vglmff",
    blurb=c("AB-Ab-aB-ab phenotype\n\n",
            "Links:    ", namesof("p", link, earg= earg, tag=TRUE), "\n",
           "Variance: multinomial type variance"),

    initialize=eval(substitute(expression({
        delete.zero.colns <- FALSE
        eval(process.categorical.data.vgam)
        predictors.names <- namesof("p", .link, earg= .earg, tag=FALSE)
        if(is.null(etastart)) {
            p <- if(is.numeric(.init.p)) rep(.init.p,n) else
                      c(sqrt(4*mustart[,4]))
            etastart <- cbind(theta2eta(p, .link, earg= .earg))
        }
    }), list(.link=link, .init.p=init.p,
             .earg=earg ))),
    inverse=eval(substitute(function(eta,extra=NULL){
        p <- eta2theta(eta, link=.link, earg= .earg)
        pp4 <- p*p/4
        cbind(AB=0.5+pp4, Ab=0.25-pp4, aB=0.25-pp4, ab=pp4) 
    }, list(.link=link,
             .earg=earg ))),
    last=eval(substitute(expression({
        misc$link <- c(p = .link)
        misc$earg = list(p= .earg )
    }), list(.link=link,
             .earg=earg ))),
    link=eval(substitute(function(mu, extra=NULL){
        p <- sqrt(4* mu[,4])
        theta2eta(p, .link, earg= .earg)
    }, list(.link=link,
             .earg=earg ))),
 
    loglikelihood=function(mu,y,w,residuals=FALSE,eta,extra=NULL)
        if(residuals) stop("loglikelihood residuals not implemented yet") else
        sum((w*y)*log(mu)),
    vfamily=c("AB.Ab.aB.ab", "vgenetic"),
    deriv=eval(substitute(expression({
        pp <- eta2theta(eta, link=.link, earg= .earg)
        p2 <- pp*pp
        nAB <- w*y[,1]
        nAb <- w*y[,2]
        naB <- w*y[,3]
        nab <- w*y[,4]
        dl.dp <- 8 * pp * (nAB/(2+p2) - (nAb+naB)/(1-p2) + nab/p2)
        dp.deta <- dtheta.deta(pp, link=.link, earg= .earg)
        dl.dp * dp.deta
    }), list(.link=link,
             .earg=earg ))),
    weight=eval(substitute(expression({
        ed2l.dp2 <- 4 * w * p2 * (1/(2+p2) + 2/(1-p2) + 1/p2)
        wz <- cbind((dp.deta^2) * ed2l.dp2)
        wz
    }), list(.link=link,
             .earg=earg ))))
}



AA.Aa.aa <- function(link="logit", earg = list(), init.pA=NULL)
{
    if(mode(link) != "character" && mode(link) != "name")
        link <- as.character(substitute(link))
    if(!is.list(earg)) earg = list()

    new("vglmff",
    blurb=c("AA-Aa-aa phenotype\n\n",
            "Links:    ", namesof("pA",link, earg= earg), "\n",
            "Variance: multinomial type variance"),
    initialize=eval(substitute(expression({
        delete.zero.colns <- FALSE
        eval(process.categorical.data.vgam)
        predictors.names <- namesof("pA", .link, earg= .earg, tag=FALSE)
        if(is.null(etastart)) {
            pA <- if(is.numeric(.init.pA)) rep(.init.pA, n) else
                      c(sqrt(mustart[,1]))
            etastart <- cbind(theta2eta(pA, .link, earg= .earg))
        }
    }), list(.link=link, .init.pA=init.pA,
             .earg=earg ))),
    inverse=eval(substitute(function(eta,extra=NULL){
        pA <- eta2theta(eta, link=.link, earg= .earg)
        pp <- pA*pA
        cbind(AA=pp, Aa=2*pA*(1-pA), aa=(1-pA)^2) 
    }, list(.link=link,
             .earg=earg ))),
    last=eval(substitute(expression({
        misc$link <- c("pA" = .link)
        misc$earg = list("pA" = .earg )
    }), list(.link=link,
             .earg=earg ))),
    link=eval(substitute(function(mu, extra=NULL){
        pA <- sqrt(mu[,1])
        theta2eta(pA, .link, earg= .earg)
    }, list(.link=link,
             .earg=earg ))),
    loglikelihood=function(mu,y,w,residuals=FALSE,eta,extra=NULL)
        if(residuals) stop("loglikelihood residuals not implemented yet") else
        sum((w*y)*log(mu)),
    vfamily=c("AA.Aa.aa", "vgenetic"),
    deriv=eval(substitute(expression({
        pA <- eta2theta(eta, link=.link, earg= .earg)
        nAA <- w*y[,1]
        nAa <- w*y[,2]
        naa <- w*y[,3]
        dl.dpA <- (2*nAA+nAa)/pA - (nAa+2*naa)/(1-pA)
        dpA.deta <- dtheta.deta(pA, link=.link, earg= .earg)
        dl.dpA * dpA.deta
    }), list(.link=link,
             .earg=earg ))),
    weight=eval(substitute(expression({
        d2l.dp2 <- (2*nAA+nAa)/pA^2 + (nAa+2*naa)/(1-pA)^2
        wz <- cbind((dpA.deta^2) * d2l.dp2)
        wz
    }), list(.link=link,
             .earg=earg ))))
}




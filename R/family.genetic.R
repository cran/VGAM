# These functions are
# Copyright (C) 1998-2013 T.W. Yee, University of Auckland.
# All rights reserved.












 G1G2G3 <- function(link = "logit",
                    ip1 = NULL, ip2 = NULL, iF = NULL) {

  link <- as.list(substitute(link))
  earg <- link2list(link)
  link <- attr(earg, "function.name")



  new("vglmff",
  blurb = c("G1-G2-G3 phenotype\n\n",
            "Links:    ",
            namesof("p1", link, earg = earg), ", ", 
            namesof("p2", link, earg = earg), ", ", 
            namesof("f",  link, earg = earg, tag = FALSE)),
  deviance = Deviance.categorical.data.vgam,
  initialize = eval(substitute(expression({
    mustart.orig <- mustart

    delete.zero.colns <- FALSE
    eval(process.categorical.data.vgam)

    if (length(mustart.orig))
      mustart <- mustart.orig

    ok.col.ny <- c("G1G1","G1G2","G1G3","G2G2","G2G3","G3G3")
    if (length(col.ny <- colnames(y)) == length(ok.col.ny) &&
       setequal(ok.col.ny, col.ny)) {
        if (!all(ok.col.ny == col.ny))
            stop("the columns of the response matrix should have ",
                 "names (output of colnames()) ordered as ",
                 "c('G1G1','G1G2','G1G3','G2G2','G2G3','G3G3')")
    }

    predictors.names <-
     c(namesof("p1", .link , earg = .earg , tag = FALSE),
       namesof("p2", .link , earg = .earg , tag = FALSE),
       namesof("f",  .link , earg = .earg , tag = FALSE))

    if (is.null(etastart)) {



      mydeterminant <- mustart[, 2] * mustart[, 3] +
                       mustart[, 2] * mustart[, 5] +
                       mustart[, 3] * mustart[, 5]
      p1 <- if (is.numeric( .ip1 )) rep( .ip1 , len = n) else
            mustart[, 2] * mustart[, 3] / mydeterminant
      p2 <- if (is.numeric( .ip2 )) rep( .ip2 , len = n) else
            mustart[, 2] * mustart[, 5] / mydeterminant
      ff <- if (is.numeric( .iF  )) rep( .iF  , len = n) else
            abs(1 - mustart[, 2] / (2 * p1 * p2))

      if (any(p1 <= 0) || any(p1 >= 1))
        stop("bad initial value for 'p1'")
      if (any(p2 <= 0) || any(p2 >= 1))
        stop("bad initial value for 'p2'")

      etastart <-
        cbind(theta2eta(p1, .link , earg = .earg ),
              theta2eta(p2, .link , earg = .earg ),
              theta2eta(ff, .link , earg = .earg ))
      mustart <- NULL  # Since etastart has been computed.

    }
  }), list( .link = link, .ip1 = ip1, .ip2 = ip2, .iF = iF,
            .earg = earg))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    p1 <- eta2theta(eta[, 1], link = .link , earg = .earg )
    p2 <- eta2theta(eta[, 2], link = .link , earg = .earg )
    f  <- eta2theta(eta[, 3], link = .link , earg = .earg )
    p3 <- abs(1 - p1 - p2)
      cbind("G1G1" = f*p1+(1-f)*p1^2,
            "G1G2" = 2*p1*p2*(1-f),
            "G1G3" = 2*p1*p3*(1-f),
            "G2G2" = f*p2+(1-f)*p2^2,
            "G2G3" = 2*p2*p3*(1-f),
            "G3G3" = f*p3+(1-f)*p3^2)
  }, list( .link = link, .earg = earg))),

  last = eval(substitute(expression({
    misc$link <-    c(p1 = .link , p2 = .link , f = .link )

    misc$earg <- list(p1 = .earg , p2 = .earg , f = .earg )

    misc$expected <- TRUE
  }), list( .link = link, .earg = earg))),

  loglikelihood = function(mu, y, w, residuals = FALSE, eta, extra = NULL)
    if (residuals)
      stop("loglikelihood residuals not implemented yet") else {
          sum(dmultinomial(x = w * y, size = w, prob = mu,
                           log = TRUE, dochecking = FALSE))
      },
  vfamily = c("G1G2G3", "vgenetic"),
  deriv = eval(substitute(expression({
    p1 <- eta2theta(eta[, 1], link = .link , earg = .earg )
    p2 <- eta2theta(eta[, 2], link = .link , earg = .earg )
    p3 <- 1-p1-p2
    f  <- eta2theta(eta[, 3], link = .link , earg = .earg )
    dP1 <- cbind(f + 2*p1*(1-f), 2*(1-f)*p2, 2*(1-f)*(1-p2-2*p1),
                0, -2*(1-f)*p2, -f - 2*p3*(1-f))
    dP2 <- cbind(0, 2*p1*(1-f), -2*(1-f)*p1, f+2*p2*(1-f),
                 2*(1-f)*(1-p1-2*p2), -f - 2*p3*(1-f))
    dP3 <- cbind(p1*(1-p1), -2*p1*p2, -2*p1*p3, p2*(1-p2), -2*p2*p3, 
                 p3*(1-p3))
    dl1 <- rowSums(y * dP1 / mu)
    dl2 <- rowSums(y * dP2 / mu)
    dl3 <- rowSums(y * dP3 / mu)
    dPP.deta <- dtheta.deta(cbind(p1, p2, f),
                            link = .link , earg = .earg )
    c(w) * cbind(dPP.deta[, 1] * dl1,
                 dPP.deta[, 2] * dl2, 
                 dPP.deta[, 3] * dl3)
  }), list( .link = link, .earg = earg))),
  weight = eval(substitute(expression({
    dPP <- array(c(dP1, dP2, dP3), c(n, 6, 3))

    wz <- matrix(as.numeric(NA), n, dimm(M)) # dimm(M)==6 because M==3
    for(i1 in 1:M)
      for(i2 in i1:M) {
        index <- iam(i1,i2, M)
        wz[,index] <- rowSums(dPP[, , i1, drop = TRUE] *
                              dPP[, , i2, drop = TRUE] / mu) *
                              dPP.deta[, i1] * dPP.deta[, i2]
    }
    c(w) * wz
  }), list( .link = link, .earg = earg))))
}



 AAaa.nohw <- function(link = "logit", ipA = NULL, iF = NULL) {

  link <- as.list(substitute(link))
  earg <- link2list(link)
  link <- attr(earg, "function.name")


  new("vglmff",
  blurb = c("AA-Aa-aa phenotype (without Hardy-Weinberg assumption)\n\n",
            "Links:    ",
            namesof("pA", link, earg = earg), ", ", 
            namesof("f",  "identity", tag = FALSE)),
  deviance = Deviance.categorical.data.vgam,
  initialize = eval(substitute(expression({
    mustart.orig <- mustart

    delete.zero.colns <- FALSE
    eval(process.categorical.data.vgam)

    if (length(mustart.orig))
      mustart <- mustart.orig

    ok.col.ny <- c("AA","Aa","aa")
    if (length(col.ny <- colnames(y)) == length(ok.col.ny) &&
        setequal(ok.col.ny, col.ny)) {
        if (!all(ok.col.ny == col.ny))
          stop("the columns of the response matrix should have names ",
               "(output of colnames()) ordered as c('AA','Aa','aa')")
    }

    predictors.names <-
        c(namesof("pA", .link , earg = .earg , tag = FALSE),
          namesof("f",  "identity", earg = list(), tag = FALSE))

    if (is.null(etastart)) {
      pA <- if (is.numeric( .ipA )) rep( .ipA , len = n) else
            c(sqrt(mustart[, 1] - mustart[, 2] / 2))
      f <- if (is.numeric( .iF )) rep( .iF , len = n) else
           rep(0.01, len = n) # 1- mustart[, 2]/(2*pA*(1-pA))
      if (any(pA <= 0) || any(pA >= 1))
        stop("bad initial value for 'pA'")
      etastart <- cbind(theta2eta(pA, .link , earg = .earg ),
                        theta2eta(f,  "identity"))
      mustart <- NULL  # Since etastart has been computed.
    }
  }), list( .link = link, .ipA = ipA, .iF = iF, .earg = earg))),

  linkinv = eval(substitute(function(eta, extra = NULL) {
    pA <- eta2theta(eta[, 1], link = .link , earg = .earg )
    f  <- eta2theta(eta[, 2], link = "identity", earg = list())
    cbind(AA = pA^2+pA*(1-pA)*f,
          Aa = 2*pA*(1-pA)*(1-f),
          aa = (1-pA)^2 + pA*(1-pA)*f)
  }, list( .link = link, .earg = earg))),

  last = eval(substitute(expression({
    misc$link <-    c(pA = .link , f = "identity")

    misc$earg <- list(pA = .earg , f = list() )

    misc$expected <- TRUE
  }), list( .link = link, .earg = earg))),


  loglikelihood = function(mu, y, w, residuals = FALSE, eta, extra = NULL)
    if (residuals)
      stop("loglikelihood residuals not implemented yet") else {
      sum(dmultinomial(x = w * y, size = w, prob = mu,
                       log = TRUE, dochecking = FALSE))
    },
  vfamily = c("AAaa.nohw", "vgenetic"),
  deriv = eval(substitute(expression({
    pA <- eta2theta(eta[, 1], link = .link , earg = .earg )
    f  <- eta2theta(eta[, 2], link = "identity")
    dP1 <- cbind(f + 2*pA*(1-f),
                 2*(1-f)*(1-2*pA),
                 -2*(1-pA) +f*(1-2*pA))
    dP2 <- cbind(pA*(1-pA),
                 -2*pA*(1-pA),
                 pA*(1-pA))
    dl1 <- rowSums(y * dP1 / mu)
    dl2 <- rowSums(y * dP2 / mu)

    dPP.deta <- dtheta.deta(pA, link = .link , earg = .earg )

    c(w) * cbind(dPP.deta * dl1,
                 dl2)
  }), list( .link = link, .earg = earg))),
  weight = eval(substitute(expression({
    dPP <- array(c(dP1, dP2), c(n, 3, 2))
    dPP.deta <- cbind(dtheta.deta(pA, link = .link , earg = .earg ),
                      dtheta.deta(f,  link = "identity"))
    wz <- matrix(as.numeric(NA), n, dimm(M)) # dimm(M)==3 because M==2
    for(i1 in 1:M)
      for(i2 in i1:M) {
        index <- iam(i1, i2, M)
        wz[,index] <- rowSums(dPP[,,i1,drop = TRUE] *
                              dPP[,,i2,drop = TRUE] / mu) *
                              dPP.deta[, i1] * dPP.deta[, i2]
      }
    c(w) * wz
  }), list( .link = link, .earg = earg))))
}




 AB.Ab.aB.ab2 <- function(link = "logit", init.p = NULL) {

  link <- as.list(substitute(link))
  earg <- link2list(link)
  link <- attr(earg, "function.name")


  new("vglmff",
  blurb = c("AB-Ab-aB-ab2 phenotype\n\n",
            "Links:    ",
            namesof("p", link, earg = earg)),
  deviance = Deviance.categorical.data.vgam,
  initialize = eval(substitute(expression({
    mustart.orig <- mustart

    delete.zero.colns <- FALSE
    eval(process.categorical.data.vgam)
    predictors.names <- namesof("p", .link , earg = .earg , tag = FALSE)

    if (length(mustart.orig))
      mustart <- mustart.orig

        ok.col.ny <- c("AB","Ab","aB","ab")
        if (length(col.ny <- colnames(y)) == length(ok.col.ny) &&
           setequal(ok.col.ny, col.ny)) {
            if (!all(ok.col.ny == col.ny))
                stop("the columns of the response matrix should have names ",
                     "(output of colnames()) ordered as c('AB','Ab','aB','ab')")
        }

        if (is.null(etastart)) {
            p.init <- if (is.numeric(.init.p)) rep(.init.p, n) else
                     c(1 - 2 * sqrt(mustart[, 4]))
            etastart <- theta2eta(p.init, .link , earg = .earg )
            mustart <- NULL  # Since etastart has been computed.
        }
    }), list( .link = link, .init.p=init.p, .earg = earg))),
    linkinv = eval(substitute(function(eta,extra = NULL) {
        p <- eta2theta(eta, link = .link , earg = .earg )
        cbind("AB" = (2+(1-p)^2),
              "Ab" = (1-(1-p)^2),
              "aB" = (1-(1-p)^2),
              "ab" = (1-p)^2) / 4
    }, list( .link = link, .earg = earg) )),

  last = eval(substitute(expression({
    misc$link <-    c(p = .link )

    misc$earg <- list(p = .earg )

    misc$expected <- TRUE
  }), list( .link = link, .earg = earg) )),


  loglikelihood = function(mu, y, w, residuals = FALSE, eta, extra = NULL)
    if (residuals)
      stop("loglikelihood residuals not implemented yet") else {
      sum(dmultinomial(x = w * y, size = w, prob = mu,
                       log = TRUE, dochecking = FALSE))
    },
  vfamily = c("AB.Ab.aB.ab2", "vgenetic"),
  deriv = eval(substitute(expression({
    pp <- eta2theta(eta, link = .link , earg = .earg )
    dP1 <- cbind(-0.5*(1-pp),
                  0.5*(1-pp),
                  0.5*(1-pp),
                 -0.5*(1-pp))
    dl1 <- rowSums(y * dP1 / mu)
    dPP.deta <- dtheta.deta(pp, link = .link , earg = .earg )
    c(w) * dPP.deta * dl1
  }), list( .link = link, .earg = earg) )),
  weight = eval(substitute(expression({
    wz <- rowSums(dP1 * dP1 / mu) * dPP.deta^2
    c(w) * wz
  }), list( .link = link, .earg = earg) )))
}



 A1A2A3 <- function(link = "logit", ip1 = NULL, ip2 = NULL) {
  link <- as.list(substitute(link))
  earg <- link2list(link)
  link <- attr(earg, "function.name")


  new("vglmff",
  blurb = c("A1A2A3 Allele System ",
            "(A1A1, A1A2, A2A2, A1A3, A2A3, A3A3)\n\n",
            "Links:    ",
            namesof("p1", link, earg = earg), ", ", 
            namesof("p2", link, earg = earg, tag = FALSE)),
  deviance = Deviance.categorical.data.vgam,
  initialize = eval(substitute(expression({
    mustart.orig <- mustart

    delete.zero.colns <- FALSE
    eval(process.categorical.data.vgam)

    if (length(mustart.orig))
      mustart <- mustart.orig

        ok.col.ny <- c("A1A1","A1A2","A2A2","A1A3","A2A3","A3A3")
        if (length(col.ny <- colnames(y)) == length(ok.col.ny) &&
           setequal(ok.col.ny, col.ny)) {
            if (!all(ok.col.ny == col.ny))
                stop("the columns of the response matrix should have names ",
                     "(output of colnames()) ordered as ",
                     "c('A1A1','A1A2','A2A2','A1A3','A2A3','A3A3')")
        }

        predictors.names <-
            c(namesof("pA", .link , earg = .earg , tag = FALSE),
              namesof("pB", .link , earg = .earg , tag = FALSE))

        if (is.null(etastart)) {
            p1 <- if (is.numeric(.ip1)) rep(.ip1, n) else
                       c(sqrt(mustart[, 1]))
            p2 <- if (is.numeric(.ip2)) rep(.ip2, n) else
                       c(sqrt(mustart[, 3]))
            etastart <- cbind(theta2eta(p1, .link , earg = .earg ),
                              theta2eta(p2, .link , earg = .earg ))
            mustart <- NULL  # Since etastart has been computed.
    }
  }), list( .link = link, .ip1 = ip1, .ip2 = ip2, .earg = earg))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    p1 <- eta2theta(eta[, 1], link = .link , earg = .earg )
    p2 <- eta2theta(eta[, 2], link = .link , earg = .earg )
    qq <- abs(1 - p1 - p2)
    cbind(A1A1 = p1*p1,
          A1A2 = 2*p1*p2,
          A2A2 = p2*p2,
          A1A3 = 2*p1*qq,
          A2A3 = 2*p2*qq,
          A3A3 = qq*qq)
  }, list( .link = link, .earg = earg))),

  last = eval(substitute(expression({
    misc$link <-    c(p1 = .link , p2 = .link )

    misc$earg <- list(p1 = .earg , p2 = .earg )

    misc$expected <- TRUE
  }), list( .link = link, .earg = earg))),


  loglikelihood = function(mu, y, w, residuals = FALSE, eta, extra = NULL)
    if (residuals)
      stop("loglikelihood residuals not implemented yet") else {
      sum(dmultinomial(x = w * y, size = w, prob = mu,
                       log = TRUE, dochecking = FALSE))
    },
  vfamily = c("A1A2A3", "vgenetic"),
  deriv = eval(substitute(expression({
    p1 <- eta2theta(eta[, 1], link = .link , earg = .earg )
    p2 <- eta2theta(eta[, 2], link = .link , earg = .earg )

    dl.dp1 <- (2*y[, 1]+y[, 2]+y[, 4])/p1 - (2*y[,6]+y[, 4]+y[,5])/(1-p1-p2)
    dl.dp2 <- (2*y[, 3]+y[, 2]+y[,5])/p2 - (2*y[,6]+y[, 4]+y[,5])/(1-p1-p2)

    dp1.deta <- dtheta.deta(p1, link = .link , earg = .earg )
    dp2.deta <- dtheta.deta(p2, link = .link , earg = .earg )

    c(w) * cbind(dl.dp1 * dp1.deta,
              dl.dp2 * dp2.deta)
  }), list( .link = link, .earg = earg))),
  weight = eval(substitute(expression({
    qq <- 1-p1-p2
    wz <- matrix(as.numeric(NA), n, dimm(M)) # dimm(M)==3 because M==2
    ned2l.dp12  <-  2 * (1/p1 + 1/qq)
    ned2l.dp22  <-  2 * (1/p2 + 1/qq)
    ned2l.dp1dp2 <-  2 / qq
    wz[, iam(1, 1, M)] <- ned2l.dp12 * dp1.deta^2
    wz[, iam(2, 2, M)] <- ned2l.dp22 * dp2.deta^2
    wz[, iam(1, 2, M)] <- ned2l.dp1dp2 * dp1.deta * dp2.deta
    c(w) * wz
  }), list( .link = link, .earg = earg))))
}




 MNSs <- function(link = "logit",
                  imS = NULL, ims = NULL, inS = NULL) {

  link <- as.list(substitute(link))
  earg <- link2list(link)
  link <- attr(earg, "function.name")


  new("vglmff",
  blurb = c("MNSs Blood Group System (MS-Ms-MNS-MNs-NS-Ns phenotype)\n\n",
            "Links:    ",
            namesof("mS", link, earg = earg), ", ", 
            namesof("ms", link, earg = earg), ", ", 
            namesof("nS", link, earg = earg, tag = FALSE)),
  deviance = Deviance.categorical.data.vgam,
  initialize = eval(substitute(expression({
    mustart.orig <- mustart

    delete.zero.colns <- FALSE
    eval(process.categorical.data.vgam)

    if (length(mustart.orig))
      mustart <- mustart.orig

        ok.col.ny <- c("MS","Ms","MNS","MNs","NS","Ns")
        if (length(col.ny <- colnames(y)) == length(ok.col.ny) &&
            setequal(ok.col.ny, col.ny)) {
        if (!all(ok.col.ny == col.ny))
            stop("the columns of the response matrix should have ",
                 "names (output of colnames()) ordered as ",
                 "c('MS','Ms','MNS','MNs','NS','Ns')")
    }

    predictors.names <-
       c(namesof("mS", .link , earg = .earg , tag = FALSE),
         namesof("ms", .link , earg = .earg , tag = FALSE),
         namesof("nS", .link , earg = .earg , tag = FALSE))

    if (is.null(etastart)) {
      ms <- if (is.numeric(.ims)) rep(.ims, n) else
                 c(sqrt(mustart[, 2]))
      ns <- c(sqrt(mustart[,6]))
      nS <- if (is.numeric(.inS)) rep(.inS, n) else
          c(-ns + sqrt(ns^2 + mustart[,5]))  # Solve a quadratic eqn
      mS <- if (is.numeric(.imS)) rep(.imS, n) else
              1-ns-ms-nS
      etastart <- cbind(theta2eta(mS, .link , earg = .earg ),
                       theta2eta(ms, .link , earg = .earg ),
                       theta2eta(nS, .link , earg = .earg ))
      mustart <- NULL  # Since etastart has been computed.
    }
  }), list( .link = link, .imS = imS, .ims = ims, .inS = inS, .earg = earg))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    mS <- eta2theta(eta[, 1], link = .link , earg = .earg )
    ms <- eta2theta(eta[, 2], link = .link , earg = .earg )
    nS <- eta2theta(eta[, 3], link = .link , earg = .earg )
    ns <- abs(1 - mS - ms - nS)
    cbind(MS  = mS^2 + 2*mS*ms,
          Ms  = ms^2,
          MNS = 2*(mS*nS + ms*nS + mS*ns),
          MNs = 2*ms*ns,
          NS  = nS^2 + 2*nS*ns,
          Ns  = ns^2)
  }, list( .link = link, .earg = earg))),

  last = eval(substitute(expression({
    misc$link <-    c(mS = .link , ms = .link , nS = .link )

    misc$earg <- list(mS = .earg , ms = .earg , nS = .earg )

    misc$expected <- TRUE
  }), list( .link = link, .earg = earg))),


  loglikelihood = function(mu, y, w, residuals = FALSE, eta, extra = NULL)
    if (residuals)
      stop("loglikelihood residuals not implemented yet") else {
      sum(dmultinomial(x = w * y, size = w, prob = mu,
                       log = TRUE, dochecking = FALSE))
    },
  vfamily = c("MNSs", "vgenetic"),
  deriv = eval(substitute(expression({
    mS <- eta2theta(eta[, 1], link = .link , earg = .earg )
    ms <- eta2theta(eta[, 2], link = .link , earg = .earg )
    nS <- eta2theta(eta[, 3], link = .link , earg = .earg )
    ns <- 1-mS-ms-nS
    dP1 <- cbind(2*(mS+ms), 0, 2*(nS+ns-mS), -2*ms, -2*nS, -2*ns)
    dP2 <- cbind(2*mS, 2*ms, 2*(nS-mS), 2*(ns-ms), -2*nS, -2*ns)
    dP3 <- cbind(0, 0, 2*ms, -2*ms,  2*ns, -2*ns) # n x 6
    dl1 <- rowSums(y * dP1 / mu)
    dl2 <- rowSums(y * dP2 / mu)
    dl3 <- rowSums(y * dP3 / mu)
    dPP.deta <- dtheta.deta(cbind(mS, ms, nS), link = .link , earg = .earg )
    c(w) * dPP.deta * cbind(dl1, dl2, dl3)
  }), list( .link = link, .earg = earg))),
  weight = eval(substitute(expression({
    dPP <- array(c(dP1,dP2,dP3), c(n,6, 3))
    wz <- matrix(as.numeric(NA), n, dimm(M)) # dimm(M)==6 because M==3
    for(i1 in 1:M)
      for(i2 in i1:M) {
        index <- iam(i1,i2, M)
        wz[,index] <- rowSums(dPP[,,i1,drop = TRUE] *
                              dPP[,,i2,drop = TRUE] / mu) *
                              dPP.deta[,i1] * dPP.deta[,i2]
    }
    c(w) * wz
  }), list( .link = link, .earg = earg))))
}






 ABO <- function(link = "logit", ipA = NULL, ipO = NULL) {
  link <- as.list(substitute(link))
  earg <- link2list(link)
  link <- attr(earg, "function.name")


  new("vglmff",
  blurb = c("ABO Blood Group System (A-B-AB-O phenotype)\n\n",
            "Links:    ",
            namesof("pA", link, earg = earg), ", ", 
            namesof("pB", link, earg = earg, tag = FALSE)),
  deviance = Deviance.categorical.data.vgam,

  initialize = eval(substitute(expression({
    mustart.orig <- mustart

    delete.zero.colns <- FALSE
    eval(process.categorical.data.vgam)

    if (length(mustart.orig))
      mustart <- mustart.orig

    ok.col.ny <- c("A","B","AB","O")
    if (length(col.ny <- colnames(y)) == length(ok.col.ny) &&
        setequal(ok.col.ny, col.ny)) {
      if (!all(ok.col.ny == col.ny))
        stop("the columns of the response matrix should have names ",
             "(output of colnames()) ordered as c('A','B','AB','O')")
    }


    predictors.names <-
      c(namesof("pA", .link , earg = .earg , tag = FALSE),
        namesof("pB", .link , earg = .earg , tag = FALSE))

    if (!length(etastart)) {
      pO <- if (is.Numeric( .ipO )) rep( .ipO , len = n) else
           c(sqrt(mustart[, 4]))
      pA <- if (is.Numeric( .ipA )) rep( .ipA , len = n) else
          c(1 - sqrt(mustart[, 2] + mustart[, 4]))
      pB <- abs(1 - pA - pO)
      etastart <- cbind(theta2eta(pA, .link , earg = .earg ),
                       theta2eta(pB, .link , earg = .earg ))
      mustart <- NULL  # Since etastart has been computed.
    }
  }), list( .link = link, .ipO = ipO, .ipA = ipA, .earg = earg))),


  linkinv = eval(substitute(function(eta, extra = NULL) {
      pA <- eta2theta(eta[, 1], link = .link , earg = .earg )
      pB <- eta2theta(eta[, 2], link = .link , earg = .earg )
      pO <- abs(1 - pA - pB)
      cbind(A  = pA*(pA+2*pO),
            B  = pB*(pB+2*pO),
            AB = 2*pA*pB,
            O  = pO*pO) 
  }, list( .link = link, .earg = earg))),

  last = eval(substitute(expression({
    misc$link <-    c(pA = .link , pB = .link )

    misc$earg <- list(pA = .earg , pB = .earg )

    misc$expected <- TRUE
  }), list( .link = link, .earg = earg))),


  loglikelihood =
  function(mu, y, w, residuals = FALSE, eta, extra = NULL)
    if (residuals)
      stop("loglikelihood residuals not implemented yet") else {
        sum(dmultinomial(x = w * y, size = w, prob = mu, log = TRUE,
                         dochecking = FALSE))
    },

  vfamily = c("ABO", "vgenetic"),

  deriv = eval(substitute(expression({
    ppp <- eta2theta(eta[, 1], link = .link , earg = .earg )
    qqq <- eta2theta(eta[, 2], link = .link , earg = .earg )
    rrr <- abs(1 - ppp - qqq)


    pbar <- 2*rrr + ppp
    qbar <- 2*rrr + qqq
    naa <- y[, 1]
    nbb <- y[, 2]
    nab <- y[, 3]
    noo <- y[, 4]

    dl.dp <- (naa+nab)/ppp -   naa/pbar - 2*nbb/qbar - 2*noo/rrr
    dl.dq <- (nbb+nab)/qqq - 2*naa/pbar -   nbb/qbar - 2*noo/rrr
    dp.deta <- dtheta.deta(ppp, link = .link , earg = .earg )
    dq.deta <- dtheta.deta(qqq, link = .link , earg = .earg )

    c(w) * cbind(dl.dp * dp.deta,
                 dl.dq * dq.deta)
  }), list( .link = link, .earg = earg))),

  weight = eval(substitute(expression({
    wz <- matrix(as.numeric(NA), n, dimm(M)) # dimm(M)==3 because M==2

    ned2l.dp2  <- (1 + 2/ppp + 4*qqq/qbar + ppp/pbar)
    ned2l.dq2  <- (1 + 2/qqq + 4*ppp/pbar + qqq/qbar)
    ned2l.dpdq <- 2 * (1 + qqq/qbar + ppp/pbar)

    wz[, iam(1, 1, M)] <- ned2l.dp2 * dp.deta^2
    wz[, iam(2, 2, M)] <- ned2l.dq2 * dq.deta^2
    wz[, iam(1, 2, M)] <- ned2l.dpdq * dp.deta * dq.deta
    c(w) * wz
  }), list( .link = link, .earg = earg))))
}




 AB.Ab.aB.ab <- function(link = "logit", init.p = NULL) {
  link <- as.list(substitute(link))
  earg <- link2list(link)
  link <- attr(earg, "function.name")


  new("vglmff",
  blurb = c("AB-Ab-aB-ab phenotype\n\n",
            "Links:    ", namesof("p", link, earg = earg, tag = TRUE)),
  deviance = Deviance.categorical.data.vgam,
  initialize = eval(substitute(expression({
    mustart.orig <- mustart

    delete.zero.colns <- FALSE
    eval(process.categorical.data.vgam)

    if (length(mustart.orig))
      mustart <- mustart.orig

    ok.col.ny <- c("AB","Ab","aB","ab")
    if (length(col.ny <- colnames(y)) == length(ok.col.ny) &&
       setequal(ok.col.ny, col.ny)) {
        if (!all(ok.col.ny == col.ny))
          stop("the columns of the response matrix should have ",
               "names (output of colnames()) ordered as ",
               "c('AB','Ab','aB','ab')")
    }

    predictors.names <- namesof("p", .link , earg = .earg , tag = FALSE)

    if (is.null(etastart)) {
      p <- if (is.numeric( .init.p )) rep(.init.p, len = n) else
          c(sqrt(4 * mustart[, 4]))
      etastart <- cbind(theta2eta(p, .link , earg = .earg ))
      mustart <- NULL  # Since etastart has been computed.
    }
  }), list( .link = link, .init.p=init.p, .earg = earg))),
  linkinv = eval(substitute(function(eta,extra = NULL) {
    p <- eta2theta(eta, link = .link , earg = .earg )
    pp4 <- p * p / 4
    cbind(AB = 0.5 + pp4,
          Ab = 0.25 - pp4,
          aB = 0.25 - pp4,
          ab = pp4) 
  }, list( .link = link, .earg = earg))),

  last = eval(substitute(expression({
    misc$link <-    c(p = .link )

    misc$earg <- list(p = .earg )

    misc$expected <- TRUE
   }), list( .link = link, .earg = earg))),


  loglikelihood = function(mu, y, w, residuals = FALSE, eta, extra = NULL)
    if (residuals)
      stop("loglikelihood residuals not implemented yet") else {
      sum(dmultinomial(x = w * y, size = w, prob = mu,
                       log = TRUE, dochecking = FALSE))
    },
  vfamily = c("AB.Ab.aB.ab", "vgenetic"),
  deriv = eval(substitute(expression({
    pp <- eta2theta(eta, link = .link , earg = .earg )

    p2 <- pp*pp
    nAB <- w * y[, 1]
    nAb <- w * y[, 2]
    naB <- w * y[, 3]
    nab <- w * y[, 4]

    dl.dp <- 8 * pp * (nAB/(2+p2) - (nAb+naB)/(1-p2) + nab/p2)

    dp.deta <- dtheta.deta(pp, link = .link , earg = .earg )

    dl.dp * dp.deta
  }), list( .link = link, .earg = earg))),
  weight = eval(substitute(expression({
    ned2l.dp2 <- 4 * p2 * (1/(2+p2) + 2/(1-p2) + 1/p2)
    wz <- cbind((dp.deta^2) * ned2l.dp2)
    c(w) * wz
  }), list( .link = link, .earg = earg))))
}



 AA.Aa.aa <- function(link = "logit", init.pA = NULL) {
  link <- as.list(substitute(link))
  earg <- link2list(link)
  link <- attr(earg, "function.name")


  new("vglmff",
  blurb = c("AA-Aa-aa phenotype\n\n",
            "Links:    ", namesof("pA", link, earg = earg)),
  deviance = Deviance.categorical.data.vgam,
  initialize = eval(substitute(expression({
    mustart.orig <- mustart

    delete.zero.colns <- FALSE
    eval(process.categorical.data.vgam)

    if (length(mustart.orig))
      mustart <- mustart.orig

    ok.col.ny <- c("AA","Aa","aa")
    if (length(col.ny <- colnames(y)) == length(ok.col.ny) &&
       setequal(ok.col.ny, col.ny)) {
        if (!all(ok.col.ny == col.ny))
          stop("the columns of the response matrix ",
               "should have names ",
               "(output of colnames()) ordered as c('AA','Aa','aa')")
    }

    predictors.names <- namesof("pA", .link , earg = .earg , tag = FALSE)

    if (is.null(etastart)) {
      pA <- if (is.numeric(.init.pA)) rep(.init.pA, n) else
                c(sqrt(mustart[, 1]))
      etastart <- cbind(theta2eta(pA, .link , earg = .earg ))
      mustart <- NULL  # Since etastart has been computed.
    }
  }), list( .link = link, .init.pA=init.pA, .earg = earg))),
  linkinv = eval(substitute(function(eta,extra = NULL) {
    pA <- eta2theta(eta, link = .link , earg = .earg )
    pp <- pA*pA
    cbind(AA = pp,
          Aa = 2*pA*(1-pA),
          aa = (1-pA)^2) 
  }, list( .link = link, .earg = earg))),

  last = eval(substitute(expression({
    misc$link <-    c("pA" = .link )

    misc$earg <- list("pA" = .earg )

    misc$expected = TRUE
  }), list( .link = link, .earg = earg))),


  loglikelihood = function(mu, y, w, residuals = FALSE, eta, extra = NULL)
    if (residuals)
      stop("loglikelihood residuals not implemented yet") else {
        sum(dmultinomial(x = w * y, size = w, prob = mu,
                         log = TRUE, dochecking = FALSE))
    },
  vfamily = c("AA.Aa.aa", "vgenetic"),
  deriv = eval(substitute(expression({
    pA  <- eta2theta(eta, link = .link , earg = .earg )
    nAA <- w * y[, 1]
    nAa <- w * y[, 2]
    naa <- w * y[, 3]
    dl.dpA <- (2*nAA+nAa)/pA - (nAa+2*naa)/(1-pA)
    dpA.deta <- dtheta.deta(pA, link = .link , earg = .earg )
    dl.dpA * dpA.deta
  }), list( .link = link, .earg = earg))),
  weight = eval(substitute(expression({
    ned2l.dp2 <- (2*nAA+nAa)/pA^2 + (nAa+2*naa)/(1-pA)^2
    wz <- cbind((dpA.deta^2) * ned2l.dp2)
    wz
  }), list( .link = link, .earg = earg))))
}







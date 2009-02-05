# These functions are
# Copyright (C) 1998-2009 T.W. Yee, University of Auckland. All rights reserved.






vnonlinear.control <- function(save.weight=TRUE, ...)
{



    list(save.weight=as.logical(save.weight)[1])
}


micmen <- function(rpar=0.001, divisor=10,
                   init1=NULL, init2=NULL,
                   link1="identity",
                   link2="identity",
                   earg1=list(), 
                   earg2=list(), 
                   dispersion=0,
                   zero=NULL)
{



    estimated.dispersion <- dispersion==0

    if(mode(link1) != "character" && mode(link1) != "name")
        link1 <- as.character(substitute(link1))
    if(mode(link2) != "character" && mode(link2) != "name")
        link2 <- as.character(substitute(link2))
    if(!is.list(earg1)) earg1 = list()
    if(!is.list(earg2)) earg2 = list()

    new("vglmff",
    blurb=c("Michaelis-Menton regression model\n",
           "Y_i = theta1 * u_i / (theta2 + u_i) + e_i\n\n",
           "Links:    ",
           namesof("theta1", link1, earg=earg1), ", ",
           namesof("theta2", link2, earg=earg2),
           "\n",
           "Variance: constant"),
    constraints=eval(substitute(expression({
        constraints <- cm.zero.vgam(constraints, x, .zero, M=2)
    }), list(.zero=zero))),
    deviance=function(mu, y, w, residuals=FALSE, eta, extra=NULL) {
        M <- if(is.matrix(y)) ncol(y) else 1
        if(residuals) {
            if(M>1) NULL else (y-mu) * sqrt(w)
        } else
            rss.vgam(y-mu, w, M=M)
    },
    initialize=eval(substitute(expression({
        if(ncol(cbind(y)) != 1)
            stop("response must be a vector or a one-column matrix")

        if(!length(Xm2))
            stop("regressor not found")
        if(ncol(as.matrix(Xm2)) != 1)
            stop(paste("regressor not found or is not a vector. Use the",
                 "'form2' argument without an intercept"))
        Xm2 = as.vector(Xm2) # Make sure
        extra$Xm2 = Xm2          # Needed for @inverse

        predictors.names <-
          c(namesof("theta1", .link1, earg= .earg1, tag=FALSE),
            namesof("theta2", .link2, earg= .earg2, tag=FALSE))

        if(length(mustart) || length(coefstart))
            stop("can't handle mustart or coefstart")
        if(!length(etastart)) {
            index <- (1:n)[Xm2>quantile(Xm2, prob=.85)]
            init1 <- median(y[index])
            init2 <- median(init1*Xm2/y - Xm2)

            if(length(.init1)) init1 = .init1
            if(length(.init2)) init2 = .init2

            etastart = cbind(
                rep(theta2eta(init1, .link1, earg= .earg1), len=n),
                rep(theta2eta(init2, .link2, earg= .earg2), len=n))
        } else {
            stop("cannot handle etastart or mustart")
        }

    }), list(.init1=init1, .init2=init2,
             .earg1=earg1, .earg2=earg2,
             .link1=link1, .link2=link2))),
    inverse=eval(substitute(function(eta, extra=NULL) {
        theta1 <- eta2theta(eta[,1], .link1, earg= .earg1)
        theta2 <- eta2theta(eta[,2], .link2, earg= .earg2)
        theta1 * extra$Xm2 / (theta2 + extra$Xm2)
    }, list(.link1=link1, .link2=link2,
            .earg1=earg1, .earg2=earg2 ))),
    last=eval(substitute(expression({
        misc$link <- c(theta1= .link1, theta2= .link2)
        misc$earg = list(theta1= .earg1, theta2= .earg2 )
        misc$rpar <- rpar
        fit$df.residual <- n - rank   # Not n.big - rank
        fit$df.total <- n             # Not n.big

        extra$Xm2 = NULL             # Regressor is in control$regressor 
        dpar <- .dispersion
        if(!dpar) {
            dpar <- sum(w * (y-mu)^2) / (n - p.big)
        }
        misc$dispersion <- dpar
        misc$default.dispersion <- 0
        misc$estimated.dispersion <- .estimated.dispersion
    }), list(.link1=link1, .link2=link2, .dispersion=dispersion,
             .earg1=earg1, .earg2=earg2,
             .estimated.dispersion=estimated.dispersion))),
    summary.dispersion=FALSE,
    vfamily=c("micmen","vnonlinear"),
    deriv=eval(substitute(expression({
        if(iter>1) { 
            rpar = max(rpar / .divisor, 1000 * .Machine$double.eps)
        } else {
            rpar = .rpar
            d3 = deriv3(~ theta1 * Xm2 / (theta2 + Xm2),
                        c("theta1","theta2"), hessian=FALSE)
        }

        theta1 <- eta2theta(eta[,1], .link1, earg= .earg1)
        theta2 <- eta2theta(eta[,2], .link2, earg= .earg2)

        if(TRUE) {
            dmus.dthetas  = attr(eval(d3), "gradient")
        } else {
            dmu.dtheta1 <- Xm2 / (theta2 + Xm2)
            dmu.dtheta2 <- -theta1 * Xm2 / (Xm2 + theta2)^2
            dmus.dthetas  = cbind(dmu.dtheta1, dmu.dtheta2)
        }

        dthetas.detas = cbind(dtheta.deta(theta1, .link1, earg= .earg1),
                              dtheta.deta(theta2, .link2, earg= .earg2))

        if(TRUE) {
            index = iam(NA, NA, M=M, both=TRUE)
            temp200809 = dmus.dthetas * dthetas.detas
            if(M>1)
                temp200809[,2:M] = temp200809[,2:M] + sqrt(rpar)
            w * (y-mu) * temp200809
        } else {
            w * (y-mu) *
            cbind(dmus.dthetas[,1] * dthetas.detas[,1],
                  dmus.dthetas[,2] * dthetas.detas[,2] + sqrt(rpar))
        }
    }), list( .link1=link1, .link2=link2, .rpar=rpar,
              .earg1=earg1, .earg2=earg2,
              .divisor=divisor))),
    weight=eval(substitute(expression({
        if(TRUE) {
            wz = dmus.dthetas[,index$row] * dmus.dthetas[,index$col] *
                 dthetas.detas[,index$row] * dthetas.detas[,index$col]
            if(M>1)
                wz[,2:M] = wz[,2:M] + rpar
        } else {
            wz = cbind((dmus.dthetas[,1] * dthetas.detas[,1])^2,
                       (dmus.dthetas[,2] * dthetas.detas[,2])^2 + rpar,
                        dmus.dthetas[,1] * dmus.dthetas[,2] * 
                        dthetas.detas[,1] * dthetas.detas[,2])
        }
        w * wz
    }), list(.link1=link1, .link2=link2))))
}



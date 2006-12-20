# These functions are
# Copyright (C) 1998-2006 T.W. Yee, University of Auckland. All rights reserved.




loglinb2 <- function(exchangeable=FALSE, zero=NULL)
{

    new("vglmff",
    blurb=c("Log-linear model for binary data\n\n",
           "Links:    ",
           "Identity: u_1, u_2, u_{12}",
           "\n"),
    constraints=eval(substitute(expression({
        constraints <- cm.vgam(matrix(c(1,1,0, 0,0,1), 3, 2), x,
                               .exchangeable, constraints, intercept.apply=TRUE)
        constraints <- cm.zero.vgam(constraints, x, .zero, M)
    }), list(.exchangeable=exchangeable, .zero=zero))),
    initialize=expression({

        y <- as.matrix(y)
        predictors.names <- c("u1", "u2", "u12")
        if(ncol(y) != 2)
            stop("ncol(y) must be = 2")

        if(is.null(mustart)) {
            mustart <- matrix(as.numeric(NA), nrow(y), 4)
            mustart[,1] <- weighted.mean((1-y[,1])*(1-y[,2]), w)
            mustart[,2] <- weighted.mean((1-y[,1])*y[,2], w)
            mustart[,3] <- weighted.mean(y[,1]*(1-y[,2]), w)
            mustart[,4] <- weighted.mean(y[,1]*y[,2], w)
            if(any(mustart==0)) 
                stop("some combinations of the response not realized") 
        }
    }),
    inverse= function(eta, extra=NULL) {
        u1 <-  eta[,1]
        u2 <-  eta[,2]
        u12 <- eta[,3]
        denom <- 1 + exp(u1) + exp(u2) + exp(u1 + u2 + u12)
        cbind("00"=1/denom,
              "01"=exp(u2) / denom,
              "10"=exp(u1) / denom,
              "11"=exp(u1+u2+u12) / denom)
    },
    last=expression({
        misc$link = c("u1" = "identity", "u2" = "identity", "u12" = "identity")
        misc$earg = list(u1=list(), u2=list(), u12=list())
    }),
    link= function(mu, extra=NULL)  {
        u0 <-  log(mu[,1]) 
        u2 <-  log(mu[,2]) - u0
        u1 <-  log(mu[,3]) - u0
        u12 <- log(mu[,4]) - u0 - u1 - u2 
        cbind(u1, u2, u12)
    },
    loglikelihood=function(mu,y,w,residuals=FALSE,eta,extra=NULL) {
        u1 <-  eta[,1]
        u2 <-  eta[,2]
        u12 <- eta[,3]
        denom <- 1 + exp(u1) + exp(u2) + exp(u1 + u2 + u12)
        u0 <- -log(denom)
        if(residuals) stop("loglikelihood residuals not implemented yet") else
        sum(w*(u0 + u1*y[,1] + u2*y[,2] + u12*y[,1]*y[,2]))
    },
    vfamily=c("loglinb2"),
    deriv=expression({
        u1 <-  eta[,1]
        u2 <-  eta[,2]
        u12 <- eta[,3]
        denom <- 1 + exp(u1) + exp(u2) + exp(u1 + u2 + u12)
        du0.du1 <- -(exp(u1) + exp(u1 + u2 + u12)) / denom 
        du0.du2 <- -(exp(u2) + exp(u1 + u2 + u12)) / denom 
        du0.du12 <- -exp(u1 + u2 + u12) / denom 
        w * cbind(du0.du1 + y[,1], 
                  du0.du2 + y[,2],
                  du0.du12 + y[,1]*y[,2]) 
    }),
    weight=expression({
        d2u0.du1.2 <- -(exp(u1) + exp(u1 + u2 + u12)) * (1+exp(u2)) / denom^2 
        d2u0.du22 <-  -(exp(u2) + exp(u1 + u2 + u12)) * (1+exp(u1)) / denom^2 
        d2u0.du122 <- -exp(u1 + u2 + u12) * (1+exp(u1)+exp(u2)) / denom^2 
        d2u0.du1u2 <- -(exp(u1 + u2 + u12) - exp(u1 + u2)) / denom^2 
        d2u0.du1u3 <- -(1 + exp(u2)) * exp(u1 + u2 + u12) / denom^2 
        d2u0.du2u3 <- -(1 + exp(u1)) * exp(u1 + u2 + u12) / denom^2 

        wz <- matrix(as.numeric(NA), n, dimm(M)) 
        wz[,iam(1,1,M)] <- -d2u0.du1.2 
        wz[,iam(2,2,M)] <- -d2u0.du22
        wz[,iam(3,3,M)] <- -d2u0.du122 
        wz[,iam(1,2,M)] <- -d2u0.du1u2
        wz[,iam(1,3,M)] <- -d2u0.du1u3
        wz[,iam(2,3,M)] <- -d2u0.du2u3
        w * wz
    }))
}


loglinb3 <- function(exchangeable=FALSE, zero=NULL)
{

    new("vglmff",
    blurb=c("Log-linear model for trivariate binary data\n\n",
           "Links:    ",
           "Identity: u1, u2, u3, u12, u13, u23",
           "\n"),
    constraints=eval(substitute(expression({
        constraints <- cm.vgam(matrix(c(1,1,1,0,0,0, 0,0,0,1,1,1), 6, 2), x,
                               .exchangeable, constraints, intercept.apply=TRUE)
        constraints <- cm.zero.vgam(constraints, x, .zero, M)
    }), list(.exchangeable=exchangeable, .zero=zero))),
    initialize=expression({
        y <- as.matrix(y)
        predictors.names <- c("u1", "u2", "u3", "u12", "u13", "u23")
        if(ncol(y) != 3)
            stop("ncol(y) must be = 3")
        extra$my.expression <- expression({
            u1 <-  eta[,1]
            u2 <-  eta[,2]
            u3 <-  eta[,3]
            u12 <- eta[,4]
            u13 <- eta[,5]
            u23 <- eta[,6]
            denom <- 1 + exp(u1) + exp(u2) + exp(u3) + exp(u1 + u2 + u12) +
                     exp(u1 + u3 + u13) + exp(u2 + u3 + u23) +
                     exp(u1 + u2 + u3 + u12 + u13 + u23)
        })
        extra$deriv.expression <- expression({
            allterms <- exp(u1+u2+u3+u12+u13+u23)
            A1 <- exp(u1) + exp(u1 + u2 + u12) + exp(u1 + u3 + u13) + allterms
            A2 <- exp(u2) + exp(u1 + u2 + u12) + exp(u2 + u3 + u23) + allterms
            A3 <- exp(u3) + exp(u3 + u2 + u23) + exp(u1 + u3 + u13) + allterms
            A12 <- exp(u1 + u2 + u12) + allterms
            A13 <- exp(u1 + u3 + u13) + allterms
            A23 <- exp(u2 + u3 + u23) + allterms
        })
        if(!length(mustart)) {
            mustart <- matrix(as.numeric(NA), nrow(y), 2^3)
            mustart[,1] <- weighted.mean((1-y[,1])*(1-y[,2])*(1-y[,3]), w)
            mustart[,2] <- weighted.mean((1-y[,1])*(1-y[,2])*y[,3], w)
            mustart[,3] <- weighted.mean((1-y[,1])*y[,2]*(1-y[,3]), w)
            mustart[,4] <- weighted.mean((1-y[,1])*y[,2]*y[,3], w)
            mustart[,5] <- weighted.mean(y[,1]*(1-y[,2])*(1-y[,3]), w)
            mustart[,6] <- weighted.mean(y[,1]*(1-y[,2])*y[,3], w)
            mustart[,7] <- weighted.mean(y[,1]*y[,2]*(1-y[,3]), w)
            mustart[,8] <- weighted.mean(y[,1]*y[,2]*y[,3], w)
            if(any(mustart==0)) 
                stop("some combinations of the response not realized") 
        }
    }),
    inverse= function(eta, extra=NULL) {
        eval(extra$my.expression)
        cbind("000"=1,
              "001"=exp(u3),
              "010"=exp(u2),
              "011"=exp(u2+u3+u23),
              "100"=exp(u1),
              "101"=exp(u1+u3+u13),
              "110"=exp(u1+u2+u12),
              "111"=exp(u1+u2+u3+u12+u13+u23)) / denom
    },
    last=expression({
        misc$link = rep("identity", length=M)
        names(misc$link) = predictors.names
        misc$earg = list(u1=list(), u2=list(), u3=list(),
                         u12=list(), u13=list(), u23=list())
    }),
    link= function(mu, extra=NULL)  {
        u0 <-  log(mu[,1])
        u3 <-  log(mu[,2]) - u0
        u2 <-  log(mu[,3]) - u0
        u23 <- log(mu[,4]) - u0 - u2 - u3
        u1 <-  log(mu[,5]) - u0
        u13 <- log(mu[,6]) - u0 - u1 - u3
        u12 <- log(mu[,7]) - u0 - u1 - u2
        cbind(u1, u2, u3, u12, u13, u23)
    },
    loglikelihood=function(mu,y,w,residuals=FALSE,eta,extra=NULL) {
        eval(extra$my.expression)
        u0 <- -log(denom)
        if(residuals) stop("loglikelihood residuals not implemented yet") else
        sum(w*(u0 + u1*y[,1] + u2*y[,2] + u3*y[,3] +u12*y[,1]*y[,2] +
               u13*y[,1]*y[,3] + u23*y[,2]*y[,3]))
    },
    vfamily=c("loglinb3"),
    deriv=expression({
        eval(extra$my.expression)
        eval(extra$deriv.expression)
        w * cbind(-A1/denom + y[,1], 
                  -A2/denom + y[,2],
                  -A3/denom + y[,3],
                  -A12/denom + y[,1]*y[,2],
                  -A13/denom + y[,1]*y[,3],
                  -A23/denom + y[,2]*y[,3])
    }),
    weight=expression({
        u0 <- -log(denom)
        dA2.du1 <- exp(u1 + u2 + u12) + allterms
        dA3.du1 <- exp(u1 + u3 + u13) + allterms
        dA3.du2 <- exp(u2 + u3 + u23) + allterms
        wz <- matrix(as.numeric(NA), n, dimm(6)) 
        expu0 <- exp(u0)
        wz[,iam(1,1,M)] <- A1 * (1 - expu0 * A1)
        wz[,iam(2,2,M)] <- A2 * (1 - expu0 * A2)
        wz[,iam(3,3,M)] <- A3 * (1 - expu0 * A3)
        wz[,iam(1,2,M)] <- (dA2.du1 - expu0 * A1 * A2)
        wz[,iam(1,3,M)] <- (dA3.du1 - expu0 * A1 * A3)
        wz[,iam(2,3,M)] <- (dA3.du2 - expu0 * A2 * A3)
        wz[,iam(4,4,M)] <- A12 * (1 - expu0 * A12)
        wz[,iam(5,5,M)] <- A13 * (1 - expu0 * A13)
        wz[,iam(6,6,M)] <- A23 * (1 - expu0 * A23)
        wz[,iam(4,6,M)] <- (allterms - expu0 * A12 * A23)
        wz[,iam(5,6,M)] <- (allterms - expu0 * A12 * A23)
        wz[,iam(4,5,M)] <- (allterms - expu0 * A12 * A13)
        wz[,iam(1,4,M)] <- A12 * (1 - expu0 * A1)
        wz[,iam(1,5,M)] <- A13 * (1 - expu0 * A1)
        wz[,iam(1,6,M)] <- (allterms - expu0 * A1 * A23)
        wz[,iam(2,4,M)] <- A12 * (1 - expu0 * A2)
        wz[,iam(2,5,M)] <- (allterms - expu0 * A2 * A13)
        wz[,iam(2,6,M)] <- A23 * (1 - expu0 * A2)
        wz[,iam(3,4,M)] <- (allterms - expu0 * A3 * A12)
        wz[,iam(3,5,M)] <- A13 * (1 - expu0 * A3)
        wz[,iam(3,6,M)] <- A23 * (1 - expu0 * A3)
        wz <- w * expu0 * wz 
        wz
    }))
}


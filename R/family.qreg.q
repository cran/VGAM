# These functions are
# Copyright (C) 1998-2006 T.W. Yee, University of Auckland. All rights reserved.











lms.bcn.control <-
lms.bcg.control <-
lms.yjn.control <- function(trace=TRUE, ...)
   list(trace=trace) 





lms.bcn <- function(percentiles=c(25,50,75),
                    zero=NULL,
                    link.sigma="loge",
                    link.mu="identity",
                    dfmu.init=4,
                    dfsigma.init=2,
                    init.lambda=1,
                    init.sigma=NULL)
{
    if(mode(link.sigma) != "character" && mode(link.sigma) != "name")
        link.sigma <- as.character(substitute(link.sigma))
    if(mode(link.mu) != "character" && mode(link.mu) != "name")
        link.mu <- as.character(substitute(link.mu))

    new("vglmff",
    blurb=c("LMS Quantile Regression (Box-Cox transformation to normality)\n",
            "Links:    ",
            "lambda",
            ", ",
            namesof("mu", link=link.mu),
            ", ",
            namesof("sigma", link=link.sigma)),
    constraints=eval(substitute(expression({
        constraints <- cm.zero.vgam(constraints, x, .zero, M)
    }), list(.zero=zero))),
    initialize=eval(substitute(expression({
        if(any(y<0, na.rm = TRUE))
            stop("negative responses not allowed")

        predictors.names <- c("lambda",
                     namesof("mu",  .link.mu,  short= TRUE),
                     namesof("sigma",  .link.sigma,  short= TRUE))
 
        if(!length(etastart)) {

            fit500=vsmooth.spline(x=x[,min(ncol(x),2)],y=y,w=w, df= .dfmu.init)
            fv.init = c(predict(fit500, x=x[,min(ncol(x),2)])$y)

            lambda.init = if(is.Numeric( .init.lambda)) .init.lambda else 1.0
            sigma.init <- if(is.null(.init.sigma)) {
                myratio <- ((y/fv.init)^lambda.init - 1) / lambda.init
                if(is.Numeric( .dfsigma.init)) {
                    fit600 = vsmooth.spline(x=x[,min(ncol(x),2)], y=myratio^2,
                                            w=w, df= .dfsigma.init)
                    sqrt(c(abs(predict(fit600, x=x[,min(ncol(x),2)])$y)))
                } else 
                    sqrt(var(myratio))
            } else .init.sigma
 
            etastart <- cbind(lambda.init,
                              theta2eta(fv.init, .link.mu),
                              theta2eta(sigma.init,  .link.sigma))
        }
    }), list(.link.sigma=link.sigma,
             .link.mu=link.mu,
             .dfmu.init=dfmu.init,
             .dfsigma.init=dfsigma.init,
             .init.lambda=init.lambda,
             .init.sigma=init.sigma))),
    inverse=eval(substitute(function(eta, extra=NULL) {
        eta[,2] <- eta2theta(eta[,2], .link.mu)
        eta[,3] <- eta2theta(eta[,3], .link.sigma)
        qtplot.lms.bcn(percentiles= .percentiles, eta=eta)
    }, list(.percentiles=percentiles,
            .link.mu=link.mu,
            .link.sigma=link.sigma))),
    last=eval(substitute(expression({
        misc$percentiles <- .percentiles
        misc$links <- c(lambda = "identity", mu = .link.mu, sigma = .link.sigma)
        misc$true.mu <- FALSE    # $fitted is not a true mu
        if(control$cdf) {
            post$cdf = cdf.lms.bcn(y, eta0=matrix(c(lambda,mymu,sigma), 
                ncol=3, dimnames=list(dimnames(x)[[1]], NULL)))
        }
    }), list(.percentiles=percentiles,
            .link.mu=link.mu,
            .link.sigma=link.sigma))),
    loglikelihood=eval(substitute(
        function(mu,y,w, residuals= FALSE, eta, extra=NULL) {
            lambda <- eta[,1]
            mu <- eta2theta(eta[,2], .link.mu)
            sigma <- eta2theta(eta[,3], .link.sigma)
            z <- ((y/mu)^lambda - 1) / (lambda * sigma)
         if(residuals) stop("loglikelihood residuals not implemented yet") else
            sum(w * (lambda * log(y/mu) - log(sigma) - 0.5*z^2))
        }, list(.link.sigma=link.sigma, .link.mu=link.mu))),
    vfamily=c("lms.bcn", "lmscreg"),
    deriv=eval(substitute(expression({
        lambda <- eta[,1]
        mymu <- eta2theta(eta[,2], .link.mu)
        sigma <- eta2theta(eta[,3], .link.sigma)
        z <- ((y/mymu)^lambda - 1) / (lambda * sigma)
        z2m1 <- z * z - 1
        d1 <- z*(z - log(y/mymu) / sigma) / lambda - z2m1 * log(y/mymu)
        d2 <- z / (mymu * sigma) + z2m1 * lambda / mymu
        d2 <- d2 * dtheta.deta(mymu, .link.mu)
        d3 <- z2m1 / sigma
        d3 <- d3 * dtheta.deta(sigma, .link.sigma)
        w * cbind(d1, d2, d3)
    }), list(.link.sigma=link.sigma, .link.mu=link.mu))),
    weight=eval(substitute(expression({
        wz <- matrix(as.numeric(NA), n, 6)
        wz[,iam(1,1,M)] <- (7 * sigma^2 / 4)
        wz[,iam(2,2,M)] <- (1 + 2*(lambda * sigma)^2) / (mymu*sigma)^2 *
                           dtheta.deta(mymu, .link.mu)^2
        wz[,iam(3,3,M)] <- (2 / sigma^2) *
                           dtheta.deta(sigma, .link.sigma)^2
        wz[,iam(1,2,M)] <- (-1 / (2 * mymu)) *
                           dtheta.deta(mymu, .link.mu)
        wz[,iam(1,3,M)] <- (lambda * sigma) *
                           dtheta.deta(sigma, .link.sigma)
        wz[,iam(2,3,M)] <- (2 * lambda / (mymu * sigma)) *
                           dtheta.deta(sigma, .link.sigma) *
                           dtheta.deta(mymu, .link.mu)
        wz * w
    }), list(.link.sigma=link.sigma, .link.mu=link.mu))))
}



lms.bcg <- function(percentiles=c(25,50,75),
                          zero=NULL,
                          link.sigma="loge",
                          link.mu="identity",
                          dfmu.init=4,
                          dfsigma.init=2,
                          init.lambda=1,
                          init.sigma=NULL)
{
    if(mode(link.sigma) != "character" && mode(link.sigma) != "name")
        link.sigma <- as.character(substitute(link.sigma))
    if(mode(link.mu) != "character" && mode(link.mu) != "name")
        link.mu <- as.character(substitute(link.mu))

    new("vglmff",
    blurb=c("LMS Quantile Regression (Box-Cox transformation to a Gamma distribution)\n",
            "Links:    ",
            "lambda",
            ", ",
            namesof("mu", link=link.mu),
            ", ",
            namesof("sigma", link=link.sigma)),
    constraints=eval(substitute(expression({
        constraints <- cm.zero.vgam(constraints, x, .zero, M)
    }), list(.zero=zero))),
    initialize=eval(substitute(expression({
        if(any(y<0, na.rm = TRUE))
            stop("negative responses not allowed")

        predictors.names <- c("lambda",
                     namesof("mu",  .link.mu,  short= TRUE),
                     namesof("sigma",  .link.sigma,  short= TRUE))

        if(!length(etastart)) {

            fit500=vsmooth.spline(x=x[,min(ncol(x),2)],y=y,w=w, df= .dfmu.init)
            fv.init = c(predict(fit500, x=x[,min(ncol(x),2)])$y)

            lambda.init = if(is.Numeric( .init.lambda)) .init.lambda else 1.0

            sigma.init <- if(is.null(.init.sigma)) {
               myratio=((y/fv.init)^lambda.init-1)/lambda.init #~(0,var=sigma^2)
                if(is.numeric( .dfsigma.init) && is.finite( .dfsigma.init)) {
                    fit600 = vsmooth.spline(x=x[,min(ncol(x),2)],
                                            y=(myratio)^2,
                                            w=w, df= .dfsigma.init)
                    sqrt(c(abs(predict(fit600, x=x[,min(ncol(x),2)])$y)))
                } else 
                    sqrt(var(myratio))
            } else .init.sigma

            etastart <- cbind(lambda.init,
                              theta2eta(fv.init,  .link.mu),
                              theta2eta(sigma.init,  .link.sigma))
        }
    }), list(.link.sigma=link.sigma,
             .link.mu=link.mu,
             .dfmu.init=dfmu.init,
             .dfsigma.init=dfsigma.init,
             .init.lambda=init.lambda,
             .init.sigma=init.sigma))),
    inverse=eval(substitute(function(eta, extra=NULL) {
        eta[,2] <- eta2theta(eta[,2], .link.mu)
        eta[,3] <- eta2theta(eta[,3], .link.sigma)
        qtplot.lms.bcg(percentiles= .percentiles, eta=eta)
    }, list(.percentiles=percentiles,
            .link.mu=link.mu,
            .link.sigma=link.sigma))),
    last=eval(substitute(expression({
        misc$percentiles <- .percentiles
        misc$links = c(lambda = "identity", mu = .link.mu, sigma = .link.sigma)
        misc$true.mu <- FALSE    # $fitted is not a true mu
        if(control$cdf) {
            post$cdf = cdf.lms.bcg(y, eta0=matrix(c(lambda,mymu,sigma), 
                ncol=3, dimnames=list(dimnames(x)[[1]], NULL)))
        }
    }), list(.percentiles=percentiles,
            .link.mu=link.mu,
            .link.sigma=link.sigma))),

    loglikelihood=eval(substitute(
        function(mu,y,w, residuals= FALSE, eta, extra=NULL) {
            lambda <- eta[,1]
            mu <- eta2theta(eta[,2], .link.mu)
            sigma <- eta2theta(eta[,3], .link.sigma)
            g <- (y/mu)^lambda
            theta <- 1 / (sigma * lambda)^2
         if(residuals) stop("loglikelihood residuals not implemented yet") else
            sum(w * (log(abs(lambda)) + theta*(log(theta)+log(g)-g) - 
                     lgamma(theta) - log(y)))
        }, list(.link.sigma=link.sigma, .link.mu=link.mu))),
    vfamily=c("lms.bcg", "lmscreg"),
    deriv=eval(substitute(expression({
        lambda <- eta[,1]
        mymu <- eta2theta(eta[,2], .link.mu)
        sigma <- eta2theta(eta[,3], .link.sigma)

        g <- (y/mymu)^lambda
        theta <- 1 / (sigma * lambda)^2
        dd <- digamma(theta)

        dl.dlambda <- (1 + 2*theta*(dd+g-1-log(theta) -
                      0.5 * (g+1)*log(g))) / lambda
        dl.dmu <- lambda * theta * (g-1) / mymu
        dl.dsigma <- 2*theta*(dd+g-log(theta * g)-1) / sigma
        dsigma.deta <- dtheta.deta(sigma, link=.link.sigma)

        cbind(dl.dlambda,
              dl.dmu * dtheta.deta(mymu, link= .link.mu),
              dl.dsigma * dsigma.deta) * w
    }), list(.link.sigma=link.sigma, .link.mu=link.mu))),
    weight=eval(substitute(expression({
        tt <- trigamma(theta)
 
        wz <- matrix(0, n, 6)

        if(TRUE) {
            part2 <- dd + 2/theta - 2*log(theta)
            wz[,iam(1,1,M)] <- (1 + theta*(tt*(1+4*theta) - 4*(1+1/theta) -
                log(theta)*(2/theta - log(theta)) + dd*part2)) / lambda^2
        } else {
            temp <- mean( g*(log(g))^2 )
            wz[,iam(1,1,M)] <- (4*theta*(theta*tt-1) -1+ theta*temp)/lambda^2
        }

        wz[,iam(2,2,M)] <- 1 / (mymu*sigma)^2  *
                           dtheta.deta(mymu, .link.mu)^2
        wz[,iam(3,3,M)] <- (4*theta*(theta*tt-1) / sigma^2) *
                           dtheta.deta(sigma, .link.sigma)^2
        wz[,iam(1,2,M)] <- -theta * (dd + 1/theta - log(theta)) / mymu
        wz[,iam(1,2,M)] <- wz[,iam(1,2,M)] * 
                           dtheta.deta(mymu, .link.mu)
        wz[,iam(1,3,M)] <- 2 * theta^1.5 * (2 * theta * tt - 2 -
                           1/theta) * dtheta.deta(sigma, .link.sigma)
        wz * w
    }), list(.link.sigma=link.sigma, .link.mu=link.mu))))
}






dy.dyj <- function(psi, lambda, epsilon=sqrt(.Machine$double.eps)) {

    L = max(length(psi), length(lambda))
    psi = rep(psi, len=L); lambda = rep(lambda, len=L);
    ifelse(psi>0, (1 + psi * lambda)^(1/lambda - 1),
                  (1 - (2-lambda) * psi)^((lambda - 1)/(2-lambda)))
}

dyj.dy <- function(y, lambda) {
    L = max(length(y), length(lambda))
    y = rep(y, len=L); lambda = rep(lambda, len=L);
    ifelse(y>0, (1 + y)^(lambda - 1), (1 - y)^(1 - lambda))
}

yeo.johnson <- function(y, lambda, derivative=0,
                        epsilon=sqrt(.Machine$double.eps), inverse= FALSE)
{

    if(length(derivative)!=1 || derivative<0 || round(derivative)!=derivative)
        stop("derivative must be a non-negative integer")
    ans <- y

    if((length(epsilon) != 1) || (epsilon <= 0))
        stop("epsilon must be a single positive number")
    lambda <- rep(lambda, len=length(y))  # lambda may be of length 1

    if(inverse) {
        if(derivative!=0)
            stop("derivative must 0 when inverse=TRUE")
        index <- y >= 0 & abs(lambda) > epsilon
        ans[index] <- (y[index]*lambda[index] + 1)^(1/lambda[index]) - 1
        index <- y >= 0 & abs(lambda) <= epsilon 
        ans[index] <- exp(y[index]) - 1
        index <- y <  0 & abs(lambda-2) > epsilon
        ans[index] <- 1-(-(2-lambda[index])*y[index]+1)^(1/(2-lambda[index]))
        index <- y <  0 & abs(lambda-2) <= epsilon
        ans[index] <- 1 - exp(-y[index])
        return(ans)
    }
    if(derivative==0) {
        index <- y >= 0 & abs(lambda) > epsilon
        ans[index] <- ((y[index]+1)^(lambda[index]) - 1) / lambda[index]
        index <- y >= 0 & abs(lambda) <= epsilon 
        ans[index] <- log(y[index]+1)
        index <- y <  0 & abs(lambda-2) > epsilon
        ans[index] <- -((-y[index]+1)^(2-lambda[index]) - 1) / (2-lambda[index])
        index <- y <  0 & abs(lambda-2) <= epsilon
        ans[index] <- -log(-y[index]+1)
    } else {
        psi <- yeo.johnson(y, lambda, derivative=derivative-1,
                           epsilon=epsilon, inverse=inverse)
        index <- y >= 0 & abs(lambda) > epsilon
        ans[index] <- ( (y[index]+1)^(lambda[index]) *
                      (log(y[index]+1))^(derivative) - derivative *
                      psi[index] ) / lambda[index]
        index <- y >= 0 & abs(lambda) <= epsilon 
        ans[index] <- (log(y[index]+1))^(derivative + 1) / (derivative + 1)
        index <- y <  0 & abs(lambda-2) > epsilon
        ans[index] <- -( (-y[index]+1)^(2-lambda[index]) *
                      (-log(-y[index]+1))^(derivative) - derivative *
                      psi[index] ) / (2-lambda[index])
        index <- y <  0 & abs(lambda-2) <= epsilon
        ans[index] <- (-log(-y[index]+1))^(derivative + 1) / (derivative + 1)
    }
    ans
}


dpsi.dlambda.yjn = function(psi, lambda, mymu, sigma, derivative=0, small=1e-8) {

    answer = matrix(as.numeric(NA), length(mymu), derivative+1) 

    if(length(psi) < length(lambda))
        psi = rep(psi, length=length(lambda))

    CC = psi >= 0
    BB = ifelse(CC, lambda, -2+lambda)
    AA = 1 + psi * BB 

    if(derivative>0) {
        answer[,1:derivative] =
            Recall(psi=psi, lambda=lambda, mymu=mymu, sigma=sigma,
                   derivative=derivative-1, small=small) 
        temp8 = answer[,derivative] * derivative
    } else { 
        temp8 = 0
    }

    answer[,1+derivative] =
           (AA * (log(AA)/BB)^derivative -
           temp8) / BB



    pos = (CC & abs(lambda) <= small) | (!CC & abs(lambda-2) <= small)
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


lms.yjn <- function(percentiles=c(25,50,75),
                    zero=NULL,
                    link.lambda="identity",
                    link.sigma="loge",
                    dfmu.init=4,
                    dfsigma.init=2,
                    init.lambda=1.0,
                    init.sigma=NULL,
                    rule=c(10,5),
                    yoffset=NULL,
                    diagW=FALSE, iters.diagW=6)
{



    if(mode(link.sigma) != "character" && mode(link.sigma) != "name")
        link.sigma <- as.character(substitute(link.sigma))
    if(mode(link.lambda) != "character" && mode(link.lambda) != "name")
        link.lambda <- as.character(substitute(link.lambda))

    rule = rule[1] # Number of points (common) for all the quadrature schemes
    if(rule != 5 && rule != 10)
        stop("only rule=5 or 10 is supported")

    new("vglmff",
    blurb=c("LMS Quantile Regression (Yeo-Johnson transformation to normality)\n",
            "Links:    ",
            namesof("lambda", link=link.lambda),
            ", mu, ",
            namesof("sigma", link=link.sigma)),
    constraints=eval(substitute(expression({
        constraints <- cm.zero.vgam(constraints, x, .zero, M)
    }), list(.zero=zero))),
    initialize=eval(substitute(expression({
        predictors.names <- c(namesof("lambda", .link.lambda, short= TRUE),
                "mu", namesof("sigma",  .link.sigma,  short= TRUE))

        y.save <- y
        yoff = if(is.Numeric( .yoffset)) .yoffset else -median(y) 
        extra$yoffset = yoff
        y <- y + yoff

        if(!length(etastart)) {

            lambda.init = if(is.Numeric( .init.lambda)) .init.lambda else 1.0

            y.tx = yeo.johnson(y, lambda.init)
            fit700=vsmooth.spline(x=x[,min(ncol(x),2)],y=y.tx,w=w, df= .dfmu.init)
            fv.init = c(predict(fit700, x=x[,min(ncol(x),2)])$y)

            sigma.init = if(!is.Numeric(.init.sigma)) {
                              if(is.Numeric( .dfsigma.init)) {
                                   fit710 = vsmooth.spline(x=x[,min(ncol(x),2)],
                                            y=(y.tx - fv.init)^2,
                                            w=w, df= .dfsigma.init)
                                   sqrt(c(abs(predict(fit710,
                                        x=x[,min(ncol(x),2)])$y)))
                              } else {
                                   sqrt( sum( w * (y.tx - fv.init)^2 ) / sum(w) )
                              }
                          } else
                              .init.sigma

            etastart <- cbind(theta2eta(lambda.init,  .link.lambda),
                              fv.init,
                              theta2eta(sigma.init,  .link.sigma))

        }
    }), list(.link.sigma=link.sigma,
             .link.lambda=link.lambda,
             .dfmu.init=dfmu.init,
             .dfsigma.init=dfsigma.init,
             .init.lambda=init.lambda,
             .yoffset=yoffset,
             .init.sigma=init.sigma))),
    inverse=eval(substitute(function(eta, extra=NULL) {
        eta[,1] = eta2theta(eta[,1], .link.lambda)
        eta[,3] = eta2theta(eta[,3], .link.sigma)
        qtplot.lms.yjn(percentiles= .percentiles, eta=eta, yoffset= extra$yoff)
    }, list(.percentiles=percentiles,
            .link.lambda=link.lambda,
            .link.sigma=link.sigma))),
    last=eval(substitute(expression({
        misc$percentiles <- .percentiles
        misc$links = c(lambda= .link.lambda, mu= "identity", sigma= .link.sigma)
        misc$true.mu <- FALSE    # $fitted is not a true mu
        misc[["yoffset"]] = extra$yoff   # zz Splus6.0 bug: sometimes the name is lost

        y <- y.save   # Restore back the value; to be attached to object

        if(control$cdf) {
            post$cdf = cdf.lms.yjn(y + misc$yoffset,
                eta0=matrix(c(lambda,mymu,sigma), 
                ncol=3, dimnames=list(dimnames(x)[[1]], NULL)))
        }
    }), list(.percentiles=percentiles,
            .link.lambda=link.lambda,
            .link.sigma=link.sigma))),
    loglikelihood=eval(substitute(
        function(mu,y,w, residuals= FALSE, eta, extra=NULL) {
            lambda <- eta2theta(eta[,1], .link.lambda)
            mu <- eta[,2]
            sigma <- eta2theta(eta[,3], .link.sigma)
            psi <- yeo.johnson(y, lambda)
         if(residuals) stop("loglikelihood residuals not implemented yet") else
            sum(w * (-log(sigma) - 0.5 * ((psi-mu)/sigma)^2 +
                     (lambda-1) * sign(y) * log(abs(y)+1)))
        }, list(.link.sigma=link.sigma, .link.lambda=link.lambda))),
    vfamily=c("lms.yjn", "lmscreg"),
    deriv=eval(substitute(expression({
        lambda <- eta2theta(eta[,1], .link.lambda)
        mymu <- eta[,2]
        sigma <- eta2theta(eta[,3], .link.sigma)

        psi <- yeo.johnson(y, lambda)
        d1 <- yeo.johnson(y, lambda, deriv=1)
        AA = (psi - mymu) / sigma 

        dl.dlambda <- -AA * d1 /sigma + sign(y) * log(abs(y)+1) 
        dl.dmu <- AA / sigma 
        dl.dsigma <- (AA^2 -1) / sigma
        dlambda.deta <- dtheta.deta(lambda, link=.link.lambda)
        dsigma.deta <- dtheta.deta(sigma, link=.link.sigma)

        cbind(dl.dlambda * dlambda.deta,
              dl.dmu,
              dl.dsigma * dsigma.deta) * w
    }), list(.link.sigma=link.sigma, .link.lambda=link.lambda))),
    weight=eval(substitute(expression({
        wz <- matrix(0, n, 6)


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
                     eps=as.double(1.0e-5))$ans #zz adjust eps for more accuracy
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

        wz[,iam(1,1,M)] <- wz[,iam(1,1,M)] * dlambda.deta^2
        wz[,iam(1,2,M)] <- wz[,iam(1,2,M)] * dlambda.deta
        wz[,iam(1,3,M)] <- wz[,iam(1,3,M)] * dsigma.deta * dlambda.deta
        if( .diagW && iter <= .iters.diagW) {
            wz[,iam(1,2,M)] = wz[,iam(1,3,M)] = 0
        }
        wz[,iam(2,3,M)] <- wz[,iam(2,3,M)] * dsigma.deta
        wz[,iam(3,3,M)] <- wz[,iam(3,3,M)] * dsigma.deta^2

        wz = wz * w
        wz
    }), list(.link.sigma=link.sigma,
             .rule=rule,
             .diagW=diagW,
             .iters.diagW=iters.diagW,
             .link.lambda=link.lambda))))
}



lmscreg.control <- function(cdf= TRUE, at.arg=NULL, x0=NULL, ...)
{

    if(!is.logical(cdf)) {
        warning("\"cdf\" is not logical; using TRUE instead")
        cdf = T
    }
    list(cdf=cdf, at.arg=at.arg, x0=x0)
}





if(FALSE)
lms.yjn1 = function(percentiles=c(25,50,75),
                    zero=NULL,
                    link.lambda="identity",
                    dfmu.init=4,
                    dfsigma.init=2,
                    init.lambda=1.0,
                    yoffset=NULL)
{
    if(mode(link.lambda) != "character" && mode(link.lambda) != "name")
        link.lambda <- as.character(substitute(link.lambda))

    new("vglmff",
    blurb=c("LMS Quantile Regression (Yeo-Johnson transformation to normality)\n",
            "Links:    ",
            namesof("lambda", link=link.lambda)),
    constraints=eval(substitute(expression({
        constraints <- cm.zero.vgam(constraints, x, .zero, M)
    }), list(.zero=zero))),
    initialize=eval(substitute(expression({
        predictors.names <- c(namesof("lambda", .link.lambda, short= TRUE))

        y.save <- y
        yoff = if(is.Numeric( .yoffset)) .yoffset else -median(y) 
        extra$yoffset = yoff
        y <- y + yoff

        if(!length(etastart)) {

            lambda.init = if(is.Numeric( .init.lambda)) .init.lambda else 1.0

            y.tx = yeo.johnson(y, lambda.init)
            fit700=vsmooth.spline(x=x[,min(ncol(x),2)],y=y.tx,w=w, df= .dfmu.init)
            fv.init = c(predict(fit700, x=x[,min(ncol(x),2)])$y)
            extra$mymu = fv.init

            sigma.init = if(TRUE) {
                              if(is.Numeric( .dfsigma.init)) {
                                   fit710 = vsmooth.spline(x=x[,min(ncol(x),2)],
                                            y=(y.tx - fv.init)^2,
                                            w=w, df= .dfsigma.init)
                                   sqrt(c(abs(predict(fit710,
                                        x=x[,min(ncol(x),2)])$y)))
                              } else {
                                   sqrt( sum( w * (y.tx - fv.init)^2 ) / sum(w) )
                              }
                          } else
                              1
            extra$sigma = sigma.init

            etastart <- cbind(theta2eta(lambda.init,  .link.lambda))
        }
    }), list(.link.lambda=link.lambda,
             .dfmu.init=dfmu.init,
             .dfsigma.init=dfsigma.init,
             .init.lambda=init.lambda,
             .yoffset=yoffset,
             ))),
    inverse=eval(substitute(function(eta, extra=NULL) {
        eta = eta2theta(eta, .link.lambda)
        eta2 = extra$mymu
        eta3 = extra$sigma
        qtplot.lms.yjn(percentiles= .percentiles,
            eta=cbind(c(eta),eta2,eta3), yoffset= extra$yoff)
    }, list(.percentiles=percentiles,
            .link.lambda=link.lambda))),
    last=eval(substitute(expression({
        misc$percentiles <- .percentiles
        misc$links <- c(lambda = .link.lambda)
        misc$true.mu <- FALSE    # $fitted is not a true mu
        misc[["yoffset"]] = extra$yoff   # zz Splus6.0 bug: sometimes the name is lost

        y <- y.save   # Restore back the value; to be attached to object

        if(control$cdf) {
            post$cdf = cdf.lms.yjn(y + misc$yoffset,
                eta0=matrix(c(lambda,mymu,sigma), 
                ncol=3, dimnames=list(dimnames(x)[[1]], NULL)))
        }
    }), list(.percentiles=percentiles,
            .link.lambda=link.lambda))),
    loglikelihood=eval(substitute(
        function(mu,y,w, residuals= FALSE, eta, extra=NULL) {
            lambda <- eta2theta(eta, .link.lambda)
            mu <- extra$mymu
            sigma <- extra$sigma
            psi <- yeo.johnson(y, lambda)
         if(residuals) stop("loglikelihood residuals not implemented yet") else
            sum(w * (-log(sigma) - 0.5 * ((psi-mu)/sigma)^2))
        }, list(.link.lambda=link.lambda))),
    vfamily=c("lms.yjn", "lmscreg"),
    deriv=eval(substitute(expression({
        lambda <- eta2theta(eta, .link.lambda)
        psi <- yeo.johnson(y, lambda)

        fit8 <- vsmooth.spline(x=x[,min(ncol(x),2)],y=psi,w=w, df= .dfmu.init)
        mymu = c(predict(fit8, x=x[,min(ncol(x),2)])$y)
        extra$mymu = mymu
        fit8 <- vsmooth.spline(x=x[,min(ncol(x),2)],y=(psi-mymu)^2,w=w,
                               df= .dfsigma.init)
        sigma = sqrt(c(predict(fit8, x=x[,min(ncol(x),2)])$y))
        extra$sigma = sigma

        d1 <- yeo.johnson(y, lambda, deriv=1)
        AA = (psi - mymu) / sigma 

        dl.dlambda = -AA * d1 / sigma
 warning("dl.dlambda is wrong")
        dlambda.deta <- dtheta.deta(lambda, link=.link.lambda)
        cbind(dl.dlambda * dlambda.deta) * w
    }), list(.dfmu.init=dfmu.init, 
             .dfsigma.init=dfsigma.init, 
             .link.lambda=link.lambda))),
    weight=eval(substitute(expression({
        wz = (d1 / sigma)^2    # Approximate
        wz = ifelse(wz <= 0, 1.0e-9, wz)
        wz = wz * w
        wz
    }), list(.link.lambda=link.lambda))))
}




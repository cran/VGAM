# These functions are
# Copyright (C) 1998-2009 T.W. Yee, University of Auckland. All rights reserved.



vgam.control <- function(all.knots=FALSE,
                         bf.epsilon=1e-7,
                         bf.maxit=30, 
                         checkwz=TRUE,
                         criterion = names(.min.criterion.VGAM), 
                         epsilon=1e-7,
                         maxit=30,
                         na.action=na.fail,
                         nk=NULL,
                         save.weight=FALSE,
                         se.fit=TRUE,
                         trace=FALSE,
                         wzepsilon = .Machine$double.eps^0.75,
                         ...)
{





    if(mode(criterion) != "character" && mode(criterion) != "name")
        criterion <- as.character(substitute(criterion))
    criterion <- pmatch(criterion[1], names(.min.criterion.VGAM), nomatch=1)
    criterion <- names(.min.criterion.VGAM)[criterion]

    if(!is.logical(checkwz) || length(checkwz) != 1)
        stop("bad input for 'checkwz'")
    if(!is.Numeric(wzepsilon, allow=1, positive=TRUE))
        stop("bad input for 'wzepsilon'")

    if(length(all.knots) > 1)
        warning("all.knots should be of length 1; using first value only")
    if(!is.Numeric(bf.epsilon, allow=1, posit=TRUE)) {
        warning("bad input for 'bf.epsilon'; using 0.00001 instead")
        bf.epsilon <- 0.00001
    }
    if(!is.Numeric(bf.maxit, allow=1, posit=TRUE, integ=TRUE)) {
        warning("bad input for 'bf.maxit'; using 20 instead")
        bf.maxit <- 20
    }
    if(!is.Numeric(epsilon, allow=1, posit=TRUE)) {
        warning("bad input for 'epsilon'; using 0.0001 instead")
        epsilon <- 0.0001
    }
    if(!is.Numeric(maxit, allow=1, posit=TRUE, integ=TRUE)) {
        warning("bad input for 'maxit'; using 20 instead")
        maxit <- 20
    }

    convergence <- expression({
        switch(criterion,
        coefficients=if(iter==1) iter<maxit else (iter<maxit &&
        max(abs(new.coeffs - old.coeffs)/(abs(old.coeffs)+epsilon)) > epsilon),
        abs(old.crit-new.crit)/(abs(old.crit)+epsilon) > epsilon && iter<maxit)
    })

    list(all.knots=as.logical(all.knots)[1],
         bf.epsilon=bf.epsilon, 
         bf.maxit=bf.maxit, 
         checkwz=checkwz,
         convergence=convergence,
         criterion=criterion,
         epsilon=epsilon, 
         maxit=maxit, 
         nk=nk,
         min.criterion = .min.criterion.VGAM,
         save.weight=as.logical(save.weight)[1],
         se.fit=as.logical(se.fit)[1],
         trace=as.logical(trace)[1],
         wzepsilon = wzepsilon)
}


vgam.nlchisq <- function(qr, resid, wz, s, deriv, U, smooth.labels,
                         assign, M, n, constraints)
{

        attr(qr, "class") = "qr" 
        class(qr) <- "qr"

    if(!is.matrix(s)) s <- as.matrix(s)
    if(!is.matrix(wz)) wz <- as.matrix(wz)
    if(!is.matrix(deriv)) deriv <- as.matrix(deriv)
    if(!is.matrix(resid)) resid <- as.matrix(resid)

    trivc <- trivial.constraints(constraints)

    ans <- rep(as.numeric(NA), length=ncol(s))
    Uderiv <- vbacksub(U, t(deriv), M=M, n=n)    # \bU_i^{-1} \biu_i
    ptr <- 0
    for(i in 1:length(smooth.labels)) {
        cmat <- constraints[[ smooth.labels[i] ]]
        index <- (ptr+1):(ptr+ncol(cmat))

        for(j in index) {
            yy <- t(cmat[,j-ptr,drop=FALSE])
            yy <- kronecker(s[,j,drop=FALSE], yy)  # n x M
            Us <- mux22(U, yy, M=M, upper=TRUE, as.matrix=TRUE)  # n * M

            Uss <- matrix(c(t(Us)), nrow=n*M, ncol=1)

            Rsw <- qr.resid(qr, Uss)

            vRsw <- matrix(Rsw, nrow=n, ncol=M, byrow=TRUE)
            newans <- vbacksub(U, t(vRsw), M=M, n=n)

            ans[j] <- sum(vRsw^2 + 2 * newans * deriv)

        }
        ptr <- ptr + ncol(cmat)
    }

    names(ans) <- dimnames(s)[[2]]
    ans
}
    






# These functions are Copyright (C) 1998-2013 T. W. Yee  All rights reserved

# 26/6/98; family.sur.q
# 20110406; renamed to family.sur.R

# zz; does or doesn't handle? : 
# vglm(Sur(mydataframe), sur, ...), i.e., 1st coln
# of mydataframe is the response. 



# History
# 20110406; editing it to bring it up to scratch.
# 20130125; trying to get SUR() going.



# --------------------------------------------------------------------
# Maybe should call this surff()??:


 SUR <- function(
                 mle.normal = FALSE,
                 divisor = c("n", "n-max(pj,pk)", "sqrt((n-pj)*(n-pk))"),
#                estimator = c("classical", "iterative"),
                 parallel = FALSE, 
                 apply.parint = TRUE,
#                zero = NULL,
                 Varcov = NULL,
                 matrix.arg = FALSE) {
# Notes:
# 1. Varcov may be assigned a solve(wz) (=solve(\bSigma)),
#    and matrix.arg tells what format it is in.
# 2. Based a little on normal1().
# 3. Set maxit = 1   for Zellner's estimator (2-stage).
#    Set maxit = 111 for iterative GLS === IGLS.


# Wrong:
# 1. "2stage"     == Zellners estimator.
#    "iterative"  == iterative GLS === IGLS.
#    "MLE.normal" == not yet done.
#    Or "maxit.sur = 2"?


# Last modified:
# 20130125; trying to get SUR() going.
# 20130126; seems to work basically but not the above arguments.
#   A lot more work needed.
# 20130130; seems to work.
#   Removed 'zero' argument.


# Yettodo:
# 2013013 ; argument 'mle.normal' is logical.


#print("20130129; in SUR()")


  lmean <- "identity"
  lsdev <- "loge"
  emean <- list()
  esdev <- list()


  if (!is.logical(mle.normal) ||
      length(mle.normal) != 1)
    stop("argument 'mle.normal' must be a single logical")

  if (!is.logical(apply.parint) ||
      length(apply.parint) != 1)
    stop("argument 'apply.parint' must be a single logical")


# if(mode(estimator) != "character" && mode(estimator) != "name")
#   estimator <- as.character(substitute(estimator))
# estimator <- match.arg(estimator,
#                      c("classical", "iterative"))[1]
#print(paste('estimator =', estimator))


  divisor <- match.arg(divisor,
      c("n", "n-max(pj,pk)", "sqrt((n-pj)*(n-pk))"))[1]
#print("divisor")
#print( divisor )

  if (mle.normal && divisor != "n")
    warning("MLE requires 'n' as the value of argument 'divisor'. ",
            "The solution will probably not be the MLE")


  ret.ff <-
  new("vglmff",
  blurb = c("Seemingly unrelated regressions"),
  constraints = eval(substitute(expression({
    constraints <- cm.vgam(matrix(1, M, 1), x,
                           .parallel , constraints,
                           apply.int = .apply.parint )
#   constraints <- cm.zero.vgam(constraints, x, .zero , M)
  }), list( .parallel = parallel,
#           .zero = zero, 
            .apply.parint = apply.parint ))),

# deviance = function(y, mu, w, residuals = FALSE,
#                     eta = NULL, extra = NULL) {
# Returns the residual sum of squares
# Nb. extra$wz is wz

#print("head(y - mu)")
#print( head(y - mu) )
#print("head(extra$wz)")
#print( head(extra$wz) )

#   M <- if (length(extra$M)) extra$M else ifelse(is.matrix(y), ncol(y), 1)
#   if (residuals) {
#     if (M > 1) NULL else (y-mu) * sqrt(extra$wz)
#   } else {
#     ResSS.vgam(y - mu, extra$wz, M = M)
#   }
# },

  infos = eval(substitute(function(...) {
    list(Musual = 1,  # zz???
#        zero = .zero ,
#        link = .link ,
         parallel = .parallel ,
         multipleResponses = TRUE )
  }, list( .parallel = parallel ))),

  initialize = eval(substitute(expression({

    if (!is.matrix(y) || ncol(y) == 1)
      stop("response must be a matrix with at least 2 columns")
    ncoly <- ncol(y)

#print("constraints")
#print( constraints )
   if (is.logical( .parallel ) &&
       .parallel &&
       !all(as.logical(trivial.constraints(constraints))))
     warning("setting 'parallel = TRUE' with nontrivial constraints may not ",
             "make sense")

   temp5 <-
    w.y.check(w = w, y = y,
#             Is.positive.y = TRUE,
              ncol.w.min = 1,
              ncol.w.max = 1,
              ncol.y.max = Inf,
              Is.integer.y = FALSE,
              Is.positive.y = FALSE,
              out.wy = TRUE,
              colsyperw = ncoly,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y
#print("head(w)")
#print( head(w) )

    if (!all(w[1, 1] == w))
      stop("all prior 'weights' must currently have equal values")


    ncoly <- ncol(y)
    Musual <- 1
    extra$ncoly <- ncoly
    extra$Musual <- Musual
    M <- Musual * ncoly


    predictors.names <- if (!length(ddd <- dimnames(y)[[2]]))
        paste("Y", 1:M, sep = "") else ddd


#   if ( .estimator == "classical")
#       maxit <- 1


# Iteration may lead to an increase in RSS 
#   if ( .estimator == "iterative")
#     half.stepsizing <- FALSE


# Assign "extra$wz" something corresponding to the M x M identity matrix.
    extra$wz <- matrix(1, nrow(x), M)


    if (!length(etastart)) {
# Note: it is a good idea to start with the OLS estimators here first.
      etastart <- matrix(0, n, M)


      Blist.early <- process.constraints(constraints, x, M,
                                         specialCM = specialCM)
#print("Blist.early")
#print( Blist.early )
      X_vlm.early  <- lm2vlm.model.matrix(x, Blist.early, xij = control$xij,
                                          Xm2 = Xm2)
#print("head(X_vlm.early)")
#print( head(X_vlm.early) )

      Hmatrices <- matrix(c(unlist(Blist.early)), nrow = M)
      jay.index <- 1:ncol(Hmatrices)


      extra$ncols_X_lm <- numeric(ncoly)
      for (jay in 1:ncoly) {
# model.matrix(fit, lapred.index = 1, type = "lm")
#print("Hmatrices")
#print( Hmatrices )
# 20121231; this code adapted from model.matrixvlm():
#       lapred.index <- jay.index[jay]
#       index0 <- Hmatrices[jay, ] != 0  # Orig.
#       Index0 <- Hmatrices[lapred.index, ] != 0
#       X_lm_jay <- X_vlm[(0:(n_lm - 1)) * M + lapred.index, Index0,
#                         drop = FALSE]

        X_lm_jay <- vlm2lm.model.matrix(x_vlm = X_vlm.early,
                                        Blist = Blist.early,
                                        which.lp = jay, M = M)
#print("head(X_lm_jay)")
#print( head(X_lm_jay) )

# This is useful, e.g,. for changing the denominator
        extra$ncols_X_lm[jay] <- ncol(X_lm_jay)

        etastart[, jay] <- y[, jay] -
                           lsfit(x = X_lm_jay, y = y[, jay],
                                 wt = c(w), intercept = FALSE)$residuals
      }  # jay
    }  # !length(etastart)
  }), list(
#           .estimator = estimator,
            .parallel = parallel 
          ))),
  linkinv = function(eta, extra = NULL) eta, 
  last = eval(substitute(expression({

    Musual <- extra$Musual
    misc$link <- c(rep( .lmean , length = ncoly))
    temp.names <- predictors.names
#   temp.names <- temp.names[interleave.VGAM(Musual * ncoly, M = Musual)]
    names(misc$link) <- temp.names
#print("head(w)")
#print( head(w) )

    misc$earg <- vector("list", Musual * ncoly)
    names(misc$earg) <- temp.names
    for(ii in 1:ncoly) {
      misc$earg[[Musual*ii]] <- .emean
    }
    names(misc$earg) <- temp.names

    misc$Musual <- Musual
    misc$expected <- TRUE
    misc$divisor <- .divisor
    misc$values.divisor <- round(n / ratio.df)

  }), list( .lmean = lmean, .lsdev = lsdev,
            .emean = emean, .esdev = esdev,
            .divisor = divisor
          ))),

# linkfun = function(mu, extra = NULL) mu,
  vfamily = "SUR",


  deriv = eval(substitute(expression({
#print("in @deriv of SUR()")
#print(paste("iter =", iter))
    mymu <- eta
    iam.indices <- iam(NA, NA, M = M, both = TRUE)
#print("iam.indices")
#print( iam.indices )
#print("y")
#print( y )
#print("mu")
#print( mu )
    resmat <- y - mymu
    Sigma.elts <- colMeans(resmat[, iam.indices$row.index] *
                           resmat[, iam.indices$col.index])

    if ( .divisor != "n") {
# Make an adjustment for the denominator (above assumes "n")
# Here, ratio.df >= 1.
      ratio.df <- n / switch( .divisor ,
        "n-max(pj,pk)" = n - pmax(extra$ncols_X_lm[iam.indices$row.index],
                                  extra$ncols_X_lm[iam.indices$col.index]),
        "sqrt((n-pj)*(n-pk))" =
        sqrt((n - extra$ncols_X_lm[iam.indices$row.index]) *
             (n - extra$ncols_X_lm[iam.indices$col.index])),
        stop("argument 'divisor' unmatched"))
#print("ratio.df")
#print( ratio.df )
      Sigma.elts <- Sigma.elts * ratio.df
    } else {
      ratio.df <- rep(1, length = M*(M+1)/2)
    }

#print("Sigma.elts")
#print( Sigma.elts )
    Sigma.mat <- matrix(0, M, M)
    Sigma.mat[cbind(iam.indices$row.index,
                    iam.indices$col.index)] <- Sigma.elts
    Sigma.mat[cbind(iam.indices$col.index,
                    iam.indices$row.index)] <- Sigma.elts

#print("Sigma.mat")
#print( Sigma.mat )
# Cholesky is more efficient than solve()
    invSigma.mat <- chol2inv(chol(Sigma.mat))
#   invSigma.mat <- solve(Sigma.mat)  # Inefficient
#print("invSigma.mat")
#print( invSigma.mat )


# dl.dmu returns \bW_i (\biy_i - \bmu_i)
    temp3 <- matrix(invSigma.mat[cbind(iam.indices$row.index,
                                       iam.indices$col.index)],
                    M*(M+1)/2, n)
    dl.dmu <- mux22(temp3, y - mymu, M = M,
                    upper = FALSE, as.matrix = TRUE)
#print("dim(dl.dmu)")
#print( dim(dl.dmu) )
#print("head(dl.dmu)")
#print( head(dl.dmu) )
#   dl.dmu <- (y - mymu) / sdev^2  # For normal1()
    dmu.deta <- dtheta.deta(mymu,   .lmean , earg = .emean )
#print("head(dmu.deta)")
#print( head(dmu.deta) )

    c(w) * dl.dmu * dmu.deta
  }), list( .lmean = lmean,
            .emean = emean,
            .divisor = divisor ))),


  weight = eval(substitute(expression({
#print("in @weight of SUR()")


# Overwrite invSigma.mat with the inverse variance, if given.
    if (length( .Varcov )) {
      Sigma.mat <- if ( .matrix.arg ) .Varcov else {
                     temp.vec <- rep( .Varcov , len = M*(M+1)/2)
                     temp.mat <- matrix(0, M, M)
                     temp.mat[cbind(iam.indices$col.index,
                                    iam.indices$row.index)] <- temp.vec
                     temp.mat[cbind(iam.indices$row.index,
                                    iam.indices$col.index)] <- temp.vec
#                    temp.mat <- chol2inv(chol(temp.mat))
                     temp.mat
                   }
      invSigma.mat <- chol2inv(chol(Sigma.mat))
    }


    wz <-
    extra$wz <- c(w) * matrix(invSigma.mat[cbind(iam.indices$col.index,
                                                 iam.indices$row.index)],
                              n, M*(M+1)/2, byrow = TRUE)
    extra$Sigma.mat <- Sigma.mat
    extra$invSigma.mat <- invSigma.mat

#print("head(wz)")
#print( head(wz) )
    wz
  }), list( .divisor = divisor,
#           .estimator = estimator,
            .Varcov = Varcov,
            .matrix.arg = matrix.arg ))))



  if (mle.normal) {
# Add a 'loglikelihood' slot to the object.
# This code based on normal1().
# Note wz is retrieved from 'extra', and 'wz' has only
# one general symmetric pos-definite matrix that is replicated
# a lot.

# Yettodo: if "all prior 'weights' must currently have equal values" is
# relaxed then have to do some code changes??

    ret.ff@loglikelihood <-
      function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
      M <- if (is.matrix(y)) ncol(y) else 1
      n <- if (is.matrix(y)) nrow(y) else length(y)

# Orig:
#     wz <- VGAM.weights.function(w = w, M = M, n = n)
# Now:
      wz <- extra$wz

      temp1 <- ResSS.vgam(y-mu, wz = wz, M = M)
# Each row of wz is the same (or should be!!)
      onewz <- if (length(extra$invSigma.mat))
                 extra$invSigma.mat else
                 (m2adefault(wz[1, , drop = FALSE], M = M))[,, 1]  # M x M
#print("onewz")
#print( onewz )
#print("extra$invSigma.mat - onewz")
#print( extra$invSigma.mat - onewz )


# 20130131; done: use det() or determinant():
      logdet <- determinant(onewz)$modulus
#print("logdet")
#print( logdet )
#       logdet <- sum(log(eigen(onewz, symmetric = TRUE,
#                               only.values = TRUE)$values))
#print("logdet2")
#print( logdet )
      logretval <- -0.5 * temp1 + 0.5 * n * logdet -
                   n * (M / 2) * log(2*pi)
#     logretval <- -(ncol(onewz) * log(2 * pi) + logdet + distval)/2
      logretval
    }
  }

  ret.ff
}




# 20130125; Below here is old stuff... i will leave this alone
# --------------------------------------------------------------------
# 20110407; Below here is old stuff... i will leave this alone
# --------------------------------------------------------------------
# --------------------------------------------------------------------
# --------------------------------------------------------------------



# Sur <- function...






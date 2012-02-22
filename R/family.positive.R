# These functions are
# Copyright (C) 1998-2012 T.W. Yee, University of Auckland.
# All rights reserved.







rhuggins91 =
  function(n, nTimePts = 5, pvars = length(xcoeff),
           xcoeff = c(-2, 1, 2),
           capeffect = -1,
           double.ch = FALSE,
           link = "logit", earg = list()
           ) {


  use.n <- if ((length.n <- length(n)) > 1) length.n else
           if (!is.Numeric(n, integer.valued = TRUE,
                           allowable.length = 1, positive = TRUE))
               stop("bad input for argument 'n'") else n
  orig.n <- use.n
  use.n <- 1.50 * use.n + 100  # Bigger due to rejections

  if (pvars == 0)
    stop("argument 'pvars' must be at least one")
  if (pvars > length(xcoeff))
    stop("argument 'pvars' is too high")
  
  if (mode(link) != "character" && mode(link) != "name")
    link = as.character(substitute(link))
  if (!is.list(earg)) earg = list()


  Ymatrix = matrix(0, use.n, nTimePts, dimnames =
                   list(as.character(1:use.n),
                        paste("y", 1:nTimePts, sep = "")))

  CHmatrix = matrix(0, use.n, nTimePts, dimnames =
                    list(as.character(1:use.n),
                         paste("ch", 0:(nTimePts-1), sep = "")))

  Xmatrix = cbind(x1 = rep(1.0, len = use.n))
  if (pvars > 1)
    Xmatrix = cbind(Xmatrix,
                    matrix(runif(n = use.n * (pvars-1)), use.n, pvars - 1,
                           dimnames = list(as.character(1:use.n),
                                           paste("x", 2:pvars, sep = ""))))


  linpred.baseline = xcoeff[1]
  if (pvars > 1)
    linpred.baseline = linpred.baseline +
                       Xmatrix[, 2:pvars, drop = FALSE] %*% xcoeff[2:pvars]
  sumrowy = rep(0, length = use.n)
  for (jlocal in 1:nTimePts) {

    CHmatrix[, jlocal] = as.numeric(sumrowy > 0) *
                         (1 + double.ch)

    linpred = linpred.baseline + (CHmatrix[, jlocal] >  0) * capeffect

    Ymatrix[, jlocal] = rbinom(use.n, size = 1,
             prob = eta2theta(linpred, link = link, earg = earg))

    sumrowy = sumrowy + Ymatrix[, jlocal]
  }


  # Strip off rows where the animals were never caught
  # Bug: problem if all values of sumrowy are zero.
  index0 = (sumrowy == 0)
  if (all(!index0))
    stop("bug in this code: cannot handle no animals being caught")
  Ymatrix = Ymatrix[!index0, , drop = FALSE]
  Xmatrix = Xmatrix[!index0, , drop = FALSE]
  CHmatrix = CHmatrix[!index0, , drop = FALSE]

  # Bug: problem if all values of sumrowy are zero:
  zCHmatrix = matrix(0, nrow(CHmatrix), ncol(CHmatrix),
                     dimnames = list(as.character(1:nrow(CHmatrix)),
                     paste("zch", 0:(ncol(CHmatrix)-1), sep = "")))


  ans = data.frame(Ymatrix, Xmatrix, CHmatrix, zCHmatrix,
                   Chistory = rep(0, length = nrow(Ymatrix)))


  ans = if (nrow(ans) >= orig.n) ans[1:orig.n, ] else {
        rbind(ans,
              Recall(n = orig.n - nrow(ans),
                     nTimePts = nTimePts, pvars = pvars,
                     xcoeff = xcoeff,
                     capeffect = capeffect,
                     link = link, earg = earg))
        }

  rownames(ans) = as.character(1:orig.n)

  attr(ans, "pvars") = pvars
  attr(ans, "nTimePts") = nTimePts
  attr(ans, "capeffect") = capeffect

  ans
}





  

dhuggins91 = function(x, prob, prob0 = prob, log = FALSE) {


  x     = as.matrix(x)
  prob  = as.matrix(prob)
  prob0 = as.matrix(prob0)
  log.arg = log
  rm(log)



  logAA0 = rowSums(log1p(-prob0))
  AA0 = exp(logAA0)

  ell1 = rowSums(x * log(prob) + (1 - x) * log1p(-prob)) - log1p(-AA0)
  if (log.arg) ell1 else exp(ell1)
}







 huggins91 = function(link = "logit", earg = list(),
                      parallel = TRUE,
                      iprob = NULL,
                      eim.not.oim = TRUE) {





  if (mode(link) != "character" && mode(link) != "name")
    link = as.character(substitute(link))
  if (!is.list(earg)) earg = list()

  if (length(iprob))
    if (!is.Numeric(iprob, positive = TRUE) ||
        max(iprob) >= 1)
      stop("argument 'iprob' should have values in (0,1)")

  if (!is.logical(eim.not.oim) ||
      length(eim.not.oim) != 1)
    stop("argument 'eim.not.oim' should be 'TRUE' or 'FALSE' only")


  new("vglmff",
  blurb = c("Huggins (1991) capture-recapture model\n\n",
            "Links:    ",
            namesof("prob1",   link, earg = earg, tag = FALSE), ", ",
            namesof("prob1.0", link, earg = earg, tag = FALSE), ", ",
            namesof("prob2",   link, earg = earg, tag = FALSE), ", ",
            namesof("prob2.0", link, earg = earg, tag = FALSE), ", ..., ",
            namesof("probT.0", link, earg = earg, tag = FALSE),
            "\n"),
  constraints = eval(substitute(expression({
    constraints <- cm.vgam(matrix(1, M, 1), x, .parallel, constraints,
                           intercept.apply = TRUE)
  }), list( .parallel = parallel ))),
  infos = eval(substitute(function(...) {
    list(Musual = 2,
         parallel = .parallel)
  }, list( .parallel = parallel ))),

  initialize = eval(substitute(expression({
    Musual = 2
    mustart.orig = mustart
    y = as.matrix(y)
    Mdiv2 = ncoly = ncol(y)
    M = Musual * ncoly

    w = matrix(w, n, ncoly)
    mustart = matrix(colSums(y) / colSums(w),
                    n, ncol(y), byrow = TRUE)
    mustart[mustart == 0] = 0.05
    mustart[mustart == 1] = 0.95

    if (ncoly == 1)
      stop("the response is univariate, therefore use posbinomial()")


    if (!all(y == 0 | y == 1))
      stop("response must contain 0s and 1s only")
    if (!all(w == 1))
      stop("argument 'weight' must contain 1s only")

    dn2 = if (is.matrix(y)) dimnames(y)[[2]] else NULL
    dn2 = if (length(dn2)) {
      paste("E[", dn2, "]", sep = "")
    } else {
      paste("prob", 1:Mdiv2, sep = "")
    }
    dn2 = c(dn2, paste(dn2, ".0", sep = ""))
    dn2 = dn2[interleave.VGAM(M, M = Musual)]
    predictors.names = namesof(dn2, .link, earg = .earg, short = TRUE)


    if (!length(etastart)) {
      mustart.use = if (length(mustart.orig)) {
        mustart.orig
      } else
      if (length( .iprob )) {
        matrix( .iprob, nrow(mustart), ncol(mustart), byrow = TRUE)
      } else {
        mustart
      }
      etastart = cbind(theta2eta(mustart.use, .link, earg = .earg ))
      etastart = kronecker(etastart, cbind(1, 1))
    }
    mustart = NULL
  }), list( .link = link, .earg = earg, .iprob = iprob ))),

  linkinv = eval(substitute(function(eta, extra = NULL) {
    Musual = 2
    Mdiv2  =  ncol(eta) / Musual
    index1 =  Musual * (1:Mdiv2) - 1
    index2 =  Musual * (1:Mdiv2) - 0

    probs.numer = eta2theta(eta[, index1], # + extra$moffset[, index1],
                            .link, earg = .earg )

    probs.denom = eta2theta(eta[, index1], .link, earg = .earg )

    logAA0 = rowSums(log1p(-probs.denom))


    AA0 = exp(logAA0)
    AAA = exp(log1p(-AA0))  # 1 - AA0
    probs.numer / AAA
  }, list( .link = link, .earg = earg ))),
  last = eval(substitute(expression({

    misc$link = rep( .link, length = M)
    names(misc$link) = dn2

    misc$earg = vector("list", M)
    names(misc$earg) = names(misc$link)
    for(ii in 1:M)
      misc$earg[[ii]] = .earg

    misc$expected = .eim.not.oim
    misc$mv       = TRUE
    misc$iprob    = .iprob
    misc$eim.not.oim = .eim.not.oim

    misc$parallel   = .parallel
  }), list( .link = link, .earg = earg,
            .parallel = parallel,
            .eim.not.oim = eim.not.oim, .iprob = iprob ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL) {

      ycounts = y
      Musual = 2
      Mdiv2  =  ncol(eta) / Musual
      index1 =  Musual * (1:Mdiv2) - 1
      index2 =  Musual * (1:Mdiv2) - 0

      probs.numer = eta2theta(eta[, index1], # + extra$moffset[, index1],
                              .link, earg = .earg )

      probs.denom = eta2theta(eta[, index1], .link, earg = .earg )


      if (residuals) stop("loglikelihood residuals ",
                          "not implemented yet") else {
        sum(dhuggins91(x = ycounts, # size = 1, # Bernoulli trials
                       prob  = probs.numer,
                       prob0 = probs.denom, # zz choose this??
                       log = TRUE))
      }
  }, list( .link = link, .earg = earg ))),
  vfamily = c("huggins91"),
  deriv = eval(substitute(expression({
    Musual = 2
    Mdiv2  =  ncol(eta) / Musual
    index1 =  Musual * (1:Mdiv2) - 1
    index2 =  Musual * (1:Mdiv2) - 0
    probs.numer = eta2theta(eta[, index1], .link, earg = .earg )


    probs.denom = eta2theta(eta[, index1], .link, earg = .earg )

    logAA0 = rowSums(log1p(-probs.denom))


    AA0 = exp(logAA0)
    AAA = exp(log1p(-AA0))  # 1 - AA0

    B_s = AA0 / (1 - probs.denom)
    B_st = array(0, c(n, Mdiv2, Mdiv2))
    for(slocal in 1:(Mdiv2-1))
      for(tlocal in (slocal+1):Mdiv2)
        B_st[, slocal, tlocal] =
        B_st[, tlocal, slocal] = B_s[, slocal] / (1 - probs.denom[, tlocal])


    Temp2 =     (1 - probs.numer)^2
    temp2 =     (1 - probs.denom)^2

    dprob1.deta1 = dtheta.deta(probs.numer, .link , earg = .earg ) # trivial
    dprob1.deta2 = dtheta.deta(probs.numer, .link , earg = .earg ) # trivial
    dprob2.deta1 = dtheta.deta(probs.denom, .link , earg = .earg ) # trivial
    dprob2.deta2 = dtheta.deta(probs.denom, .link , earg = .earg ) # trivial

    dl.dprob1 =  y / probs.numer  - (1 - y) / (1 - probs.numer)
    dl.dprob2 =  -B_s / AAA
    dl.deta1  =  dl.dprob1 * dprob1.deta1
    dl.deta2  =  dl.dprob2 * dprob2.deta1
    dl.deta2  =  dl.dprob2 * dprob2.deta2 # zz

    deriv.ans = cbind(dl.deta1 + dl.deta2,
                      dl.deta1 + dl.deta2)
    deriv.ans = deriv.ans[, interleave.VGAM(M, M = Musual)]
    deriv.ans = deriv.ans / Musual   # Matches with CCCC

    deriv.ans
  }), list( .link = link, .earg = earg ))),

  weight = eval(substitute(expression({
    ed2l.dprob1.2 = 1 / (probs.numer * AAA) + 1 / Temp2 -
                    probs.numer / (AAA * Temp2) - (B_s / AAA)^2

    od2l.dprob1.2 =  y / probs.numer^2  + (1 - y) / (1 - probs.numer)^2 -
                     (B_s / AAA)^2




    d2prob1.deta1.2 = d2theta.deta2(probs.numer, .link , earg = .earg )
    d2prob1.deta2.2 = d2theta.deta2(probs.numer, .link , earg = .earg )
    d2prob1.deta12  = d2theta.deta2(probs.numer, .link , earg = .earg )
    d2prob2.deta1.2 = d2theta.deta2(probs.denom, .link , earg = .earg )
    d2prob2.deta12  = d2theta.deta2(probs.denom, .link , earg = .earg )


    wz = matrix(0, n, dimm(M))
    wz[, index1] <-
    wz[, index2] <-
    if ( .eim.not.oim ) {
       ed2l.dprob1.2 * (dprob1.deta1^2) # +
    } else {
      od2l.dprob1.2 * (dprob1.deta1^2) -
      (dl.dprob1 + dl.dprob2) * d2prob1.deta1.2
    }

    for(slocal in 1:(Mdiv2-1))
      for(tlocal in (slocal+1):Mdiv2)
        wz[, iam(Musual*slocal - 1,
                 Musual*tlocal - 1, M = M)] =
        wz[, iam(Musual*slocal    ,
                 Musual*tlocal    , M = M)] =
              dprob2.deta1[, slocal] *
              dprob2.deta1[, tlocal] *
            (B_st[, slocal, tlocal] +
             B_s [, slocal] *
             B_s [, tlocal] / AAA) / (-AAA)


    wz = wz / Musual   # Matches with CCCC


    wz
  }), list( .link = link, .earg = earg, .eim.not.oim = eim.not.oim ))))
}











dposnegbin = function(x, size, prob = NULL, munb = NULL, log = FALSE) {
  if (length(munb)) {
    if (length(prob))
      stop("'prob' and 'munb' both specified")
    prob <- size/(size + munb)
  }
  if (!is.logical(log.arg <- log)) stop("bad input for 'log'")
  rm(log)

  L = max(length(x), length(prob), length(size))
  x = rep(x, len = L); prob = rep(prob, len = L);
  size = rep(size, len = L);

  ans = dnbinom(x = x, size = size, prob = prob, log = log.arg)
  index0 = (x == 0)

  if (log.arg) {
    ans[ index0] = log(0.0)
    ans[!index0] = ans[!index0] - log1p(-dnbinom(x = 0 * x[!index0],
                   size = size[!index0], prob = prob[!index0]))
  } else {
    ans[ index0] = 0.0
    ans[!index0] = ans[!index0] / pnbinom(q = 0 * x[!index0],
                   size = size[!index0], prob = prob[!index0],
                   lower.tail = FALSE)
  }
  ans
}


pposnegbin = function(q, size, prob = NULL, munb = NULL) {
  if (length(munb)) {
    if (length(prob))
      stop("'prob' and 'munb' both specified")
    prob <- size/(size + munb)
  }
  L = max(length(q), length(prob), length(size))
  q = rep(q, len = L); prob = rep(prob, len = L); size = rep(size, len = L)

  ifelse(q < 1, 0, (pnbinom(q,   size = size, prob = prob) -
                    dnbinom(q*0, size = size, prob = prob))
                  / pnbinom(q*0, size = size, prob = prob,
                            lower.tail = FALSE))
}


qposnegbin = function(p, size, prob = NULL, munb = NULL) {
  if (length(munb)) {
    if (length(prob))
      stop("'prob' and 'munb' both specified")
    prob <- size / (size + munb)
  }
  if (!is.Numeric(p, positive = TRUE) ||
      any(p >= 1))
    stop("bad input for argument 'p'")
  qnbinom(p * pnbinom(q = p*0, size = size, prob = prob,
                      lower.tail = FALSE) +
              dnbinom(x = p*0, size = size, prob = prob),
          size = size, prob = prob)
}


rposnegbin = function(n, size, prob = NULL, munb = NULL) {


  use.n = if ((length.n <- length(n)) > 1) length.n else
          if (!is.Numeric(n, integer.valued = TRUE, allowable.length = 1, positive = TRUE))
              stop("bad input for argument 'n'") else n

  if (length(munb)) {
    if (length(prob))
      stop("'prob' and 'munb' both specified")
    prob <- size / (size + munb)
  }
  ans = rnbinom(use.n, size = size, prob = prob)
  index = (ans == 0)
  size = rep(size, len = use.n)
  prob = rep(prob, len = use.n)
  while(any(index, na.rm = TRUE)) {
    more = rnbinom(n = sum(index), size = size[index],
                   prob = prob[index])
    ans[index] = more
    index = (ans == 0)
  }
  ans
}








posnegbinomial.control <- function(save.weight = TRUE, ...)
{
  list(save.weight = save.weight)
}



 posnegbinomial = function(lmunb = "loge", lsize = "loge",
                           emunb = list(), esize = list(),
                           isize = NULL, zero = -2,
                           nsimEIM = 250,
                           shrinkage.init = 0.95, imethod = 1)
{

  if (!is.Numeric(imethod, allowable.length = 1, integer.valued = TRUE, positive = TRUE) ||
     imethod > 2)
    stop("argument 'imethod' must be 1 or 2")
  if (length(isize) && !is.Numeric(isize, positive = TRUE))
      stop("bad input for argument 'isize'")
  if (!is.Numeric(shrinkage.init, allowable.length = 1) ||
     shrinkage.init < 0 ||
     shrinkage.init > 1)
    stop("bad input for argument 'shrinkage.init'")

  if (mode(lmunb) != "character" && mode(lmunb) != "name")
      lmunb = as.character(substitute(lmunb))
  if (mode(lsize) != "character" && mode(lsize) != "name")
      lsize = as.character(substitute(lsize))

  if (!is.list(emunb)) emunb = list()
  if (!is.list(esize)) esize = list()


  if (!is.Numeric(nsimEIM, allowable.length = 1, positive = TRUE, integer.valued = TRUE))
    stop("argument 'nsimEIM' must be a positive integer")
  if (nsimEIM <= 30)
    warning("argument 'nsimEIM' should be greater than 30, say")


  new("vglmff",
  blurb = c("Positive-negative binomial distribution\n\n",
              "Links:    ",
              namesof("munb", lmunb, earg = emunb ), ", ",
              namesof("size", lsize, earg = esize ), "\n",
              "Mean:     munb / (1 - (size / (size + munb))^size)"),
  constraints = eval(substitute(expression({

    dotzero = .zero
    Musual = 2
    eval(negzero.expression)
  }), list( .zero = zero ))),
  initialize = eval(substitute(expression({
    Musual = 2

    if (any(y == 0))
      stop("there are zero values in the response")
    y = as.matrix(y) 
    M = 2 * ncol(y) 
    extra$NOS = NOS = ncoly = ncol(y)  # Number of species

    predictors.names = c(
      namesof(if (NOS == 1) "munb" else
              paste("munb", 1:NOS, sep = ""),
              .lmunb, earg = .emunb, tag = FALSE),
      namesof(if (NOS == 1) "size" else
              paste("size", 1:NOS, sep = ""),
              .lsize, earg = .esize, tag = FALSE))
    predictors.names = predictors.names[interleave.VGAM(M, M = Musual)]

    if (!length(etastart)) {
      mu.init = y
      for(iii in 1:ncol(y)) {
        use.this = if ( .imethod == 1) {
          weighted.mean(y[, iii], w)
        } else {
          median(y[,iii])
        }
        mu.init[, iii] = (1 - .sinit) * y[, iii] + .sinit * use.this
      }

      if ( is.Numeric( .isize )) {
        kmat0 = matrix( .isize , nrow = n, ncol = NOS, byrow = TRUE)
      } else {
        posnegbinomial.Loglikfun =
            function(kmat, y, x, w, extraargs) {
            munb = extraargs
              sum(w * dposnegbin(x = y, size = kmat, munb = munb,
                                 log = TRUE))
              }
            k.grid = 2^((-6):6)
            kmat0 = matrix(0, nrow = n, ncol = NOS)
            for(spp. in 1:NOS) {
              kmat0[, spp.] = getMaxMin(k.grid,
                                objfun = posnegbinomial.Loglikfun,
                                y = y[, spp.], x = x, w = w,
                                extraargs = mu.init[, spp.])
            }
      }
      p00 = (kmat0 / (kmat0 + mu.init))^kmat0
      etastart =
        cbind(
              theta2eta(mu.init * (1 - p00), .lmunb, earg = .emunb ),
              theta2eta(kmat0,               .lsize, earg = .esize ))
      etastart = etastart[,interleave.VGAM(M, M = Musual), drop = FALSE]
    }
  }), list( .lmunb = lmunb, .lsize = lsize, .isize = isize,
            .emunb = emunb, .esize = esize,
            .sinit = shrinkage.init,
            .imethod = imethod ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    Musual = 2
    NOS = ncol(eta) / Musual
    munb = eta2theta(eta[, Musual*(1:NOS)-1, drop = FALSE],
                     .lmunb, earg = .emunb )
    kmat = eta2theta(eta[, Musual*(1:NOS),   drop = FALSE],
                     .lsize, earg = .esize )
    po0 = (kmat / (kmat + munb))^kmat
    munb / (1 - po0)
  }, list( .lsize = lsize, .lmunb = lmunb,
           .esize = esize, .emunb = emunb ))),
  last = eval(substitute(expression({
    temp0303 = c(rep( .lmunb , length = NOS),
                 rep( .lsize , length = NOS))
    names(temp0303) =
       c(if (NOS == 1) "munb" else paste("munb", 1:NOS, sep = ""),
         if (NOS == 1) "size" else paste("size", 1:NOS, sep = ""))
    temp0303  = temp0303[interleave.VGAM(M, M = Musual)]
    misc$link = temp0303  # Already named

    misc$earg = vector("list", Musual*NOS)
    names(misc$earg) = names(misc$link)
    for(ii in 1:NOS) {
      misc$earg[[Musual*ii-1]] = .emunb
      misc$earg[[Musual*ii  ]] = .esize
    }

    misc$nsimEIM = .nsimEIM
    misc$imethod = .imethod
  }), list( .lmunb = lmunb, .lsize = lsize,
            .emunb = emunb, .esize = esize,
            .nsimEIM = nsimEIM, .imethod = imethod ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
    Musual = 2
    NOS = ncol(eta) / Musual
    munb = eta2theta(eta[, Musual*(1:NOS)-1, drop = FALSE],
                     .lmunb, earg = .emunb )
    kmat = eta2theta(eta[, Musual*(1:NOS)  , drop = FALSE],
                     .lsize, earg = .esize )
    if (residuals)
      stop("loglikelihood residuals not implemented yet") else {
      sum(w * dposnegbin(x = y, size = kmat, munb = munb, log = TRUE))
    }
  }, list( .lmunb = lmunb, .lsize = lsize,
           .emunb = emunb, .esize = esize ))),

  vfamily = c("posnegbinomial"),
  deriv = eval(substitute(expression({
    Musual = 2
    NOS = extra$NOS

    munb = eta2theta(eta[, Musual*(1:NOS)-1, drop = FALSE],
                     .lmunb , earg = .emunb )
    kmat = eta2theta(eta[, Musual*(1:NOS)  , drop = FALSE],
                     .lsize , earg = .esize )

    dmunb.deta = dtheta.deta(munb, .lmunb, earg = .emunb )
    dsize.deta = dtheta.deta(kmat, .lsize, earg = .esize )
    NOS = ncol(eta) / Musual


    tempk = kmat / (kmat + munb)
    tempm = munb / (kmat + munb)
    prob0  = tempk^kmat
    oneminusf0  = 1 - prob0
    df0.dmunb   = -tempk * prob0
    df0.dkmat   = prob0 * (tempm + log(tempk))
    df02.dmunb2 = prob0 * tempk / (kmat + munb) - tempk * df0.dmunb
    df02.dkmat2 = (prob0 / kmat) * tempm^2
    df02.dkmat.dmunb = prob0 * (-tempk) * (tempm + log(tempk)) -
                       tempm * prob0 / (kmat + munb)


    dl.dmunb = y / munb - (y + kmat) / (munb + kmat) +
               df0.dmunb / oneminusf0
    dl.dsize = digamma(y + kmat) - digamma(kmat) -
               (y + kmat)/(munb + kmat) + 1 + log(tempk) +
               df0.dkmat / oneminusf0

    myderiv = c(w) * cbind(dl.dmunb * dmunb.deta,
                           dl.dsize * dsize.deta)
    myderiv[, interleave.VGAM(M, M = Musual)]
  }), list( .lmunb = lmunb, .lsize = lsize,
            .emunb = emunb, .esize = esize ))),
  weight = eval(substitute(expression({
    run.varcov =
    wz = matrix(0.0, n, 4*NOS-1)




    if (FALSE) {
    usualmeanY =  munb
    meanY = usualmeanY / oneminusf0
    ed2l.dmu2 = meanY / munb^2 -
                (meanY + kmat) / (munb + kmat)^2 -
                df02.dmunb2 / oneminusf0 -
                (df0.dmunb / oneminusf0)^2
    }





    {
      ind2 = iam(NA, NA, M = Musual, both = TRUE, diag = TRUE)
      for(ii in 1:( .nsimEIM )) {
        ysim = rposnegbin(n = n*NOS, mu = c(munb), size = c(kmat))
        dim(ysim) = c(n, NOS)

        dl.dmunb = ysim / munb - (ysim + kmat) / (munb + kmat) +
                   df0.dmunb / oneminusf0
        dl.dsize = digamma(ysim + kmat) - digamma(kmat) -
                   (ysim + kmat) / (munb + kmat) + 1 + log(tempk) +
                   df0.dkmat / oneminusf0

        for(kk in 1:NOS) {
          temp2 = cbind(dl.dmunb[, kk],
                        dl.dsize[, kk]) *
                  cbind(dmunb.deta[, kk],
                        dsize.deta[, kk])
          small.varcov = temp2[, ind2$row.index] *
                         temp2[, ind2$col.index]

          run.varcov[, ((kk-1)*Musual+1):(kk*Musual)] =
          run.varcov[, ((kk-1)*Musual+1):(kk*Musual)] +
            c(small.varcov[, 1:Musual])
          run.varcov[, M + (kk-1)*Musual + 1] =
          run.varcov[, M + (kk-1)*Musual + 1] +
            c(small.varcov[, Musual + 1])
        }
      } # ii

      run.varcov = cbind(run.varcov / .nsimEIM )
      wz = if (intercept.only)
          matrix(colMeans(run.varcov),
                 n, ncol(run.varcov), byrow = TRUE) else run.varcov

    }


    c(w) * wz
  }), list( .nsimEIM = nsimEIM ))))
}





dposgeom = function(x, prob, log = FALSE) {
  dgeom(x - 1, prob = prob, log = log)
}


pposgeom = function(q, prob) {
  if (!is.Numeric(prob, positive = TRUE))
    stop("bad input for argument 'prob'")
  L = max(length(q), length(prob))
  if (length(q)    != L) q    = rep(q,    len = L);
  if (length(prob) != L) prob = rep(prob, len = L);
  ifelse(q < 1, 0, (pgeom(q, prob) - prob) / (1 - prob))
}


qposgeom = function(p, prob) {
  if (!is.Numeric(prob, positive = TRUE))
    stop("bad input for argument 'prob'")
  if (!is.Numeric(p, positive = TRUE) || any(p >= 1))
    stop("bad input for argument 'p'")
  qgeom(p * (1 - prob) + prob, prob = prob)
}


rposgeom = function(n, prob) {
  use.n = if ((length.n <- length(n)) > 1) length.n else
          if (!is.Numeric(n, integer.valued = TRUE, allowable.length = 1, positive = TRUE))
              stop("bad input for argument 'n'") else n

  ans = rgeom(use.n, prob = prob)
  prob = rep(prob, len = use.n)
  index = (ans == 0)
  while(any(index)) {
    more = rgeom(n = sum(index), prob[index])
    ans[index] = more
    index = (ans == 0)
  }
  ans
}







dpospois = function(x, lambda, log = FALSE) {
  if (!is.logical(log.arg <- log)) stop("bad input for 'log'")
  rm(log)

  if (!is.Numeric(lambda, positive = TRUE))
    stop("bad input for argument 'lambda'")
  L = max(length(x), length(lambda))
  x = rep(x, len = L); lambda = rep(lambda, len = L); 
  ans = if (log.arg) {
    ifelse(x == 0, log(0.0), dpois(x, lambda, log = TRUE) -
           log1p(-exp(-lambda)))
  } else {
    ifelse(x == 0, 0, -dpois(x, lambda) / expm1(-lambda))
  }
  ans
}


ppospois = function(q, lambda) {
  if (!is.Numeric(lambda, positive = TRUE))
    stop("bad input for argument 'lambda'")
  L = max(length(q), length(lambda))
  q = rep(q, len = L); lambda = rep(lambda, len = L); 
  ifelse(q < 1, 0, (ppois(q, lambda) - exp(-lambda)) / (-expm1(-lambda)))
}


qpospois = function(p, lambda) {
  if (!is.Numeric(lambda, positive = TRUE))
    stop("bad input for argument 'lambda'")
  if (!is.Numeric(p, positive = TRUE) || any(p >= 1))
    stop("bad input for argument 'p'")
  qpois(p * (-expm1(-lambda)) + exp(-lambda), lambda)
}


rpospois = function(n, lambda) {
  use.n = if ((length.n <- length(n)) > 1) length.n else
          if (!is.Numeric(n, integer.valued = TRUE, allowable.length = 1, positive = TRUE))
              stop("bad input for argument 'n'") else n

  if (any(lambda == 0))
    stop("no zero values allowed for argument 'lambda'")
  ans = rpois(use.n, lambda)
  lambda = rep(lambda, len = use.n)
  index = (ans == 0)
  while(any(index)) {
    more = rpois(n = sum(index), lambda[index])
    ans[index] = more
    index = (ans == 0)
  }
  ans
}




 pospoisson = function(link = "loge", earg = list(), expected = TRUE,
                       ilambda = NULL, imethod = 1)
{

  if (mode(link) != "character" && mode(link) != "name")
    link <- as.character(substitute(link))
  if (!is.list(earg)) earg <- list()

  if (!is.logical(expected) || length(expected) != 1)
    stop("bad input for argument 'expected'")
  if (length( ilambda) && !is.Numeric(ilambda, positive = TRUE))
    stop("bad input for argument 'ilambda'")

  if (!is.Numeric(imethod, allowable.length = 1, integer.valued = TRUE, positive = TRUE) ||
    imethod > 3)
    stop("argument 'imethod' must be 1 or 2 or 3")

  new("vglmff",
  blurb = c("Positive-Poisson distribution\n\n",
            "Links:    ",
            namesof("lambda", link, earg = earg, tag = FALSE)),
  infos = eval(substitute(function(...) {
    list(Musual = 1,
         link = .link,
         earg = .earg)
  }, list( .link = link, .earg = earg ))),

  initialize = eval(substitute(expression({
    y <- as.matrix(y)

    if (any(y < 1))
        stop("all y values must be in 1,2,3,...")
    if (any(y != round(y )))
        stop("the response must be integer-valued")

    predictors.names <-
      namesof(paste("lambda", if (ncol(y) > 1) 1:ncol(y) else "", sep = ""),
              .link, earg = .earg, tag = FALSE)

    if ( .imethod == 1) {
      lambda.init <- apply(y, 2, median) + 1/8
      lambda.init <- matrix(lambda.init, n, ncol(y), byrow = TRUE)
    } else if ( .imethod == 2) {
      lambda.init <- apply(y, 2, weighted.mean, w = w) + 1/8
      lambda.init <- matrix(lambda.init, n, ncol(y), byrow = TRUE)
    } else {
      lambda.init <- -y / expm1(-y)
    }
    if (length( .ilambda))
      lambda.init <- lambda.init * 0 + .ilambda

    if (!length(etastart))
        etastart <- theta2eta(lambda.init, .link, earg = .earg)
  }), list( .link = link, .earg = earg,
            .ilambda = ilambda, .imethod = imethod ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    lambda <- eta2theta(eta, .link, earg = .earg )
    -lambda / expm1(-lambda)
  }, list( .link = link, .earg = earg ))),
  last = eval(substitute(expression({
    misc$expected <- .expected

    misc$link <- rep( .link, len = M)
    names(misc$link) <- if (M == 1) "lambda" else
                       paste("lambda", 1:M, sep = "")

    misc$earg <- vector("list", M)
    names(misc$earg) <- names(misc$link)
    for(ii in 1:M)
      misc$earg[[ii]] <- .earg
  }), list( .link = link, .earg = earg, .expected = expected ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
    lambda <- eta2theta(eta, .link, earg = .earg ) 
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      sum(w * dpospois(x = y, lambda = lambda, log = TRUE))
    }
  }, list( .link = link, .earg = earg ))),
  vfamily = c("pospoisson"),
  deriv = eval(substitute(expression({
    lambda <- eta2theta(eta, .link, earg = .earg ) 
    temp6 <- expm1(lambda)
    dl.dlambda <- y / lambda - 1 - 1 / temp6
    dlambda.deta <- dtheta.deta(lambda, .link, earg = .earg )
    w * dl.dlambda * dlambda.deta
  }), list( .link = link, .earg = earg ))),
  weight = eval(substitute(expression({
    if ( .expected ) {
      ed2l.dlambda2 <- (temp6 + 1) * (1/lambda - 1/temp6) / temp6
      wz <- (dlambda.deta^2) * ed2l.dlambda2
    } else {
      d2l.dlambda2 <- y / lambda^2 - (temp6 + 1) / temp6^2
      d2lambda.deta2 <- d2theta.deta2(lambda, .link, earg = .earg)
      wz <- (dlambda.deta^2) * d2l.dlambda2 - dl.dlambda * d2lambda.deta2
    }
    c(w) * wz
  }), list( .link = link, .earg = earg, .expected = expected ))))
}








pposbinom = function(q, size, prob, lower.tail = TRUE, log.p = FALSE) {
  if (!is.Numeric(prob, positive = TRUE)) 
    stop("no zero or non-numeric values allowed for argument 'prob'")
  L = max(length(q), length(size), length(prob))
  q = rep(q, len = L); size = rep(size, len = L); prob = rep(prob, len = L)
  ifelse(q < 1, 0,
        (pbinom(q = q, size = size, prob = prob, lower.tail = lower.tail,
                log.p = log.p) - (1-prob)^size) / (1 - (1-prob)^size))
}


qposbinom = function(p, size, prob, lower.tail = TRUE, log.p = FALSE) {
  if (!is.Numeric(prob, positive = TRUE)) 
      stop("no zero or non-numeric values allowed for argument 'prob'")
  if (!is.Numeric(p, positive = TRUE) || any(p >= 1))
      stop("bad input for argument 'p'")
  qbinom(p = p * (1 - (1-prob)^size) + (1-prob)^size, size = size,
         prob = prob, lower.tail = lower.tail, log.p = log.p)
}


rposbinom = function(n, size, prob) {
  use.n = if ((length.n <- length(n)) > 1) length.n else
          if (!is.Numeric(n, integer.valued = TRUE, allowable.length = 1, positive = TRUE))
              stop("bad input for argument 'n'") else n

  if (any(prob == 0))
    stop("no zero values allowed for argument 'prob'")
  ans = rbinom(n = use.n, size = size, prob = prob)
  index = (ans == 0)
  size = rep(size, len=length(ans))
  prob = rep(prob, len=length(ans))
  while(any(index)) {
    more = rbinom(n = sum(index), size[index], prob = prob[index])
    ans[index] = more
    index = (ans == 0)
  }
  ans
}


dposbinom = function(x, size, prob, log = FALSE) {
  log.arg = log
  rm(log)
  L = max(length(x), length(size), length(prob))
  x = rep(x, len = L); size = rep(size, len = L);
                       prob = rep(prob, len = L);

  answer = NaN * x
  is0 <- (x == 0)
  ok2 <- (prob > 0) & (prob <= 1) &
         (size == round(size)) & (size > 0)
  answer = dbinom(x = x, size = size, prob = prob, log = TRUE) -
           log1p(-dbinom(x = 0*x, size = size, prob = prob))
  answer[!ok2] = NaN
  if (log.arg) {
    answer[is0 & ok2]  = log(0.0)
  } else {
    answer = exp(answer)
    answer[is0 & ok2] = 0.0
  }
  answer
}







 posbinomial = function(link = "logit", earg = list(),
                        mv = FALSE, parallel = FALSE, zero = NULL) {


  if (mode(link) != "character" && mode(link) != "name")
    link = as.character(substitute(link))
  if (!is.list(earg)) earg = list()


  new("vglmff",
  blurb = c("Positive-binomial distribution\n\n",
            "Links:    ",
            if (mv)
            c(namesof("prob1", link, earg = earg, tag = FALSE),
              ",...,",
              namesof("probM", link, earg = earg, tag = FALSE)) else
            namesof("prob", link, earg = earg, tag = FALSE),
            "\n"),
  constraints = eval(substitute(expression({
    constraints <- cm.vgam(matrix(1, M, 1), x, .parallel, constraints)
    constraints <- cm.zero.vgam(constraints, x, .zero, M)
  }), list( .parallel = parallel, .zero = zero ))),
  infos = eval(substitute(function(...) {
    list(Musual = 1,
         zero = .zero)
  }, list( .zero = zero ))),

  initialize = eval(substitute(expression({

    mustart.orig = mustart
    if ( .mv ) {
      y = as.matrix(y)
      M = ncoly = ncol(y)
      extra$orig.w = w
      w = as.matrix(w)  # Added 20110308
      mustart = matrix(colSums(y) / colSums(w),
                       n, ncol(y), byrow = TRUE)

    } else {
      eval(binomialff(link = .link, earg = .earg )@initialize)
    }


    if ( .mv ) {

      dn2 = if (is.matrix(y)) dimnames(y)[[2]] else NULL
      dn2 = if (length(dn2)) {
        paste("E[", dn2, "]", sep = "")
      } else {
        paste("prob", 1:M, sep = "")
      }
      predictors.names = namesof(if (M > 1) dn2 else
        "prob", .link, earg = .earg, short = TRUE)

      w = matrix(w, n, ncoly)
      y = y / w  # Now sample proportion
    } else {
      predictors.names =
        namesof("prob", .link, earg = .earg , tag = FALSE)
    }

    if (length(extra)) extra$w = w else extra = list(w = w)

    if (!length(etastart)) {
      mustart.use = if (length(mustart.orig)) mustart.orig else mustart
      etastart = cbind(theta2eta(mustart.use, .link, earg = .earg ))
    }
    mustart = NULL
  }), list( .link = link, .earg = earg, .mv = mv ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    w = extra$w
    mymu = eta2theta(eta, .link, earg = .earg )
    nvec = if ( .mv ) {
             w
           } else {
             if (is.numeric(extra$orig.w)) round(w / extra$orig.w) else
               round(w)
           }
    mymu / (1.0 - (1.0 - mymu)^(nvec))
  },
  list( .link = link, .earg = earg, .mv = mv ))),
  last = eval(substitute(expression({
    extra$w   = NULL   # Kill it off 


    misc$link = rep( .link, length = M)
    names(misc$link) = if (M > 1) dn2 else "prob"

    misc$earg = vector("list", M)
    names(misc$earg) = names(misc$link)
    for(ii in 1:M) misc$earg[[ii]] = .earg

    misc$expected = TRUE

    misc$mv   = .mv
    w = as.numeric(w)
  }), list( .link = link, .earg = earg, .mv = mv ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL) {

      ycounts = if ( .mv ) {
                  round(y * extra$orig.w)
                } else {
                  if (is.numeric(extra$orig.w)) y * w / extra$orig.w else
                  y * w # Convert proportions to counts
                }
      nvec = if ( .mv ) {
               w
             } else {
               if (is.numeric(extra$orig.w)) round(w / extra$orig.w) else
                 round(w)
             }
      use.orig.w = if (is.numeric(extra$orig.w)) extra$orig.w else 1
    mymu = eta2theta(eta, .link, earg = .earg )

    if (residuals) stop("loglikelihood residuals ",
                        "not implemented yet") else {
      sum(use.orig.w * dposbinom(x = ycounts, size = nvec,
                                 prob = mymu, log = TRUE))
    }
  }, list( .link = link, .earg = earg, .mv = mv ))),
  vfamily = c("posbinomial"),
  deriv = eval(substitute(expression({
    use.orig.w = if (is.numeric(extra$orig.w)) extra$orig.w else
                 rep(1, n)

    nvec = if ( .mv ) {
             w
           } else {
             if (is.numeric(extra$orig.w)) round(w / extra$orig.w) else
               round(w)
           }
    mymu = eta2theta(eta, .link, earg = .earg )
    dmu.deta = dtheta.deta(mymu, .link, earg = .earg )

    temp1 = 1 - (1 - mymu)^nvec
    temp2 =     (1 - mymu)^2
    temp3 =     (1 - mymu)^(nvec-2)

    dl.dmu = y / mymu - (1 - y) / (1 - mymu) -
             (1 - mymu) * temp3 / temp1

    c(w) * dl.dmu * dmu.deta
  }), list( .link = link, .earg = earg, .mv = mv ))),
  weight = eval(substitute(expression({
    ed2l.dmu2 = 1 / (mymu * temp1) + 1 / temp2 -
                mymu / (temp1 * temp2) -
                (nvec-1) * temp3 / temp1 -
                nvec * (temp2^(nvec-1)) / temp1^2
    wz = c(w) * ed2l.dmu2 * dmu.deta^2
    wz
  }), list( .link = link, .earg = earg, .mv = mv ))))
}







 rasch = function(lability = "identity",    eability = list(),
                  ldifficulty = "identity", edifficulty = list(),
                  iability = NULL,
                  idifficulty = NULL,
                  parallel = TRUE) {



  if (mode(labil) != "character" && mode(labil) != "name")
    labil = as.character(substitute(labil))

  if (!is.list(eabil)) eabil = list()
  if (!is.list(ediff)) ediff = list()

  if (length(iability))
    if (!is.Numeric(iability))
      stop("bad input in argument 'iability'")
  if (length(idifficulty))
    if (!is.Numeric(idifficulty))
      stop("bad input in argument 'idifficulty'")

  labil = lability
  eabil = eability
  ldiff = ldifficulty
  ediff = edifficulty


  new("vglmff",
  blurb = c("Rasch model\n\n",
            "Links:    ",
            namesof("ability",    labil, earg = eabil, tag = FALSE), ", ",
            namesof("difficulty", ldiff, earg = ediff, tag = FALSE),
            "\n"),

  initialize = eval(substitute(expression({
    mustart.orig = mustart
    y = as.matrix(y)
    extra$ncoly = ncoly = ncol(y)
    M = n + ncoly  # number of ability and number of item parameters


    mustart = matrix(apply(y, 2, weighted.mean, w = w),
                     n, ncol(y), byrow = TRUE)
    mustart[mustart == 0] = 0.05
    mustart[mustart == 1] = 0.95

    if (ncoly == 1)
      stop("the response is univariate, therefore use binomialff()")


    if (!all(y == 0 | y == 1))
      stop("response must contain 0s and 1s only")
    if (any(w <= 0))
      stop("argument 'weights' must contain positive values only")


    dn2 = if (is.matrix(y)) dimnames(y)[[2]] else NULL
    dn2 = as.character(1:ncoly)
    dn2 = as.character(1:nrow(y))
    dn2 = if (length(dn2)) {
      paste("ability", dn2, sep = "")
    } else {
      paste("zz", 1:Mdiv2, sep = "")
    }
    dn2 = c(dn2, paste("item", as.character(1:nrow(y)), sep = ""))
    predictors.names =
      namesof(dn2, .labil, earg = .eability, short = TRUE)




    if (!length(etastart)) {

      init.abil = runif(n) / (1 + colSums(y) - (1:n))
      init.diff = -logit(apply(y, 2, weighted.mean, w = w), inverse = TRUE)

      etastart =
        cbind(matrix(init.abil, n, n), #   byrow = TRUE ,
              matrix(init.diff, n, ncoly, byrow = TRUE))

    }
  }), list( .labil = labil, .eabil = eabil,
            .ldiff = ldiff, .ediff = ediff,
            .iability = iability,
            .idifficulty = idifficulty ))),

  linkinv = eval(substitute(function(eta, extra = NULL) {
    myprobs = eta2theta(eta, "logit", earg = list())
    myprobs
  }, list( .labil = labil, .eabil = eabil,
           .ldiff = ldiff, .ediff = ediff ))),
  last = eval(substitute(expression({

    misc$link = c(rep( .labil, length = n),
                  rep( .ldiff, length = ncoly))

    names(misc$link) = dn2

    misc$earg = vector("list", M)
    names(misc$earg) = names(misc$link)
    for(ii in 1:n)
      misc$earg[[ii]] = .eabil
    for(ii in 1:ncoly)
      misc$earg[[n + ii]] = .ediff

    misc$expected = TRUE
    misc$iability    = .iability
    misc$idifficulty = .idifficulty

  }), list( .labil = labil, .eabil = eabil,
            .ldiff = ldiff, .ediff = ediff,
            .iability = iability,
            .idifficulty = idifficulty ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
      if (residuals) stop("loglikelihood residuals ",
                          "not implemented yet") else {
        sum(w * (y * log(mu) + (1 - y) * log1p(-mu)))
      }
  }, list( .labil = labil, .eabil = eabil,
           .ldiff = ldiff, .ediff = ediff ))),
  vfamily = c("rasch"),
  deriv = eval(substitute(expression({
 print("head(mu)")
 print( head(mu) )
    dabil.deta = 1
    ddiff.deta = 1

    dl.dabil =   matrix(colSums(y - mu), n, n)
    dl.ddiff =  -cbind(y - mu)

    deriv.ans = cbind(dl.dabil * dabil.deta,
                      dl.ddiff * ddiff.deta)
 print("head(deriv.ans)")
 print( head(deriv.ans) )

    deriv.ans
  }), list( .labil = labil, .eabil = eabil,
            .ldiff = ldiff, .ediff = ediff ))),

  weight = eval(substitute(expression({

    wz = matrix(0, n, dimm(M))
    wz[, 1:M] = sqrt( .Machine$double.eps )


    tmp1 = colSums(mu * (1 - mu))
    for (ii in 1:n)
      wz[ii, ii] = tmp1[ii]



    wz[, n + (1:ncoly)] = mu * (1 - mu)


    for (ii in 1:n)
      for (jay in 1:ncoly)
        wz[ii, iam(ii, jay, M = M)] = -mu[ii, jay] * (1 - mu[ii, jay])

 print("head(wz)")
 print( head(wz) )

    wz = wz * w
    wz
  }), list( .labil = labil, .eabil = eabil ))))
}








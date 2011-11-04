# These functions are
# Copyright (C) 1998-2011 T.W. Yee, University of Auckland.
# All rights reserved.













 abbott = function(link0 = "logit", earg0 = list(),
                   link1 = "logit", earg1 = list(),
                   iprob0 = NULL, iprob1 = NULL,
                   fitted.type = c("observed", "treatment", "control"),
                   mux.offdiagonal = 0.98,
                   zero = 1) {


  fitted.type <- match.arg(fitted.type,
                           c("observed", "treatment", "control"),
                           several.ok =TRUE)


  if (mode(link0) !=  "character" && mode(link0) !=  "name")
    link0 <- as.character(substitute(link0))
  if (!is.list(earg0)) earg0 = list()

  if (mode(link1) !=  "character" && mode(link1) !=  "name")
    link1 <- as.character(substitute(link1))
  if (!is.list(earg1)) earg1 = list()

  if (!is.Numeric(mux.offdiagonal, allow = 1) ||
      mux.offdiagonal >= 1 ||
      mux.offdiagonal < 0)
    stop("argument 'mux.offdiagonal' must be in the interval [0, 1)")


  new("vglmff",
  blurb = c("Abbott's model for binary responses\n",
            "mu = prob0 + (1 - prob0) * prob1\n",
            "where 'prob0' is the 'control' mortality and\n",
            "'prob1' is the 'treatment' mortality and\n",
            "'mu' is the 'observed' mortality\n\n",
            "Links: ",
            namesof("prob0", link0, earg = earg0), ",  ",
            namesof("prob1", link1, earg = earg1)),
  constraints = eval(substitute(expression({
      constraints <- cm.zero.vgam(constraints, x, .zero, M)
  }), list( .zero = zero # ,
           ))),
  initialize = eval(substitute(expression({

    eval(binomialff(link = .link0)@initialize) # w, y, mustart are assigned


    predictors.names <-
      c(namesof("prob0", .link0, earg = .earg0, short =  TRUE),
        namesof("prob1", .link1, earg = .earg1, short =  TRUE))


    if (is.null(etastart)) {
      prob0.init <- if (length( .iprob0 )) {
          rep( .iprob0, len = n)
        } else {
          mustart / 2
        }

      prob1.init <- if (length( .iprob1 )) {
          rep( .iprob1, len = n)
        } else {
          mustart / 2
        }


      mustart <- NULL


        etastart <-
         cbind(theta2eta(prob0.init, link = .link0 , earg = .earg0 ),
               theta2eta(prob1.init, link = .link1 , earg = .earg1 ))
    }
  }), list( .link0 = link0, .earg0 = earg0,
            .link1 = link1, .earg1 = earg1,
            .iprob0 = iprob0, .iprob1 = iprob1 ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    prob0 <- eta2theta(eta[, 1], .link0 , earg = .earg0 )
    prob1 <- eta2theta(eta[, 2], .link1 , earg = .earg1 )

    con.fv = prob0
    trt.fv = prob1
    obs.fv = prob0 + (1 - prob0) * prob1



    ans = cbind("observed"  = obs.fv,
                "treatment" = trt.fv,
                "control"   = con.fv)

                   
    ans[, .fitted.type , drop = FALSE]
  }, list( .link0 = link0, .earg0 = earg0,
           .link1 = link1, .earg1 = earg1,
           .fitted.type = fitted.type ))),

  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL) {


      prob0 <- eta2theta(eta[, 1], .link0, earg = .earg0 )
      prob1 <- eta2theta(eta[, 2], .link1, earg = .earg1 )
      mymu = prob0 + (1 - prob0) * prob1


      if (residuals) {
        w * (y / mymu - (1-y) / (1 - mymu))
      } else {
        ycounts = if (is.numeric(extra$orig.w)) y * w / extra$orig.w else
                    y * w # Convert proportions to counts
        nvec = if (is.numeric(extra$orig.w)) round(w / extra$orig.w) else
                    round(w)
        smallno = 1.0e6 * .Machine$double.eps
        smallno = sqrt(.Machine$double.eps)
        if (max(abs(ycounts - round(ycounts))) > smallno)
          warning("converting 'ycounts' to integer in @loglikelihood")
        ycounts = round(ycounts)
        sum((if (is.numeric(extra$orig.w)) extra$orig.w else 1) *
            dbinom(x = ycounts, size = nvec, prob = mymu, log = TRUE))
      }
  }, list( .link0 = link0, .earg0 = earg0,
           .link1 = link1, .earg1 = earg1 ))),


  last = eval(substitute(expression({
    misc$link <-    c(prob0 = .link0 , prob1 = .link1 )
    misc$earg <- list(prob0 = .earg0 , prob1 = .earg1 )
    misc$mux.offdiagonal = .mux.offdiagonal
    misc$fitted.type = .fitted.type
    misc$true.mu = ( .fitted.type == "observed")


  }), list( .link0 = link0, .earg0 = earg0,
            .link1 = link1, .earg1 = earg1,
            .mux.offdiagonal = mux.offdiagonal,
            .fitted.type = fitted.type
          ))),
  vfamily = c("abbott", "vquantal"),
  deriv = eval(substitute(expression({
    prob0 <- eta2theta(eta[, 1], .link0, earg = .earg0 )
    prob1 <- eta2theta(eta[, 2], .link1, earg = .earg1 )
    dprob0.deta <- dtheta.deta(prob0, .link0 , earg = .earg0 )
    dprob1.deta <- dtheta.deta(prob1, .link1 , earg = .earg1 )


    mymu = prob0 + (1 - prob0) * prob1


    dl.dmu <- y / mymu - (1 - y) / (1 - mymu)
    dmu.dprob0 <- 1 - prob1
    dmu.dprob1 <- 1 - prob0
    dl.dprob0 <- dl.dmu * dmu.dprob0 
    dl.dprob1 <- dl.dmu * dmu.dprob1 


    c(w) * cbind(dl.dprob0 * dprob0.deta,
                 dl.dprob1 * dprob1.deta)
  }), list( .link0 = link0, .earg0 = earg0,
            .link1 = link1, .earg1 = earg1 ))),
    weight = eval(substitute(expression({


    ed2l.dmu2 <- 1 / (mymu * (1-mymu))
    ed2l.dprob02     <- ed2l.dmu2 * dmu.dprob0^2
    ed2l.dprob12     <- ed2l.dmu2 * dmu.dprob1^2
    ed2l.dprob1prob2 <-             ( 1)  # seems sort of ok but slow cvgc
    ed2l.dprob1prob2 <-             ( 0)  # kill it
    ed2l.dprob1prob2 <- ed2l.dmu2 * ( 1)  # dont seem to work

    ed2l.dprob1prob2 <- ed2l.dmu2 * dmu.dprob1 * dmu.dprob0 *
                        .mux.offdiagonal

    od2l.dmu2 <- y / mymu^2 + (1 - y) / (1 - mymu)^2
    od2l.dprob02     <- od2l.dmu2 * dmu.dprob0^2
    od2l.dprob12     <- od2l.dmu2 * dmu.dprob1^2
    od2l.dprob1prob2 <- od2l.dmu2 * dmu.dprob1 * dmu.dprob0 + dl.dmu


    wz <- cbind(ed2l.dprob02 * dprob0.deta^2,
                ed2l.dprob12 * dprob1.deta^2,
                ed2l.dprob1prob2 * dprob1.deta * dprob0.deta)


 if (FALSE)
    wz <- cbind(od2l.dprob02 * dprob0.deta^2,
                od2l.dprob12 * dprob1.deta^2,
                od2l.dprob1prob2 * dprob1.deta * dprob0.deta)


    c(w) * wz
  }), list( .link0 = link0, .earg0 = earg0,
            .link1 = link1, .earg1 = earg1,
             .mux.offdiagonal = mux.offdiagonal ))))
}










if (FALSE)
 Abbott = function(lprob1 = "elogit",
                   eprob1 = list(min = 0, max = 1), # For now, that is
                   lprob0 = "logit", eprob0 = list(),
                   iprob0 = NULL, iprob1 = NULL,
                   nointercept = 2, # NULL,
                   zero = 1) {




 stop("does not work")



  if (mode(lprob1) !=  "character" && mode(lprob1) !=  "name")
    lprob1 <- as.character(substitute(lprob1))
  if (!is.list(eprob1)) eprob1 = list()

  if (mode(lprob0) !=  "character" && mode(lprob0) !=  "name")
    lprob0 <- as.character(substitute(lprob0))
  if (!is.list(eprob0)) eprob0 = list()


  new("vglmff",
  blurb = c("Abbott's model for binary response\n",
            "mu = prob0 + prob1\n",
            "where 'prob0' is the control mortality and\n",
            "'prob1' is the treatment mortality\n\n",
            "Links: ",
            namesof("prob0", lprob0, earg = eprob0), ",  ",
            namesof("prob1", lprob1, earg = eprob1)),
  constraints = eval(substitute(expression({
      constraints <- cm.zero.vgam(constraints, x, .zero, M)
      constraints = cm.nointercept.vgam(constraints, x, .nointercept, M)
  }), list( .zero = zero,
            .nointercept = nointercept ))),


  initialize = eval(substitute(expression({
 print("here1")

    eval(binomialff(link = .lprob1)@initialize) # w, y, mustart are assigned

 print("here2")
 print("summary(mustart)")
 print( summary(mustart) )

    predictors.names <-
      c(namesof("prob0", .lprob0, earg = .eprob0, short =  TRUE),
        namesof("prob1", .lprob1, earg = .eprob1, short =  TRUE))


    if (is.null(etastart)) {
      prob0.init <- if (length( .iprob0 )) {
          rep( .iprob0, len = n)
        } else {
          mustart / 2
        }

      prob1.init <- if (length( .iprob1 )) {
          rep( .iprob1, len = n)
        } else {
          mustart * 1 / 4
        }


      mustart <- NULL


 print("prob0.init ")
 print( sort(prob0.init) )
 print("prob1.init ")
 print( sort(prob1.init) )


        eprob1 = list(min = prob0.init, max = 1)
        etastart <-
         cbind(theta2eta(prob0.init, link = .lprob0 , earg = .eprob0 ),
               theta2eta(prob1.init, link = .lprob1 , earg =  eprob1 ))
 print("head(etastart)")
 print( head(etastart) )
    }
  }), list( .lprob1 = lprob1, .eprob1 = eprob1,
            .lprob0 = lprob0, .eprob0 = eprob0,
            .iprob0 = iprob0, .iprob1 = iprob1 ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    prob0 <- eta2theta(eta[, 1], .lprob0 , earg = .eprob0)

    eprob1 = list(min = prob0, max = 1)

    prob1 <- eta2theta(eta[, 2], .lprob1 , earg =  eprob1)
    prob0 + prob1
  }, list( .lprob1 = lprob1, .eprob1 = eprob1,
           .lprob0 = lprob0, .eprob0 = eprob0 ))),
  last = eval(substitute(expression({
    eprob1 = list(min = prob0, max = 1)
    misc$link <-    c(prob0 = .lprob0, prob1 = .lprob1)
    misc$earg <- list(prob0 = .eprob0, prob1 =  eprob1)

    misc$nointercept = .nointercept
  }), list( .lprob1 = lprob1, .eprob1 = eprob1,
            .lprob0 = lprob0, .eprob0 = eprob0,
            .nointercept = nointercept ))),
  vfamily = c("Abbott", "vquantal"),
  deriv = eval(substitute(expression({
    prob0 <- eta2theta(eta[,1], .lprob0, earg = .eprob0)

    eprob1 = list(min = prob0, max = 1)
    prob1 <- eta2theta(eta[,2], .lprob1, earg =  eprob1)
    dprob0.deta <- dtheta.deta(prob0, .lprob0 , earg = .eprob0 )
    dprob1.deta <- dtheta.deta(prob1, .lprob1 , earg =  eprob1 )

    dl.dmu <- y / mu - (1 - y) / (1 - mu)
    dmu.dprob0 <- 1 # - prob1
    dmu.dprob1 <- 1 # - prob0
    dl.dprob0 <- dl.dmu * dmu.dprob0 
    dl.dprob1 <- dl.dmu * dmu.dprob1 

    c(w) * cbind(dl.dmu * dmu.dprob0 * dprob0.deta,
                 dl.dmu * dmu.dprob1 * dprob1.deta)
  }), list( .lprob1 = lprob1, .eprob1 = eprob1,
            .lprob0 = lprob0, .eprob0 = eprob0 ))),
    weight = eval(substitute(expression({


    ed2l.dmu2 <- 1 / (mu * (1-mu))
    ed2l.dprob02 <- ed2l.dmu2 * dmu.dprob0^2
    ed2l.dprob12 <- ed2l.dmu2 * dmu.dprob1^2

    wz <- cbind(ed2l.dprob02 * dprob0.deta^2,
                ed2l.dprob12 * dprob1.deta^2)

 print("head(wz)")
 print( head(wz) )
    c(w) * wz
  }), list( .lprob1 = lprob1, .eprob1 = eprob1,
            .lprob0 = lprob0, .eprob0 = eprob0 ))))
}

















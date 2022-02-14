# These functions are
# Copyright (C) 1998-2022 T.W. Yee, University of Auckland.
# All rights reserved.






















amazon.col <- "#3b7a57"
avocado.col <-  "#568203"
indigo.col <- "#6D5ACF"  # Indigo http://hexcolor16.com/6d5acf
iris.col <- "#5a4fcf"
turquoise.col <- "#30d5c8"  # for truncation
dirt.col <- "#9b7653"
deer.col <- "#ba8759"
desire.col <- "#ea3c53"   # A dark red colour



peach.col <- "#ffe5b4"
azure.col <- "#007fff"  # Unneeded
artichoke.col <- "#8f9779"
asparagus.col <- "#87a96b"




























 pd.damlm <-
  function(ind.num, ind.den,
           prob.num, prob.den = prob.num,
           is.dipped = c(FALSE, FALSE)) {
  if (is.dipped[1] && !is.dipped[2])
    is.dipped <- c(FALSE, FALSE)

  if (!any(is.dipped)) {
    return(prob.num * ((ind.num == ind.den) - prob.den))
  }

  if (!is.dipped[1] && is.dipped[2])  # Negated.
    return( prob.num * prob.den)

  if (all(is.dipped))  # Negated.
    return(-((ind.num == ind.den) - prob.num) * prob.den)
  stop("am confused... should never reach here")
}  # pd.damlm




 if (FALSE)
 pd.damlm.old <-
  function(num.a = NULL, num.d = NULL,
           den.a = NULL, den.d = NULL,
           prob.num, prob.den = prob.num,
           Denom = NA, eta.d.max = 0) {

  if (length(num.d) && length(den.d)) {
    if (num.a == num.d)
      return((1 - prob.num) * (prob.num - exp(eta.d.max) / Denom)) else
      return(-prob.num * (prob.den - exp(eta.d.max) / Denom))
  }
  if (length(num.a) && length(den.d)) {
    return(-prob.num * (prob.den - exp(eta.d.max) / Denom))
  }
  if (length(num.a) && length(den.a)) {
    if (num.a == num.d) return(prob.num * (1 - prob.num))
  }
  return(-prob.num * prob.den)
}  # pd.damlm.old








 get.indices.gaitd <- function(ind, indeta) {
  ind.b <- indeta[ind, 'launch']
  ind.b <- c(na.omit(ind.b))
  ind.e <- indeta[ind, 'finish']
  ind.e <- c(na.omit(ind.e))
  ind.z <- NULL
  if (length(ind.e) > 0)
    for (jay in seq(length(ind.e)))
      ind.z <- c(ind.z, seq(from = ind.b[jay], to = ind.e[jay]))
  ind.z
}  # get.ind.gaitd









 meangaitd <-
   function(theta.p,
            fam = c("pois", "log", "zeta"),  # "genpois0",
            a.mix = NULL, i.mix = NULL, d.mix = NULL,  # 
            a.mlm = NULL, i.mlm = NULL, d.mlm = NULL,  # 
            truncate = NULL, max.support = Inf,
            pobs.mix = 0,  # scalar
            pobs.mlm = 0,  # vector of length a.mlm
            pstr.mix = 0,  # scalar
            pstr.mlm = 0,  # vector of length a.mlm
            pdip.mix = 0,  # scalar
            pdip.mlm = 0,  # vector of length d.mlm
            byrow.aid = FALSE,  # Applies to 'pobs.mlm' & 'pstr.mlm'
            theta.a = theta.p,  # scalar, else a 2-vector
            theta.i = theta.p,  # scalar, else a 2-vector
            theta.d = theta.p,  # scalar, else a 2-vector
            ...
           ) {  # ... ignored currently.




  fam.choices <- c("pois", "log", "zeta")
  fam <- match.arg(fam[1], fam.choices)[1]
  baseparams.argnames <-
    switch(fam, "pois" = "lambda", "log" = "shape", "zeta" = "shape")


  gaitd.errorcheck(a.mix, a.mlm, i.mix, i.mlm,
                   d.mix, d.mlm, truncate, max.support)

  MM <- switch(fam, "pois" = 1, "log" = 1, "zeta" = 1)
  if (MM != 1 && MM !=  2)
    stop("can only handle 1 or 2 parameters")
  
  if (!length(a.mix)) pobs.mix <- 0  # Make sure for all
  if (!length(a.mlm)) pobs.mlm <- 0
  if (!length(i.mix)) pstr.mix <- 0
  if (!length(i.mlm)) pstr.mlm <- 0
  if (!length(d.mix)) pdip.mix <- 0
  if (!length(d.mlm)) pdip.mlm <- 0

  if (length(pobs.mix) != 1) stop("bad input for argument 'pobs.mix'")
  if (length(pstr.mix) != 1) stop("bad input for argument 'pstr.mix'")
  if (length(pdip.mix) != 1) stop("bad input for argument 'pdip.mix'")
  if (length(a.mlm) && length(pobs.mlm) > length(a.mlm))
    warning("bad input for argument 'pobs.mlm'?")
  if (length(i.mlm) && length(pstr.mlm) > length(i.mlm))
    warning("bad input for argument 'pstr.mlm'?")
  if (length(d.mlm) && length(pdip.mlm) > length(d.mlm))
    warning("bad input for argument 'pdip.mlm'?")

  if (length(a.mlm))
    pobs.mlm <- matrix(pobs.mlm, 1,  # length(xx),
                       length(a.mlm), byrow = byrow.aid)
  if (length(i.mlm))
    pstr.mlm <- matrix(pstr.mlm, 1,  # length(xx),
                       length(i.mlm), byrow = byrow.aid)
  if (length(d.mlm))
    pdip.mlm <- matrix(pdip.mlm, 1,  # length(xx),
                       length(d.mlm), byrow = byrow.aid)

      
  alist <- list(  # x = xx,  # theta.p,
            a.mix = a.mix, a.mlm = a.mlm,
            i.mix = i.mix, i.mlm = i.mlm,
            d.mix = d.mix, d.mlm = d.mlm,
            truncate = truncate, max.support = max.support,
            pobs.mix = pobs.mix, pobs.mlm = pobs.mlm,
            pstr.mix = pstr.mix, pstr.mlm = pstr.mlm,
            pdip.mix = pdip.mix, pdip.mlm = pdip.mlm,
            byrow.aid = byrow.aid)


    alist[[paste0(baseparams.argnames[1], ".p")]] <- theta.p[1]
    alist[[paste0(baseparams.argnames[1], ".a")]] <- theta.a[1]
    alist[[paste0(baseparams.argnames[1], ".i")]] <- theta.i[1]
    alist[[paste0(baseparams.argnames[1], ".d")]] <- theta.d[1]
    if (MM == 2) {
      alist[[paste0(baseparams.argnames[2], ".p")]] <- theta.p[2]
      alist[[paste0(baseparams.argnames[2], ".a")]] <- theta.a[2]
      alist[[paste0(baseparams.argnames[2], ".i")]] <- theta.i[2]
      alist[[paste0(baseparams.argnames[2], ".d")]] <- theta.d[2]
    }

  mlist <- alist
  mlist$type.fitted <- "All"
  mlist$moments2 <- TRUE
  mom.fun <- paste0("moments.gaitdcombo.", fam)
  Bits <- do.call(mom.fun, mlist)


  if (length(c(a.mix, a.mlm, i.mix, i.mlm, d.mix, d.mlm))) {

    Denom.p <- as.vector(Bits[["cdf.max.s"]] - Bits[["SumT0.p"]] -
                         Bits[["SumA0.mix.p"]] - Bits[["SumA0.mlm.p"]])
    if (any(Denom.p <= 0))
      stop("0s found in the denominator (variable 'Denom.p')")
    Numer <- as.vector(1 -
               (if (length(a.mix)) pobs.mix else 0) -
               (if (length(i.mix)) pstr.mix else 0) +
               (if (length(d.mix)) pdip.mix else 0) -
               (if (length(a.mlm)) rowSums(rbind(pobs.mlm)) else 0) -
               (if (length(i.mlm)) rowSums(rbind(pstr.mlm)) else 0) +
               (if (length(d.mlm)) rowSums(rbind(pdip.mlm)) else 0))
    if (!all(is.finite(Numer)))
      warning("variable 'Numer' contains non-finite values")
    if (min(Numer, na.rm = TRUE) < 0)
      warning("variable 'Numer' has negative values")
  }  # Inflation or deflation

  c(Bits$mean)
}  # meangaitd











 Trunc <- function(Range, mux = 2, location = 0, omits = TRUE) {
  if (!is.finite(mux) || length(mux) != 1 || round(mux) != mux ||
      mux < 1)
    stop("argument 'mux' must be a positive integer")
  if (any(!is.finite(Range)) || length(Range) < 2 ||
      !all(round(Range) == Range))
    stop("bad input in argument 'Range'")
  if (length(Range) > 2)
    Range <- range(Range, na.rm = TRUE)  # For vglm() and na.omit()
  Min <- Range[1]
  Max <- Range[2]

  allx <- location + (mux * Min):(mux * Max)
  multiples <- location + mux * (Min:Max)
  if (omits) {
    ans <- setdiff(allx, multiples)
    if (length(ans)) ans else NULL
  } else {
    multiples
  }
}  # Trunc














 gaitdzeta <-
  function(a.mix = NULL, i.mix = NULL, 
           d.mix = NULL,
           a.mlm = NULL, i.mlm = NULL,  # Unstructured probs are
           d.mlm = NULL,                # contiguous
           truncate = NULL, max.support = Inf,
           zero = c("pobs", "pstr", "pdip"),  # Pruned, handles all 6
           eq.ap = TRUE, eq.ip = TRUE, eq.dp = TRUE,
           parallel.a = FALSE, parallel.i = FALSE, parallel.d = FALSE,
           lshape.p = "loglink",
           lshape.a = lshape.p,  # "logitlink", 20201117
           lshape.i = lshape.p,  # "logitlink", 20201117
           lshape.d = lshape.p,  # "logitlink", 20211011
           type.fitted = c("mean", "shapes",
                           "pobs.mlm", "pstr.mlm", "pdip.mlm",
                           "pobs.mix", "pstr.mix", "pdip.mix",
                           "Pobs.mix", "Pstr.mix", "Pdip.mix",
                           "nonspecial", "Numer", "Denom.p",
                           "sum.mlm.i", "sum.mix.i",
                           "sum.mlm.d", "sum.mix.d",
                           "ptrunc.p", "cdf.max.s"),
           gshape.p = -expm1(-ppoints(7)),
           gpstr.mix = ppoints(7) / 3,  # ppoints(9) / 2,
           gpstr.mlm = ppoints(7) / (3 + length(i.mlm)),
           imethod = 1,
           mux.init = c(0.75, 0.5, 0.75),   # Order is A, I, D.
           ishape.p = NULL, ishape.a = ishape.p,
           ishape.i = ishape.p, ishape.d = ishape.p,
           ipobs.mix = NULL, ipstr.mix = NULL,  # 0.25, 
           ipdip.mix = NULL,   # 0.01,   # Easy but inflexible 0.01
           ipobs.mlm = NULL, ipstr.mlm = NULL,  # 0.25, 
           ipdip.mlm = NULL,   # 0.01,   # NULL, Easy but inflexible
           byrow.aid = FALSE,
           ishrinkage = 0.95,
           probs.y = 0.35) {

  mux.init <- rep_len(mux.init, 3)
  if (length(a.mix) == 0) a.mix <- NULL
  if (length(i.mix) == 0) i.mix <- NULL
  if (length(d.mix) == 0) d.mix <- NULL
  if (length(a.mlm) == 0) a.mlm <- NULL
  if (length(i.mlm) == 0) i.mlm <- NULL
  if (length(d.mlm) == 0) d.mlm <- NULL
  if (length(truncate) == 0) truncate <- NULL



  lowsup <- 1
  gaitd.errorcheck(a.mix, a.mlm, i.mix, i.mlm,
                   d.mix, d.mlm, truncate, max.support,
                   min.support = lowsup)
  la.mix <- length(a.mix <- sort(a.mix))
  li.mix <- length(i.mix <- sort(i.mix))
  ld.mix <- length(d.mix <- sort(d.mix))
  la.mlm <- length(a.mlm)
  li.mlm <- length(i.mlm)
  ld.mlm <- length(d.mlm)
  ltruncat <- length(truncate <- sort(truncate))
  ltrunc.use <- ltruncat > 0 || !is.infinite(max.support) 


  if (is.finite(max.support))
    stop("argument 'max.support' must be 'Inf'.")


  lshape.p <- as.list(substitute(lshape.p))
  eshape.p <- link2list(lshape.p)
  lshape.p <- attr(eshape.p, "function.name")
  lshape.p.save <- lshape.p

  lpobs.mix <- "multilogitlink"  # \omega_p
  epobs.mix <- list()  # zz NULL for now 20200907 coz 'multilogitlink'
  eshape.a <- link2list(lshape.a)
  lshape.a <- attr(eshape.a, "function.name")

  lpstr.mix <- "multilogitlink"  # \phi_p
  epstr.mix <- list()  # zz NULL for now 20200907 coz 'multilogitlink'
  lpdip.mix <- "multilogitlink"  # zz unsure 20211002
  epdip.mix <- list()  # zz unsure 20211002
  eshape.i <- link2list(lshape.i)
  lshape.i <- attr(eshape.i, "function.name")
  eshape.d <- link2list(lshape.d)
  lshape.d <- attr(eshape.d, "function.name")
  lshape.p.save <- lshape.p
  gshape.p.save <- gshape.p


  if (is.vector(zero) && is.character(zero) && length(zero) == 3) {
    if (li.mix + li.mlm == 0)
      zero <- setdiff(zero, "pstr")
    if (la.mix + la.mlm == 0)
      zero <- setdiff(zero, "pobs")
    if (ld.mix + ld.mlm == 0)
      zero <- setdiff(zero, "pdip")
    if (length(zero) == 0)
      zero <- NULL  # Better than character(0)
  }


  lall.len <- la.mix + li.mix + ld.mix + la.mlm + li.mlm + ld.mlm
  if (lall.len + ltruncat == 0 && is.infinite(max.support))
    return(eval(substitute(
       zetaff(lshape = .lshape.p.save ,
              gshape = .gshape.p.save ,
              zero = NULL),
           list( .lshape.p.save = lshape.p.save,
                 .gshape.p.save = gshape.p.save ))))

  if (!is.logical(eq.ap) || length(eq.ap) != 1)
    stop("argument 'eq.ap' must be a single logical")
  if (!is.logical(eq.ip) || length(eq.ip) != 1)
    stop("argument 'eq.ip' must be a single logical")
  if (!is.logical(parallel.a) || length(parallel.a) != 1)
    stop("argument 'parallel.a' must be a single logical")
  if (!is.logical(parallel.i) || length(parallel.i) != 1)
    stop("argument 'parallel.i' must be a single logical")
  if (!is.logical(parallel.d) || length(parallel.d) != 1)
    stop("argument 'parallel.d' must be a single logical")


  if (FALSE) {  # Comment this out to allow default eq.ap = TRUE, etc.
  if (la.mix <= 1 && eq.ap)
    stop("<= one unstructured altered value (no 'shape.a')",
         ", so setting 'eq.ap = TRUE' is meaningless")
  if (li.mix <= 1 && eq.ip)
    stop("<= one unstructured inflated value (no 'shape.i')",
            ", so setting 'eq.ip = TRUE' is meaningless")
  if (ld.mix <= 1 && eq.dp)
    stop("<= one unstructured deflated value (no 'shape.d')",
         ", so setting 'eq.dp = TRUE' is meaningless")
  if (la.mlm <= 1 && parallel.a)  # Only \omega_1
    stop("<= one altered mixture probability, 'pobs", a.mlm,
            "', so setting 'parallel.a = TRUE' is meaningless")
  if (li.mlm <= 1 && parallel.i)  # Only \phi_1
    stop("<= one inflated mixture probability, 'pstr", i.mlm,
            "', so setting 'parallel.i = TRUE' is meaningless")
  if (ld.mlm <= 1 && parallel.d)  # Only \psi_1
    stop("<= one deflated mixture probability, 'pdip", d.mlm,
         "', so setting 'parallel.d = TRUE' is meaningless")
  }  # FALSE


  type.fitted.choices <-
            c("mean", "shapes",
              "pobs.mlm", "pstr.mlm", "pdip.mlm",
              "pobs.mix", "pstr.mix", "pdip.mix",
              "Pobs.mix", "Pstr.mix", "Pdip.mix",
              "nonspecial", "Numer", "Denom.p",
              "sum.mlm.i", "sum.mix.i",
              "sum.mlm.d", "sum.mix.d",
              "ptrunc.p", "cdf.max.s")
  type.fitted <- match.arg(type.fitted[1], type.fitted.choices)[1]

  tmp7a <- if (la.mlm) paste0("pobs.mlm", a.mlm) else NULL
  tmp7b <- if (li.mlm) paste0("pstr.mlm", i.mlm) else NULL
  tmp7c <- if (ld.mlm) paste0("pdip.mlm", d.mlm) else NULL
  tmp3 <- c(shape.p = lshape.p,
            pobs.mix = if (la.mix) "multilogitlink" else NULL,
            shape.a = if (la.mix > 1) lshape.a else NULL,
            pstr.mix = if (li.mix) "multilogitlink" else NULL,
            shape.i = if (li.mix > 1) lshape.i else NULL,
            pdip.mix = if (ld.mix) "multilogitlink" else NULL,
            shape.d = if (ld.mix > 1) lshape.d else NULL,
            if (la.mlm) rep("multilogitlink", la.mlm) else NULL,
            if (li.mlm) rep("multilogitlink", li.mlm) else NULL,
            if (ld.mlm) rep("multilogitlink", ld.mlm) else NULL)
  Ltmp3 <- length(tmp3) 
  if (la.mlm + li.mlm + ld.mlm)
    names(tmp3)[(Ltmp3 - la.mlm - li.mlm - ld.mlm + 1):Ltmp3] <-
      c(tmp7a, tmp7b, tmp7c)
  par1or2 <- 1  # 2
  tmp3.TF <- c(TRUE, la.mix > 0, la.mix > 1,
                     li.mix > 0, li.mix > 1,
                     ld.mix > 0, ld.mix > 1,
                     la.mlm > 0, li.mlm > 0, ld.mlm > 0)
  indeta.finish <- cumsum(c(par1or2, 1, par1or2,
                                     1, par1or2,
                                     1, par1or2,
                            la.mlm, li.mlm, ld.mlm,
                            ld.mlm + 1) * c(tmp3.TF, 1))
  indeta.launch <- c(1, 1 + head(indeta.finish, -1))

  indeta.launch <- head(indeta.launch, -1)
  indeta.finish <- head(indeta.finish, -1)
  indeta.launch[!tmp3.TF] <- NA  # Not to be accessed
  indeta.finish[!tmp3.TF] <- NA  # Not to be accessed
  indeta <- cbind(launch = indeta.launch,
                  finish = indeta.finish)
  rownames(indeta) <- c("shape.p",
                        "pobs.mix", "shape.a",
                        "pstr.mix", "shape.i",
                        "pdip.mix", "shape.d",
                        "pobs.mlm", "pstr.mlm", "pdip.mlm")
  M1 <- max(indeta, na.rm = TRUE)
  predictors.names <- tmp3  # Passed into @infos and @initialize.
      

  blurb1 <- "Z"   # zz1
  if (la.mlm + la.mix) blurb1 <- "Generally-altered z"
  if (li.mlm + li.mix) blurb1 <- "Generally-inflated z"
  if (ltrunc.use) blurb1 <- "Generally-truncated z"
  if ( (la.mlm + la.mix) &&  (li.mlm + li.mix) && !ltrunc.use)
    blurb1 <- "Generally-altered and -inflated z"
  if ( (la.mlm + la.mix) && !(li.mlm + li.mix) &&  ltrunc.use)
    blurb1 <- "Generally-altered and -truncated z"
  if (!(la.mlm + la.mix) &&  (li.mlm + li.mix) &&  ltrunc.use)
    blurb1 <- "Generally-inflated and -truncated z"
  if ( (la.mlm + la.mix) &&  (li.mlm + li.mix) &&  ltrunc.use)
    blurb1 <- "Generally-altered, -inflated and -truncated z"

  if (ld.mlm + ld.mix) blurb1 <-
    c(blurb1,
      if (la.mlm + la.mix + li.mlm + li.mix) "and " else "Generally",
      "-deflated ")



      
  new("vglmff",
  blurb = c(blurb1, "eta regression\n",
            "(GAITD-zeta(shape.p)-",
                   "zeta(shape.a)-MLM-",
                   "zeta(shape.i)-MLM-\n",
                   "zeta(shape.d)-MLM generally)\n\n",
            "Links: ",
            namesof("shape.p", lshape.p, earg = eshape.p,
                    tag = FALSE),
            if (la.mix > 0) c(", ", "multilogit(pobs.mix)"),
            if (la.mix > 1) c(", ",
            namesof("shape.a",  lshape.a, eshape.a, tag = FALSE)),
            if (la.mix && li.mix) ", \n       ",
            if (li.mix > 0) c(  if (la.mix) "" else ", ",
            "multilogit(pstr.mix)"),
            if (li.mix > 1) c(", ",
            namesof("shape.i",  lshape.i, eshape.i, tag = FALSE)),
            if (li.mix && ld.mix) ", \n       ",
            if (ld.mix > 0) c(  if (li.mix) "" else ", ",
            "multilogit(pdip.mix)"),
            if (ld.mix > 1) c(", ",
            namesof("shape.d",  lshape.d, eshape.d, tag = FALSE)),
            if (la.mlm) paste0(",\n",
              paste0("       multilogit(", tmp7a, collapse = "),\n"),
            ")") else NULL,
            if (li.mlm) paste0(",\n",
              paste0("       multilogit(", tmp7b, collapse = "),\n"),
              ")") else NULL,
            if (ld.mlm) paste0(",\n",
              paste0("       multilogit(", tmp7c, collapse = "),\n"),
            ")") else NULL),
  constraints = eval(substitute(expression({
    M1 <- max(extra$indeta, na.rm = TRUE)
    la.mix <- ( .la.mix )
    li.mix <- ( .li.mix )
    ld.mix <- ( .ld.mix )
    la.mlm <- ( .la.mlm )
    li.mlm <- ( .li.mlm )
    ld.mlm <- ( .ld.mlm )

    use.mat.mlm.a <- if (la.mlm) {
      if ( .parallel.a ) matrix(1, la.mlm, 1) else diag(la.mlm)
    } else {
      NULL
    }
    use.mat.mlm.i <- if (li.mlm) {
       if ( .parallel.i ) matrix(1, li.mlm, 1) else diag(li.mlm)
    } else {
      NULL
    }
    use.mat.mlm.d <- if (ld.mlm) {
       if ( .parallel.d ) matrix(1, ld.mlm, 1) else diag(ld.mlm)
    } else {
      NULL
    }

    if (la.mlm + li.mlm + ld.mlm == 0) {
      Use.mat <- use.mat.mlm <- cbind(M)  # shape.p only
    }
    if (la.mlm + li.mlm + ld.mlm) {
      nc1 <- if (length(use.mat.mlm.a)) ncol(use.mat.mlm.a) else 0
      nc2 <- if (length(use.mat.mlm.i)) ncol(use.mat.mlm.i) else 0
      nc3 <- if (length(use.mat.mlm.d)) ncol(use.mat.mlm.d) else 0
      use.mat.mlm <- cbind(1, matrix(0, 1, nc1 + nc2 + nc3))
      if (la.mlm)
        use.mat.mlm <- rbind(use.mat.mlm,
                             cbind(matrix(0, la.mlm, 1),
                                   use.mat.mlm.a,
                                   if (length(use.mat.mlm.i) == 0)
                                   NULL else matrix(0, la.mlm, nc2),
                                   if (length(use.mat.mlm.d) == 0)
                                   NULL else matrix(0, la.mlm, nc3)))
      if (li.mlm )
       use.mat.mlm <-
         rbind(use.mat.mlm,
               cbind(matrix(0, li.mlm, 1 + nc1),
                     use.mat.mlm.i,
                     matrix(0, li.mlm, nc3)))
      if (ld.mlm)
        use.mat.mlm <-
          rbind(use.mat.mlm,  # zz1 next line:
                cbind(matrix(0, ld.mlm, 1 + nc1 + nc2),
                      use.mat.mlm.d))
    }  # la.mlm + li.mlm






    
    tmp3.TF <- ( .tmp3.TF )   # Logical of length 10.


  use.mat.mix <- cm3gaitd( .eq.ap , .eq.ip , .eq.dp , npar = 1)



    tmp3.subset <- tmp3.TF[-(8:10)]
    use.mat.mix <- use.mat.mix[tmp3.subset, , drop = FALSE]
    notall0 <- function(x) !all(x == 0)
    use.mat.mix <- use.mat.mix[, apply(use.mat.mix, 2, notall0),
                               drop = FALSE]


    if (la.mix + li.mix + ld.mix > 0)
      Use.mat <- use.mat.mix





    if (la.mlm + li.mlm + ld.mlm > 0) {
      Use.mat <- rbind(use.mat.mix,
                       matrix(0, nrow(use.mat.mlm) - 1,  # bottom
                                 ncol(use.mat.mix)))
      Use.mat <- cbind(Use.mat,
                       matrix(0, nrow(Use.mat),  # RHS
                                 ncol(use.mat.mlm) - 1))
      Use.mat[row(Use.mat) > nrow(use.mat.mix) &
              col(Use.mat) > ncol(use.mat.mix)] <- use.mat.mlm[-1, -1]
    }  # la.mlm + li.mlm + ld.mlm > 0



    if (is.null(constraints)) {
      constraints <-
        cm.VGAM(Use.mat, x = x, apply.int = TRUE,  # FALSE
                bool = .eq.ap || .eq.ip || .eq.dp ||
                       .parallel.a || .parallel.i || .parallel.d ,
                constraints = constraints)  # FALSE
    }  # is.null(constraints)



    if (la.mix + li.mix + ld.mix + la.mlm + li.mlm + ld.mlm)
      constraints <-
        cm.zero.VGAM(constraints, x = x, .zero , M = M, M1 = M1,
                     predictors.names = paste0(predictors.names,
                                        names(predictors.names)))
  }),
  list( .zero = zero, .tmp3.TF = tmp3.TF,  # .tmp3 = tmp3,
        .eq.ap = eq.ap, .eq.ip = eq.ip, .eq.dp = eq.dp,
        .parallel.a = parallel.a, .parallel.i = parallel.i,
        .parallel.d = parallel.d,
        .la.mlm = la.mlm, .li.mlm = li.mlm, .ld.mlm = ld.mlm,
        .la.mix = la.mix, .li.mix = li.mix, .ld.mix = ld.mix ))),
  infos = eval(substitute(function(...) {
    list(M1 = .M1 ,
         Q1 = 1,
         dpqrfun = "gaitdzeta",
         link = .predictors.names ,  # ...strips... from above
         link1parameter = as.logical( .lall.len <= 2),  # <= 1 safer
         mixture.links  = any(c( .la.mlm , .li.mlm , .ld.mlm ,
                                 .la.mix , .li.mix ,
                                 .ld.mix ) > 1),  # FALSE if NULL
         a.mix = as.vector( .a.mix ),  # Handles NULL
         a.mlm = as.vector( .a.mlm ),
         i.mix = as.vector( .i.mix ),
         i.mlm = as.vector( .i.mlm ),
         d.mix = as.vector( .d.mix ),
         d.mlm = as.vector( .d.mlm ),
         truncate = as.vector( .truncate ),
         max.support = as.vector( .max.support ),
         Support  = c( .lowsup , Inf, 1),  # a(b)c format as a,c,b.
         expected = TRUE,
         multipleResponses = FALSE,  # zetaff can be called if TRUE
         parameters.names = names( .predictors.names ),
         parent.name = c("zetaff", "zeta"),
         type.fitted  = as.vector( .type.fitted ),
         type.fitted.choices = ( .type.fitted.choices ),
         baseparams.argnames  = "shape",
         MM1 = 1,  # One parameter for 1 response (shape). Needed.
         zero = .zero )
  }, list( .zero = zero, .lowsup = lowsup,
           .type.fitted = type.fitted,
           .type.fitted.choices = type.fitted.choices,
           .lshape.p = lshape.p, .eshape.p = eshape.p,
           .lshape.a = lshape.a, .eshape.a = eshape.a,
           .lshape.i = lshape.i, .eshape.i = eshape.i,
           .a.mlm = a.mlm, .i.mlm = i.mlm, .d.mlm = d.mlm,
           .a.mix = a.mix, .i.mix = i.mix, .d.mix = d.mix,
           .la.mlm = la.mlm, .li.mlm = li.mlm, .ld.mlm = ld.mlm,
           .la.mix = la.mix, .li.mix = li.mix, .ld.mix = ld.mix,
           .truncate = truncate, .max.support = max.support,
           .predictors.names = predictors.names,
           .M1 = M1, .lall.len = lall.len
         ))),
  initialize = eval(substitute(expression({
    extra$indeta <- ( .indeta )  # Avoids recomputing it several times
    la.mix <- length((a.mix <- as.vector( .a.mix )))
    li.mix <- length((i.mix <- as.vector( .i.mix )))
    ld.mix <- length((d.mix <- as.vector( .d.mix )))
    la.mlm <- length((a.mlm <- as.vector( .a.mlm )))
    li.mlm <- length((i.mlm <- as.vector( .i.mlm )))
    ld.mlm <- length((d.mlm <- as.vector( .d.mlm )))
    lall.len <- la.mix + li.mix + ld.mix +
                la.mlm + li.mlm + ld.mlm
    truncate <- as.vector( .truncate )
    ltruncat <- length(truncate)
    M1 <- max(extra$indeta, na.rm = TRUE)
    NOS <- NCOL(y)
    M <- NOS * M1
    tmp3.TF <- ( .tmp3.TF )

    temp5 <-
    w.y.check(w = w, y = y,
              Is.nonnegative.y = TRUE,
              Is.integer.y = TRUE,
              ncol.w.max = 1,  # Since max.support = 9 is possible
              ncol.y.max = 1,
              out.wy = TRUE, colsyperw = 1, maximize = TRUE)
    w <- temp5$w
    y <- temp5$y
    glist <- y.gaitcombo.check(y, truncate = truncate,
                               a.mlm = a.mlm, a.mix = a.mix,
                               i.mlm = i.mlm, i.mix = i.mix,
                               d.mlm = d.mlm, d.mix = d.mix,
                               max.support = .max.support ,
                               min.support = .min.support )
    extra$skip.mix.a <- glist$skip.mix.a
    extra$skip.mix.i <- glist$skip.mix.i
    extra$skip.mix.d <- glist$skip.mix.d
    extra$skip.mlm.a <- glist$skip.mlm.a
    extra$skip.mlm.i <- glist$skip.mlm.i
    extra$skip.mlm.d <- glist$skip.mlm.d


    extra$NOS <- NOS <- ncoly <- ncol(y)  # Number of species
    extra$type.fitted <- as.vector( .type.fitted )
    extra$mux.init <- as.vector( .mux.init )
    extra$colnames.y  <- colnames(y)
    extra$M1 <- M1
    extra$index.M <- iam(NA, NA, M, both = TRUE)  # Used in @weight



    predictors.names <- ( .predictors.names )  # Got it, named


    if (!length(etastart)) {
      shape.p.init <- if (length( .ishape.p )) .ishape.p else {
        zetaff.Loglikfun <- function(shapeval, y, x, w, extraargs) {
          sum(c(w) * dzeta(x = y, shape = shapeval, log = TRUE))
        }
        shape.p.grid <- ( .gshape.p )
        grid.search(shape.p.grid, objfun = zetaff.Loglikfun,
                    y = y, w = w)
      }
      shape.p.init <- rep(shape.p.init, length = n)

      shape.d.init <-
      shape.a.init <- shape.i.init <- shape.p.init  # Needed


      etastart <- matrix(nrow = n, ncol = M,
        theta2eta(shape.p.init, .lshape.p , earg = .eshape.p ))




      mux.more.a <- extra$mux.init[1]   # 0.75 Err to slightly smaller
      init.pobs.mix <- numeric(n)
      if (tmp3.TF[ 2]) {  # la.mix > 0
        init.pobs.mix <- if (length( .ipobs.mix )) {
        rep_len( .ipobs.mix , n)
      } else {
        is.a.mix1 <- rowSums(extra$skip.mix.a) > 0
        rep(mux.more.a * sum(w[is.a.mix1]) / sum(w), n)
      }  
      }  # la.mix > 0


      if (tmp3.TF[ 3]) {  # Assign coln 3; la.mix > 1
        shape.a.init <- if (length( .ishape.a ))
          rep_len( .ishape.a , n) else shape.p.init  # A vector
        etastart[, 3] <-
          theta2eta(shape.a.init, .lshape.a , earg = .eshape.a )
      }

      init.pstr.mix <- init.pdip.mix <- numeric(n)
      try.gridsearch.pstr.mix <- FALSE
      if (tmp3.TF[ 4]) {  # li.mix > 0
        init.pstr.mix <- if (length( .ipstr.mix )) {
          rep_len( .ipstr.mix , n)
        } else {
          try.gridsearch.pstr.mix <- TRUE
          numeric(n)  # Overwritten by gridsearch
        }
      }  # li.mix > 0


      if (tmp3.TF[ 5]) {  # li.mix > 1
        shape.i.init <- if (length( .ishape.i ))
          rep_len( .ishape.i , n) else shape.p.init  # A vector
        etastart[, (extra$indeta[5, 'launch'])] <-
          theta2eta(shape.i.init, .lshape.i , earg = .eshape.i )
      }  # li.mix > 1


      if (tmp3.TF[ 8]) {  #  la.mlm
        init.pobs.mlm <- if (length( .ipobs.mlm )) {
          matrix( .ipobs.mlm , n, la.mlm, byrow = .byrow.aid )
        } else {
          mux.more.a <- extra$mux.init[1]                            
          init.pobs.mlm <- colSums(c(w) * extra$skip.mlm.a) / colSums(w)
          init.pobs.mlm <- init.pobs.mlm * as.vector( mux.more.a )
          matrix(init.pobs.mlm, n, la.mlm, byrow = TRUE)
        }
      } else {
        init.pobs.mlm <- matrix(0, n, 1)
      }


      try.gridsearch.pstr.mlm <- FALSE
      if (tmp3.TF[ 9]) {  #  li.mlm
        try.gridsearch.pstr.mlm <- !(length( .ipstr.mlm ))
        init.pstr.mlm <- 0  # Might be overwritten by gridsearch

        if (length( .ipstr.mlm ))
          init.pstr.mlm <- as.vector( .ipstr.mlm )

        init.pstr.mlm <- matrix(init.pstr.mlm, n, li.mlm,
                                byrow = .byrow.aid )
      } else {
        init.pstr.mlm <- matrix(0, n, 1)
      }

        

      init.pdip.mlm <- matrix(0, n, 2)  # rowSums() needs > 1 colns.





      gaitdzeta.Loglikfun1.mix <-
        function(pstr.mix.val, y, x, w, extraargs) {
    sum(c(w) *
        dgaitdzeta(y, pstr.mix = pstr.mix.val,
                  pstr.mlm    = extraargs$pstr.mlm,  # Differs here
                  shape.p     = extraargs$shape.p,
                  shape.a     = extraargs$shape.a,
                  shape.i     = extraargs$shape.i,
                  shape.d     = extraargs$shape.d,
                  a.mix       = extraargs$a.mix,
                  a.mlm       = extraargs$a.mlm,
                  i.mix       = extraargs$i.mix,
                  i.mlm       = extraargs$i.mlm,
                  d.mix       = extraargs$d.mix,
                  d.mlm       = extraargs$d.mlm,
                  max.support = extraargs$max.support,
                  truncate    = extraargs$truncate,
                  pobs.mix    = extraargs$pobs.mix,
                  pobs.mlm    = extraargs$pobs.mlm,
                  pdip.mix    = extraargs$pdip.mix,
                  pdip.mlm    = extraargs$pdip.mlm, log = TRUE))
  }

 gaitdzeta.Loglikfun1.mlm <-
     function(pstr.mlm.val, y, x, w, extraargs) {
    sum(c(w) *
        dgaitdzeta(y, pstr.mlm  = pstr.mlm.val,
                   pstr.mix    = extraargs$pstr.mix,  # Differs here
                   shape.p     = extraargs$shape.p,
                   shape.a     = extraargs$shape.a,
                   shape.i     = extraargs$shape.i,
                   shape.d     = extraargs$shape.d,
                   a.mix       = extraargs$a.mix,
                   a.mlm       = extraargs$a.mlm,
                   i.mix       = extraargs$i.mix,
                   i.mlm       = extraargs$i.mlm,
                   d.mix       = extraargs$d.mix,
                   d.mlm       = extraargs$d.mlm,
                   max.support = extraargs$max.support,
                   truncate    = extraargs$truncate,
                   pobs.mix    = extraargs$pobs.mix,
                   pobs.mlm    = extraargs$pobs.mlm,
                   pdip.mix    = extraargs$pdip.mix,
                   pdip.mlm    = extraargs$pdip.mlm, log = TRUE))
  }

 gaitdzeta.Loglikfun2 <-
     function(pstr.mix.val, pstr.mlm.val, y, x, w, extraargs) {
    sum(c(w) *
        dgaitdzeta(y, pstr.mix = pstr.mix.val,
                   pstr.mlm    = pstr.mlm.val,
                   shape.p     = extraargs$shape.p,
                   shape.a     = extraargs$shape.a,
                   shape.i     = extraargs$shape.i,
                   shape.d     = extraargs$shape.d,
                   a.mix       = extraargs$a.mix,
                   a.mlm       = extraargs$a.mlm,
                   i.mix       = extraargs$i.mix,
                   i.mlm       = extraargs$i.mlm,
                   d.mix       = extraargs$d.mix,
                   d.mlm       = extraargs$d.mlm,
                   max.support = extraargs$max.support,
                   truncate    = extraargs$truncate,
                   pobs.mix    = extraargs$pobs.mix,
                   pobs.mlm    = extraargs$pobs.mlm,
                   pdip.mix    = extraargs$pdip.mix,
                   pdip.mlm    = extraargs$pdip.mlm, log = TRUE))
  }



      if (li.mix + li.mlm) {
        extraargs <- list(
              shape.p     = shape.p.init,
              shape.a     = shape.a.init,
              shape.i     = shape.i.init,
              shape.d     = shape.d.init,
              a.mix       = a.mix,
              a.mlm       = a.mlm,
              i.mix       = i.mix,
              i.mlm       = i.mlm,
              d.mix       = d.mix,
              d.mlm       = d.mlm,
              truncate    = truncate,
              max.support = as.vector( .max.support ),
              pobs.mix    = init.pobs.mix ,
              pobs.mlm    = init.pobs.mlm ,
              pdip.mix    = init.pdip.mix ,
              pdip.mlm    = init.pdip.mlm )
        pre.warn <- options()$warn
        options(warn = -1)  # Ignore warnings during gridsearch 

        try.this <-
          if (try.gridsearch.pstr.mix && try.gridsearch.pstr.mlm) {
            grid.search2( .gpstr.mix ,  .gpstr.mlm ,
                         objfun = gaitdzeta.Loglikfun2,
                         y = y, w = w, extraargs = extraargs,
                         ret.objfun = TRUE)
          } else if (try.gridsearch.pstr.mix) {
            extraargs$pstr.mlm <- init.pstr.mlm
            grid.search ( .gpstr.mix ,
                         objfun = gaitdzeta.Loglikfun1.mix,
                         y = y, w = w, extraargs = extraargs,
                         ret.objfun = TRUE)
          } else if (try.gridsearch.pstr.mlm) {
            extraargs$pstr.mix <- init.pstr.mix
            grid.search ( .gpstr.mlm ,
                         objfun = gaitdzeta.Loglikfun1.mlm,
                         y = y, w = w, extraargs = extraargs,
                         ret.objfun = TRUE)
          }

        options(warn = pre.warn)  # Restore warnings 
        if (any(is.na(try.this)))
          warning("gridsearch returned NAs. It's going to crash.",
                  immediate. = TRUE)
        if (try.gridsearch.pstr.mix && try.gridsearch.pstr.mlm) {
          init.pstr.mix <- rep_len(try.this["Value1"], n)
          init.pstr.mlm <- matrix(try.this["Value2"], n, li.mlm)
          if (any(is.na(try.this)))
            stop("Crashing. ",
                 "Try something like 'gpstr.mix = seq(5) / 100'",
                 " and/or 'gpstr.mlm = seq(5) / 100'.")
        } else if (try.gridsearch.pstr.mix) {
          init.pstr.mix <- rep_len(try.this["Value"], n)
          if (any(is.na(try.this)))
            stop("Crashing. ",
                 "Try something like 'gpstr.mix = seq(5) / 100'.")
        } else if (try.gridsearch.pstr.mlm) {
          init.pstr.mlm <- matrix(try.this["Value"], n, li.mlm)
          if (any(is.na(try.this)))
            stop("Crashing. ",
                 "Try something like 'gpstr.mlm = seq(5) / 100'.")
        }
      }  # la.mix + lnf.mix



      mux.more.d <- extra$mux.init[3]   
      if (ld.mix) {
        init.pdip.mix <- if (length( .ipdip.mix ))
          rep_len( .ipdip.mix, n) else {
          is.d.mix1 <- rowSums(extra$skip.mix.d) > 0
          rep(mux.more.d * sum(w[is.d.mix1]) / sum(w), n) 
        }
      }  # ld.mix

      if (ld.mlm) {
        init.pdip.mlm <- if (length( .ipdip.mlm ))
          matrix( .ipdip.mlm, n, ld.mlm, byrow = TRUE) else {
          is.d.mlm1 <- rowSums(extra$skip.mlm.d) > 0
          matrix(mux.more.d * (sum(w[is.d.mlm1]) / sum(w)) / ld.mlm,
                 n, ld.mlm)
        }
      }  # ld.mlm




        


      while (any((vecTF <- init.pobs.mix + init.pstr.mix +  # -
                           init.pdip.mix +
                           rowSums(init.pobs.mlm) +
                           rowSums(init.pstr.mlm) +  # -
                           rowSums(init.pdip.mlm) > 0.96875))) {
        init.pobs.mix[vecTF]   <- 0.875 * init.pobs.mix[vecTF]
        init.pstr.mix[vecTF]   <- 0.875 * init.pstr.mix[vecTF]
        init.pdip.mix[vecTF]   <- 0.875 * init.pdip.mix[vecTF]
        init.pobs.mlm[vecTF, ] <- 0.875 * init.pobs.mlm[vecTF, ]
        init.pstr.mlm[vecTF, ] <- 0.875 * init.pstr.mlm[vecTF, ]
        init.pdip.mlm[vecTF, ] <- 0.875 * init.pdip.mlm[vecTF, ]
      }  # while

      Numer.init1 <- 1 - rowSums(init.pobs.mlm) -
                         rowSums(init.pstr.mlm) -  # +
                         rowSums(init.pdip.mlm) -
                         init.pobs.mix - init.pstr.mix -  # +
                         init.pdip.mix  # Differs from 'Numer'.
        
      etastart.z <- if (lall.len == 0) NULL else {
        tmp.mat <- cbind(if (tmp3.TF[ 2]) init.pobs.mix else NULL,
                         if (tmp3.TF[ 4]) init.pstr.mix else NULL,
                         if (tmp3.TF[ 6]) init.pdip.mix else NULL,
                         if (tmp3.TF[ 8]) init.pobs.mlm else NULL,
                         if (tmp3.TF[ 9]) init.pstr.mlm else NULL,
                         if (tmp3.TF[10]) init.pdip.mlm else NULL,
                         Numer.init1)
        multilogitlink(tmp.mat)
      }  # etastart.z
      if (!is.matrix(etastart.z)) etastart.z <- cbind(etastart.z)

      nextone <- 1  # Might not be used actually
      if (tmp3.TF[ 2]) {
        etastart[, 2] <- etastart.z[, nextone]
        nextone <- nextone + 1
      }
      if (tmp3.TF[ 4]) {  # Coln 2 or 4
        etastart[, (extra$indeta[4, 'launch'])] <- etastart.z[, nextone]
        nextone <- nextone + 1
      }
      if (tmp3.TF[ 6]) {  # Coln 2 or 4 or 6
        etastart[, (extra$indeta[6, 'launch'])] <- etastart.z[, nextone]
        nextone <- nextone + 1
      }
      if (tmp3.TF[ 8]) {
        ind8 <- (extra$indeta[8, 'launch']):(extra$indeta[8, 'finish'])
        etastart[, ind8] <- etastart.z[, nextone:(nextone+la.mlm - 1)]
        nextone <- nextone + la.mlm
      }
      if (tmp3.TF[ 9]) {
        ind9 <- (extra$indeta[9, 'launch']):(extra$indeta[9, 'finish'])
        etastart[, ind9] <- etastart.z[, nextone:(nextone+li.mlm - 1)]
        nextone <- nextone + li.mlm
      }
      if (tmp3.TF[10]) {
        ind0 <- (extra$indeta[10, 'launch']):(extra$indeta[10, 'finish'])
        etastart[, ind0] <- etastart.z[, nextone:(nextone + ld.mlm - 1)]
        if (ncol(etastart.z) != nextone + ld.mlm - 1)
          stop("miscalculation")
      }
    }
  }), list(
    .lshape.p = lshape.p, .eshape.p = eshape.p,
    .lshape.a = lshape.a, .eshape.a = eshape.a,
    .lshape.i = lshape.i, .eshape.i = eshape.i,
    .lshape.d = lshape.d, .eshape.d = eshape.d,
    .ishape.p = ishape.p,
    .ishape.a = ishape.a,
    .ishape.i = ishape.i,
    .ishape.d = ishape.d,
    .lpstr.mix = lpstr.mix, .lpobs.mix = lpobs.mix,
    .epstr.mix = epstr.mix, .epobs.mix = epobs.mix,
    .lpdip.mix = lpdip.mix,
    .epdip.mix = epdip.mix,
    .ipstr.mix = ipstr.mix, .ipobs.mix = ipobs.mix,
    .ipstr.mlm = ipstr.mlm, .ipobs.mlm = ipobs.mlm,
    .ipdip.mix = ipdip.mix,
    .ipdip.mlm = ipdip.mlm,
    .byrow.aid = byrow.aid,
    .a.mix = a.mix, .i.mix = i.mix, .d.mix = d.mix,
    .a.mlm = a.mlm, .i.mlm = i.mlm, .d.mlm = d.mlm,
    .truncate = truncate, .max.support = max.support,
    .min.support = lowsup,
    .tmp3.TF = tmp3.TF,  # .tmp3 = tmp3,
    .predictors.names = predictors.names,
    .mux.init = mux.init,
    .gshape.p = gshape.p,
    .gpstr.mix = gpstr.mix,   # .gpdip.mix = gpdip.mix,
    .gpstr.mlm = gpstr.mlm,   # .gpdip.mlm = gpdip.mlm,
    .ishrinkage = ishrinkage, .probs.y = probs.y,
    .indeta = indeta,
    .imethod = imethod, .type.fitted = type.fitted ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    n.obs <- NROW(eta)
    type.fitted <-
      if (length(extra$type.fitted)) extra$type.fitted else {
        warning("cannot find 'type.fitted'. Returning the 'mean'.")
        "mean"
      }
    type.fitted <-
      match.arg(type.fitted[1],
                c("mean", "shapes",
                  "pobs.mlm", "pstr.mlm", "pdip.mlm",
                  "pobs.mix", "pstr.mix", "pdip.mix",
                  "Pobs.mix", "Pstr.mix", "Pdip.mix",
                  "nonspecial", "Numer", "Denom.p",
                  "sum.mlm.i", "sum.mix.i",
                  "sum.mlm.d", "sum.mix.d",
                  "ptrunc.p", "cdf.max.s"))[1]

    if (!is.matrix(eta)) eta <- as.matrix(eta)
    la.mix <- length((a.mix <- as.vector( .a.mix )))
    li.mix <- length((i.mix <- as.vector( .i.mix )))
    ld.mix <- length((d.mix <- as.vector( .d.mix )))
    la.mlm <- length((a.mlm <- as.vector( .a.mlm )))
    li.mlm <- length((i.mlm <- as.vector( .i.mlm )))
    ld.mlm <- length((d.mlm <- as.vector( .d.mlm )))
    truncate <- as.vector( .truncate )
    max.support <- as.vector( .max.support )
    morework <- type.fitted != "mean"  # For efficiency

    lall.len <- la.mix + li.mix + ld.mix + la.mlm + li.mlm + ld.mlm
    pobs.mix <- pstr.mix <- pdip.mix <- 0  # 4 rowSums()
    pobs.mlm <- pstr.mlm <- pdip.mlm <- 0  # matrix(0, NROW(eta), 1)
    shape.p <- cbind(eta2theta(eta[, 1], .lshape.p , .eshape.p ))
    ind.shape.z <- 1  # Points to shape.p only.
    shape.a <- shape.i <-
    shape.d <- shape.p  # Needed; and answer not corrupted
    tmp3.TF <- ( .tmp3.TF )   # Logical of length 10.

    if (any(tmp3.TF[c(3, 5, 7)])) {  # At least one shape.[aid]
      ind.shape.z <- extra$indeta[c(1, 3, 5, 7), 'launch']  # Vectors
      ind.shape.z <- c(na.omit(ind.shape.z))  # At least one value
      shape.a <- if (!tmp3.TF[ 3]) shape.p else
        eta2theta(eta[, extra$indeta[3, 1]], .lshape.a , .eshape.a )
      shape.i <- if (!tmp3.TF[ 5]) shape.p else
        eta2theta(eta[, extra$indeta[5, 1]], .lshape.i , .eshape.i )
      shape.d <- if (!tmp3.TF[ 7]) shape.p else
        eta2theta(eta[, extra$indeta[7, 1]], .lshape.d , .eshape.d )
    }  # la.mix + li.mix + ld.mix > 0






    if (lall.len) {  # An MLM was fitted
      allprobs <-
        multilogitlink(eta[, -ind.shape.z, drop = FALSE],
                       inverse = TRUE)  # rowSums == 1
      if (anyNA(allprobs))
        warning("there are NAs here in slot linkinv")
      if (min(allprobs) == 0 || max(allprobs) == 1)
        warning("fitted probabilities numerically 0 or 1 occurred")


      Nextone <- 0  # Might not be used actually
      if (tmp3.TF[ 2])
        pobs.mix <- allprobs[, (Nextone <- Nextone + 1)]
      if (tmp3.TF[ 4])
        pstr.mix <- allprobs[, (Nextone <- Nextone + 1)]
      if (tmp3.TF[ 6])
        pdip.mix <- allprobs[, (Nextone <- Nextone + 1)]

      if (tmp3.TF[ 8]) {
        ind8 <- (Nextone + 1):(Nextone + la.mlm)
        pobs.mlm <- allprobs[, ind8, drop = FALSE]
        dimnames(pobs.mlm) <- list(rownames(eta), as.character(a.mlm))
        Nextone <- Nextone + la.mlm
      }
      if (tmp3.TF[ 9]) {
        ind9 <- (Nextone + 1):(Nextone + li.mlm)
        pstr.mlm <- allprobs[, ind9, drop = FALSE]
        dimnames(pstr.mlm) <- list(rownames(eta), as.character(i.mlm))
        Nextone <- Nextone + li.mlm
      }
      if (tmp3.TF[10]) {
        ind10 <- (Nextone + 1):(Nextone + ld.mlm)
        pdip.mlm <- allprobs[, ind10, drop = FALSE]
        dimnames(pdip.mlm) <- list(rownames(eta), as.character(d.mlm))
        Nextone <- Nextone + ld.mlm  # Not needed
      }
    }  # lall.len

    M1 <- max(extra$indeta, na.rm = TRUE)
    NOS <- NCOL(eta) / M1

    Bits <- moments.gaitdcombo.zeta(shape.p,
              pobs.mix = pobs.mix, pobs.mlm = pobs.mlm,
              pstr.mix = pstr.mix, pstr.mlm = pstr.mlm,
              pdip.mix = pdip.mix, pdip.mlm = pdip.mlm,
              a.mix = a.mix, i.mix = i.mix, d.mix = d.mix,
              a.mlm = a.mlm, i.mlm = i.mlm, d.mlm = d.mlm,
              shape.a = shape.a, shape.i = shape.i,
              shape.d = shape.d,
              truncate = truncate, max.support = max.support)


  if (morework) {
    Denom.p <- c(Bits[["cdf.max.s"]]   - Bits[["SumT0.p"]] -
                 Bits[["SumA0.mix.p"]] - Bits[["SumA0.mlm.p"]])

    if (any(Denom.p == 0)) {
      smallval <- min(Denom.p[Denom.p > 0])
      Denom.p[Denom.p == 0] <- 1e-09  # smallval
      warning("0s found in variable 'Denom.p'. Trying to fix it.")
    }
    Numer <- c(1 - pobs.mix - pstr.mix + pdip.mix -
               (if (la.mlm) rowSums(pobs.mlm) else 0) -
               (if (li.mlm) rowSums(pstr.mlm) else 0) +
               (if (ld.mlm) rowSums(pdip.mlm) else 0))
    probns <- Numer * (1 -
      (c(Bits[["SumI0.mix.p"]] + Bits[["SumI0.mlm.p"]]) +
       c(Bits[["SumD0.mix.p"]] + Bits[["SumD0.mlm.p"]])) / Denom.p)
  }  # morework

    if (!la.mlm && type.fitted %in% c("pobs.mlm")) {
      warning("No altered MLM values; returning an NA")
      return(NA)
    }
    if (!li.mlm && type.fitted %in% c("sum.mlm.i", "pstr.mlm")) {
      warning("No inflated MLM values; returning an NA")
      return(NA)
    }
    if (!ld.mlm && type.fitted %in% c("sum.mlm.d", "pdip.mlm")) {
      warning("No deflated MLM values; returning an NA")
      return(NA)
    }
    if (!la.mix && type.fitted %in% c("Pobs.mix")) {
      warning("No altered mixture values; returning an NA")
      return(NA)
    }
    if (!li.mix && type.fitted %in% c("sum.mix.i", "Pstr.mix")) {
      warning("No inflated mixture values; returning an NA")
      return(NA)
    }
    if (!ld.mix && type.fitted %in% c("sum.mix.d", "Pdip.mix")) {
      warning("No deflated mixture values; returning an NA")
      return(NA)
    }

    if (la.mix && morework) {
      tmp13 <-  # dpois() does not retain the matrix format
        dzeta(matrix(a.mix,   n.obs, la.mix, byrow = TRUE),
              matrix(shape.a, n.obs, la.mix)) / (
        c(Bits[["SumA0.mix.a"]]))
      dim(tmp13) <- c(n.obs, la.mix)
      dimnames(tmp13) <- list(rownames(eta), as.character(a.mix))
      propn.mat.a <- tmp13
    }  # la.mix

    if (li.mix && morework) {
      tmp55 <-  # dpois() does not retain the matrix format
        dzeta(matrix(i.mix,   n.obs, li.mix, byrow = TRUE),
              matrix(shape.i, n.obs, li.mix)) / (
        c(Bits[["SumI0.mix.i"]]))
      dim(tmp55) <- c(n.obs, li.mix)
      dimnames(tmp55) <- list(rownames(eta), as.character(i.mix))
      propn.mat.i <- tmp55  # Correct dimension
    }  # li.mix


    if (ld.mix && morework) {
      tmp55 <-  # dpois() does not retain the matrix format
        dzeta(matrix(d.mix,   n.obs, ld.mix, byrow = TRUE),
              matrix(shape.d, n.obs, ld.mix)) / (
        c(Bits[["SumD0.mix.d"]]))
      dim(tmp55) <- c(n.obs, ld.mix)
      dimnames(tmp55) <- list(rownames(eta), as.character(d.mix))
      propn.mat.d <- tmp55  # Correct dimension
    }  # ld.mix

    ans <- switch(type.fitted,
      "mean"       = Bits[["mean"]],  # Unconditional mean
      "shapes"    = cbind(shape.p,
                           if (tmp3.TF[ 3]) shape.a else NULL,
                           if (tmp3.TF[ 5]) shape.i else NULL,
                           if (tmp3.TF[ 7]) shape.d else NULL),
      "pobs.mlm"   = pobs.mlm,  # aka omegamat, n x la.mlm
      "pstr.mlm"   = pstr.mlm,  # aka phimat,   n x li.mlm
      "pdip.mlm"   = pdip.mlm,  # aka psimat,   n x ld.mlm
      "pobs.mix"   = pobs.mix,  # n-vector
      "pstr.mix"   = pstr.mix,  # n-vector
      "pdip.mix"   = pdip.mix,  # n-vector
      "Pobs.mix"   = c(pobs.mix) * propn.mat.a,  # matrix
      "Pstr.mix"   = c(pstr.mix) * propn.mat.i,
      "Pdip.mix"   = c(pdip.mix) * propn.mat.d,
      "nonspecial" = probns,
      "Numer"      = Numer,
      "Denom.p"    = Denom.p,
      "sum.mlm.i"  =  pstr.mlm + Numer *
             dzeta(matrix(i.mlm,   n.obs, li.mlm, byrow = TRUE),
                   matrix(shape.p, n.obs, li.mlm)) / Denom.p,
      "sum.mlm.d"  = -pdip.mlm + Numer *
             dzeta(matrix(d.mlm,   n.obs, ld.mlm, byrow = TRUE),
                   matrix(shape.p, n.obs, ld.mlm)) / Denom.p,
      "sum.mix.i"  =  c(pstr.mix) * propn.mat.i + Numer *
             dzeta(matrix(i.mix,   n.obs, li.mix, byrow = TRUE),
                   matrix(shape.p, n.obs, li.mix)) / Denom.p,
      "sum.mix.d"  = -c(pdip.mix) * propn.mat.d + Numer *
             dzeta(matrix(d.mix,   n.obs, ld.mix, byrow = TRUE),
                   matrix(shape.p, n.obs, ld.mix)) / Denom.p,
      "ptrunc.p"   = Bits[["SumT0.p"]] + 1 - Bits[["cdf.max.s"]],
      "cdf.max.s"  = Bits[["cdf.max.s"]])  # Pr(y <= max.support)



    ynames.pobs.mlm <- as.character(a.mlm)  # Works with NULLs
    ynames.pstr.mlm <- as.character(i.mlm)  # Works with NULLs
    ynames.pdip.mlm <- as.character(d.mlm)  # Works with NULLs
    if (length(ans))
      label.cols.y(ans, NOS = NOS, colnames.y =
      switch(type.fitted,
             "shapes"   = c("shape.p", "shape.a",  # Some colns NA
                 "shape.i", "shape.d")[(tmp3.TF[c(1, 3, 5, 7)])],
             "Pobs.mix"  = as.character(a.mix),
             "sum.mix.i" = ,  #
             "Pstr.mix"  = as.character(i.mix),
             "sum.mix.d" = ,  #
             "Pdip.mix"  = as.character(d.mix),
             "pobs.mlm"  = ynames.pobs.mlm,
             "sum.mlm.i" = ,  #
             "pstr.mlm"  = ynames.pstr.mlm,
             "sum.mlm.d" = ,  #
             "pdip.mlm"  = ynames.pdip.mlm,
             extra$colnames.y)) else ans
  }, list(
    .lshape.p = lshape.p, .eshape.p = eshape.p,
    .lshape.a = lshape.a, .eshape.a = eshape.a,
    .lshape.i = lshape.i, .eshape.i = eshape.i,
    .lshape.d = lshape.d, .eshape.d = eshape.d,
    .lpstr.mix = lpstr.mix, .lpobs.mix = lpobs.mix,
    .lpdip.mix = lpdip.mix,
    .epstr.mix = epstr.mix, .epobs.mix = epobs.mix,
    .epdip.mix = epdip.mix,
    .tmp3.TF = tmp3.TF,
    .a.mix = a.mix, .i.mix = i.mix, .d.mix = d.mix,
    .a.mlm = a.mlm, .i.mlm = i.mlm, .d.mlm = d.mlm,
    .truncate = truncate, .max.support = max.support ))),
  last = eval(substitute(expression({
    pred.names <- c( .predictors.names )  # Save it
    link.names <- as.vector( .predictors.names )
    parameter.names <- names(pred.names)
    predictors.names <- NULL
    for (jay in seq(M))
      predictors.names <- c(predictors.names,
        namesof(parameter.names[jay], link.names[jay], tag = FALSE,
                earg = list()))  # This line isnt perfect; info is lost
    misc$predictors.names <- predictors.names  # Useful for coef()
    misc$link <- link.names  # 
    names(misc$link) <- parameter.names  # 


    misc$earg <- vector("list", M1)
    names(misc$earg) <- names(misc$link)
    misc$earg[[1]] <- ( .eshape.p )  # First one always there
    iptr <- 1
    if (tmp3.TF[ 2])
      misc$earg[[(iptr <- iptr + 1)]] <- list()  # multilogitlink
    if (tmp3.TF[ 3])
      misc$earg[[(iptr <- iptr + 1)]] <- ( .eshape.a )
    if (tmp3.TF[ 4])
      misc$earg[[(iptr <- iptr + 1)]] <- list()  # See below
    if (tmp3.TF[ 5])
      misc$earg[[(iptr <- iptr + 1)]] <- ( .eshape.i )
    if (tmp3.TF[ 6])
      misc$earg[[(iptr <- iptr + 1)]] <- list()  # See below
    if (tmp3.TF[ 7])
      misc$earg[[(iptr <- iptr + 1)]] <- ( .eshape.d )
    if (tmp3.TF[ 8]) {  # la.mlm
      for (ii in seq(la.mlm)) {
        misc$earg[[(iptr <- iptr + 1)]] <- list()
      }  # ii
    }  # la.mlm
    if (tmp3.TF[ 9]) {  # li.mlm
      for (ii in seq(li.mlm)) {
        misc$earg[[(iptr <- iptr + 1)]] <- list()
      }  # ii
    }  # li.mlm
    if (tmp3.TF[10]) { # ld.mlm
      for (ii in seq(ld.mlm)) {
        misc$earg[[(iptr <- iptr + 1)]] <- list()
      }  # ii
    }  # ld.mlm
  }), list(
    .lshape.p = lshape.p, .eshape.p = eshape.p,
    .lshape.a = lshape.a, .eshape.a = eshape.a,
    .lshape.i = lshape.i, .eshape.i = eshape.i,
    .lshape.d = lshape.d, .eshape.d = eshape.d,
    .lpstr.mix = lpstr.mix, .lpobs.mix = lpobs.mix,
    .lpdip.mix = lpdip.mix,
    .epstr.mix = epstr.mix, .epobs.mix = epobs.mix,
    .epdip.mix = epdip.mix,
    .predictors.names = predictors.names,
    .tmp3.TF = tmp3.TF,  # .tmp3 = tmp3,
    .a.mix = a.mix, .i.mix = i.mix, .d.mix = d.mix,
    .a.mlm = a.mlm, .i.mlm = i.mlm, .d.mlm = d.mlm,
    .truncate = truncate, .max.support = max.support ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL, summation = TRUE) {
    if (!is.matrix(eta)) eta <- as.matrix(eta)
    la.mix <- length((a.mix <- as.vector( .a.mix )))
    li.mix <- length((i.mix <- as.vector( .i.mix )))
    ld.mix <- length((d.mix <- as.vector( .d.mix )))
    la.mlm <- length((a.mlm <- as.vector( .a.mlm )))
    li.mlm <- length((i.mlm <- as.vector( .i.mlm )))
    ld.mlm <- length((d.mlm <- as.vector( .d.mlm )))
    truncate <- as.vector( .truncate )

    lall.len <- la.mix + li.mix + ld.mix + la.mlm + li.mlm + ld.mlm
    pobs.mix <- pstr.mix <- pdip.mix <- 0  # 4 rowSums()
    pobs.mlm <- pstr.mlm <- pdip.mlm <- 0  # matrix(0, NROW(eta), 1)
    shape.p <- cbind(eta2theta(eta[, 1], .lshape.p , .eshape.p ))
    ind.shape.z <- 1  # Points to shape.p only.
    shape.a <- shape.i <-
    shape.d <- shape.p  # Needed and doesnt corrupt the answer

    if (any(tmp3.TF[c(3, 5, 7)])) {  # At least one shape.[aid]
      ind.shape.z <- extra$indeta[c(1, 3, 5, 7), 'launch']  # Vectors
      ind.shape.z <- c(na.omit(ind.shape.z))  # At least one value

      shape.a <- if (!tmp3.TF[ 3]) shape.p else
        eta2theta(eta[, extra$indeta[3, 1]], .lshape.a , .eshape.a )
      shape.i <- if (!tmp3.TF[ 5]) shape.p else
        eta2theta(eta[, extra$indeta[5, 1]], .lshape.i , .eshape.i )
      shape.d <- if (!tmp3.TF[ 7]) shape.p else
        eta2theta(eta[, extra$indeta[7, 1]], .lshape.d , .eshape.d )
    }  # la.mix + li.mix + ld.mix > 0

    if (lall.len) {  # An MLM was fitted
      allprobs <-
        multilogitlink(eta[, -ind.shape.z, drop = FALSE],
                       refLevel = "(Last)",  # Make sure
                       inverse = TRUE)  # rowSums == 1
      if (anyNA(allprobs))
        warning("there are NAs here in slot linkinv")
      if (min(allprobs) == 0 || max(allprobs) == 1)
        warning("fitted probabilities numerically 0 or 1 occurred")

      Nextone <- 0  # Might not be used actually; 0, not 1
      if (tmp3.TF[ 2])
        pobs.mix <- allprobs[, (Nextone <- Nextone + 1)]
      if (tmp3.TF[ 4])
        pstr.mix <- allprobs[, (Nextone <- Nextone + 1)]
      if (tmp3.TF[ 6])
        pdip.mix <- allprobs[, (Nextone <- Nextone + 1)]

      if (tmp3.TF[ 8]) {
        ind8 <- (Nextone + 1):(Nextone + la.mlm)
        pobs.mlm <- allprobs[, ind8, drop = FALSE]
        dimnames(pobs.mlm) <- list(rownames(eta), as.character(a.mlm))
        Nextone <- Nextone + la.mlm
      }
      if (tmp3.TF[ 9]) {
        ind9 <- (Nextone + 1):(Nextone + li.mlm)
        pstr.mlm <- allprobs[, ind9, drop = FALSE]
        dimnames(pstr.mlm) <- list(rownames(eta), as.character(i.mlm))
        Nextone <- Nextone + li.mlm
      }
      if (tmp3.TF[10]) {
        ind10 <- (Nextone + 1):(Nextone + ld.mlm)
        pdip.mlm <- allprobs[, ind10, drop = FALSE]
        dimnames(pdip.mlm) <- list(rownames(eta), as.character(d.mlm))
        Nextone <- Nextone + ld.mlm  # Not needed
      }
    }  # lall.len


    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) *
        dgaitdzeta(y, shape.p, log = TRUE,  # byrow.aid = F,
                  a.mix = a.mix, i.mix = i.mix, d.mix = d.mix,
                  a.mlm = a.mlm, i.mlm = i.mlm, d.mlm = d.mlm,
                  truncate = truncate,
                  max.support = as.vector( .max.support ),
                  shape.a = shape.a, shape.i = shape.i, 
                  shape.d = shape.d,
                  pobs.mix = pobs.mix, pstr.mix = pstr.mix,
                  pdip.mix = pdip.mix,
                  pobs.mlm = pobs.mlm, pstr.mlm = pstr.mlm,
                  pdip.mlm = pdip.mlm)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list(
    .lshape.p = lshape.p, .eshape.p = eshape.p,
    .lshape.a = lshape.a, .eshape.a = eshape.a,
    .lshape.i = lshape.i, .eshape.i = eshape.i,
    .lshape.d = lshape.d, .eshape.d = eshape.d,
    .lpstr.mix = lpstr.mix, .lpobs.mix = lpobs.mix,
    .lpdip.mix = lpdip.mix,
    .epstr.mix = epstr.mix, .epobs.mix = epobs.mix,
    .epdip.mix = epdip.mix,
    .a.mix = a.mix, .i.mix = i.mix, .d.mix = d.mix,
    .a.mlm = a.mlm, .i.mlm = i.mlm, .d.mlm = d.mlm,
    .truncate = truncate, .max.support = max.support ))),
  vfamily = c("gaitdzeta"),
  validparams = eval(substitute(function(eta, y, extra = NULL) {
    la.mix <- length((a.mix <- as.vector( .a.mix )))
    li.mix <- length((i.mix <- as.vector( .i.mix )))
    ld.mix <- length((d.mix <- as.vector( .d.mix )))
    la.mlm <- length((a.mlm <- as.vector( .a.mlm )))
    li.mlm <- length((i.mlm <- as.vector( .i.mlm )))
    ld.mlm <- length((d.mlm <- as.vector( .d.mlm )))
    lall.len <- la.mix + li.mix + ld.mix + la.mlm + li.mlm + ld.mlm
    small. <- 1e-14
    pobs.mix <- pstr.mix <- pdip.mix <- small.   # 4 rowSums():
    pobs.mlm <- pstr.mlm <- pdip.mlm <- matrix(small., NROW(eta), 1) 
    shape.a <- shape.i <- shape.d <- 0.5  # Needed

    if (!is.matrix(eta)) eta <- as.matrix(eta)
    shape.p <-
      cbind(eta2theta(eta[, 1], .lshape.p , earg = .eshape.p ))
    ind.shape.z <- 1  # Points to shape.p only.
    if (any(tmp3.TF[c(3, 5, 7)])) {  # At least one shape.[aid]
      ind.shape.z <- extra$indeta[c(1, 3, 5, 7), 1]  # Vectors
      ind.shape.z <- c(na.omit(ind.shape.z))  # At least one value

      shape.a <- if (!tmp3.TF[ 3]) shape.p else
        eta2theta(eta[, extra$indeta[3, 1]], .lshape.a , .eshape.a )
      shape.i <- if (!tmp3.TF[ 5]) shape.p else
        eta2theta(eta[, extra$indeta[5, 1]], .lshape.i , .eshape.i )
      shape.d <- if (!tmp3.TF[ 7]) shape.p else
        eta2theta(eta[, extra$indeta[7, 1]], .lshape.d , .eshape.d )
    }  # la.mix + li.mix + ld.mix > 0

    
    if (lall.len) {  # A MLM was fitted
      allprobs <-
        multilogitlink(eta[, -ind.shape.z, drop = FALSE],
                       inverse = TRUE)  # rowSums == 1

      Nextone <- 0  # Might not be used actually; 0, not 1
      if (tmp3.TF[ 2])
        pobs.mix <- allprobs[, (Nextone <- Nextone + 1)]
      if (tmp3.TF[ 4])
        pstr.mix <- allprobs[, (Nextone <- Nextone + 1)]
      if (tmp3.TF[ 6])
        pdip.mix <- allprobs[, (Nextone <- Nextone + 1)]

      if (tmp3.TF[ 8]) {
        ind8 <- (Nextone + 1):(Nextone + la.mlm)
        pobs.mlm <- allprobs[, ind8, drop = FALSE]
        dimnames(pobs.mlm) <- list(rownames(eta), as.character(a.mlm))
        Nextone <- Nextone + la.mlm
      }
      if (tmp3.TF[ 9]) {
        ind9 <- (Nextone + 1):(Nextone + li.mlm)
        pstr.mlm <- allprobs[, ind9, drop = FALSE]
        dimnames(pstr.mlm) <- list(rownames(eta), as.character(i.mlm))
        Nextone <- Nextone + li.mlm
      }
      if (tmp3.TF[10]) {
        ind10 <- (Nextone + 1):(Nextone + ld.mlm)
        pdip.mlm <- allprobs[, ind10, drop = FALSE]
        dimnames(pdip.mlm) <- list(rownames(eta), as.character(d.mlm))
        Nextone <- Nextone + ld.mlm  # Not needed
      }
    }  # lall.len

    okay.mlm <-
      all(is.finite(pobs.mlm)) && all(0 < pobs.mlm) &&
      all(is.finite(pstr.mlm)) && all(0 < pstr.mlm) &&
      all(is.finite(pdip.mlm)) && all(0 < pdip.mlm)
    okay.mix <-
      all(is.finite(shape.p)) && all(0 < shape.p) &&
      all(is.finite(shape.a)) && all(0 < shape.a) &&
      all(is.finite(shape.i)) && all(0 < shape.i) &&
      all(is.finite(shape.d)) && all(0 < shape.d) &&
      all(is.finite(pobs.mix)) && all(0 < pobs.mix) &&
      all(is.finite(pstr.mix)) && all(0 < pstr.mix) &&
      all(is.finite(pdip.mix)) && all(0 < pdip.mix) &&
      all(pobs.mix + pstr.mix + pdip.mix +
          rowSums(pobs.mlm) + rowSums(pstr.mlm) +
          rowSums(pdip.mlm) < 1)  # Combined
    okay.mlm && okay.mix
  }, list(
    .lshape.p = lshape.p, .eshape.p = eshape.p,
    .lshape.a = lshape.a, .eshape.a = eshape.a,
    .lshape.i = lshape.i, .eshape.i = eshape.i,
    .lshape.d = lshape.d, .eshape.d = eshape.d,
    .lpstr.mix = lpstr.mix, .lpobs.mix = lpobs.mix,
    .lpdip.mix = lpdip.mix,
    .epstr.mix = epstr.mix, .epobs.mix = epobs.mix,
    .epdip.mix = epdip.mix,
    .a.mix = a.mix, .i.mix = i.mix, .d.mix = d.mix,
    .a.mlm = a.mlm, .i.mlm = i.mlm, .d.mlm = d.mlm,
    .truncate = truncate, .max.support = max.support ))),
  simslot = eval(substitute(
  function(object, nsim) {
    pwts <- if (length(pwts <- object@prior.weights) > 0)
              pwts else weights(object, type = "prior")
    if (any(pwts != 1))
      warning("ignoring prior weights")
    eta <- predict(object)
    la.mix <- length((a.mix <- as.vector( .a.mix )))
    li.mix <- length((i.mix <- as.vector( .i.mix )))
    ld.mix <- length((d.mix <- as.vector( .d.mix )))
    la.mlm <- length((a.mlm <- as.vector( .a.mlm )))
    li.mlm <- length((i.mlm <- as.vector( .i.mlm )))
    ld.mlm <- length((d.mlm <- as.vector( .d.mlm )))
    truncate <- as.vector( .truncate )

    extra <- object@extra
    lall.len <- la.mix + li.mix + ld.mix + la.mlm + li.mlm + ld.mlm
    pobs.mix <- pstr.mix <- pdip.mix <- 0  # 4 rowSums()
    pobs.mlm <- pstr.mlm <- pdip.mlm <- 0  # matrix(0, NROW(eta), 1)
    shape.p <- cbind(eta2theta(eta[, 1], .lshape.p , .eshape.p ))
    ind.shape.z <- 1  # Points to shape.p only.
    shape.a <- shape.i <-
    shape.d <- shape.p  # Needed; and answer not corrupted
    tmp3.TF <- ( .tmp3.TF )


    if (any(tmp3.TF[c(3, 5, 7)])) {  # At least one shape.[aid]
      ind.shape.z <- extra$indeta[c(1, 3, 5, 7), 'launch']  # Vectors
      ind.shape.z <- c(na.omit(ind.shape.z))  # At least one value

      shape.a <- if (!tmp3.TF[ 3]) shape.p else
        eta2theta(eta[, extra$indeta[3, 1]], .lshape.a , .eshape.a )
      shape.i <- if (!tmp3.TF[ 5]) shape.p else
        eta2theta(eta[, extra$indeta[5, 1]], .lshape.i , .eshape.i )
      shape.d <- if (!tmp3.TF[ 7]) shape.p else
        eta2theta(eta[, extra$indeta[7, 1]], .lshape.d , .eshape.d )
    }  # la.mix + li.mix + ld.mix > 0


    if (lall.len) {  # A AMLM was fitted
      allprobs <-
        multilogitlink(eta[, -ind.shape.z, drop = FALSE],
                       inverse = TRUE)  # rowSums == 1
      Nextone <- 0  # Might not be used actually; 0, not 1
      if (tmp3.TF[ 2])
        pobs.mix <- allprobs[, (Nextone <- Nextone + 1)]
      if (tmp3.TF[ 4])
        pstr.mix <- allprobs[, (Nextone <- Nextone + 1)]
      if (tmp3.TF[ 6])
        pdip.mix <- allprobs[, (Nextone <- Nextone + 1)]

      if (tmp3.TF[ 8]) {
        ind8 <- (Nextone + 1):(Nextone + la.mlm)
        pobs.mlm <- allprobs[, ind8, drop = FALSE]
        dimnames(pobs.mlm) <- list(rownames(eta), as.character(a.mlm))
        Nextone <- Nextone + la.mlm
      }
      if (tmp3.TF[ 9]) {
        ind9 <- (Nextone + 1):(Nextone + li.mlm)
        pstr.mlm <- allprobs[, ind9, drop = FALSE]
        dimnames(pstr.mlm) <- list(rownames(eta), as.character(i.mlm))
        Nextone <- Nextone + li.mlm
      }
      if (tmp3.TF[10]) {
        ind10 <- (Nextone + 1):(Nextone + ld.mlm)
        pdip.mlm <- allprobs[, ind10, drop = FALSE]
        dimnames(pdip.mlm) <- list(rownames(eta), as.character(d.mlm))
        Nextone <- Nextone + ld.mlm  # Not needed
      }
    }  # lall.len

    rgaitdzeta(nsim * length(shape.p), shape.p,
               pobs.mlm = pobs.mlm, pstr.mlm = pstr.mlm,
               pobs.mix = pobs.mix, pstr.mix = pstr.mix,
               pdip.mix = pdip.mix, pdip.mlm = pdip.mlm,
               shape.a = shape.a, shape.i = shape.i, 
               shape.d = shape.d,
               a.mix = a.mix, i.mix = i.mix, d.mix = d.mix,
               a.mlm = a.mlm, i.mlm = i.mlm, d.mlm = d.mlm,
               truncate = .truncate , max.support = .max.support )
  }, list(
    .lshape.p = lshape.p, .eshape.p = eshape.p,
    .lshape.a = lshape.a, .eshape.a = eshape.a,
    .lshape.i = lshape.i, .eshape.i = eshape.i,
    .lshape.d = lshape.d, .eshape.d = eshape.d,
    .lpstr.mix = lpstr.mix, .lpobs.mix = lpobs.mix,
    .lpdip.mix = lpdip.mix,
    .epstr.mix = epstr.mix, .epobs.mix = epobs.mix,
    .epdip.mix = epdip.mix,
    .tmp3.TF = tmp3.TF,
    .a.mix = a.mix, .i.mix = i.mix, .d.mix = d.mix,
    .a.mlm = a.mlm, .i.mlm = i.mlm, .d.mlm = d.mlm,
    .truncate = truncate, .max.support = max.support ))),
  deriv = eval(substitute(expression({

    tmp3.TF <- ( .tmp3.TF )
    calA.p  <- tmp3.TF[ 2]
    calI.p  <- tmp3.TF[ 4]
    calD.p  <- tmp3.TF[ 6]
    calA.np <- tmp3.TF[ 8]
    calI.np <- tmp3.TF[ 9]
    calD.np <- tmp3.TF[10]

    Denom1.a <- Denom1.i <- Denom1.d <-
                Denom2.i <- Denom2.d <- 0  # Denom2.a is unneeded

    if (!is.matrix(eta)) eta <- as.matrix(eta)
    la.mix <- length((a.mix <- as.vector( .a.mix )))
    li.mix <- length((i.mix <- as.vector( .i.mix )))
    ld.mix <- length((d.mix <- as.vector( .d.mix )))
    la.mlm <- length((a.mlm <- as.vector( .a.mlm )))
    li.mlm <- length((i.mlm <- as.vector( .i.mlm )))
    ld.mlm <- length((d.mlm <- as.vector( .d.mlm )))
    truncate <- as.vector( .truncate )
    max.support <- as.vector( .max.support )



    
    lall.len <- la.mix + li.mix + ld.mix + la.mlm + li.mlm + ld.mlm
    pobs.mix <- pstr.mix <- pdip.mix <- 0  # 4 rowSums()
    pobs.mlm <- pstr.mlm <- pdip.mlm <- 0  # matrix(0, NROW(eta), 1)
    shape.p <- cbind(eta2theta(eta[, 1], .lshape.p , .eshape.p ))
    ind.shape.z <- 1  # Points to shape.p only.
    shape.a <- shape.i <-
    shape.d <- shape.p  # Needed; and answer not corrupted

    if (any(tmp3.TF[c(3, 5, 7)])) {  # At least one shape.[aid]
      ind.shape.z <- extra$indeta[c(1, 3, 5, 7), 'launch']  # Vectors
      ind.shape.z <- c(na.omit(ind.shape.z))  # At least one value

      shape.a <- if (!tmp3.TF[ 3]) shape.p else
        eta2theta(eta[, extra$indeta[3, 1]], .lshape.a , .eshape.a )
      shape.i <- if (!tmp3.TF[ 5]) shape.p else
        eta2theta(eta[, extra$indeta[5, 1]], .lshape.i , .eshape.i )
      shape.d <- if (!tmp3.TF[ 7]) shape.p else
        eta2theta(eta[, extra$indeta[7, 1]], .lshape.d , .eshape.d )
    }  # la.mix + li.mix + ld.mix > 0



    if (lall.len) {  # A MLM was fitted.
      allprobs <-
        multilogitlink(eta[, -ind.shape.z, drop = FALSE],
                       refLevel = "(Last)",  # Make sure
                       inverse = TRUE)  # rowSums == 1
      minprob.baseline <- min(allprobs[, ncol(allprobs)], na.rm = TRUE)
      if (anyNA(allprobs))
        warning("there are NAs here in slot linkinv")
      if (min(allprobs) == 0 || max(allprobs) == 1) {
        warning("fitted probabilities numerically 0 or 1 occurred")
      } else
      if (minprob.baseline < 0.10)
        warning("Minimum baseline (reserve) probability close to 0")
      if (control$trace)
        cat("Minimum baseline (reserve) probability = ",
            format(minprob.baseline, digits = 3), "\n")


      Nextone <- 0  # Might not be used actually; 0, not 1
      if (tmp3.TF[ 2])
        pobs.mix <- allprobs[, (Nextone <- Nextone + 1)]
      if (tmp3.TF[ 4])
        pstr.mix <- allprobs[, (Nextone <- Nextone + 1)]
      if (tmp3.TF[ 6])
        pdip.mix <- allprobs[, (Nextone <- Nextone + 1)]

      if (tmp3.TF[ 8]) {
        ind8 <- (Nextone + 1):(Nextone + la.mlm)
        pobs.mlm <- allprobs[, ind8, drop = FALSE]
        dimnames(pobs.mlm) <- list(rownames(eta), as.character(a.mlm))
        Nextone <- Nextone + la.mlm
      }
      if (tmp3.TF[ 9]) {
        ind9 <- (Nextone + 1):(Nextone + li.mlm)
        pstr.mlm <- allprobs[, ind9, drop = FALSE]
        dimnames(pstr.mlm) <- list(rownames(eta), as.character(i.mlm))
        Nextone <- Nextone + li.mlm
      }
      if (tmp3.TF[10]) {
        ind10 <- (Nextone + 1):(Nextone + ld.mlm)
        pdip.mlm <- allprobs[, ind10, drop = FALSE]
        dimnames(pdip.mlm) <- list(rownames(eta), as.character(d.mlm))
        Nextone <- Nextone + ld.mlm  # Not needed
      }
    }  # lall.len


    ltruncat <- length(truncate)
    M1 <- max(extra$indeta, na.rm = TRUE)
    NOS <- ncol(eta) / M1  # extra$NOS
    if (NOS != 1) stop("can only handle 1 response")



    is.a.mixed <- if (tmp3.TF[ 2])
      rowSums(extra$skip.mix.a) > 0 else rep(FALSE, n)
    is.i.mixed <- if (tmp3.TF[ 4])
      rowSums(extra$skip.mix.i) > 0 else rep(FALSE, n)
    is.d.mixed <- if (tmp3.TF[ 6])
      rowSums(extra$skip.mix.d) > 0 else rep(FALSE, n)
    is.a.mlmed <- if (tmp3.TF[ 8])
      rowSums(extra$skip.mlm.a) > 0 else rep(FALSE, n)
    is.i.mlmed <- if (tmp3.TF[ 9])
      rowSums(extra$skip.mlm.i) > 0 else rep(FALSE, n)
    is.d.mlmed <- if (tmp3.TF[10])
      rowSums(extra$skip.mlm.d) > 0 else rep(FALSE, n)

    is.ns <- !is.a.mlmed & !is.i.mlmed  & !is.d.mlmed  &
             !is.a.mixed & !is.i.mixed  & !is.d.mixed  # & !is.truncd


    prob.mlm.a <- if (la.mlm) rowSums(pobs.mlm) else 0  # scalar okay
    prob.mlm.i <- if (li.mlm) rowSums(pstr.mlm) else 0  # scalar okay
    prob.mlm.d <- if (ld.mlm) rowSums(pdip.mlm) else 0  # scalar okay


   pmf.deriv1 <- function(y, shape) {
      -dzeta(y, shape) * (log(y) +
        zeta(shape + 1, deriv = 1) / zeta(shape + 1))
    }
    pmf.deriv2 <- function(y, shape) {
      tmp2 <- zeta(shape + 1, deriv = 1) / zeta(shape + 1)
      dzeta(y, shape) * ((log(y) + tmp2)^2 -
        zeta(shape + 1, deriv = 2) / zeta(shape + 1) + tmp2^2)
    }
    

    sumD.mix.1a.p <- sumD.mix.2a.p <- matrix(0, n, NOS)
    if (la.mix > 0) {  # \calA_p
      DA.mix.0mat.a <-  # Matches naming convention further below
      DA.mix.1mat.a <- matrix(0, n, la.mix)
      for (jay in seq(la.mix)) {
        aval <- a.mix[jay]
        sumD.mix.1a.p <- sumD.mix.1a.p + pmf.deriv1(aval, shape.p)
        sumD.mix.2a.p <- sumD.mix.2a.p + pmf.deriv2(aval, shape.p)
        pmf.a <- dzeta(aval, shape.a)
        DA.mix.0mat.a[, jay] <- pmf.a
        DA.mix.1mat.a[, jay] <- pmf.deriv1(aval, shape.a)
      }
      Denom1.a <- rowSums(DA.mix.1mat.a)  # aka sumD.mix.1a.a
    }  # la.mix > 0





    if (li.mix) {
      DI.mix.0mat.i <-  # wrt inflated distribution
      DI.mix.1mat.i <- DI.mix.2mat.i <- matrix(0, n, li.mix)
      DP.mix.0mat.i <-  # wrt parent distribution
      DP.mix.1mat.i <- DP.mix.2mat.i <- matrix(0, n, li.mix)
        for (jay in seq(li.mix)) {
          ival <- i.mix[jay]
          pmf.i <- dzeta(ival, shape.i)
          DI.mix.0mat.i[, jay] <- pmf.i
          DI.mix.1mat.i[, jay] <- pmf.deriv1(ival, shape.i)
          DI.mix.2mat.i[, jay] <- pmf.deriv2(ival, shape.i)
          pmf.p <- dzeta(ival, shape.p)
          DP.mix.0mat.i[, jay] <- pmf.p
          DP.mix.1mat.i[, jay] <- pmf.deriv1(ival, shape.p)
          DP.mix.2mat.i[, jay] <- pmf.deriv2(ival, shape.p)
        }  # jay
      Denom1.i <- rowSums(DI.mix.1mat.i)
      Denom2.i <- rowSums(DI.mix.2mat.i)
    }  # li.mix


    if (ld.mix) {
      DD.mix.0mat.d <-  # wrt deflated distribution
      DD.mix.1mat.d <- DD.mix.2mat.d <- matrix(0, n, ld.mix)
      DP.mix.0mat.d <-  # wrt parent distribution
      DP.mix.1mat.d <- DP.mix.2mat.d <- matrix(0, n, ld.mix)
        for (jay in seq(ld.mix)) {
          dval <- d.mix[jay]
          pmf.d <- dzeta(dval, shape.d)
          DD.mix.0mat.d[, jay] <- pmf.d
          DD.mix.1mat.d[, jay] <- pmf.deriv1(dval, shape.d)
          DD.mix.2mat.d[, jay] <- pmf.deriv2(dval, shape.d)
          pmf.p <- dzeta(dval, shape.p)
          DP.mix.0mat.d[, jay] <- pmf.p
          DP.mix.1mat.d[, jay] <- pmf.deriv1(dval, shape.p)
          DP.mix.2mat.d[, jay] <- pmf.deriv2(dval, shape.p)
        }  # jay
      Denom1.d <- rowSums(DD.mix.1mat.d)
      Denom2.d <- rowSums(DD.mix.2mat.d)
    }  # ld.mix

    Bits <- moments.gaitdcombo.zeta(shape.p,
              pobs.mix = pobs.mix, pobs.mlm = pobs.mlm,
              pstr.mix = pstr.mix, pstr.mlm = pstr.mlm,
              pdip.mix = pdip.mix, pdip.mlm = pdip.mlm,
              a.mix = a.mix, i.mix = i.mix, d.mix = d.mix,
              a.mlm = a.mlm, i.mlm = i.mlm, d.mlm = d.mlm,
              shape.a = shape.a, shape.i = shape.i,
              shape.d = shape.d,
              truncate = truncate, max.support = max.support)


    sumD.mlm.1a.p <- sumD.mlm.2a.p <- matrix(0, n, NOS)
    if (la.mlm)
      for (aval in a.mlm) {
        sumD.mlm.1a.p <- sumD.mlm.1a.p + pmf.deriv1(aval, shape.p)
        sumD.mlm.2a.p <- sumD.mlm.2a.p + pmf.deriv2(aval, shape.p)
      }


    Denom0.p <- c(Bits[["cdf.max.s"]]   - Bits[["SumT0.p"]] -
                  Bits[["SumA0.mix.p"]] - Bits[["SumA0.mlm.p"]])
    Numer <- 1 - pobs.mix - pstr.mix - prob.mlm.a - prob.mlm.i +
                 pdip.mix +            prob.mlm.d
    Denom0.a <- c(Bits[["SumA0.mix.a"]])  # Not .p
    Denom0.i <- c(Bits[["SumI0.mix.i"]])
    Denom0.d <- c(Bits[["SumD0.mix.d"]])


 

    Dp.mlm.0Mat.i <-  # wrt parent distribution
    Dp.mlm.1Mat.i <- Dp.mlm.2Mat.i <- matrix(0, n, NOS)
    if (li.mlm > 0) {
      Dp.mlm.0Mat.i <-  # wrt parent distribution
      Dp.mlm.1Mat.i <- Dp.mlm.2Mat.i <- matrix(0, n, li.mlm)
      for (jay in seq(li.mlm)) {
        ival <- i.mlm[jay]
        pmf.p <- dzeta(ival, shape.p)
        Dp.mlm.0Mat.i[, jay] <- pmf.p
        Dp.mlm.1Mat.i[, jay] <- pmf.deriv1(ival, shape.p)
        Dp.mlm.2Mat.i[, jay] <- pmf.deriv2(ival, shape.p)
      }  # jay
    }  # li.mlm



    Dp.mlm.0Mat.d <-  # wrt parent distribution
    Dp.mlm.1Mat.d <- Dp.mlm.2Mat.d <- matrix(0, n, NOS)
    if (ld.mlm > 0) {
      Dp.mlm.0Mat.d <-  # wrt parent distribution
      Dp.mlm.1Mat.d <- Dp.mlm.2Mat.d <- matrix(0, n, ld.mlm)
      for (jay in seq(ld.mlm)) {
        dval <- d.mlm[jay]
        pmf.p <- dzeta(dval, shape.p)
        Dp.mlm.0Mat.d[, jay] <- pmf.p
        Dp.mlm.1Mat.d[, jay] <- pmf.deriv1(dval, shape.p)
        Dp.mlm.2Mat.d[, jay] <- pmf.deriv2(dval, shape.p)
      }  # jay
    }  # ld.mlm





    

    sumD.1t.p <- sumD.2t.p <-
    sumD.1t.a <- sumD.2t.a <-
    sumD.1t.i <- sumD.2t.i <-
    sumD.1t.d <- sumD.2t.d <- matrix(0, n, NOS)
    if (ltruncat)
      for (tval in truncate) {
        sumD.1t.p <- sumD.1t.p + pmf.deriv1(tval, shape.p)
        sumD.2t.p <- sumD.2t.p + pmf.deriv2(tval, shape.p)
        sumD.1t.a <- sumD.1t.a + pmf.deriv1(tval, shape.a)
        sumD.2t.a <- sumD.2t.a + pmf.deriv2(tval, shape.a)
        sumD.1t.i <- sumD.1t.i + pmf.deriv1(tval, shape.i)
        sumD.2t.i <- sumD.2t.i + pmf.deriv2(tval, shape.i)
        sumD.1t.d <- sumD.1t.d + pmf.deriv1(tval, shape.d)
        sumD.2t.d <- sumD.2t.d + pmf.deriv2(tval, shape.d)
      }




      
    if (is.finite(max.support)) {
      stop("argument 'max.support' must be 'Inf'.")
    }  # is.finite(max.support)

     



    Denom1.p <- c(-sumD.1t.p - sumD.mlm.1a.p - sumD.mix.1a.p)
    Denom2.p <- c(-sumD.2t.p - sumD.mlm.2a.p - sumD.mix.2a.p)




    d0B.PI.mlm <- Dp.mlm.0Mat.i / Denom0.p
    d1B.PI.mlm <- Dp.mlm.1Mat.i / Denom0.p -  # This is most general
                  Dp.mlm.0Mat.i * Denom1.p / Denom0.p^2
    d2B.PI.mlm <-     Dp.mlm.2Mat.i / Denom0.p -
                  2 * Dp.mlm.1Mat.i * Denom1.p / Denom0.p^2 -
                      Dp.mlm.0Mat.i * Denom2.p / Denom0.p^2 +
                  2 * Dp.mlm.0Mat.i * (Denom1.p^2) / Denom0.p^3





    d0B.PD.mlm <- Dp.mlm.0Mat.d / Denom0.p
    d1B.PD.mlm <- Dp.mlm.1Mat.d / Denom0.p -  # This is most general
                  Dp.mlm.0Mat.d * Denom1.p / Denom0.p^2
    d2B.PD.mlm <-     Dp.mlm.2Mat.d / Denom0.p -
                  2 * Dp.mlm.1Mat.d * Denom1.p / Denom0.p^2 -
                      Dp.mlm.0Mat.d * Denom2.p / Denom0.p^2 +
                  2 * Dp.mlm.0Mat.d * (Denom1.p^2) / Denom0.p^3


    DELTA.i.mlm <- if (li.mlm > 0) {
      Numer * d0B.PI.mlm + pstr.mlm  # n x li.mlm.
    } else {
      matrix(0, n, 1)  # If li.mlm == 0, for rowSums().
    }

    DELTA.d.mlm <- if (ld.mlm > 0) {
      Numer * d0B.PD.mlm - pdip.mlm  # n x ld.mlm.
    } else {
      matrix(0, n, 1)  # If ld.mlm == 0, for rowSums().
    }


    if (li.mix > 0) {
      d0A.i <- DI.mix.0mat.i / Denom0.i
      d0B.PI.mix <- DP.mix.0mat.i / Denom0.p
      DELTA.i.mix <- Numer * d0B.PI.mix + pstr.mix * d0A.i

      d1A.i <- (DI.mix.1mat.i - DI.mix.0mat.i *
                Denom1.i / Denom0.i) / Denom0.i
      d2A.i <- (DI.mix.2mat.i - (2 * DI.mix.1mat.i * Denom1.i +
                DI.mix.0mat.i * Denom2.i) / Denom0.i +
            2 * DI.mix.0mat.i * (Denom1.i / Denom0.i)^2) / Denom0.i

      d1B.PI.mix <- DP.mix.1mat.i / Denom0.p -
                    DP.mix.0mat.i * Denom1.p / Denom0.p^2
      d2B.PI.mix <-     DP.mix.2mat.i /  Denom0.p -
                    2 * DP.mix.1mat.i *  Denom1.p    / Denom0.p^2 -
                        DP.mix.0mat.i *  Denom2.p    / Denom0.p^2 +
                    2 * DP.mix.0mat.i * (Denom1.p^2) / Denom0.p^3
    }  # li.mix > 0



    if (ld.mix > 0) {
      d0A.d <- DD.mix.0mat.d / Denom0.d
      d0B.PD.mix <- DP.mix.0mat.d / Denom0.p
      DELTA.d.mix <- Numer * d0B.PD.mix - pdip.mix * d0A.d

      d1A.d <- (DD.mix.1mat.d - DD.mix.0mat.d *
                Denom1.d / Denom0.d) / Denom0.d
      d2A.d <- (DD.mix.2mat.d - (2 * DD.mix.1mat.d * Denom1.d +
                DD.mix.0mat.d * Denom2.d) / Denom0.d +
            2 * DD.mix.0mat.d * (Denom1.d / Denom0.d)^2) / Denom0.d

      d1B.PD.mix <- DP.mix.1mat.d / Denom0.p -
                    DP.mix.0mat.d * Denom1.p / Denom0.p^2
      d2B.PD.mix <-     DP.mix.2mat.d /  Denom0.p -
                    2 * DP.mix.1mat.d *  Denom1.p    / Denom0.p^2 -
                        DP.mix.0mat.d *  Denom2.p    / Denom0.p^2 +
                    2 * DP.mix.0mat.d * (Denom1.p^2) / Denom0.p^3
    }  # ld.mix > 0




    if (la.mix) {
      d0A.a <- DA.mix.0mat.a / Denom0.a
      d1A.a <- DA.mix.1mat.a / Denom0.a -
               DA.mix.0mat.a * Denom1.a / Denom0.a^2
    }  # la.mix







    fred0.p <- zeta(shape.p + 1)
    fred1.p <- zeta(shape.p + 1, deriv = 1)
    dl.dshape.p <- -log(y) - fred1.p / fred0.p  # Usual formula
    dl.dshape.p[!is.ns] <- 0  # For is.a.mixed & is.a.mlmed

    dl.dshape.a <- dl.dshape.i <- dl.dshape.d <- numeric(n)
    dl.dpstr.mix <- (-1) / Numer  # \notin A, I, T, D
    dl.dpstr.mix[is.a.mixed] <- 0
    dl.dpstr.mix[is.a.mlmed] <- 0
    dl.dpdip.mix <- (+1) / Numer  # \notin A, I, T, D
    dl.dpdip.mix[is.a.mixed] <- 0
    dl.dpdip.mix[is.a.mlmed] <- 0
    dl.dpobs.mix <- numeric(n)  # 0 for \calA_{np}
    dl.dpobs.mix[is.ns] <- (-1) / Numer[is.ns]
    dl.dpobs.mlm <-
    dl.dpstr.mlm <- matrix(0, n, 1)  # May be unneeded
    dl.dpdip.mlm <- matrix(0, n, max(1, ld.mlm))  # Initzed if used.
    dl.dpdip.mlm[is.ns, ] <- 1 / Numer[is.ns]



    if (tmp3.TF[ 8] && la.mlm) {  # aka \calA_{np}
      dl.dpobs.mlm <- matrix(-1 / Numer, n, la.mlm)  # \notin calS
      dl.dpobs.mlm[!is.ns, ] <- 0  # For a.mix only really
      for (jay in seq(la.mlm)) {
        aval <- a.mlm[jay]
        is.alt.j.mlm <- extra$skip.mlm.a[, jay]  # Logical vector
        tmp7a <- 1 / pobs.mlm[is.alt.j.mlm, jay]
        dl.dpobs.mlm[is.alt.j.mlm, jay] <- tmp7a
      }  # jay
    }  # la.mlm



    dl.dshape.p[is.ns] <- dl.dshape.p[is.ns] -
                           (Denom1.p / Denom0.p)[is.ns]


    
    if (tmp3.TF[ 9] && li.mlm > 0) {  # aka \calI_{np}
      dl.dpstr.mlm <- matrix(-1 / Numer, n, li.mlm)
      dl.dpstr.mlm[!is.ns, ] <- 0  # For a.mlm and a.mix

      for (jay in seq(li.mlm)) {
        is.inf.j.mlm <- extra$skip.mlm.i[, jay]  # Logical vector
        tmp7i <- Numer * d1B.PI.mlm[, jay] / DELTA.i.mlm[, jay]
        dl.dshape.p[is.inf.j.mlm] <- tmp7i[is.inf.j.mlm]


        tmp9i <- d0B.PI.mlm[, jay] / DELTA.i.mlm[, jay]
        n.tmp <- -tmp9i[is.inf.j.mlm]
        p.tmp <- +tmp9i[is.inf.j.mlm]
        if (tmp3.TF[ 8] && la.mlm) dl.dpobs.mlm[is.inf.j.mlm, ] <- n.tmp
        if (tmp3.TF[ 2] && la.mix) dl.dpobs.mix[is.inf.j.mlm  ] <- n.tmp
        if (tmp3.TF[ 4] && li.mix) dl.dpstr.mix[is.inf.j.mlm  ] <- n.tmp
        if (tmp3.TF[10] && ld.mlm) dl.dpdip.mlm[is.inf.j.mlm, ] <- p.tmp
        if (tmp3.TF[ 6] && ld.mix) dl.dpdip.mix[is.inf.j.mlm  ] <- p.tmp


        tmp8 <- (1 - d0B.PI.mlm[, jay]) / DELTA.i.mlm[, jay]
        dl.dpstr.mlm[is.inf.j.mlm, ] <- n.tmp  # tmp9[is.inf.j.mlm]
        dl.dpstr.mlm[is.inf.j.mlm, jay] <- tmp8[is.inf.j.mlm]
      }  # jay
    }  # li.mlm > 0






    if (tmp3.TF[10] && ld.mlm > 0) {  # aka \calD_{np}

      for (jay in seq(ld.mlm)) {
        is.def.j.mlm <- extra$skip.mlm.d[, jay]  # Logical vector
        tmp7d <- Numer * d1B.PD.mlm[, jay] / DELTA.d.mlm[, jay]
        dl.dshape.p[is.def.j.mlm] <- tmp7d[is.def.j.mlm]  # 20211020

 
        tmp9d <- d0B.PD.mlm[, jay] / DELTA.d.mlm[, jay]
        p.tmp <- +tmp9d[is.def.j.mlm]
        n.tmp <- -tmp9d[is.def.j.mlm]
        if (tmp3.TF[ 9] && li.mlm) dl.dpstr.mlm[is.def.j.mlm, ] <- n.tmp
        if (tmp3.TF[ 4] && li.mix) dl.dpstr.mix[is.def.j.mlm  ] <- n.tmp
        if (tmp3.TF[ 8] && la.mlm) dl.dpobs.mlm[is.def.j.mlm, ] <- n.tmp
        if (tmp3.TF[ 2] && la.mix) dl.dpobs.mix[is.def.j.mlm  ] <- n.tmp
        if (tmp3.TF[ 6] && ld.mix) dl.dpdip.mix[is.def.j.mlm  ] <- p.tmp
                                   dl.dpdip.mlm[is.def.j.mlm, ] <- p.tmp
        dl.dpdip.mlm[is.def.j.mlm, jay] <-
        dl.dpdip.mlm[is.def.j.mlm, jay] -
        1 / DELTA.d.mlm[is.def.j.mlm, jay]

      }  # jay
    }  # ld.mlm > 0




    



    if (tmp3.TF[ 2] && la.mix) {  # aka \calA_{p}
      dl.dpobs.mix[is.a.mixed] <- 1 / pobs.mix[is.a.mixed]

      if (tmp3.TF[ 3] && la.mix > 1)
        for (jay in seq(la.mix)) {
          is.alt.j.mix <- extra$skip.mix.a[, jay]  # Logical vector
          tmp2 <- d1A.a[, jay] / d0A.a[, jay]
          dl.dshape.a[is.alt.j.mix] <- tmp2[is.alt.j.mix]  # ccc.
        }  # jay
    }  # la.mix




    if (tmp3.TF[ 4] && li.mix > 0) {  # aka \calI_{p}
      for (jay in seq(li.mix)) {
        ival <- i.mix[jay]
        is.inf.j.mix <- extra$skip.mix.i[, jay]  # Logical vector
        tmp7b <- Numer * d1B.PI.mix[, jay] / DELTA.i.mix[, jay]
        dl.dshape.p[is.inf.j.mix] <- tmp7b[is.inf.j.mix]
        tmp8 <- (d0A.i[, jay] - d0B.PI.mix[, jay]) / DELTA.i.mix[, jay]
        dl.dpstr.mix[is.inf.j.mix] <- tmp8[is.inf.j.mix]
        if (li.mix > 1) {
          tmp2 <- pstr.mix * d1A.i[, jay] / DELTA.i.mix[, jay]
          dl.dshape.i[is.inf.j.mix] <- tmp2[is.inf.j.mix]
        }





        tmp9i <- d0B.PI.mix[, jay] / DELTA.i.mix[, jay]
        n.tmp <- -tmp9i[is.inf.j.mix]
        p.tmp <- +tmp9i[is.inf.j.mix]
        if (tmp3.TF[ 2] && la.mix) dl.dpobs.mix[is.inf.j.mix  ] <- n.tmp
        if (tmp3.TF[ 8] && la.mlm) dl.dpobs.mlm[is.inf.j.mix, ] <- n.tmp
        if (tmp3.TF[ 9] && li.mlm) dl.dpstr.mlm[is.inf.j.mix, ] <- n.tmp
        if (tmp3.TF[10] && ld.mlm) dl.dpdip.mlm[is.inf.j.mix, ] <- p.tmp
        if (tmp3.TF[ 6] && ld.mix) dl.dpdip.mix[is.inf.j.mix  ] <- p.tmp

      }  # jay
    }  # li.mix > 0




    if (tmp3.TF[ 6] && ld.mix > 0) {  # aka \calD_{p}
      for (jay in seq(ld.mix)) {
        dval <- d.mix[jay]
        is.def.j.mix <- extra$skip.mix.d[, jay]  # Logical vector
        tmp7b <- Numer * d1B.PD.mix[, jay] / DELTA.d.mix[, jay]
        dl.dshape.p[is.def.j.mix] <- tmp7b[is.def.j.mix]
        tmp8 <- (d0B.PD.mix[, jay] - d0A.d[, jay]) / DELTA.d.mix[, jay]
        dl.dpdip.mix[is.def.j.mix] <- tmp8[is.def.j.mix]

        if (ld.mix > 1) {
          tmp2 <- (-pdip.mix) * d1A.d[, jay] / DELTA.d.mix[, jay]
          dl.dshape.d[is.def.j.mix] <- tmp2[is.def.j.mix]
        }


        tmp9d <- d0B.PD.mix[, jay] / DELTA.d.mix[, jay]
        n.tmp <- -tmp9d[is.def.j.mix]
        p.tmp <- +tmp9d[is.def.j.mix]
        if (tmp3.TF[ 9] && li.mlm) dl.dpstr.mlm[is.def.j.mix, ] <- n.tmp
        if (tmp3.TF[ 4] && li.mix) dl.dpstr.mix[is.def.j.mix  ] <- n.tmp
        if (tmp3.TF[ 8] && la.mlm) dl.dpobs.mlm[is.def.j.mix, ] <- n.tmp
        if (tmp3.TF[ 2] && la.mix) dl.dpobs.mix[is.def.j.mix  ] <- n.tmp
        if (tmp3.TF[10] && ld.mlm) dl.dpdip.mlm[is.def.j.mix, ] <- p.tmp

      }  # jay
   }  # ld.mix > 0







    new.ansd <- matrix(0, n, M)  # Same dimension as eta
    tmp3.TF <- !is.na(rowSums(extra$indeta))



    if (lall.len) {  # An MLM fitted
      all6.dldp <- cbind(if (tmp3.TF[ 2]) dl.dpobs.mix else NULL,
                         if (tmp3.TF[ 4]) dl.dpstr.mix else NULL,
                         if (tmp3.TF[ 6]) dl.dpdip.mix else NULL,
                         if (tmp3.TF[ 8]) dl.dpobs.mlm else NULL,
                         if (tmp3.TF[ 9]) dl.dpstr.mlm else NULL,
                         if (tmp3.TF[10]) dl.dpdip.mlm else NULL)




      rSs.tmp <- rowSums(allprobs[, -ncol(allprobs), drop = FALSE] *
                         all6.dldp)
      new.ansd[, -ind.shape.z] <- allprobs[, -ncol(allprobs)] *
                                   (all6.dldp - rSs.tmp)
    }  # lall.len


    

    dshape.p.deta <- dtheta.deta(shape.p, .lshape.p , .eshape.p )
    if (tmp3.TF[ 3])
      dshape.a.deta <- dtheta.deta(shape.a, .lshape.a , .eshape.a )
    if (tmp3.TF[ 5])
      dshape.i.deta <- dtheta.deta(shape.i, .lshape.i , .eshape.i )
    if (tmp3.TF[ 7])
      dshape.d.deta <- dtheta.deta(shape.d, .lshape.d , .eshape.d )
    new.ansd[, 1] <- dl.dshape.p * dshape.p.deta
    if (tmp3.TF[ 3])
      new.ansd[, extra$indeta[3, 1]] <- dl.dshape.a * dshape.a.deta
    if (tmp3.TF[ 5])
      new.ansd[, extra$indeta[5, 1]] <- dl.dshape.i * dshape.i.deta
    if (tmp3.TF[ 7])
      new.ansd[, extra$indeta[7, 1]] <- dl.dshape.d * dshape.d.deta
    onecoln.indeta <- extra$indeta[1:7, ]  # One coln params only
    onecoln.indeta <- na.omit(onecoln.indeta)  # Only those present
    allcnames <- c(rownames(onecoln.indeta),
                   as.character(c(a.mlm, i.mlm, d.mlm)))
    colnames(new.ansd) <- allcnames

 


    c(w) * new.ansd
  }), list(
    .lshape.p = lshape.p, .eshape.p = eshape.p,
    .lshape.a = lshape.a, .eshape.a = eshape.a,
    .lshape.i = lshape.i, .eshape.i = eshape.i,
    .lshape.d = lshape.d, .eshape.d = eshape.d,
    .lpstr.mix = lpstr.mix, .lpobs.mix = lpobs.mix,
    .lpdip.mix = lpdip.mix,
    .epstr.mix = epstr.mix, .epobs.mix = epobs.mix,
    .epdip.mix = epdip.mix,
    .a.mix = a.mix, .i.mix = i.mix, .d.mix = d.mix,
    .a.mlm = a.mlm, .i.mlm = i.mlm, .d.mlm = d.mlm,
    .tmp3.TF = tmp3.TF,  # .tmp3 = tmp3,
    .truncate = truncate, .max.support = max.support ))),

  weight = eval(substitute(expression({  # gaitdzeta



    wz <- matrix(0, n, M * (M + 1) / 2)  # The complete size
    probns <- Numer * (1 -
      (c(Bits[["SumI0.mix.p"]] + Bits[["SumI0.mlm.p"]]) +
       c(Bits[["SumD0.mix.p"]] + Bits[["SumD0.mlm.p"]])) / Denom0.p)



    if (min(probns) < 0 || 1 < max(probns))
      stop("variable 'probns' for P(nonspecial) is out of range")







    zero0n <- numeric(n)
    ned2l.dpobs.mix.shape.p <- zero0n  # mB overwritten below [4279]
    ned2l.dpobs.mix.shape.a <- zero0n  # Fini; (2, 3) element
    ned2l.dpobs.mix.shape.i <- zero0n  # mB overwritten below
    ned2l.dpobs.mix.shape.d <- zero0n  # mB overwritten below
    ned2l.dpstr.mix.shape.p <- zero0n  # Optional (1, 4) element
    ned2l.dpstr.mix.shape.a <- zero0n  # Final; nothing to do
    ned2l.dpstr.mix.shape.i <- zero0n  # mB overwritten below
    ned2l.dpstr.mix.shape.d <- zero0n  # mB overwritten below
    ned2l.dpdip.mix.shape.p <- zero0n  # Optional (1, 6) element




      posn.pobs.mix <- as.vector(extra$indeta[ 2, 'launch'])
      posn.shape.a <- as.vector(extra$indeta[ 3, 'launch'])
      posn.pstr.mix <- as.vector(extra$indeta[ 4, 'launch'])
      posn.shape.i <- as.vector(extra$indeta[ 5, 'launch'])
      posn.pdip.mix <- as.vector(extra$indeta[ 6, 'launch'])
      posn.shape.d <- as.vector(extra$indeta[ 7, 'launch'])
      posn.pobs.mlm <- as.vector(extra$indeta[ 8, 'launch'])
      posn.pstr.mlm <- as.vector(extra$indeta[ 9, 'launch'])
      posn.pdip.mlm <- as.vector(extra$indeta[10, 'launch'])





    ned2l.dpdip.mix2         <-  # Elt (6, 6)
    ned2l.dpstr.mix2         <-  # Elt (4, 4). Unchanged by deflation.
    ned2l.dpobs.mlm.pstr.mix <-  # Elts (4, >=8). (((09)))
    ned2l.dpobs.mix.pstr.mix <- +probns / Numer^2  # ccc Elt (2, 4)
    if (all(c(la.mix, li.mlm) > 0))  # (((08)))
    ned2l.dpobs.mix.pstr.mlm <- matrix( probns / Numer^2, n, li.mlm)
    if (all(c(li.mix, li.mlm) > 0))  # (((10)))
    ned2l.dpstr.mix.pstr.mlm <- matrix( probns / Numer^2, n, li.mlm)
    if (all(c(ld.mix, ld.mlm) > 0))  # (((21)))
    ned2l.dpdip.mix.pdip.mlm <- matrix( probns / Numer^2, n, ld.mlm)


    ned2l.dpobs.mlm.pdip.mix <-  # Elts (6, >=8). (((19)))
    ned2l.dpstr.mix.pdip.mix <-  # Elt (4, 6)
    ned2l.dpobs.mix.pdip.mix <- -probns / Numer^2  # ccc Elt (2, 6)
    if (all(c(la.mix, ld.mlm) > 0))  # (((17)))
    ned2l.dpobs.mix.pdip.mlm <- matrix(-probns / Numer^2, n, ld.mlm)
    if (all(c(li.mix, ld.mlm) > 0))  # (((18)))
    ned2l.dpstr.mix.pdip.mlm <- matrix(-probns / Numer^2, n, ld.mlm)
    if (all(c(ld.mix, li.mlm) > 0))  # (((20)))
    ned2l.dpdip.mix.pstr.mlm <- matrix(-probns / Numer^2, n, li.mlm)


    


    ned2l.dshape.p2 <- probns * (
      zeta(shape.p + 1, deriv = 2) / fred0.p -
      (fred1.p / fred0.p)^2 +  # ccc
      Denom2.p / Denom0.p - (Denom1.p / Denom0.p)^2) + 
    (if (tmp3.TF[ 4] && li.mix) Numer *
    rowSums(Numer * (d1B.PI.mix^2) / DELTA.i.mix - d2B.PI.mix) else 0) +
    (if (tmp3.TF[ 9] && li.mlm) Numer *
    rowSums(Numer * (d1B.PI.mlm^2) / DELTA.i.mlm - d2B.PI.mlm) else 0) +
    (if (tmp3.TF[ 6] && ld.mix) Numer *
    rowSums(Numer * (d1B.PD.mix^2) / DELTA.d.mix - d2B.PD.mix) else 0) +
    (if (tmp3.TF[10] && ld.mlm) Numer *  # nnn.
    rowSums(Numer * (d1B.PD.mlm^2) / DELTA.d.mlm - d2B.PD.mlm) else 0)


    wz[, iam(1, 1, M)] <- ned2l.dshape.p2 * dshape.p.deta^2


    ned2l.dpobs.mix2 <- 1 / pobs.mix + probns / Numer^2
    if (tmp3.TF[ 4] && li.mix > 0) {
      ned2l.dpobs.mix2 <-  # More just below, ccc
      ned2l.dpobs.mix2 + rowSums(d0B.PI.mix^2 / DELTA.i.mix)
    }
    if (tmp3.TF[ 9] && li.mlm > 0) {
      ned2l.dpobs.mix2 <-  # ccc.
      ned2l.dpobs.mix2 + rowSums(d0B.PI.mlm^2 / DELTA.i.mlm)
    }
    if (tmp3.TF[ 6] && ld.mix > 0) {
      ned2l.dpobs.mix2 <-  # nnn
      ned2l.dpobs.mix2 + rowSums(d0B.PD.mix^2 / DELTA.d.mix)
    }
    if (tmp3.TF[10] && ld.mlm > 0) {
      ned2l.dpobs.mix2 <-  # nnn
      ned2l.dpobs.mix2 + rowSums(d0B.PD.mlm^2 / DELTA.d.mlm)
    }
    if (tmp3.TF[ 2] && la.mix > 0)
      wz[, iam(2, 2, M)] <- ned2l.dpobs.mix2  # Link done later


    if (tmp3.TF[ 3] && la.mix > 1) {
      ned2l.dshape.a2 <- pobs.mix * (
        rowSums((DA.mix.1mat.a^2) / DA.mix.0mat.a) / Denom0.a -
        (Denom1.a / Denom0.a)^2)  # ccc.
      wz[, iam(3, 3, M)] <- ned2l.dshape.a2 * dshape.a.deta^2
    }




    if (tmp3.TF[ 4] && li.mix > 0) {
      ned2l.dpstr.mix2 <-
      ned2l.dpstr.mix2 +
        rowSums((d0A.i - d0B.PI.mix)^2 / DELTA.i.mix)

      if (tmp3.TF[ 2] && la.mix > 0)
        ned2l.dpobs.mix.shape.p <-
        ned2l.dpobs.mix.shape.p +
          rowSums(d1B.PI.mix * (1 - Numer * d0B.PI.mix / DELTA.i.mix))

      ned2l.dpstr.mix.shape.p <-
      ned2l.dpstr.mix.shape.p + rowSums(
        d1B.PI.mix * (1 + Numer * (d0A.i - d0B.PI.mix) / DELTA.i.mix))

      if (tmp3.TF[ 6])
        ned2l.dpdip.mix.shape.p <-
        ned2l.dpdip.mix.shape.p - rowSums(
          d1B.PI.mix * (1 - Numer * d0B.PI.mix / DELTA.i.mix))

      if (all(tmp3.TF[c(2, 4)]))
        ned2l.dpobs.mix.pstr.mix <-  # ccc
        ned2l.dpobs.mix.pstr.mix +
          rowSums(-d0B.PI.mix * (d0A.i - d0B.PI.mix) / DELTA.i.mix)

      if (all(tmp3.TF[c(4, 6)]))
        ned2l.dpstr.mix.pdip.mix <-
        ned2l.dpstr.mix.pdip.mix + rowSums(
          d0B.PI.mix * (d0A.i - d0B.PI.mix) / DELTA.i.mix)

      if (!is.na(posn.pdip.mix)) {
        ned2l.dpdip.mix2 <-
        ned2l.dpdip.mix2 + rowSums(d0B.PI.mix^2 / DELTA.i.mix)
      }
    }  # (tmp3.TF[ 4] && li.mix > 0)







    if (all(tmp3.TF[c(2, 4,  9)])) {  # was la.mix > 0 & DELTA.i.mix
      ned2l.dpobs.mix.pstr.mix <-  # ccc
      ned2l.dpobs.mix.pstr.mix + rowSums(d0B.PI.mlm^2 / DELTA.i.mlm)
    }
    if (all(tmp3.TF[c(2, 4,  6)])) {  # ==  ld.mix > 0 & DELTA.d.mix
      ned2l.dpobs.mix.pstr.mix <-  # nnn
      ned2l.dpobs.mix.pstr.mix + rowSums(d0B.PD.mix^2 / DELTA.d.mix)
    }
    if (all(tmp3.TF[c(2, 4, 10)])) {  # ==  ld.mlm > 0 & DELTA.d.mlm
      ned2l.dpobs.mix.pstr.mix <-  # nnn.
      ned2l.dpobs.mix.pstr.mix + rowSums(d0B.PD.mlm^2 / DELTA.d.mlm)
    }
    if (!is.na(posn.pobs.mix) && !is.na(posn.pstr.mix))
      wz[, iam(posn.pobs.mix, posn.pstr.mix, M)] <-
        ned2l.dpobs.mix.pstr.mix  # Link done later



    if (all(tmp3.TF[c(2, 6)]))
      ned2l.dpobs.mix.pdip.mix <-  # nnn
      ned2l.dpobs.mix.pdip.mix +
        rowSums( d0B.PD.mix * (d0A.d - d0B.PD.mix) / DELTA.d.mix)

    if (all(tmp3.TF[c(2, 6,  9)])) {  # ==  li.mlm > 0 & DELTA.i.mix
      ned2l.dpobs.mix.pdip.mix <-  # nnn
      ned2l.dpobs.mix.pdip.mix - rowSums(d0B.PI.mlm^2 / DELTA.i.mlm)
    }
    if (all(tmp3.TF[c(2, 6,  4)])) {  # ==  li.mix > 0 & DELTA.i.mix
      ned2l.dpobs.mix.pdip.mix <-  # nnn
      ned2l.dpobs.mix.pdip.mix - rowSums(d0B.PI.mix^2 / DELTA.i.mix)
    }
    if (all(tmp3.TF[c(2, 6, 10)])) {  # ==  ld.mlm > 0 & DELTA.d.mlm
      ned2l.dpobs.mix.pdip.mix <-  # nnn.
      ned2l.dpobs.mix.pdip.mix - rowSums(d0B.PD.mlm^2 / DELTA.d.mlm)
    }
    if (!is.na(posn.pobs.mix) && !is.na(posn.pdip.mix))
      wz[, iam(posn.pobs.mix, posn.pdip.mix, M)] <-
        ned2l.dpobs.mix.pdip.mix  # Link done later





    if (tmp3.TF[ 5] && li.mix > 1) {  # \calI_{p}, includes \theta_i.

      ned2l.dshape.p.shape.i <- pstr.mix * Numer *
        rowSums(d1A.i * d1B.PI.mix / DELTA.i.mix)  # ccc.
      wz[, iam(1, posn.shape.i, M)] <- ned2l.dshape.p.shape.i *
        dshape.p.deta * dshape.i.deta  # All links done here

      ned2l.dshape.i2 <- pstr.mix *
        rowSums(pstr.mix * (d1A.i^2) / DELTA.i.mix - d2A.i)  # ccc.
      wz[, iam(posn.shape.i, posn.shape.i, M)] <-
        ned2l.dshape.i2 * dshape.i.deta^2

      if (tmp3.TF[ 2]) {  # tmp3.TF[ 4] is TRUE, given tmp3.TF[ 5]
        ned2l.dpobs.mix.shape.i <-
          rowSums(-pstr.mix * d1A.i * d0B.PI.mix / DELTA.i.mix)  # ccc.
        wz[, iam(posn.pobs.mix, posn.shape.i, M)] <-
          ned2l.dpobs.mix.shape.i  # * dshape.i.deta done later
      }

      if (tmp3.TF[ 4]) {
        ned2l.dpstr.mix.shape.i <- rowSums(  # ccc.
          d1A.i * (pstr.mix * (d0A.i - d0B.PI.mix) / DELTA.i.mix - 1))
        wz[, iam(posn.pstr.mix, posn.shape.i, M)] <-
          ned2l.dpstr.mix.shape.i  # * dshape.i.deta done later
      }

      if (all(tmp3.TF[c(5, 6)])) {
        ned2l.dpdip.mix.shape.i <- rowSums(
          (-pstr.mix) * d0B.PI.mix * d1A.i / DELTA.i.mix)
        wz[, iam(posn.pdip.mix, posn.shape.i, M)] <-
          ned2l.dpdip.mix.shape.i  # link done later
      }

      if (tmp3.TF[ 8]) {
        ned2l.dpobs.mlm.shape.i <- rowSums(
          -pstr.mix * d0B.PI.mix * d1A.i / DELTA.i.mix)  # ccc.
        for (uuu in seq(la.mlm))
          wz[, iam(posn.pobs.mlm - 1 + uuu, posn.shape.i, M)] <-
            ned2l.dpobs.mlm.shape.i  # * dshape.i.deta done later
      }


    }  # (tmp3.TF[ 5] && li.mix > 1)





    if (tmp3.TF[ 6] && ld.mix > 0) {  # \calD_{p}, maybe w. \theta_d

      if (tmp3.TF[ 2] && la.mix > 0)
        ned2l.dpobs.mix.shape.p <-
        ned2l.dpobs.mix.shape.p +
          rowSums(d1B.PD.mix * (1 - Numer * d0B.PD.mix / DELTA.d.mix))

      ned2l.dpstr.mix.shape.p <-
      ned2l.dpstr.mix.shape.p + rowSums(
        d1B.PD.mix * (1 - Numer * d0B.PD.mix / DELTA.d.mix))

      ned2l.dpdip.mix.shape.p <-
      ned2l.dpdip.mix.shape.p - rowSums(
        d1B.PD.mix * (1 + Numer * (d0A.d - d0B.PD.mix) / DELTA.d.mix))

      if (!is.na(posn.pstr.mix)) {
        ned2l.dpstr.mix2 <-
        ned2l.dpstr.mix2 + rowSums(d0B.PD.mix^2 / DELTA.d.mix)
      }

      if (all(tmp3.TF[c(4, 6)]))
        ned2l.dpstr.mix.pdip.mix <-
        ned2l.dpstr.mix.pdip.mix + rowSums(
          d0B.PD.mix * (d0A.d - d0B.PD.mix) / DELTA.d.mix)

      ned2l.dpdip.mix2 <-
      ned2l.dpdip.mix2 +
        rowSums((d0A.d - d0B.PD.mix)^2 / DELTA.d.mix)

    }  # (tmp3.TF[ 6] && ld.mix > 0)





    if (tmp3.TF[ 7] && ld.mix > 1) {  # \calD_{p}, includes \theta_d


      ned2l.dshape.p.shape.d <- (-pdip.mix) * Numer *
        rowSums(d1A.d * d1B.PD.mix / DELTA.d.mix)  # nnn.
      wz[, iam(1, posn.shape.d, M)] <- ned2l.dshape.p.shape.d *
        dshape.p.deta * dshape.d.deta  # All links done here

      if (tmp3.TF[ 2]) {  # tmp3.TF[ 6] is TRUE, given tmp3.TF[ 7]
        ned2l.dpobs.mix.shape.d <-
          rowSums(pdip.mix * d1A.d * d0B.PD.mix / DELTA.d.mix)  # nnn.
        wz[, iam(posn.pobs.mix, posn.shape.d, M)] <-
          ned2l.dpobs.mix.shape.d  # link done later
      }

      if (tmp3.TF[ 4]) {
        ned2l.dpstr.mix.shape.d <- rowSums(
          pdip.mix * d1A.d * d0B.PD.mix / DELTA.d.mix)
        wz[, iam(posn.pstr.mix, posn.shape.d, M)] <-
          ned2l.dpstr.mix.shape.d  # * dshape.i.deta done later
      }

        ned2l.dpdip.mix.shape.d <- rowSums(
          d1A.d * (1 + pdip.mix * (d0A.d - d0B.PD.mix) / DELTA.d.mix))
        wz[, iam(posn.pdip.mix, posn.shape.d, M)] <-
          ned2l.dpdip.mix.shape.d  # * dshape.d.deta done later

      ned2l.dshape.d2 <- pdip.mix *
        rowSums(pdip.mix * (d1A.d^2) / DELTA.d.mix + d2A.d)  # nnn.
      wz[, iam(posn.shape.d, posn.shape.d, M)] <-
        ned2l.dshape.d2 * dshape.d.deta^2

      if (tmp3.TF[ 8]) {
        ned2l.dpobs.mlm.shape.d <- rowSums(
           pdip.mix * d0B.PD.mix * d1A.d / DELTA.d.mix)  # nnn.
        for (uuu in seq(la.mlm))
          wz[, iam(posn.pobs.mlm - 1 + uuu, posn.shape.d, M)] <-
            ned2l.dpobs.mlm.shape.d  # * dshape.d.deta done later
      }

    }  # (tmp3.TF[ 7] && ld.mix > 1)




        
    if (tmp3.TF[ 9] && li.mlm > 0) {  # \calI_{np}, includes \phi_s.

      if (la.mix && tmp3.TF[ 2])
        ned2l.dpobs.mix.shape.p <-  # ccc
        ned2l.dpobs.mix.shape.p +
          rowSums(d1B.PI.mlm * (1 - Numer * d0B.PI.mlm / DELTA.i.mlm))

      ned2l.dpstr.mix.shape.p <-  # ccc.
      ned2l.dpstr.mix.shape.p + rowSums(
        d1B.PI.mlm * (1 - Numer * d0B.PI.mlm / DELTA.i.mlm))

      if (tmp3.TF[ 6])
        ned2l.dpdip.mix.shape.p <-
        ned2l.dpdip.mix.shape.p - rowSums(
          d1B.PI.mlm * (1 - Numer * d0B.PI.mlm / DELTA.i.mlm))

      if (!is.na(posn.pstr.mix)) {
        ned2l.dpstr.mix2 <-
        ned2l.dpstr.mix2 + rowSums(d0B.PI.mlm^2 / DELTA.i.mlm)
      }

      if (all(tmp3.TF[c(4, 6)]))
        ned2l.dpstr.mix.pdip.mix <-
        ned2l.dpstr.mix.pdip.mix - rowSums(d0B.PI.mlm^2 / DELTA.i.mlm)

      if (!is.na(posn.pdip.mix)) {
        ned2l.dpdip.mix2 <-
        ned2l.dpdip.mix2 + rowSums(d0B.PI.mlm^2 / DELTA.i.mlm)
      }

    }  # tmp3.TF[ 9] && li.mlm > 0



        
    if (tmp3.TF[10] && ld.mlm > 0) {  # \calD_{np}, includes \psi_s.


      if (la.mix && tmp3.TF[ 2])
        ned2l.dpobs.mix.shape.p <-  # nnn.
        ned2l.dpobs.mix.shape.p +
          rowSums(d1B.PD.mlm * (1 - Numer * d0B.PD.mlm / DELTA.d.mlm))

      ned2l.dpstr.mix.shape.p <-  # nnn.
      ned2l.dpstr.mix.shape.p + rowSums(
        d1B.PD.mlm * (1 - Numer * d0B.PD.mlm / DELTA.d.mlm))

      if (tmp3.TF[ 6])
        ned2l.dpdip.mix.shape.p <-
        ned2l.dpdip.mix.shape.p - rowSums(
          d1B.PD.mlm * (1 - Numer * d0B.PD.mlm / DELTA.d.mlm))

      if (!is.na(posn.pstr.mix)) {
        ned2l.dpstr.mix2 <-
        ned2l.dpstr.mix2 + rowSums(d0B.PD.mlm^2 / DELTA.d.mlm)
      }

      if (all(tmp3.TF[c(4, 6)]))
        ned2l.dpstr.mix.pdip.mix <-
        ned2l.dpstr.mix.pdip.mix - rowSums(d0B.PD.mlm^2 / DELTA.d.mlm)

      if (!is.na(posn.pdip.mix)) {
        ned2l.dpdip.mix2 <-
        ned2l.dpdip.mix2 + rowSums(d0B.PD.mlm^2 / DELTA.d.mlm)
      }



    }  # tmp3.TF[10] && ld.mlm > 0






      




    if (!is.na(posn.pobs.mix))  # Optional (1, 2) element:
      wz[, iam(1, posn.pobs.mix, M)] <-
        ned2l.dpobs.mix.shape.p  # One link done later
    if (!is.na(posn.pstr.mix))  # Optional (1, 4) element
      wz[, iam(1, posn.pstr.mix, M)] <-
        ned2l.dpstr.mix.shape.p  # One link done later
    if (!is.na(posn.pdip.mix))  # Optional (1, 6) element
      wz[, iam(1, posn.pdip.mix, M)] <-
        ned2l.dpdip.mix.shape.p  # One link done later
    if (!is.na(posn.pstr.mix) &&
        !is.na(posn.pdip.mix))  # Optional (4, 6) element
      wz[, iam(posn.pstr.mix, posn.pdip.mix, M)] <-
        ned2l.dpstr.mix.pdip.mix  # Links done later zz1

    if (!is.na(posn.pstr.mix))  # Optional (4, 4) element
      wz[, iam(posn.pstr.mix,  # Link done later
               posn.pstr.mix, M)] <- ned2l.dpstr.mix2

    if (!is.na(posn.pdip.mix))  # Optional (6, 6) element
      wz[, iam(posn.pdip.mix,  # Link done later
               posn.pdip.mix, M)] <- ned2l.dpdip.mix2







    if (tmp3.TF[ 8] && la.mlm) {  # \calA_{np}, includes \omega_s
      ofset <- posn.pobs.mlm - 1  # 7 for GAITD combo
      for (uuu in seq(la.mlm)) {  # Diagonal elts only
        wz[, iam(ofset + uuu,
                 ofset + uuu, M)] <- 1 / pobs.mlm[, uuu]
      }  # uuu

      tmp8a <- probns / Numer^2
      if (tmp3.TF[ 4] && li.mix)
        tmp8a <- tmp8a + rowSums((d0B.PI.mix^2) / DELTA.i.mix)
      if (tmp3.TF[ 9] && li.mlm)
        tmp8a <- tmp8a + rowSums((d0B.PI.mlm^2) / DELTA.i.mlm)
      if (tmp3.TF[ 6] && ld.mix)
        tmp8a <- tmp8a + rowSums((d0B.PD.mix^2) / DELTA.d.mix)
      if (tmp3.TF[10] && ld.mlm)
        tmp8a <- tmp8a + rowSums((d0B.PD.mlm^2) / DELTA.d.mlm)
      for (uuu in seq(la.mlm))  # All elts
        for (vvv in uuu:la.mlm)
          wz[, iam(ofset + uuu, ofset + vvv, M)] <-
          wz[, iam(ofset + uuu, ofset + vvv, M)] + tmp8a  # All elts
    }  # la.mlm


 

    if (tmp3.TF[ 8] && la.mlm) {

      init0.i.val <- init0.d.val <- 0
      if (tmp3.TF[ 9] && li.mlm) init0.i.val <-
        rowSums(d1B.PI.mlm * (1 - Numer * d0B.PI.mlm / DELTA.i.mlm))
      if (tmp3.TF[10] && ld.mlm) init0.d.val <-
        rowSums(d1B.PD.mlm * (1 - Numer * d0B.PD.mlm / DELTA.d.mlm))
      ned2l.dpobs.mlm.shape.p <- init0.i.val + init0.d.val  # Vector

      if (tmp3.TF[ 4] && li.mix)
        ned2l.dpobs.mlm.shape.p <-
        ned2l.dpobs.mlm.shape.p + rowSums(
          d1B.PI.mix * (1 - Numer * d0B.PI.mix / DELTA.i.mix))
      if (tmp3.TF[ 6] && ld.mix)
        ned2l.dpobs.mlm.shape.p <-
        ned2l.dpobs.mlm.shape.p + rowSums(  # nnn
          d1B.PD.mix * (1 - Numer * d0B.PD.mix / DELTA.d.mix))

      ofset <- posn.pobs.mlm - 1  # 5 for combo
      for (vvv in seq(la.mlm))  # ccc.
        wz[, iam(1, ofset + vvv, M)] <- ned2l.dpobs.mlm.shape.p
    }  # la.mlm > 0





    if (tmp3.TF[ 9] && li.mlm > 0) {  # \calI_{np}, includes \phi_s
      init0.val <- probns / Numer^2
      if (li.mix)
        init0.val <- init0.val + rowSums((d0B.PI.mix^2) / DELTA.i.mix)
      if (ld.mix)  # nnn
        init0.val <- init0.val + rowSums((d0B.PD.mix^2) / DELTA.d.mix)
      if (ld.mlm)  # nnn
        init0.val <- init0.val + rowSums((d0B.PD.mlm^2) / DELTA.d.mlm)
      ned2l.dpstr.mlm2 <-
        matrix(init0.val, n, li.mlm * (li.mlm + 1) / 2)
      for (uuu in seq(li.mlm))
        for (sss in seq(li.mlm))
          ned2l.dpstr.mlm2[, iam(uuu, uuu, li.mlm)] <-
          ned2l.dpstr.mlm2[, iam(uuu, uuu, li.mlm)] +
            ((sss == uuu) - d0B.PI.mlm[, sss])^2 / DELTA.i.mlm[, sss]
      if (li.mlm > 1) {
        for (uuu in seq(li.mlm - 1))
          for (vvv in (uuu + 1):li.mlm)
            for (sss in seq(li.mlm))
              ned2l.dpstr.mlm2[, iam(uuu, vvv, li.mlm)] <-
              ned2l.dpstr.mlm2[, iam(uuu, vvv, li.mlm)] +
              ((sss == uuu) - d0B.PI.mlm[, sss]) *
              ((sss == vvv) - d0B.PI.mlm[, sss]) / DELTA.i.mlm[, sss]
      }  # if (li.mlm > 1)

      ofset <- posn.pstr.mlm - 1
      for (uuu in seq(li.mlm))
        for (vvv in uuu:li.mlm)
          wz[, iam(ofset + uuu, ofset + vvv, M)] <-
            ned2l.dpstr.mlm2[, iam(uuu, vvv, li.mlm)] 
    }  # li.mlm > 0






    if (tmp3.TF[10] && ld.mlm > 0) {  # \calD_{np}, includes \psi_s
      init0.val <- probns / Numer^2
      if (ld.mix)
        init0.val <- init0.val + rowSums((d0B.PD.mix^2) / DELTA.d.mix)
      if (li.mix)
        init0.val <- init0.val + rowSums((d0B.PI.mix^2) / DELTA.i.mix)
      if (li.mlm)
        init0.val <- init0.val + rowSums((d0B.PI.mlm^2) / DELTA.i.mlm)
      ned2l.dpdip.mlm2 <-
        matrix(init0.val, n, ld.mlm * (ld.mlm + 1) / 2)
      for (uuu in seq(ld.mlm))
        for (sss in seq(ld.mlm))
          ned2l.dpdip.mlm2[, iam(uuu, uuu, ld.mlm)] <-
          ned2l.dpdip.mlm2[, iam(uuu, uuu, ld.mlm)] +
            (d0B.PD.mlm[, sss] - (sss == uuu))^2 / DELTA.d.mlm[, sss]
      if (ld.mlm > 1) {
        for (uuu in seq(ld.mlm - 1))
          for (vvv in (uuu + 1):ld.mlm)
            for (sss in seq(ld.mlm))
              ned2l.dpdip.mlm2[, iam(uuu, vvv, ld.mlm)] <-
              ned2l.dpdip.mlm2[, iam(uuu, vvv, ld.mlm)] +
              (d0B.PD.mlm[, sss] - (sss == uuu)) *
              (d0B.PD.mlm[, sss] - (sss == vvv)) / DELTA.d.mlm[, sss]
      }  # if (ld.mlm > 1)

      ofset <- posn.pdip.mlm - 1
      for (uuu in seq(ld.mlm))
        for (vvv in uuu:ld.mlm)
          wz[, iam(ofset + uuu, ofset + vvv, M)] <-
            ned2l.dpdip.mlm2[, iam(uuu, vvv, ld.mlm)] 
    }  # ld.mlm > 0









    if (tmp3.TF[ 9] && li.mlm > 0) {
      ned2l.dpstr.mlm.theta.p <- matrix(0, n, li.mlm)
      for (vvv in seq(li.mlm))
        for (sss in seq(li.mlm))
          ned2l.dpstr.mlm.theta.p[, vvv] <-
          ned2l.dpstr.mlm.theta.p[, vvv] +
          d1B.PI.mlm[, sss] * (1 + Numer *
          (max(0, sss == vvv) - d0B.PI.mlm[, sss]) / (
          DELTA.i.mlm[, sss]))
      if (li.mix && tmp3.TF[ 4])
        ned2l.dpstr.mlm.theta.p <-
        ned2l.dpstr.mlm.theta.p +
        rowSums(d1B.PI.mix * (1 - Numer * d0B.PI.mix / DELTA.i.mix))
      if (ld.mix && tmp3.TF[ 6])
        ned2l.dpstr.mlm.theta.p <-  # nnn
        ned2l.dpstr.mlm.theta.p +
        rowSums(d1B.PD.mix * (1 - Numer * d0B.PD.mix / DELTA.d.mix))
      if (ld.mlm && tmp3.TF[10])
        ned2l.dpstr.mlm.theta.p <-  # nnn.
        ned2l.dpstr.mlm.theta.p +
        rowSums(d1B.PD.mlm * (1 - Numer * d0B.PD.mlm / DELTA.d.mlm))
      ofset <- posn.pstr.mlm - 1
      for (vvv in seq(li.mlm))  # ccc.
        wz[, iam(1, ofset + vvv, M)] <- ned2l.dpstr.mlm.theta.p[, vvv]
    }  # li.mlm > 0




    if (tmp3.TF[10] && ld.mlm > 0) {
      ned2l.dpdip.mlm.theta.p <- matrix(0, n, ld.mlm)
      for (vvv in seq(ld.mlm))
        for (sss in seq(ld.mlm))
          ned2l.dpdip.mlm.theta.p[, vvv] <-
          ned2l.dpdip.mlm.theta.p[, vvv] -  # Minus
          d1B.PD.mlm[, sss] * (1 + Numer *
          (max(0, sss == vvv) - d0B.PD.mlm[, sss]) / (
          DELTA.d.mlm[, sss]))
      if (ld.mix && tmp3.TF[ 6])
        ned2l.dpdip.mlm.theta.p <-
        ned2l.dpdip.mlm.theta.p -  # Minus
        rowSums(d1B.PD.mix * (1 - Numer * d0B.PD.mix / DELTA.d.mix))
      if (li.mix && tmp3.TF[ 4])
        ned2l.dpdip.mlm.theta.p <-
        ned2l.dpdip.mlm.theta.p -  # Minus
        rowSums(d1B.PI.mix * (1 - Numer * d0B.PI.mix / DELTA.i.mix))
      if (li.mlm && tmp3.TF[ 9])
        ned2l.dpdip.mlm.theta.p <-  # nnn.
        ned2l.dpdip.mlm.theta.p -  # Minus
        rowSums(d1B.PI.mlm * (1 - Numer * d0B.PI.mlm / DELTA.i.mlm))
      ofset <- posn.pdip.mlm - 1
      for (vvv in seq(ld.mlm))  # nnn.
        wz[, iam(1, ofset + vvv, M)] <- ned2l.dpdip.mlm.theta.p[, vvv]
    }  # ld.mlm > 0








    if (li.mlm && li.mix > 1) {

      ned2l.dpstr.mlm.theta.i <-  # Not a matrix, just a vector
        rowSums(-pstr.mix * d0B.PI.mix * d1A.i / DELTA.i.mix)

      for (vvv in seq(li.mlm))
        wz[, iam(posn.shape.i, posn.pstr.mlm - 1 + vvv, M)] <-
          ned2l.dpstr.mlm.theta.i  # ccc.
    }  # li.mlm && li.mix > 1





    if (ld.mlm && ld.mix > 1) {

      ned2l.dpdip.mlm.theta.d <-  # Not a matrix, just a vector
        rowSums(pdip.mix * d0B.PD.mix * d1A.d / DELTA.d.mix)

      for (vvv in seq(ld.mlm))
        wz[, iam(posn.shape.d, posn.pdip.mlm - 1 + vvv, M)] <-
          ned2l.dpdip.mlm.theta.d  # nnn.
    }  # ld.mlm && ld.mix > 1



    if (ld.mlm && li.mix > 1) {

      ned2l.dpdip.mlm.theta.i <-  # Not a matrix, just a vector
        rowSums(-pstr.mix * d0B.PI.mix * d1A.i / DELTA.i.mix)

      for (vvv in seq(ld.mlm))
        wz[, iam(posn.shape.i, posn.pdip.mlm - 1 + vvv, M)] <-
          ned2l.dpdip.mlm.theta.i  # nnn.
    }  # ld.mlm && li.mix > 1



   


    if (li.mlm && ld.mix > 1) {

      ned2l.dpstr.mlm.theta.d <-  # Not a matrix, just a vector
        rowSums(pdip.mix * d0B.PD.mix * d1A.d / DELTA.d.mix)

      for (vvv in seq(li.mlm))
        wz[, iam(posn.shape.d, posn.pstr.mlm - 1 + vvv, M)] <-
          ned2l.dpstr.mlm.theta.d  # nnn.
    }  # li.mlm && ld.mix > 1






    if (all(c(la.mlm, li.mlm) > 0)) {
      ned2l.dpobs.mlm.pstr.mlm <-
        array(probns / Numer^2, c(n, la.mlm, li.mlm))
      for (uuu in seq(la.mlm))
        for (vvv in seq(li.mlm))
          for (sss in seq(li.mlm))
            ned2l.dpobs.mlm.pstr.mlm[, uuu, vvv] <- 
            ned2l.dpobs.mlm.pstr.mlm[, uuu, vvv] - d0B.PI.mlm[, sss] *
              ((sss == vvv) - d0B.PI.mlm[, sss]) / DELTA.i.mlm[, sss]
      if (tmp3.TF[ 4] && li.mix)
        ned2l.dpobs.mlm.pstr.mlm <-
        ned2l.dpobs.mlm.pstr.mlm + rowSums(d0B.PI.mix^2 / DELTA.i.mix)
      if (tmp3.TF[ 6] && ld.mix)
        ned2l.dpobs.mlm.pstr.mlm <-  # nnn
        ned2l.dpobs.mlm.pstr.mlm + rowSums(d0B.PD.mix^2 / DELTA.d.mix)
      if (tmp3.TF[10] && ld.mlm)
        ned2l.dpobs.mlm.pstr.mlm <-  # nnn
        ned2l.dpobs.mlm.pstr.mlm + rowSums(d0B.PD.mlm^2 / DELTA.d.mlm)
      ofset.pobs <- posn.pobs.mlm - 1
      ofset.pstr <- posn.pstr.mlm - 1
      for (uuu in seq(la.mlm))
        for (vvv in seq(li.mlm))
          wz[, iam(ofset.pobs + uuu, ofset.pstr + vvv, M)] <-
            ned2l.dpobs.mlm.pstr.mlm[, uuu, vvv] 
    }  # all(c(la.mlm, li.mlm) > 0)







    if (all(c(li.mlm, ld.mlm) > 0)) {
      ned2l.dpstr.mlm.pdip.mlm <-
        array(-probns / Numer^2, c(n, li.mlm, ld.mlm))
      for (uuu in seq(li.mlm))
        for (vvv in seq(ld.mlm))
          for (sss in seq(li.mlm))
            ned2l.dpstr.mlm.pdip.mlm[, uuu, vvv] <- 
            ned2l.dpstr.mlm.pdip.mlm[, uuu, vvv] + d0B.PI.mlm[, sss] *
              ((sss == uuu) - d0B.PI.mlm[, sss]) / DELTA.i.mlm[, sss]
      for (uuu in seq(li.mlm))
        for (vvv in seq(ld.mlm))
          for (sss in seq(ld.mlm))
            ned2l.dpstr.mlm.pdip.mlm[, uuu, vvv] <- 
            ned2l.dpstr.mlm.pdip.mlm[, uuu, vvv] + d0B.PD.mlm[, sss] *
              ((sss == vvv) - d0B.PD.mlm[, sss]) / DELTA.d.mlm[, sss]
      if (tmp3.TF[ 4] && li.mix)
        ned2l.dpstr.mlm.pdip.mlm <-
        ned2l.dpstr.mlm.pdip.mlm - rowSums(d0B.PI.mix^2 / DELTA.i.mix)
      if (tmp3.TF[ 6] && ld.mix)
        ned2l.dpstr.mlm.pdip.mlm <-  # nnn.
        ned2l.dpstr.mlm.pdip.mlm - rowSums(d0B.PD.mix^2 / DELTA.d.mix)
      ofset.pstr <- posn.pstr.mlm - 1
      ofset.pdip <- posn.pdip.mlm - 1
      for (uuu in seq(li.mlm))
        for (vvv in seq(ld.mlm))
          wz[, iam(ofset.pstr + uuu, ofset.pdip + vvv, M)] <-
            ned2l.dpstr.mlm.pdip.mlm[, uuu, vvv] 
    }  # all(c(li.mlm, ld.mlm) > 0)







    if (all(c(la.mlm, ld.mlm) > 0)) {
      ned2l.dpobs.mlm.pdip.mlm <-
        array(-probns / Numer^2, c(n, la.mlm, ld.mlm))
      for (uuu in seq(la.mlm))
        for (vvv in seq(ld.mlm))
          for (sss in seq(ld.mlm))
            ned2l.dpobs.mlm.pdip.mlm[, uuu, vvv] <- 
            ned2l.dpobs.mlm.pdip.mlm[, uuu, vvv] + d0B.PD.mlm[, sss] *
              ((sss == vvv) - d0B.PD.mlm[, sss]) / DELTA.d.mlm[, sss]
      if (tmp3.TF[ 4] && li.mix)
        ned2l.dpobs.mlm.pdip.mlm <-
        ned2l.dpobs.mlm.pdip.mlm - rowSums(d0B.PI.mix^2 / DELTA.i.mix)
      if (tmp3.TF[ 9] && li.mlm)
        ned2l.dpobs.mlm.pdip.mlm <-
        ned2l.dpobs.mlm.pdip.mlm - rowSums(d0B.PI.mlm^2 / DELTA.i.mlm)
      if (tmp3.TF[ 6] && ld.mix)
        ned2l.dpobs.mlm.pdip.mlm <-
        ned2l.dpobs.mlm.pdip.mlm - rowSums(d0B.PD.mix^2 / DELTA.d.mix)
      ofset.pobs <- posn.pobs.mlm - 1
      ofset.pdip <- posn.pdip.mlm - 1
      for (uuu in seq(la.mlm))
        for (vvv in seq(ld.mlm))
          wz[, iam(ofset.pobs + uuu, ofset.pdip + vvv, M)] <-
            ned2l.dpobs.mlm.pdip.mlm[, uuu, vvv] 
    }  # all(c(la.mlm, li.mlm) > 0)






    if (all(c(la.mix, la.mlm) > 0)) {
      ned2l.dpobs.mix.pobs.mlm <- probns / Numer^2  # Initialize
      if (li.mix)  # tmp3.TF[ 4]
        ned2l.dpobs.mix.pobs.mlm <-
        ned2l.dpobs.mix.pobs.mlm + rowSums(d0B.PI.mix^2 / DELTA.i.mix)
      if (li.mlm)  # tmp3.TF[ 7]
        ned2l.dpobs.mix.pobs.mlm <-
        ned2l.dpobs.mix.pobs.mlm + rowSums(d0B.PI.mlm^2 / DELTA.i.mlm)
      if (ld.mix)  # tmp3.TF[ 6]   nnn
        ned2l.dpobs.mix.pobs.mlm <-
        ned2l.dpobs.mix.pobs.mlm + rowSums(d0B.PD.mix^2 / DELTA.d.mix)
      if (ld.mlm)  # tmp3.TF[10]   nnn
        ned2l.dpobs.mix.pobs.mlm <-
        ned2l.dpobs.mix.pobs.mlm + rowSums(d0B.PD.mlm^2 / DELTA.d.mlm)

      for (uuu in seq(la.mlm))  # ccc.
        wz[, iam(posn.pobs.mix, posn.pobs.mlm - 1 + uuu, M)] <-
          ned2l.dpobs.mix.pobs.mlm  # Link done later
    }



    if (all(c(la.mix, li.mlm) > 0)) {  # all(tmp3.TF[c(2, 9)])
      if (li.mix)  # tmp3.TF[ 4]
        ned2l.dpobs.mix.pstr.mlm <-
        ned2l.dpobs.mix.pstr.mlm + rowSums(d0B.PI.mix^2 / DELTA.i.mix)
      if (ld.mix)  # tmp3.TF[ 6]
        ned2l.dpobs.mix.pstr.mlm <-  # nnn
        ned2l.dpobs.mix.pstr.mlm + rowSums(d0B.PD.mix^2 / DELTA.d.mix)
      if (ld.mlm)  # tmp3.TF[10]
        ned2l.dpobs.mix.pstr.mlm <-  # nnn; + is correct, not -
        ned2l.dpobs.mix.pstr.mlm + rowSums(d0B.PD.mlm^2 / DELTA.d.mlm)
      for (uuu in seq(li.mlm))
        for (sss in seq(li.mlm))
          ned2l.dpobs.mix.pstr.mlm[, uuu] <-
          ned2l.dpobs.mix.pstr.mlm[, uuu] -
            ((sss == uuu) - d0B.PI.mlm[, sss]) *
                            d0B.PI.mlm[, sss] / DELTA.i.mlm[, sss]
      for (uuu in seq(li.mlm))  # ccc.
        wz[, iam(posn.pobs.mix,
                 posn.pstr.mlm - 1 + uuu, M)] <-
          ned2l.dpobs.mix.pstr.mlm[, uuu]  # Link done later
    }  # all(c(la.mix, li.mlm) > 0)






    if (all(c(la.mix, ld.mlm) > 0)) {  # all(tmp3.TF[c(2, 10)])
      if (li.mix)  # tmp3.TF[ 4]
        ned2l.dpobs.mix.pdip.mlm <-
        ned2l.dpobs.mix.pdip.mlm - rowSums(d0B.PI.mix^2 / DELTA.i.mix)
      if (li.mlm)  # tmp3.TF[ 9]
        ned2l.dpobs.mix.pdip.mlm <-
        ned2l.dpobs.mix.pdip.mlm - rowSums(d0B.PI.mlm^2 / DELTA.i.mlm)
      if (ld.mix)  # tmp3.TF[ 6]
        ned2l.dpobs.mix.pdip.mlm <-
        ned2l.dpobs.mix.pdip.mlm - rowSums(d0B.PD.mix^2 / DELTA.d.mix)
      for (uuu in seq(ld.mlm))
        for (sss in seq(ld.mlm))
          ned2l.dpobs.mix.pdip.mlm[, uuu] <-
          ned2l.dpobs.mix.pdip.mlm[, uuu] +
            ((sss == uuu) - d0B.PD.mlm[, sss]) *
                            d0B.PD.mlm[, sss] / DELTA.d.mlm[, sss]
      for (uuu in seq(ld.mlm))  # nnn.
        wz[, iam(posn.pobs.mix,
                 posn.pdip.mlm - 1 + uuu, M)] <-
          ned2l.dpobs.mix.pdip.mlm[, uuu]  # Link done later
    }  # all(c(la.mix, ld.mlm) > 0)





    if (all(c(li.mix, la.mlm) > 0)) {  # all(tmp3.TF[c(4, 8)])
      if (li.mlm)  # tmp3.TF[ 9]
        ned2l.dpobs.mlm.pstr.mix <-
        ned2l.dpobs.mlm.pstr.mix + rowSums(d0B.PI.mlm^2 / DELTA.i.mlm)
      if (ld.mix)  # tmp3.TF[ 6]
        ned2l.dpobs.mlm.pstr.mix <-  # nnn
        ned2l.dpobs.mlm.pstr.mix + rowSums(d0B.PD.mix^2 / DELTA.d.mix)
      if (ld.mlm)  # tmp3.TF[10]
        ned2l.dpobs.mlm.pstr.mix <-  # nnn
        ned2l.dpobs.mlm.pstr.mix + rowSums(d0B.PD.mlm^2 / DELTA.d.mlm)
        ned2l.dpobs.mlm.pstr.mix <-  # tmp3.TF[ 4] && li.mix
        ned2l.dpobs.mlm.pstr.mix -
        rowSums((d0A.i - d0B.PI.mix) * d0B.PI.mix / DELTA.i.mix)

      for (uuu in seq(la.mlm))  # ccc.
        wz[, iam(posn.pstr.mix,
                 posn.pobs.mlm - 1 + uuu, M)] <-
          ned2l.dpobs.mlm.pstr.mix  # Link done later
    }  # all(c(li.mix, la.mlm) > 0







    if (all(c(ld.mix, la.mlm) > 0)) {  # all(tmp3.TF[c(6, 8)])
      if (ld.mlm)  # tmp3.TF[10]
        ned2l.dpobs.mlm.pdip.mix <-
        ned2l.dpobs.mlm.pdip.mix - rowSums(d0B.PD.mlm^2 / DELTA.d.mlm)
      if (li.mix)  # tmp3.TF[ 4]
        ned2l.dpobs.mlm.pdip.mix <-
        ned2l.dpobs.mlm.pdip.mix - rowSums(d0B.PI.mix^2 / DELTA.i.mix)
      if (li.mlm)  # tmp3.TF[ 9]
        ned2l.dpobs.mlm.pdip.mix <-
        ned2l.dpobs.mlm.pdip.mix - rowSums(d0B.PI.mlm^2 / DELTA.i.mlm)
        ned2l.dpobs.mlm.pdip.mix <-  # all(tmp3.TF[c(6, 8)]) 
        ned2l.dpobs.mlm.pdip.mix +
        rowSums((d0A.d - d0B.PD.mix) * d0B.PD.mix / DELTA.d.mix)

      for (uuu in seq(la.mlm))  # nnn.
        wz[, iam(posn.pdip.mix,
                 posn.pobs.mlm - 1 + uuu, M)] <-
          ned2l.dpobs.mlm.pdip.mix  # Link done later
    }  # all(c(ld.mix, la.mlm) > 0






    if (all(c(li.mix, li.mlm) > 0)) {  # all(tmp3.TF[c(4, 9)])
      for (uuu in seq(li.mlm))  # tmp3.TF[ 9]
        for (sss in seq(li.mlm))
          ned2l.dpstr.mix.pstr.mlm[, uuu] <-
          ned2l.dpstr.mix.pstr.mlm[, uuu] -
            ((sss == uuu) - d0B.PI.mlm[, sss]) *
              d0B.PI.mlm[, sss] / DELTA.i.mlm[, sss]
      ned2l.dpstr.mix.pstr.mlm <-
      ned2l.dpstr.mix.pstr.mlm -
      rowSums((d0A.i - d0B.PI.mix) * d0B.PI.mix / DELTA.i.mix)
      if (ld.mix)  # tmp3.TF[ 6]
        ned2l.dpstr.mix.pstr.mlm <-  # nnn
        ned2l.dpstr.mix.pstr.mlm + rowSums(d0B.PD.mix^2 / DELTA.d.mix)
      if (ld.mlm)  # tmp3.TF[10]
        ned2l.dpstr.mix.pstr.mlm <-  # nnn
        ned2l.dpstr.mix.pstr.mlm + rowSums(d0B.PD.mlm^2 / DELTA.d.mlm)

      for (uuu in seq(li.mlm))  # Copy it. ccc.
        wz[, iam(posn.pstr.mix,
                 posn.pstr.mlm - 1 + uuu, M)] <-
          ned2l.dpstr.mix.pstr.mlm[, uuu]  # Link done later
    }  # all(c(li.mix, li.mlm) > 0



    if (all(c(ld.mix, ld.mlm) > 0)) {  # all(tmp3.TF[c(6, 10)])
      for (uuu in seq(ld.mlm))  # tmp3.TF[ 9]
        for (sss in seq(ld.mlm))
          ned2l.dpdip.mix.pdip.mlm[, uuu] <-
          ned2l.dpdip.mix.pdip.mlm[, uuu] -
            ((sss == uuu) - d0B.PD.mlm[, sss]) *
              d0B.PD.mlm[, sss] / DELTA.d.mlm[, sss]
      if (ld.mix)  # tmp3.TF[ 6]
        ned2l.dpdip.mix.pdip.mlm <-
        ned2l.dpdip.mix.pdip.mlm -
        rowSums((d0A.d - d0B.PD.mix) * d0B.PD.mix / DELTA.d.mix)
      if (li.mix)  # tmp3.TF[ 4]
        ned2l.dpdip.mix.pdip.mlm <-
        ned2l.dpdip.mix.pdip.mlm + rowSums(d0B.PI.mix^2 / DELTA.i.mix)
      if (li.mlm)  # tmp3.TF[ 9]
        ned2l.dpdip.mix.pdip.mlm <-
        ned2l.dpdip.mix.pdip.mlm + rowSums(d0B.PI.mlm^2 / DELTA.i.mlm)

      for (uuu in seq(ld.mlm))  # Copy it. ccc.
        wz[, iam(posn.pdip.mix,
                 posn.pdip.mlm - 1 + uuu, M)] <-
          ned2l.dpdip.mix.pdip.mlm[, uuu]  # Link done later
    }  # all(c(ld.mix, ld.mlm) > 0









    if (all(c(ld.mix, li.mlm) > 0)) {  # all(tmp3.TF[c(4, 9)])
      for (uuu in seq(li.mlm))  # tmp3.TF[ 9]
        for (sss in seq(li.mlm))
          ned2l.dpdip.mix.pstr.mlm[, uuu] <-
          ned2l.dpdip.mix.pstr.mlm[, uuu] +
            ((sss == uuu) - d0B.PI.mlm[, sss]) *
              d0B.PI.mlm[, sss] / DELTA.i.mlm[, sss]
      if (ld.mix)  # tmp3.TF[ 6]
        ned2l.dpdip.mix.pstr.mlm <-
        ned2l.dpdip.mix.pstr.mlm +
        rowSums((d0A.d - d0B.PD.mix) * d0B.PD.mix / DELTA.d.mix)
      if (li.mix)  # tmp3.TF[ 4]
        ned2l.dpdip.mix.pstr.mlm <-
        ned2l.dpdip.mix.pstr.mlm - rowSums(d0B.PI.mix^2 / DELTA.i.mix)
      if (ld.mlm)  # tmp3.TF[10]
        ned2l.dpdip.mix.pstr.mlm <-
        ned2l.dpdip.mix.pstr.mlm - rowSums(d0B.PD.mlm^2 / DELTA.d.mlm)

      for (uuu in seq(li.mlm))  # Copy it. ccc.
        wz[, iam(posn.pdip.mix,
                 posn.pstr.mlm - 1 + uuu, M)] <-
          ned2l.dpdip.mix.pstr.mlm[, uuu]  # Link done later
    }  # all(c(ld.mix, li.mlm) > 0



    if (all(c(li.mix, ld.mlm) > 0)) {  # all(tmp3.TF[c(4, 10)])
      for (uuu in seq(ld.mlm))  # tmp3.TF[10]
        for (sss in seq(ld.mlm))
          ned2l.dpstr.mix.pdip.mlm[, uuu] <-
          ned2l.dpstr.mix.pdip.mlm[, uuu] +
            ((sss == uuu) - d0B.PD.mlm[, sss]) *
              d0B.PD.mlm[, sss] / DELTA.d.mlm[, sss]
      if (li.mix)  # tmp3.TF[ 4]
        ned2l.dpstr.mix.pdip.mlm <-
        ned2l.dpstr.mix.pdip.mlm +
        rowSums((d0A.i - d0B.PI.mix) * d0B.PI.mix / DELTA.i.mix)
      if (ld.mix)  # tmp3.TF[ 6]
        ned2l.dpstr.mix.pdip.mlm <-
        ned2l.dpstr.mix.pdip.mlm - rowSums(d0B.PD.mix^2 / DELTA.d.mix)
      if (li.mlm)  # tmp3.TF[ 9]
        ned2l.dpstr.mix.pdip.mlm <-  # nnn.
        ned2l.dpstr.mix.pdip.mlm - rowSums(d0B.PI.mlm^2 / DELTA.i.mlm)

      for (uuu in seq(ld.mlm))  # Copy it. ccc.
        wz[, iam(posn.pstr.mix,
                 posn.pdip.mlm - 1 + uuu, M)] <-
          ned2l.dpstr.mix.pdip.mlm[, uuu]  # Link done later
    }  # all(c(li.mix, ld.mlm) > 0)




 
 
 
    if (lall.len) {
      wz.6 <- matrix(0, n, M * (M + 1) / 2)  # Or == 0 * wz
      ind.rc <- setdiff(1:M, ind.shape.z)  # Contiguous rows and
      lind.rc <- length(ind.rc)  # cols of the DAMLM



 # Copy in the thetas values: the looping is overkill.
      for (uuu in ind.shape.z)
        for (sss in seq(M))
          wz.6[, iam(uuu, sss, M)] <- wz[, iam(uuu, sss, M)]


 


 



      speed.up <- intercept.only && (
                  length(offset) == 1 || all(offset[1] == offset))



      IND.mlm <- iam(NA, NA, lind.rc, both = TRUE, diag = TRUE)
      n.use <- if (speed.up) 2 else n  # For sandwich.mlm





      if (!length(extra$ind.wz.match)) {
        Imat <- matrix(NA, lind.rc, lind.rc)
        for (jay in seq(lind.rc)) {
          iptr <- jay
          for (kay in (ind.rc[jay]):M) {
            if (!any(kay %in% ind.shape.z)) {
              Imat[jay, iptr] <-
                which(extra$index.M$row == ind.rc[jay] &
                      extra$index.M$col == kay)
              iptr <- iptr + 1
            }  # if
          }  # kay
        }  # jay
        ind.wz.match <- Imat[cbind(IND.mlm$row.ind,
                                   IND.mlm$col.ind)]
        extra$ind.wz.match <- ind.wz.match  # Assign it once
      }  # !length(extra$ind.wz.match)


      filling <- if (speed.up)
        wz[1:n.use, extra$ind.wz.match, drop = FALSE] else
        wz[, extra$ind.wz.match, drop = FALSE]









      M.mlm <- lind.rc
      if (is.null(extra$iamlist)) {
        extra$iamlist <- iamlist <- 
          iam(NA, NA, M = M.mlm, both = TRUE)
        if (M.mlm > 1) {  # Offdiagonal elts
          extra$iamlist.nod <- iamlist.nod <-
            iam(NA, NA, M.mlm, both = TRUE, diag = FALSE)
        }
      }  # is.null(extra$iamlist)
      iamlist <- extra$iamlist
      iamlist.nod <- extra$iamlist.nod
      MM12.mlm <- M.mlm * (M.mlm + 1) / 2

      Qf3 <- rowSums(filling[, 1:M.mlm, drop = FALSE] *  # Diag elts
                     (allprobs[1:n.use, 1:M.mlm, drop = FALSE])^2)
      if (M.mlm > 1)  # Offdiagonal elts
        Qf3 <- Qf3 + 2 * rowSums(allprobs[1:n.use, iamlist.nod$row] *
                     filling[, -(1:M.mlm), drop = FALSE] *  # n-vector
                                 allprobs[1:n.use, iamlist.nod$col])
      Qf3 <- matrix(Qf3, n.use, MM12.mlm)



      Qf2rowsums <- matrix(0, n.use, M.mlm)  # rowsums stored columnwise
      for (want in seq(M.mlm)) {  # Want the equivalent of rowSums(Qf2a)
        iamvec <- iam(want, 1:M.mlm, M = M.mlm)  # Diagonals included
        Qf2rowsums[, want] <- rowSums(filling[, iamvec, drop = FALSE] *
                                      allprobs[1:n.use, 1:M.mlm])
      }  # want
      Qf2a <- Qf2rowsums[, iamlist$row]
      Qf2b <- Qf2rowsums[, iamlist$col]


      Qform <- filling - Qf2a - Qf2b + Qf3  # n x MM12.mlm
      Qform <- Qform *
               allprobs[1:n.use, iamlist$row, drop = FALSE] *
               allprobs[1:n.use, iamlist$col, drop = FALSE]



      wz.6[, extra$ind.wz.match] <- if (speed.up)
        matrix(Qform[1, ], n, ncol(Qform), byrow = TRUE) else c(Qform)











 


      dstar.deta <- cbind(dshape.p.deta,
                          if (tmp3.TF[ 3]) dshape.a.deta else NULL,
                          if (tmp3.TF[ 5]) dshape.i.deta else NULL,
                          if (tmp3.TF[ 7]) dshape.d.deta else NULL)
      iptr <- 0
      if (length(ind.shape.z))
      for (uuu in ind.shape.z) {  # Could delete 3 for shape.a (orthog)
        iptr <- iptr + 1
        for (ttt in seq(lind.rc)) {
          wz.6[, iam(uuu, ind.rc[ttt], M)] <- 0  # Initialize
          for (sss in seq(lind.rc)) {
            wz.6[, iam(uuu, ind.rc[ttt], M)] <-
            wz.6[, iam(uuu, ind.rc[ttt], M)] +
              allprobs[, sss] * (max(0, sss == ttt) - allprobs[, ttt]) *
              wz[, iam(uuu, ind.rc[sss], M)] * dstar.deta[, iptr]
          }  # sss
        }  # ttt
      }  # uuu

      wz <- wz.6  # Completed
    }  # lall.len



    if (lall.len) {  # A MLM was fitted
      mytiny <- (allprobs <       sqrt(.Machine$double.eps)) |
                (allprobs > 1.0 - sqrt(.Machine$double.eps))
      atiny <- rowSums(mytiny) > 0
      if (any(atiny)) {
        ind.diags <- setdiff(1:M, ind.shape.z)  # Exclude thetas
        wz[atiny, ind.diags] <- .Machine$double.eps +
        wz[atiny, ind.diags] * (1 + .Machine$double.eps^0.5)
      }
    }  # lall.len




    c(w) * wz
  }), list( .truncate = truncate ))))
}  # gaitdzeta







































 gaitdlog <-
  function(a.mix = NULL, i.mix = NULL, 
           d.mix = NULL,
           a.mlm = NULL, i.mlm = NULL,  # Unstructured probs are
           d.mlm = NULL,                # contiguous
           truncate = NULL, max.support = Inf,
           zero = c("pobs", "pstr", "pdip"),  # Pruned, handles all 6
           eq.ap = TRUE, eq.ip = TRUE, eq.dp = TRUE,
           parallel.a = FALSE, parallel.i = FALSE, parallel.d = FALSE,
           lshape.p = "logitlink",
           lshape.a = lshape.p,  # "logitlink", 20201117
           lshape.i = lshape.p,  # "logitlink", 20201117
           lshape.d = lshape.p,  # "logitlink", 20211011
           type.fitted = c("mean", "shapes",
                           "pobs.mlm", "pstr.mlm", "pdip.mlm",
                           "pobs.mix", "pstr.mix", "pdip.mix",
                           "Pobs.mix", "Pstr.mix", "Pdip.mix",
                           "nonspecial", "Numer", "Denom.p",
                           "sum.mlm.i", "sum.mix.i",
                           "sum.mlm.d", "sum.mix.d",
                           "ptrunc.p", "cdf.max.s"),
           gshape.p = -expm1(-7 * ppoints(12)),
           gpstr.mix = ppoints(7) / 3,  # ppoints(9) / 2,
           gpstr.mlm = ppoints(7) / (3 + length(i.mlm)),
           imethod = 1,
           mux.init = c(0.75, 0.5, 0.75),   # Order is A, I, D.
           ishape.p = NULL, ishape.a = ishape.p,
           ishape.i = ishape.p, ishape.d = ishape.p,
           ipobs.mix = NULL, ipstr.mix = NULL,  # 0.25, 
           ipdip.mix = NULL,   # 0.01,   # Easy but inflexible 0.01
           ipobs.mlm = NULL, ipstr.mlm = NULL,  # 0.25, 
           ipdip.mlm = NULL,   # 0.01,   # NULL, Easy but inflexible
           byrow.aid = FALSE,
           ishrinkage = 0.95,
           probs.y = 0.35) {

  mux.init <- rep_len(mux.init, 3)
  if (length(a.mix) == 0) a.mix <- NULL
  if (length(i.mix) == 0) i.mix <- NULL
  if (length(d.mix) == 0) d.mix <- NULL
  if (length(a.mlm) == 0) a.mlm <- NULL
  if (length(i.mlm) == 0) i.mlm <- NULL
  if (length(d.mlm) == 0) d.mlm <- NULL
  if (length(truncate) == 0) truncate <- NULL



  lowsup <- 1
  gaitd.errorcheck(a.mix, a.mlm, i.mix, i.mlm,
                   d.mix, d.mlm, truncate, max.support,
                   min.support = lowsup)
  la.mix <- length(a.mix <- sort(a.mix))
  li.mix <- length(i.mix <- sort(i.mix))
  ld.mix <- length(d.mix <- sort(d.mix))
  la.mlm <- length(a.mlm)
  li.mlm <- length(i.mlm)
  ld.mlm <- length(d.mlm)
  ltruncat <- length(truncate <- sort(truncate))
  ltrunc.use <- ltruncat > 0 || !is.infinite(max.support) 

  lshape.p <- as.list(substitute(lshape.p))
  eshape.p <- link2list(lshape.p)
  lshape.p <- attr(eshape.p, "function.name")
  lshape.p.save <- lshape.p

  lpobs.mix <- "multilogitlink"  # \omega_p
  epobs.mix <- list()  # zz NULL for now 20200907 coz 'multilogitlink'
  eshape.a <- link2list(lshape.a)
  lshape.a <- attr(eshape.a, "function.name")

  lpstr.mix <- "multilogitlink"  # \phi_p
  epstr.mix <- list()  # zz NULL for now 20200907 coz 'multilogitlink'
  lpdip.mix <- "multilogitlink"  # zz unsure 20211002
  epdip.mix <- list()  # zz unsure 20211002
  eshape.i <- link2list(lshape.i)
  lshape.i <- attr(eshape.i, "function.name")
  eshape.d <- link2list(lshape.d)
  lshape.d <- attr(eshape.d, "function.name")
  lshape.p.save <- lshape.p
  gshape.p.save <- gshape.p


  if (is.vector(zero) && is.character(zero) && length(zero) == 3) {
    if (li.mix + li.mlm == 0)
      zero <- setdiff(zero, "pstr")
    if (la.mix + la.mlm == 0)
      zero <- setdiff(zero, "pobs")
    if (ld.mix + ld.mlm == 0)
      zero <- setdiff(zero, "pdip")
    if (length(zero) == 0)
      zero <- NULL  # Better than character(0)
  }


  lall.len <- la.mix + li.mix + ld.mix + la.mlm + li.mlm + ld.mlm
  if (lall.len + ltruncat == 0 && is.infinite(max.support))
    return(eval(substitute(
        logff(lshape = .lshape.p.save ,
              gshape = .gshape.p.save ,
              zero = NULL),
           list( .lshape.p.save = lshape.p.save,
                 .gshape.p.save = gshape.p.save ))))

  if (!is.logical(eq.ap) || length(eq.ap) != 1)
    stop("argument 'eq.ap' must be a single logical")
  if (!is.logical(eq.ip) || length(eq.ip) != 1)
    stop("argument 'eq.ip' must be a single logical")
  if (!is.logical(parallel.a) || length(parallel.a) != 1)
    stop("argument 'parallel.a' must be a single logical")
  if (!is.logical(parallel.i) || length(parallel.i) != 1)
    stop("argument 'parallel.i' must be a single logical")
  if (!is.logical(parallel.d) || length(parallel.d) != 1)
    stop("argument 'parallel.d' must be a single logical")


  if (FALSE) {  # Comment this out to allow default eq.ap = TRUE, etc.
  if (la.mix <= 1 && eq.ap)
    stop("<= one unstructured altered value (no 'shape.a')",
         ", so setting 'eq.ap = TRUE' is meaningless")
  if (li.mix <= 1 && eq.ip)
    stop("<= one unstructured inflated value (no 'shape.i')",
            ", so setting 'eq.ip = TRUE' is meaningless")
  if (ld.mix <= 1 && eq.dp)
    stop("<= one unstructured deflated value (no 'shape.d')",
         ", so setting 'eq.dp = TRUE' is meaningless")
  if (la.mlm <= 1 && parallel.a)  # Only \omega_1
    stop("<= one altered mixture probability, 'pobs", a.mlm,
            "', so setting 'parallel.a = TRUE' is meaningless")
  if (li.mlm <= 1 && parallel.i)  # Only \phi_1
    stop("<= one inflated mixture probability, 'pstr", i.mlm,
            "', so setting 'parallel.i = TRUE' is meaningless")
  if (ld.mlm <= 1 && parallel.d)  # Only \psi_1
    stop("<= one deflated mixture probability, 'pdip", d.mlm,
         "', so setting 'parallel.d = TRUE' is meaningless")
  }  # FALSE


  type.fitted.choices <-
            c("mean", "shapes",
              "pobs.mlm", "pstr.mlm", "pdip.mlm",
              "pobs.mix", "pstr.mix", "pdip.mix",
              "Pobs.mix", "Pstr.mix", "Pdip.mix",
              "nonspecial", "Numer", "Denom.p",
              "sum.mlm.i", "sum.mix.i",
              "sum.mlm.d", "sum.mix.d",
              "ptrunc.p", "cdf.max.s")
  type.fitted <- match.arg(type.fitted[1], type.fitted.choices)[1]

  tmp7a <- if (la.mlm) paste0("pobs.mlm", a.mlm) else NULL
  tmp7b <- if (li.mlm) paste0("pstr.mlm", i.mlm) else NULL
  tmp7c <- if (ld.mlm) paste0("pdip.mlm", d.mlm) else NULL
  tmp3 <- c(shape.p = lshape.p,
            pobs.mix = if (la.mix) "multilogitlink" else NULL,
            shape.a = if (la.mix > 1) lshape.a else NULL,
            pstr.mix = if (li.mix) "multilogitlink" else NULL,
            shape.i = if (li.mix > 1) lshape.i else NULL,
            pdip.mix = if (ld.mix) "multilogitlink" else NULL,
            shape.d = if (ld.mix > 1) lshape.d else NULL,
            if (la.mlm) rep("multilogitlink", la.mlm) else NULL,
            if (li.mlm) rep("multilogitlink", li.mlm) else NULL,
            if (ld.mlm) rep("multilogitlink", ld.mlm) else NULL)
  Ltmp3 <- length(tmp3) 
  if (la.mlm + li.mlm + ld.mlm)
    names(tmp3)[(Ltmp3 - la.mlm - li.mlm - ld.mlm + 1):Ltmp3] <-
      c(tmp7a, tmp7b, tmp7c)
  par1or2 <- 1  # 2
  tmp3.TF <- c(TRUE, la.mix > 0, la.mix > 1,
                     li.mix > 0, li.mix > 1,
                     ld.mix > 0, ld.mix > 1,
                     la.mlm > 0, li.mlm > 0, ld.mlm > 0)
  indeta.finish <- cumsum(c(par1or2, 1, par1or2,
                                     1, par1or2,
                                     1, par1or2,
                            la.mlm, li.mlm, ld.mlm,
                            ld.mlm + 1) * c(tmp3.TF, 1))
  indeta.launch <- c(1, 1 + head(indeta.finish, -1))

  indeta.launch <- head(indeta.launch, -1)
  indeta.finish <- head(indeta.finish, -1)
  indeta.launch[!tmp3.TF] <- NA  # Not to be accessed
  indeta.finish[!tmp3.TF] <- NA  # Not to be accessed
  indeta <- cbind(launch = indeta.launch,
                  finish = indeta.finish)
  rownames(indeta) <- c("shape.p",
                        "pobs.mix", "shape.a",
                        "pstr.mix", "shape.i",
                        "pdip.mix", "shape.d",
                        "pobs.mlm", "pstr.mlm", "pdip.mlm")
  M1 <- max(indeta, na.rm = TRUE)
  predictors.names <- tmp3  # Passed into @infos and @initialize.
      

  blurb1 <- "L"   # zz1
  if (la.mlm + la.mix) blurb1 <- "Generally-altered l"
  if (li.mlm + li.mix) blurb1 <- "Generally-inflated l"
  if (ltrunc.use) blurb1 <- "Generally-truncated l"
  if ( (la.mlm + la.mix) &&  (li.mlm + li.mix) && !ltrunc.use)
    blurb1 <- "Generally-altered and -inflated l"
  if ( (la.mlm + la.mix) && !(li.mlm + li.mix) &&  ltrunc.use)
    blurb1 <- "Generally-altered and -truncated l"
  if (!(la.mlm + la.mix) &&  (li.mlm + li.mix) &&  ltrunc.use)
    blurb1 <- "Generally-inflated and -truncated l"
  if ( (la.mlm + la.mix) &&  (li.mlm + li.mix) &&  ltrunc.use)
    blurb1 <- "Generally-altered, -inflated and -truncated l"

  if (ld.mlm + ld.mix) blurb1 <-
    c(blurb1,
      if (la.mlm + la.mix + li.mlm + li.mix) "and " else "Generally",
      "-deflated ")



      
  new("vglmff",
  blurb = c(blurb1, "ogarithmic regression\n",
            "(GAITD-Log(shape.p)-",
                   "Log(shape.a)-MLM-",
                   "Log(shape.i)-MLM-\n",
                   "Log(shape.d)-MLM generally)\n\n",
            "Links: ",
            namesof("shape.p", lshape.p, earg = eshape.p,
                    tag = FALSE),
            if (la.mix > 0) c(", ", "multilogit(pobs.mix)"),
            if (la.mix > 1) c(", ",
            namesof("shape.a",  lshape.a, eshape.a, tag = FALSE)),
            if (la.mix && li.mix) ", \n       ",
            if (li.mix > 0) c(  if (la.mix) "" else ", ",
            "multilogit(pstr.mix)"),
            if (li.mix > 1) c(", ",
            namesof("shape.i",  lshape.i, eshape.i, tag = FALSE)),
            if (li.mix && ld.mix) ", \n       ",
            if (ld.mix > 0) c(  if (li.mix) "" else ", ",
            "multilogit(pdip.mix)"),
            if (ld.mix > 1) c(", ",
            namesof("shape.d",  lshape.d, eshape.d, tag = FALSE)),
            if (la.mlm) paste0(",\n",
              paste0("       multilogit(", tmp7a, collapse = "),\n"),
            ")") else NULL,
            if (li.mlm) paste0(",\n",
              paste0("       multilogit(", tmp7b, collapse = "),\n"),
              ")") else NULL,
            if (ld.mlm) paste0(",\n",
              paste0("       multilogit(", tmp7c, collapse = "),\n"),
            ")") else NULL),
  constraints = eval(substitute(expression({
    M1 <- max(extra$indeta, na.rm = TRUE)
    la.mix <- ( .la.mix )
    li.mix <- ( .li.mix )
    ld.mix <- ( .ld.mix )
    la.mlm <- ( .la.mlm )
    li.mlm <- ( .li.mlm )
    ld.mlm <- ( .ld.mlm )

    use.mat.mlm.a <- if (la.mlm) {
      if ( .parallel.a ) matrix(1, la.mlm, 1) else diag(la.mlm)
    } else {
      NULL
    }
    use.mat.mlm.i <- if (li.mlm) {
       if ( .parallel.i ) matrix(1, li.mlm, 1) else diag(li.mlm)
    } else {
      NULL
    }
    use.mat.mlm.d <- if (ld.mlm) {
       if ( .parallel.d ) matrix(1, ld.mlm, 1) else diag(ld.mlm)
    } else {
      NULL
    }

    if (la.mlm + li.mlm + ld.mlm == 0) {
      Use.mat <- use.mat.mlm <- cbind(M)  # shape.p only
    }
    if (la.mlm + li.mlm + ld.mlm) {
      nc1 <- if (length(use.mat.mlm.a)) ncol(use.mat.mlm.a) else 0
      nc2 <- if (length(use.mat.mlm.i)) ncol(use.mat.mlm.i) else 0
      nc3 <- if (length(use.mat.mlm.d)) ncol(use.mat.mlm.d) else 0
      use.mat.mlm <- cbind(1, matrix(0, 1, nc1 + nc2 + nc3))
      if (la.mlm)
        use.mat.mlm <- rbind(use.mat.mlm,
                             cbind(matrix(0, la.mlm, 1),
                                   use.mat.mlm.a,
                                   if (length(use.mat.mlm.i) == 0)
                                   NULL else matrix(0, la.mlm, nc2),
                                   if (length(use.mat.mlm.d) == 0)
                                   NULL else matrix(0, la.mlm, nc3)))
      if (li.mlm )
       use.mat.mlm <-
         rbind(use.mat.mlm,
               cbind(matrix(0, li.mlm, 1 + nc1),
                     use.mat.mlm.i,
                     matrix(0, li.mlm, nc3)))
      if (ld.mlm)
        use.mat.mlm <-
          rbind(use.mat.mlm,  # zz1 next line:
                cbind(matrix(0, ld.mlm, 1 + nc1 + nc2),
                      use.mat.mlm.d))
    }  # la.mlm + li.mlm






    
    tmp3.TF <- ( .tmp3.TF )   # Logical of length 10.


  use.mat.mix <- cm3gaitd( .eq.ap , .eq.ip , .eq.dp , npar = 1)



    tmp3.subset <- tmp3.TF[-(8:10)]
    use.mat.mix <- use.mat.mix[tmp3.subset, , drop = FALSE]
    notall0 <- function(x) !all(x == 0)
    use.mat.mix <- use.mat.mix[, apply(use.mat.mix, 2, notall0),
                               drop = FALSE]


    if (la.mix + li.mix + ld.mix > 0)
      Use.mat <- use.mat.mix





    if (la.mlm + li.mlm + ld.mlm > 0) {
      Use.mat <- rbind(use.mat.mix,
                       matrix(0, nrow(use.mat.mlm) - 1,  # bottom
                                 ncol(use.mat.mix)))
      Use.mat <- cbind(Use.mat,
                       matrix(0, nrow(Use.mat),  # RHS
                                 ncol(use.mat.mlm) - 1))
      Use.mat[row(Use.mat) > nrow(use.mat.mix) &
              col(Use.mat) > ncol(use.mat.mix)] <- use.mat.mlm[-1, -1]
    }  # la.mlm + li.mlm + ld.mlm > 0



    if (is.null(constraints)) {
      constraints <-
        cm.VGAM(Use.mat, x = x, apply.int = TRUE,  # FALSE
                bool = .eq.ap || .eq.ip || .eq.dp ||
                       .parallel.a || .parallel.i || .parallel.d ,
                constraints = constraints)  # FALSE
    }  # is.null(constraints)



    if (la.mix + li.mix + ld.mix + la.mlm + li.mlm + ld.mlm)
      constraints <-
        cm.zero.VGAM(constraints, x = x, .zero , M = M, M1 = M1,
                     predictors.names = paste0(predictors.names,
                                        names(predictors.names)))
  }),
  list( .zero = zero, .tmp3.TF = tmp3.TF,  # .tmp3 = tmp3,
        .eq.ap = eq.ap, .eq.ip = eq.ip, .eq.dp = eq.dp,
        .parallel.a = parallel.a, .parallel.i = parallel.i,
        .parallel.d = parallel.d,
        .la.mlm = la.mlm, .li.mlm = li.mlm, .ld.mlm = ld.mlm,
        .la.mix = la.mix, .li.mix = li.mix, .ld.mix = ld.mix ))),
  infos = eval(substitute(function(...) {
    list(M1 = .M1 ,
         Q1 = 1,
         dpqrfun = "gaitdlog",
         link = .predictors.names ,  # ...strips... from above
         link1parameter = as.logical( .lall.len <= 2),  # <= 1 safer
         mixture.links  = any(c( .la.mlm , .li.mlm , .ld.mlm ,
                                 .la.mix , .li.mix ,
                                 .ld.mix ) > 1),  # FALSE if NULL
         a.mix = as.vector( .a.mix ),  # Handles NULL
         a.mlm = as.vector( .a.mlm ),
         i.mix = as.vector( .i.mix ),
         i.mlm = as.vector( .i.mlm ),
         d.mix = as.vector( .d.mix ),
         d.mlm = as.vector( .d.mlm ),
         truncate = as.vector( .truncate ),
         max.support = as.vector( .max.support ),
         Support  = c( .lowsup , Inf, 1),  # a(b)c format as a,c,b.
         expected = TRUE,
         multipleResponses = FALSE,  # poissonff can be called if TRUE
         parameters.names = names( .predictors.names ),
         parent.name = c("logff", "log"),
         type.fitted  = as.vector( .type.fitted ),
         type.fitted.choices = ( .type.fitted.choices ),
         baseparams.argnames  = "shape",
         MM1 = 1,  # One parameter for 1 response (shape). Needed.
         zero = .zero )
  }, list( .zero = zero, .lowsup = lowsup,
           .type.fitted = type.fitted,
           .type.fitted.choices = type.fitted.choices,
           .lshape.p = lshape.p, .eshape.p = eshape.p,
           .lshape.a = lshape.a, .eshape.a = eshape.a,
           .lshape.i = lshape.i, .eshape.i = eshape.i,
           .a.mlm = a.mlm, .i.mlm = i.mlm, .d.mlm = d.mlm,
           .a.mix = a.mix, .i.mix = i.mix, .d.mix = d.mix,
           .la.mlm = la.mlm, .li.mlm = li.mlm, .ld.mlm = ld.mlm,
           .la.mix = la.mix, .li.mix = li.mix, .ld.mix = ld.mix,
           .truncate = truncate, .max.support = max.support,
           .predictors.names = predictors.names,
           .M1 = M1, .lall.len = lall.len
         ))),
  initialize = eval(substitute(expression({
    extra$indeta <- ( .indeta )  # Avoids recomputing it several times
    la.mix <- length((a.mix <- as.vector( .a.mix )))
    li.mix <- length((i.mix <- as.vector( .i.mix )))
    ld.mix <- length((d.mix <- as.vector( .d.mix )))
    la.mlm <- length((a.mlm <- as.vector( .a.mlm )))
    li.mlm <- length((i.mlm <- as.vector( .i.mlm )))
    ld.mlm <- length((d.mlm <- as.vector( .d.mlm )))
    lall.len <- la.mix + li.mix + ld.mix +
                la.mlm + li.mlm + ld.mlm
    truncate <- as.vector( .truncate )
    ltruncat <- length(truncate)
    M1 <- max(extra$indeta, na.rm = TRUE)
    NOS <- NCOL(y)
    M <- NOS * M1
    tmp3.TF <- ( .tmp3.TF )

    temp5 <-
    w.y.check(w = w, y = y,
              Is.nonnegative.y = TRUE,
              Is.integer.y = TRUE,
              ncol.w.max = 1,  # Since max.support = 9 is possible
              ncol.y.max = 1,
              out.wy = TRUE, colsyperw = 1, maximize = TRUE)
    w <- temp5$w
    y <- temp5$y
    glist <- y.gaitcombo.check(y, truncate = truncate,
                               a.mlm = a.mlm, a.mix = a.mix,
                               i.mlm = i.mlm, i.mix = i.mix,
                               d.mlm = d.mlm, d.mix = d.mix,
                               max.support = .max.support ,
                               min.support = .min.support )
    extra$skip.mix.a <- glist$skip.mix.a
    extra$skip.mix.i <- glist$skip.mix.i
    extra$skip.mix.d <- glist$skip.mix.d
    extra$skip.mlm.a <- glist$skip.mlm.a
    extra$skip.mlm.i <- glist$skip.mlm.i
    extra$skip.mlm.d <- glist$skip.mlm.d


    extra$NOS <- NOS <- ncoly <- ncol(y)  # Number of species
    extra$type.fitted <- as.vector( .type.fitted )
    extra$mux.init <- as.vector( .mux.init )
    extra$colnames.y  <- colnames(y)
    extra$M1 <- M1
    extra$index.M <- iam(NA, NA, M, both = TRUE)  # Used in @weight



    predictors.names <- ( .predictors.names )  # Got it, named


    if (!length(etastart)) {
      shape.p.init <- if (length( .ishape.p )) .ishape.p else {
        logff.Loglikfun <- function(shapeval, y, x, w, extraargs) {
          sum(c(w) * dlog(x = y, shape = shapeval, log = TRUE))
        }
        shape.p.grid <- ( .gshape.p )
        grid.search(shape.p.grid, objfun = logff.Loglikfun, y = y, w = w)
      }
      shape.p.init <- rep(shape.p.init, length = n)

      shape.d.init <-
      shape.a.init <- shape.i.init <- shape.p.init  # Needed


      etastart <- matrix(nrow = n, ncol = M,
        theta2eta(shape.p.init, .lshape.p , earg = .eshape.p ))




      mux.more.a <- extra$mux.init[1]   # 0.75 Err to slightly smaller
      init.pobs.mix <- numeric(n)
      if (tmp3.TF[ 2]) {  # la.mix > 0
        init.pobs.mix <- if (length( .ipobs.mix )) {
        rep_len( .ipobs.mix , n)
      } else {
        is.a.mix1 <- rowSums(extra$skip.mix.a) > 0
        rep(mux.more.a * sum(w[is.a.mix1]) / sum(w), n)
      }  
      }  # la.mix > 0


      if (tmp3.TF[ 3]) {  # Assign coln 3; la.mix > 1
        shape.a.init <- if (length( .ishape.a ))
          rep_len( .ishape.a , n) else shape.p.init  # A vector
        etastart[, 3] <-
          theta2eta(shape.a.init, .lshape.a , earg = .eshape.a )
      }

      init.pstr.mix <- init.pdip.mix <- numeric(n)
      try.gridsearch.pstr.mix <- FALSE
      if (tmp3.TF[ 4]) {  # li.mix > 0
        init.pstr.mix <- if (length( .ipstr.mix )) {
          rep_len( .ipstr.mix , n)
        } else {
          try.gridsearch.pstr.mix <- TRUE
          numeric(n)  # Overwritten by gridsearch
        }
      }  # li.mix > 0


      if (tmp3.TF[ 5]) {  # li.mix > 1
        shape.i.init <- if (length( .ishape.i ))
          rep_len( .ishape.i , n) else shape.p.init  # A vector
        etastart[, (extra$indeta[5, 'launch'])] <-
          theta2eta(shape.i.init, .lshape.i , earg = .eshape.i )
      }  # li.mix > 1


      if (tmp3.TF[ 8]) {  #  la.mlm
        init.pobs.mlm <- if (length( .ipobs.mlm )) {
          matrix( .ipobs.mlm , n, la.mlm, byrow = .byrow.aid )
        } else {
          mux.more.a <- extra$mux.init[1]                            
          init.pobs.mlm <- colSums(c(w) * extra$skip.mlm.a) / colSums(w)
          init.pobs.mlm <- init.pobs.mlm * as.vector( mux.more.a )
          matrix(init.pobs.mlm, n, la.mlm, byrow = TRUE)
        }
      } else {
        init.pobs.mlm <- matrix(0, n, 1)
      }


      try.gridsearch.pstr.mlm <- FALSE
      if (tmp3.TF[ 9]) {  #  li.mlm
        try.gridsearch.pstr.mlm <- !(length( .ipstr.mlm ))
        init.pstr.mlm <- 0  # Might be overwritten by gridsearch

        if (length( .ipstr.mlm ))
          init.pstr.mlm <- as.vector( .ipstr.mlm )

        init.pstr.mlm <- matrix(init.pstr.mlm, n, li.mlm,
                                byrow = .byrow.aid )
      } else {
        init.pstr.mlm <- matrix(0, n, 1)
      }

        

      init.pdip.mlm <- matrix(0, n, 2)  # rowSums() needs > 1 colns.





      gaitdlog.Loglikfun1.mix <-
        function(pstr.mix.val, y, x, w, extraargs) {
    sum(c(w) *
        dgaitdlog(y, pstr.mix = pstr.mix.val,
                  pstr.mlm    = extraargs$pstr.mlm,  # Differs here
                  shape.p     = extraargs$shape.p,
                  shape.a     = extraargs$shape.a,
                  shape.i     = extraargs$shape.i,
                  shape.d     = extraargs$shape.d,
                  a.mix       = extraargs$a.mix,
                  a.mlm       = extraargs$a.mlm,
                  i.mix       = extraargs$i.mix,
                  i.mlm       = extraargs$i.mlm,
                  d.mix       = extraargs$d.mix,
                  d.mlm       = extraargs$d.mlm,
                  max.support = extraargs$max.support,
                  truncate    = extraargs$truncate,
                  pobs.mix    = extraargs$pobs.mix,
                  pobs.mlm    = extraargs$pobs.mlm,
                  pdip.mix    = extraargs$pdip.mix,
                  pdip.mlm    = extraargs$pdip.mlm, log = TRUE))
  }

 gaitdlog.Loglikfun1.mlm <-
     function(pstr.mlm.val, y, x, w, extraargs) {
    sum(c(w) *
        dgaitdlog(y, pstr.mlm  = pstr.mlm.val,
                   pstr.mix    = extraargs$pstr.mix,  # Differs here
                   shape.p     = extraargs$shape.p,
                   shape.a     = extraargs$shape.a,
                   shape.i     = extraargs$shape.i,
                   shape.d     = extraargs$shape.d,
                   a.mix       = extraargs$a.mix,
                   a.mlm       = extraargs$a.mlm,
                   i.mix       = extraargs$i.mix,
                   i.mlm       = extraargs$i.mlm,
                   d.mix       = extraargs$d.mix,
                   d.mlm       = extraargs$d.mlm,
                   max.support = extraargs$max.support,
                   truncate    = extraargs$truncate,
                   pobs.mix    = extraargs$pobs.mix,
                   pobs.mlm    = extraargs$pobs.mlm,
                   pdip.mix    = extraargs$pdip.mix,
                   pdip.mlm    = extraargs$pdip.mlm, log = TRUE))
  }

 gaitdlog.Loglikfun2 <-
     function(pstr.mix.val, pstr.mlm.val, y, x, w, extraargs) {
    sum(c(w) *
        dgaitdlog(y, pstr.mix = pstr.mix.val,
                   pstr.mlm    = pstr.mlm.val,
                   shape.p     = extraargs$shape.p,
                   shape.a     = extraargs$shape.a,
                   shape.i     = extraargs$shape.i,
                   shape.d     = extraargs$shape.d,
                   a.mix       = extraargs$a.mix,
                   a.mlm       = extraargs$a.mlm,
                   i.mix       = extraargs$i.mix,
                   i.mlm       = extraargs$i.mlm,
                   d.mix       = extraargs$d.mix,
                   d.mlm       = extraargs$d.mlm,
                   max.support = extraargs$max.support,
                   truncate    = extraargs$truncate,
                   pobs.mix    = extraargs$pobs.mix,
                   pobs.mlm    = extraargs$pobs.mlm,
                   pdip.mix    = extraargs$pdip.mix,
                   pdip.mlm    = extraargs$pdip.mlm, log = TRUE))
  }



      if (li.mix + li.mlm) {
        extraargs <- list(
              shape.p     = shape.p.init,
              shape.a     = shape.a.init,
              shape.i     = shape.i.init,
              shape.d     = shape.d.init,
              a.mix       = a.mix,
              a.mlm       = a.mlm,
              i.mix       = i.mix,
              i.mlm       = i.mlm,
              d.mix       = d.mix,
              d.mlm       = d.mlm,
              truncate    = truncate,
              max.support = as.vector( .max.support ),
              pobs.mix    = init.pobs.mix ,
              pobs.mlm    = init.pobs.mlm ,
              pdip.mix    = init.pdip.mix ,
              pdip.mlm    = init.pdip.mlm )
        pre.warn <- options()$warn
        options(warn = -1)  # Ignore warnings during gridsearch 

        try.this <-
          if (try.gridsearch.pstr.mix && try.gridsearch.pstr.mlm) {
            grid.search2( .gpstr.mix ,  .gpstr.mlm ,
                         objfun = gaitdlog.Loglikfun2,
                         y = y, w = w, extraargs = extraargs,
                         ret.objfun = TRUE)
          } else if (try.gridsearch.pstr.mix) {
            extraargs$pstr.mlm <- init.pstr.mlm
            grid.search ( .gpstr.mix ,
                         objfun = gaitdlog.Loglikfun1.mix,
                         y = y, w = w, extraargs = extraargs,
                         ret.objfun = TRUE)
          } else if (try.gridsearch.pstr.mlm) {
            extraargs$pstr.mix <- init.pstr.mix
            grid.search ( .gpstr.mlm ,
                         objfun = gaitdlog.Loglikfun1.mlm,
                         y = y, w = w, extraargs = extraargs,
                         ret.objfun = TRUE)
          }

        options(warn = pre.warn)  # Restore warnings 
        if (any(is.na(try.this)))
          warning("gridsearch returned NAs. It's going to crash.",
                  immediate. = TRUE)
        if (try.gridsearch.pstr.mix && try.gridsearch.pstr.mlm) {
          init.pstr.mix <- rep_len(try.this["Value1"], n)
          init.pstr.mlm <- matrix(try.this["Value2"], n, li.mlm)
          if (any(is.na(try.this)))
            stop("Crashing. ",
                 "Try something like 'gpstr.mix = seq(5) / 100'",
                 " and/or 'gpstr.mlm = seq(5) / 100'.")
        } else if (try.gridsearch.pstr.mix) {
          init.pstr.mix <- rep_len(try.this["Value"], n)
          if (any(is.na(try.this)))
            stop("Crashing. ",
                 "Try something like 'gpstr.mix = seq(5) / 100'.")
        } else if (try.gridsearch.pstr.mlm) {
          init.pstr.mlm <- matrix(try.this["Value"], n, li.mlm)
          if (any(is.na(try.this)))
            stop("Crashing. ",
                 "Try something like 'gpstr.mlm = seq(5) / 100'.")
        }
      }  # la.mix + lnf.mix



      mux.more.d <- extra$mux.init[3]   
      if (ld.mix) {
        init.pdip.mix <- if (length( .ipdip.mix ))
          rep_len( .ipdip.mix, n) else {
          is.d.mix1 <- rowSums(extra$skip.mix.d) > 0
          rep(mux.more.d * sum(w[is.d.mix1]) / sum(w), n) 
        }
      }  # ld.mix

      if (ld.mlm) {
        init.pdip.mlm <- if (length( .ipdip.mlm ))
          matrix( .ipdip.mlm, n, ld.mlm, byrow = TRUE) else {
          is.d.mlm1 <- rowSums(extra$skip.mlm.d) > 0
          matrix(mux.more.d * (sum(w[is.d.mlm1]) / sum(w)) / ld.mlm,
                 n, ld.mlm)
        }
      }  # ld.mlm




        


      while (any((vecTF <- init.pobs.mix + init.pstr.mix +  # -
                           init.pdip.mix +
                           rowSums(init.pobs.mlm) +
                           rowSums(init.pstr.mlm) +  # -
                           rowSums(init.pdip.mlm) > 0.96875))) {
        init.pobs.mix[vecTF]   <- 0.875 * init.pobs.mix[vecTF]
        init.pstr.mix[vecTF]   <- 0.875 * init.pstr.mix[vecTF]
        init.pdip.mix[vecTF]   <- 0.875 * init.pdip.mix[vecTF]
        init.pobs.mlm[vecTF, ] <- 0.875 * init.pobs.mlm[vecTF, ]
        init.pstr.mlm[vecTF, ] <- 0.875 * init.pstr.mlm[vecTF, ]
        init.pdip.mlm[vecTF, ] <- 0.875 * init.pdip.mlm[vecTF, ]
      }  # while

      Numer.init1 <- 1 - rowSums(init.pobs.mlm) -
                         rowSums(init.pstr.mlm) -  # +
                         rowSums(init.pdip.mlm) -
                         init.pobs.mix - init.pstr.mix -  # +
                         init.pdip.mix  # Differs from 'Numer'.
        
      etastart.z <- if (lall.len == 0) NULL else {
        tmp.mat <- cbind(if (tmp3.TF[ 2]) init.pobs.mix else NULL,
                         if (tmp3.TF[ 4]) init.pstr.mix else NULL,
                         if (tmp3.TF[ 6]) init.pdip.mix else NULL,
                         if (tmp3.TF[ 8]) init.pobs.mlm else NULL,
                         if (tmp3.TF[ 9]) init.pstr.mlm else NULL,
                         if (tmp3.TF[10]) init.pdip.mlm else NULL,
                         Numer.init1)
        multilogitlink(tmp.mat)
      }  # etastart.z
      if (!is.matrix(etastart.z)) etastart.z <- cbind(etastart.z)

      nextone <- 1  # Might not be used actually
      if (tmp3.TF[ 2]) {
        etastart[, 2] <- etastart.z[, nextone]
        nextone <- nextone + 1
      }
      if (tmp3.TF[ 4]) {  # Coln 2 or 4
        etastart[, (extra$indeta[4, 'launch'])] <- etastart.z[, nextone]
        nextone <- nextone + 1
      }
      if (tmp3.TF[ 6]) {  # Coln 2 or 4 or 6
        etastart[, (extra$indeta[6, 'launch'])] <- etastart.z[, nextone]
        nextone <- nextone + 1
      }
      if (tmp3.TF[ 8]) {
        ind8 <- (extra$indeta[8, 'launch']):(extra$indeta[8, 'finish'])
        etastart[, ind8] <- etastart.z[, nextone:(nextone+la.mlm - 1)]
        nextone <- nextone + la.mlm
      }
      if (tmp3.TF[ 9]) {
        ind9 <- (extra$indeta[9, 'launch']):(extra$indeta[9, 'finish'])
        etastart[, ind9] <- etastart.z[, nextone:(nextone+li.mlm - 1)]
        nextone <- nextone + li.mlm
      }
      if (tmp3.TF[10]) {
        ind0 <- (extra$indeta[10, 'launch']):(extra$indeta[10, 'finish'])
        etastart[, ind0] <- etastart.z[, nextone:(nextone + ld.mlm - 1)]
        if (ncol(etastart.z) != nextone + ld.mlm - 1)
          stop("miscalculation")
      }
    }
  }), list(
    .lshape.p = lshape.p, .eshape.p = eshape.p,
    .lshape.a = lshape.a, .eshape.a = eshape.a,
    .lshape.i = lshape.i, .eshape.i = eshape.i,
    .lshape.d = lshape.d, .eshape.d = eshape.d,
    .ishape.p = ishape.p,
    .ishape.a = ishape.a,
    .ishape.i = ishape.i,
    .ishape.d = ishape.d,
    .lpstr.mix = lpstr.mix, .lpobs.mix = lpobs.mix,
    .epstr.mix = epstr.mix, .epobs.mix = epobs.mix,
    .lpdip.mix = lpdip.mix,
    .epdip.mix = epdip.mix,
    .ipstr.mix = ipstr.mix, .ipobs.mix = ipobs.mix,
    .ipstr.mlm = ipstr.mlm, .ipobs.mlm = ipobs.mlm,
    .ipdip.mix = ipdip.mix,
    .ipdip.mlm = ipdip.mlm,
    .byrow.aid = byrow.aid,
    .a.mix = a.mix, .i.mix = i.mix, .d.mix = d.mix,
    .a.mlm = a.mlm, .i.mlm = i.mlm, .d.mlm = d.mlm,
    .truncate = truncate, .max.support = max.support,
    .min.support = lowsup,
    .tmp3.TF = tmp3.TF,  # .tmp3 = tmp3,
    .predictors.names = predictors.names,
    .mux.init = mux.init,
    .gshape.p = gshape.p,
    .gpstr.mix = gpstr.mix,   # .gpdip.mix = gpdip.mix,
    .gpstr.mlm = gpstr.mlm,   # .gpdip.mlm = gpdip.mlm,
    .ishrinkage = ishrinkage, .probs.y = probs.y,
    .indeta = indeta,
    .imethod = imethod, .type.fitted = type.fitted ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    n.obs <- NROW(eta)
    type.fitted <-
      if (length(extra$type.fitted)) extra$type.fitted else {
        warning("cannot find 'type.fitted'. Returning the 'mean'.")
        "mean"
      }
    type.fitted <-
      match.arg(type.fitted[1],
                c("mean", "shapes",
                  "pobs.mlm", "pstr.mlm", "pdip.mlm",
                  "pobs.mix", "pstr.mix", "pdip.mix",
                  "Pobs.mix", "Pstr.mix", "Pdip.mix",
                  "nonspecial", "Numer", "Denom.p",
                  "sum.mlm.i", "sum.mix.i",
                  "sum.mlm.d", "sum.mix.d",
                  "ptrunc.p", "cdf.max.s"))[1]

    if (!is.matrix(eta)) eta <- as.matrix(eta)
    la.mix <- length((a.mix <- as.vector( .a.mix )))
    li.mix <- length((i.mix <- as.vector( .i.mix )))
    ld.mix <- length((d.mix <- as.vector( .d.mix )))
    la.mlm <- length((a.mlm <- as.vector( .a.mlm )))
    li.mlm <- length((i.mlm <- as.vector( .i.mlm )))
    ld.mlm <- length((d.mlm <- as.vector( .d.mlm )))
    truncate <- as.vector( .truncate )
    max.support <- as.vector( .max.support )
    morework <- type.fitted != "mean"  # For efficiency

    lall.len <- la.mix + li.mix + ld.mix + la.mlm + li.mlm + ld.mlm
    pobs.mix <- pstr.mix <- pdip.mix <- 0  # 4 rowSums()
    pobs.mlm <- pstr.mlm <- pdip.mlm <- 0  # matrix(0, NROW(eta), 1)
    shape.p <- cbind(eta2theta(eta[, 1], .lshape.p , .eshape.p ))
    ind.shape.z <- 1  # Points to shape.p only.
    shape.a <- shape.i <-
    shape.d <- shape.p  # Needed; and answer not corrupted
    tmp3.TF <- ( .tmp3.TF )   # Logical of length 10.

    if (any(tmp3.TF[c(3, 5, 7)])) {  # At least one shape.[aid]
      ind.shape.z <- extra$indeta[c(1, 3, 5, 7), 'launch']  # Vectors
      ind.shape.z <- c(na.omit(ind.shape.z))  # At least one value
      shape.a <- if (!tmp3.TF[ 3]) shape.p else
        eta2theta(eta[, extra$indeta[3, 1]], .lshape.a , .eshape.a )
      shape.i <- if (!tmp3.TF[ 5]) shape.p else
        eta2theta(eta[, extra$indeta[5, 1]], .lshape.i , .eshape.i )
      shape.d <- if (!tmp3.TF[ 7]) shape.p else
        eta2theta(eta[, extra$indeta[7, 1]], .lshape.d , .eshape.d )
    }  # la.mix + li.mix + ld.mix > 0






    if (lall.len) {  # An MLM was fitted
      allprobs <-
        multilogitlink(eta[, -ind.shape.z, drop = FALSE],
                       inverse = TRUE)  # rowSums == 1
      if (anyNA(allprobs))
        warning("there are NAs here in slot linkinv")
      if (min(allprobs) == 0 || max(allprobs) == 1)
        warning("fitted probabilities numerically 0 or 1 occurred")


      Nextone <- 0  # Might not be used actually
      if (tmp3.TF[ 2])
        pobs.mix <- allprobs[, (Nextone <- Nextone + 1)]
      if (tmp3.TF[ 4])
        pstr.mix <- allprobs[, (Nextone <- Nextone + 1)]
      if (tmp3.TF[ 6])
        pdip.mix <- allprobs[, (Nextone <- Nextone + 1)]

      if (tmp3.TF[ 8]) {
        ind8 <- (Nextone + 1):(Nextone + la.mlm)
        pobs.mlm <- allprobs[, ind8, drop = FALSE]
        dimnames(pobs.mlm) <- list(rownames(eta), as.character(a.mlm))
        Nextone <- Nextone + la.mlm
      }
      if (tmp3.TF[ 9]) {
        ind9 <- (Nextone + 1):(Nextone + li.mlm)
        pstr.mlm <- allprobs[, ind9, drop = FALSE]
        dimnames(pstr.mlm) <- list(rownames(eta), as.character(i.mlm))
        Nextone <- Nextone + li.mlm
      }
      if (tmp3.TF[10]) {
        ind10 <- (Nextone + 1):(Nextone + ld.mlm)
        pdip.mlm <- allprobs[, ind10, drop = FALSE]
        dimnames(pdip.mlm) <- list(rownames(eta), as.character(d.mlm))
        Nextone <- Nextone + ld.mlm  # Not needed
      }
    }  # lall.len

    M1 <- max(extra$indeta, na.rm = TRUE)
    NOS <- NCOL(eta) / M1

    Bits <- moments.gaitdcombo.log(shape.p,
              pobs.mix = pobs.mix, pobs.mlm = pobs.mlm,
              pstr.mix = pstr.mix, pstr.mlm = pstr.mlm,
              pdip.mix = pdip.mix, pdip.mlm = pdip.mlm,
              a.mix = a.mix, i.mix = i.mix, d.mix = d.mix,
              a.mlm = a.mlm, i.mlm = i.mlm, d.mlm = d.mlm,
              shape.a = shape.a, shape.i = shape.i,
              shape.d = shape.d,
              truncate = truncate, max.support = max.support)


  if (morework) {
    Denom.p <- c(Bits[["cdf.max.s"]]   - Bits[["SumT0.p"]] -
                 Bits[["SumA0.mix.p"]] - Bits[["SumA0.mlm.p"]])

    if (any(Denom.p == 0)) {
      smallval <- min(Denom.p[Denom.p > 0])
      Denom.p[Denom.p == 0] <- 1e-09  # smallval
      warning("0s found in variable 'Denom.p'. Trying to fix it.")
    }
    Numer <- c(1 - pobs.mix - pstr.mix + pdip.mix -
               (if (la.mlm) rowSums(pobs.mlm) else 0) -
               (if (li.mlm) rowSums(pstr.mlm) else 0) +
               (if (ld.mlm) rowSums(pdip.mlm) else 0))
    probns <- Numer * (1 -
      (c(Bits[["SumI0.mix.p"]] + Bits[["SumI0.mlm.p"]]) +
       c(Bits[["SumD0.mix.p"]] + Bits[["SumD0.mlm.p"]])) / Denom.p)
  }  # morework

    if (!la.mlm && type.fitted %in% c("pobs.mlm")) {
      warning("No altered MLM values; returning an NA")
      return(NA)
    }
    if (!li.mlm && type.fitted %in% c("sum.mlm.i", "pstr.mlm")) {
      warning("No inflated MLM values; returning an NA")
      return(NA)
    }
    if (!ld.mlm && type.fitted %in% c("sum.mlm.d", "pdip.mlm")) {
      warning("No deflated MLM values; returning an NA")
      return(NA)
    }
    if (!la.mix && type.fitted %in% c("Pobs.mix")) {
      warning("No altered mixture values; returning an NA")
      return(NA)
    }
    if (!li.mix && type.fitted %in% c("sum.mix.i", "Pstr.mix")) {
      warning("No inflated mixture values; returning an NA")
      return(NA)
    }
    if (!ld.mix && type.fitted %in% c("sum.mix.d", "Pdip.mix")) {
      warning("No deflated mixture values; returning an NA")
      return(NA)
    }

    if (la.mix && morework) {
      tmp13 <-  # dpois() does not retain the matrix format
        dlog(matrix(a.mix,   n.obs, la.mix, byrow = TRUE),
             matrix(shape.a, n.obs, la.mix)) / (
        c(Bits[["SumA0.mix.a"]]))
      dim(tmp13) <- c(n.obs, la.mix)
      dimnames(tmp13) <- list(rownames(eta), as.character(a.mix))
      propn.mat.a <- tmp13
    }  # la.mix

    if (li.mix && morework) {
      tmp55 <-  # dpois() does not retain the matrix format
        dlog(matrix(i.mix,   n.obs, li.mix, byrow = TRUE),
             matrix(shape.i, n.obs, li.mix)) / (
        c(Bits[["SumI0.mix.i"]]))
      dim(tmp55) <- c(n.obs, li.mix)
      dimnames(tmp55) <- list(rownames(eta), as.character(i.mix))
      propn.mat.i <- tmp55  # Correct dimension
    }  # li.mix


    if (ld.mix && morework) {
      tmp55 <-  # dpois() does not retain the matrix format
        dlog(matrix(d.mix,   n.obs, ld.mix, byrow = TRUE),
             matrix(shape.d, n.obs, ld.mix)) / (
        c(Bits[["SumD0.mix.d"]]))
      dim(tmp55) <- c(n.obs, ld.mix)
      dimnames(tmp55) <- list(rownames(eta), as.character(d.mix))
      propn.mat.d <- tmp55  # Correct dimension
    }  # ld.mix

    ans <- switch(type.fitted,
      "mean"       = Bits[["mean"]],  # Unconditional mean
      "shapes"    = cbind(shape.p,
                           if (tmp3.TF[ 3]) shape.a else NULL,
                           if (tmp3.TF[ 5]) shape.i else NULL,
                           if (tmp3.TF[ 7]) shape.d else NULL),
      "pobs.mlm"   = pobs.mlm,  # aka omegamat, n x la.mlm
      "pstr.mlm"   = pstr.mlm,  # aka phimat,   n x li.mlm
      "pdip.mlm"   = pdip.mlm,  # aka psimat,   n x ld.mlm
      "pobs.mix"   = pobs.mix,  # n-vector
      "pstr.mix"   = pstr.mix,  # n-vector
      "pdip.mix"   = pdip.mix,  # n-vector
      "Pobs.mix"   = c(pobs.mix) * propn.mat.a,  # matrix
      "Pstr.mix"   = c(pstr.mix) * propn.mat.i,
      "Pdip.mix"   = c(pdip.mix) * propn.mat.d,
      "nonspecial" = probns,
      "Numer"      = Numer,
      "Denom.p"    = Denom.p,
      "sum.mlm.i"  =  pstr.mlm + Numer *
             dlog(matrix(i.mlm,   n.obs, li.mlm, byrow = TRUE),
                  matrix(shape.p, n.obs, li.mlm)) / Denom.p,
      "sum.mlm.d"  = -pdip.mlm + Numer *
             dlog(matrix(d.mlm,   n.obs, ld.mlm, byrow = TRUE),
                  matrix(shape.p, n.obs, ld.mlm)) / Denom.p,
      "sum.mix.i"  =  c(pstr.mix) * propn.mat.i + Numer *
             dlog(matrix(i.mix,   n.obs, li.mix, byrow = TRUE),
                  matrix(shape.p, n.obs, li.mix)) / Denom.p,
      "sum.mix.d"  = -c(pdip.mix) * propn.mat.d + Numer *
             dlog(matrix(d.mix,   n.obs, ld.mix, byrow = TRUE),
                  matrix(shape.p, n.obs, ld.mix)) / Denom.p,
      "ptrunc.p"   = Bits[["SumT0.p"]] + 1 - Bits[["cdf.max.s"]],
      "cdf.max.s"  = Bits[["cdf.max.s"]])  # Pr(y <= max.support)



    ynames.pobs.mlm <- as.character(a.mlm)  # Works with NULLs
    ynames.pstr.mlm <- as.character(i.mlm)  # Works with NULLs
    ynames.pdip.mlm <- as.character(d.mlm)  # Works with NULLs
    if (length(ans))
      label.cols.y(ans, NOS = NOS, colnames.y =
      switch(type.fitted,
             "shapes"   = c("shape.p", "shape.a",  # Some colns NA
                 "shape.i", "shape.d")[(tmp3.TF[c(1, 3, 5, 7)])],
             "Pobs.mix"  = as.character(a.mix),
             "sum.mix.i" = ,  #
             "Pstr.mix"  = as.character(i.mix),
             "sum.mix.d" = ,  #
             "Pdip.mix"  = as.character(d.mix),
             "pobs.mlm"  = ynames.pobs.mlm,
             "sum.mlm.i" = ,  #
             "pstr.mlm"  = ynames.pstr.mlm,
             "sum.mlm.d" = ,  #
             "pdip.mlm"  = ynames.pdip.mlm,
             extra$colnames.y)) else ans
  }, list(
    .lshape.p = lshape.p, .eshape.p = eshape.p,
    .lshape.a = lshape.a, .eshape.a = eshape.a,
    .lshape.i = lshape.i, .eshape.i = eshape.i,
    .lshape.d = lshape.d, .eshape.d = eshape.d,
    .lpstr.mix = lpstr.mix, .lpobs.mix = lpobs.mix,
    .lpdip.mix = lpdip.mix,
    .epstr.mix = epstr.mix, .epobs.mix = epobs.mix,
    .epdip.mix = epdip.mix,
    .tmp3.TF = tmp3.TF,
    .a.mix = a.mix, .i.mix = i.mix, .d.mix = d.mix,
    .a.mlm = a.mlm, .i.mlm = i.mlm, .d.mlm = d.mlm,
    .truncate = truncate, .max.support = max.support ))),
  last = eval(substitute(expression({
    pred.names <- c( .predictors.names )  # Save it
    link.names <- as.vector( .predictors.names )
    parameter.names <- names(pred.names)
    predictors.names <- NULL
    for (jay in seq(M))
      predictors.names <- c(predictors.names,
        namesof(parameter.names[jay], link.names[jay], tag = FALSE,
                earg = list()))  # This line isnt perfect; info is lost
    misc$predictors.names <- predictors.names  # Useful for coef()
    misc$link <- link.names  # 
    names(misc$link) <- parameter.names  # 


    misc$earg <- vector("list", M1)
    names(misc$earg) <- names(misc$link)
    misc$earg[[1]] <- ( .eshape.p )  # First one always there
    iptr <- 1
    if (tmp3.TF[ 2])
      misc$earg[[(iptr <- iptr + 1)]] <- list()  # multilogitlink
    if (tmp3.TF[ 3])
      misc$earg[[(iptr <- iptr + 1)]] <- ( .eshape.a )
    if (tmp3.TF[ 4])
      misc$earg[[(iptr <- iptr + 1)]] <- list()  # See below
    if (tmp3.TF[ 5])
      misc$earg[[(iptr <- iptr + 1)]] <- ( .eshape.i )
    if (tmp3.TF[ 6])
      misc$earg[[(iptr <- iptr + 1)]] <- list()  # See below
    if (tmp3.TF[ 7])
      misc$earg[[(iptr <- iptr + 1)]] <- ( .eshape.d )
    if (tmp3.TF[ 8]) {  # la.mlm
      for (ii in seq(la.mlm)) {
        misc$earg[[(iptr <- iptr + 1)]] <- list()
      }  # ii
    }  # la.mlm
    if (tmp3.TF[ 9]) {  # li.mlm
      for (ii in seq(li.mlm)) {
        misc$earg[[(iptr <- iptr + 1)]] <- list()
      }  # ii
    }  # li.mlm
    if (tmp3.TF[10]) { # ld.mlm
      for (ii in seq(ld.mlm)) {
        misc$earg[[(iptr <- iptr + 1)]] <- list()
      }  # ii
    }  # ld.mlm
  }), list(
    .lshape.p = lshape.p, .eshape.p = eshape.p,
    .lshape.a = lshape.a, .eshape.a = eshape.a,
    .lshape.i = lshape.i, .eshape.i = eshape.i,
    .lshape.d = lshape.d, .eshape.d = eshape.d,
    .lpstr.mix = lpstr.mix, .lpobs.mix = lpobs.mix,
    .lpdip.mix = lpdip.mix,
    .epstr.mix = epstr.mix, .epobs.mix = epobs.mix,
    .epdip.mix = epdip.mix,
    .predictors.names = predictors.names,
    .tmp3.TF = tmp3.TF,  # .tmp3 = tmp3,
    .a.mix = a.mix, .i.mix = i.mix, .d.mix = d.mix,
    .a.mlm = a.mlm, .i.mlm = i.mlm, .d.mlm = d.mlm,
    .truncate = truncate, .max.support = max.support ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL, summation = TRUE) {
    if (!is.matrix(eta)) eta <- as.matrix(eta)
    la.mix <- length((a.mix <- as.vector( .a.mix )))
    li.mix <- length((i.mix <- as.vector( .i.mix )))
    ld.mix <- length((d.mix <- as.vector( .d.mix )))
    la.mlm <- length((a.mlm <- as.vector( .a.mlm )))
    li.mlm <- length((i.mlm <- as.vector( .i.mlm )))
    ld.mlm <- length((d.mlm <- as.vector( .d.mlm )))
    truncate <- as.vector( .truncate )

    lall.len <- la.mix + li.mix + ld.mix + la.mlm + li.mlm + ld.mlm
    pobs.mix <- pstr.mix <- pdip.mix <- 0  # 4 rowSums()
    pobs.mlm <- pstr.mlm <- pdip.mlm <- 0  # matrix(0, NROW(eta), 1)
    shape.p <- cbind(eta2theta(eta[, 1], .lshape.p , .eshape.p ))
    ind.shape.z <- 1  # Points to shape.p only.
    shape.a <- shape.i <-
    shape.d <- shape.p  # Needed and doesnt corrupt the answer

    if (any(tmp3.TF[c(3, 5, 7)])) {  # At least one shape.[aid]
      ind.shape.z <- extra$indeta[c(1, 3, 5, 7), 'launch']  # Vectors
      ind.shape.z <- c(na.omit(ind.shape.z))  # At least one value

      shape.a <- if (!tmp3.TF[ 3]) shape.p else
        eta2theta(eta[, extra$indeta[3, 1]], .lshape.a , .eshape.a )
      shape.i <- if (!tmp3.TF[ 5]) shape.p else
        eta2theta(eta[, extra$indeta[5, 1]], .lshape.i , .eshape.i )
      shape.d <- if (!tmp3.TF[ 7]) shape.p else
        eta2theta(eta[, extra$indeta[7, 1]], .lshape.d , .eshape.d )
    }  # la.mix + li.mix + ld.mix > 0

    if (lall.len) {  # An MLM was fitted
      allprobs <-
        multilogitlink(eta[, -ind.shape.z, drop = FALSE],
                       refLevel = "(Last)",  # Make sure
                       inverse = TRUE)  # rowSums == 1
      if (anyNA(allprobs))
        warning("there are NAs here in slot linkinv")
      if (min(allprobs) == 0 || max(allprobs) == 1)
        warning("fitted probabilities numerically 0 or 1 occurred")

      Nextone <- 0  # Might not be used actually; 0, not 1
      if (tmp3.TF[ 2])
        pobs.mix <- allprobs[, (Nextone <- Nextone + 1)]
      if (tmp3.TF[ 4])
        pstr.mix <- allprobs[, (Nextone <- Nextone + 1)]
      if (tmp3.TF[ 6])
        pdip.mix <- allprobs[, (Nextone <- Nextone + 1)]

      if (tmp3.TF[ 8]) {
        ind8 <- (Nextone + 1):(Nextone + la.mlm)
        pobs.mlm <- allprobs[, ind8, drop = FALSE]
        dimnames(pobs.mlm) <- list(rownames(eta), as.character(a.mlm))
        Nextone <- Nextone + la.mlm
      }
      if (tmp3.TF[ 9]) {
        ind9 <- (Nextone + 1):(Nextone + li.mlm)
        pstr.mlm <- allprobs[, ind9, drop = FALSE]
        dimnames(pstr.mlm) <- list(rownames(eta), as.character(i.mlm))
        Nextone <- Nextone + li.mlm
      }
      if (tmp3.TF[10]) {
        ind10 <- (Nextone + 1):(Nextone + ld.mlm)
        pdip.mlm <- allprobs[, ind10, drop = FALSE]
        dimnames(pdip.mlm) <- list(rownames(eta), as.character(d.mlm))
        Nextone <- Nextone + ld.mlm  # Not needed
      }
    }  # lall.len


    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) *
        dgaitdlog(y, shape.p, log = TRUE,  # byrow.aid = F,
                  a.mix = a.mix, i.mix = i.mix, d.mix = d.mix,
                  a.mlm = a.mlm, i.mlm = i.mlm, d.mlm = d.mlm,
                  truncate = truncate,
                  max.support = as.vector( .max.support ),
                  shape.a = shape.a, shape.i = shape.i, 
                  shape.d = shape.d,
                  pobs.mix = pobs.mix, pstr.mix = pstr.mix,
                  pdip.mix = pdip.mix,
                  pobs.mlm = pobs.mlm, pstr.mlm = pstr.mlm,
                  pdip.mlm = pdip.mlm)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list(
    .lshape.p = lshape.p, .eshape.p = eshape.p,
    .lshape.a = lshape.a, .eshape.a = eshape.a,
    .lshape.i = lshape.i, .eshape.i = eshape.i,
    .lshape.d = lshape.d, .eshape.d = eshape.d,
    .lpstr.mix = lpstr.mix, .lpobs.mix = lpobs.mix,
    .lpdip.mix = lpdip.mix,
    .epstr.mix = epstr.mix, .epobs.mix = epobs.mix,
    .epdip.mix = epdip.mix,
    .a.mix = a.mix, .i.mix = i.mix, .d.mix = d.mix,
    .a.mlm = a.mlm, .i.mlm = i.mlm, .d.mlm = d.mlm,
    .truncate = truncate, .max.support = max.support ))),
  vfamily = c("gaitdlog"),
  validparams = eval(substitute(function(eta, y, extra = NULL) {
    la.mix <- length((a.mix <- as.vector( .a.mix )))
    li.mix <- length((i.mix <- as.vector( .i.mix )))
    ld.mix <- length((d.mix <- as.vector( .d.mix )))
    la.mlm <- length((a.mlm <- as.vector( .a.mlm )))
    li.mlm <- length((i.mlm <- as.vector( .i.mlm )))
    ld.mlm <- length((d.mlm <- as.vector( .d.mlm )))
    lall.len <- la.mix + li.mix + ld.mix + la.mlm + li.mlm + ld.mlm
    small. <- 1e-14
    pobs.mix <- pstr.mix <- pdip.mix <- small.   # 4 rowSums():
    pobs.mlm <- pstr.mlm <- pdip.mlm <- matrix(small., NROW(eta), 1) 
    shape.a <- shape.i <- shape.d <- 0.5  # Needed

    if (!is.matrix(eta)) eta <- as.matrix(eta)
    shape.p <-
      cbind(eta2theta(eta[, 1], .lshape.p , earg = .eshape.p ))
    ind.shape.z <- 1  # Points to shape.p only.
    if (any(tmp3.TF[c(3, 5, 7)])) {  # At least one shape.[aid]
      ind.shape.z <- extra$indeta[c(1, 3, 5, 7), 1]  # Vectors
      ind.shape.z <- c(na.omit(ind.shape.z))  # At least one value

      shape.a <- if (!tmp3.TF[ 3]) shape.p else
        eta2theta(eta[, extra$indeta[3, 1]], .lshape.a , .eshape.a )
      shape.i <- if (!tmp3.TF[ 5]) shape.p else
        eta2theta(eta[, extra$indeta[5, 1]], .lshape.i , .eshape.i )
      shape.d <- if (!tmp3.TF[ 7]) shape.p else
        eta2theta(eta[, extra$indeta[7, 1]], .lshape.d , .eshape.d )
    }  # la.mix + li.mix + ld.mix > 0

    
    if (lall.len) {  # A MLM was fitted
      allprobs <-
        multilogitlink(eta[, -ind.shape.z, drop = FALSE],
                       inverse = TRUE)  # rowSums == 1

      Nextone <- 0  # Might not be used actually; 0, not 1
      if (tmp3.TF[ 2])
        pobs.mix <- allprobs[, (Nextone <- Nextone + 1)]
      if (tmp3.TF[ 4])
        pstr.mix <- allprobs[, (Nextone <- Nextone + 1)]
      if (tmp3.TF[ 6])
        pdip.mix <- allprobs[, (Nextone <- Nextone + 1)]

      if (tmp3.TF[ 8]) {
        ind8 <- (Nextone + 1):(Nextone + la.mlm)
        pobs.mlm <- allprobs[, ind8, drop = FALSE]
        dimnames(pobs.mlm) <- list(rownames(eta), as.character(a.mlm))
        Nextone <- Nextone + la.mlm
      }
      if (tmp3.TF[ 9]) {
        ind9 <- (Nextone + 1):(Nextone + li.mlm)
        pstr.mlm <- allprobs[, ind9, drop = FALSE]
        dimnames(pstr.mlm) <- list(rownames(eta), as.character(i.mlm))
        Nextone <- Nextone + li.mlm
      }
      if (tmp3.TF[10]) {
        ind10 <- (Nextone + 1):(Nextone + ld.mlm)
        pdip.mlm <- allprobs[, ind10, drop = FALSE]
        dimnames(pdip.mlm) <- list(rownames(eta), as.character(d.mlm))
        Nextone <- Nextone + ld.mlm  # Not needed
      }
    }  # lall.len

    okay.mlm <-
      all(is.finite(pobs.mlm)) && all(0 < pobs.mlm) &&
      all(is.finite(pstr.mlm)) && all(0 < pstr.mlm) &&
      all(is.finite(pdip.mlm)) && all(0 < pdip.mlm)
    okay.mix <-
      all(is.finite(shape.p)) && all(0 < shape.p) && all(shape.p < 1) &&
      all(is.finite(shape.a)) && all(0 < shape.a) && all(shape.a < 1) &&
      all(is.finite(shape.i)) && all(0 < shape.i) && all(shape.i < 1) &&
      all(is.finite(shape.d)) && all(0 < shape.d) && all(shape.d < 1) &&
      all(is.finite(pobs.mix)) && all(0 < pobs.mix) &&
      all(is.finite(pstr.mix)) && all(0 < pstr.mix) &&
      all(is.finite(pdip.mix)) && all(0 < pdip.mix) &&
      all(pobs.mix + pstr.mix + pdip.mix +
          rowSums(pobs.mlm) + rowSums(pstr.mlm) +
          rowSums(pdip.mlm) < 1)  # Combined
    okay.mlm && okay.mix
  }, list(
    .lshape.p = lshape.p, .eshape.p = eshape.p,
    .lshape.a = lshape.a, .eshape.a = eshape.a,
    .lshape.i = lshape.i, .eshape.i = eshape.i,
    .lshape.d = lshape.d, .eshape.d = eshape.d,
    .lpstr.mix = lpstr.mix, .lpobs.mix = lpobs.mix,
    .lpdip.mix = lpdip.mix,
    .epstr.mix = epstr.mix, .epobs.mix = epobs.mix,
    .epdip.mix = epdip.mix,
    .a.mix = a.mix, .i.mix = i.mix, .d.mix = d.mix,
    .a.mlm = a.mlm, .i.mlm = i.mlm, .d.mlm = d.mlm,
    .truncate = truncate, .max.support = max.support ))),
  simslot = eval(substitute(
  function(object, nsim) {
    pwts <- if (length(pwts <- object@prior.weights) > 0)
              pwts else weights(object, type = "prior")
    if (any(pwts != 1))
      warning("ignoring prior weights")
    eta <- predict(object)
    la.mix <- length((a.mix <- as.vector( .a.mix )))
    li.mix <- length((i.mix <- as.vector( .i.mix )))
    ld.mix <- length((d.mix <- as.vector( .d.mix )))
    la.mlm <- length((a.mlm <- as.vector( .a.mlm )))
    li.mlm <- length((i.mlm <- as.vector( .i.mlm )))
    ld.mlm <- length((d.mlm <- as.vector( .d.mlm )))
    truncate <- as.vector( .truncate )

    extra <- object@extra
    lall.len <- la.mix + li.mix + ld.mix + la.mlm + li.mlm + ld.mlm
    pobs.mix <- pstr.mix <- pdip.mix <- 0  # 4 rowSums()
    pobs.mlm <- pstr.mlm <- pdip.mlm <- 0  # matrix(0, NROW(eta), 1)
    shape.p <- cbind(eta2theta(eta[, 1], .lshape.p , .eshape.p ))
    ind.shape.z <- 1  # Points to shape.p only.
    shape.a <- shape.i <-
    shape.d <- shape.p  # Needed; and answer not corrupted
    tmp3.TF <- ( .tmp3.TF )


    if (any(tmp3.TF[c(3, 5, 7)])) {  # At least one shape.[aid]
      ind.shape.z <- extra$indeta[c(1, 3, 5, 7), 'launch']  # Vectors
      ind.shape.z <- c(na.omit(ind.shape.z))  # At least one value

      shape.a <- if (!tmp3.TF[ 3]) shape.p else
        eta2theta(eta[, extra$indeta[3, 1]], .lshape.a , .eshape.a )
      shape.i <- if (!tmp3.TF[ 5]) shape.p else
        eta2theta(eta[, extra$indeta[5, 1]], .lshape.i , .eshape.i )
      shape.d <- if (!tmp3.TF[ 7]) shape.p else
        eta2theta(eta[, extra$indeta[7, 1]], .lshape.d , .eshape.d )
    }  # la.mix + li.mix + ld.mix > 0


    if (lall.len) {  # A AMLM was fitted
      allprobs <-
        multilogitlink(eta[, -ind.shape.z, drop = FALSE],
                       inverse = TRUE)  # rowSums == 1
      Nextone <- 0  # Might not be used actually; 0, not 1
      if (tmp3.TF[ 2])
        pobs.mix <- allprobs[, (Nextone <- Nextone + 1)]
      if (tmp3.TF[ 4])
        pstr.mix <- allprobs[, (Nextone <- Nextone + 1)]
      if (tmp3.TF[ 6])
        pdip.mix <- allprobs[, (Nextone <- Nextone + 1)]

      if (tmp3.TF[ 8]) {
        ind8 <- (Nextone + 1):(Nextone + la.mlm)
        pobs.mlm <- allprobs[, ind8, drop = FALSE]
        dimnames(pobs.mlm) <- list(rownames(eta), as.character(a.mlm))
        Nextone <- Nextone + la.mlm
      }
      if (tmp3.TF[ 9]) {
        ind9 <- (Nextone + 1):(Nextone + li.mlm)
        pstr.mlm <- allprobs[, ind9, drop = FALSE]
        dimnames(pstr.mlm) <- list(rownames(eta), as.character(i.mlm))
        Nextone <- Nextone + li.mlm
      }
      if (tmp3.TF[10]) {
        ind10 <- (Nextone + 1):(Nextone + ld.mlm)
        pdip.mlm <- allprobs[, ind10, drop = FALSE]
        dimnames(pdip.mlm) <- list(rownames(eta), as.character(d.mlm))
        Nextone <- Nextone + ld.mlm  # Not needed
      }
    }  # lall.len

    rgaitdlog(nsim * length(shape.p), shape.p,
               pobs.mlm = pobs.mlm, pstr.mlm = pstr.mlm,
               pobs.mix = pobs.mix, pstr.mix = pstr.mix,
               pdip.mix = pdip.mix, pdip.mlm = pdip.mlm,
               shape.a = shape.a, shape.i = shape.i, 
               shape.d = shape.d,
               a.mix = a.mix, i.mix = i.mix, d.mix = d.mix,
               a.mlm = a.mlm, i.mlm = i.mlm, d.mlm = d.mlm,
               truncate = .truncate , max.support = .max.support )
  }, list(
    .lshape.p = lshape.p, .eshape.p = eshape.p,
    .lshape.a = lshape.a, .eshape.a = eshape.a,
    .lshape.i = lshape.i, .eshape.i = eshape.i,
    .lshape.d = lshape.d, .eshape.d = eshape.d,
    .lpstr.mix = lpstr.mix, .lpobs.mix = lpobs.mix,
    .lpdip.mix = lpdip.mix,
    .epstr.mix = epstr.mix, .epobs.mix = epobs.mix,
    .epdip.mix = epdip.mix,
    .tmp3.TF = tmp3.TF,
    .a.mix = a.mix, .i.mix = i.mix, .d.mix = d.mix,
    .a.mlm = a.mlm, .i.mlm = i.mlm, .d.mlm = d.mlm,
    .truncate = truncate, .max.support = max.support ))),
  deriv = eval(substitute(expression({

    tmp3.TF <- ( .tmp3.TF )
    calA.p  <- tmp3.TF[ 2]
    calI.p  <- tmp3.TF[ 4]
    calD.p  <- tmp3.TF[ 6]
    calA.np <- tmp3.TF[ 8]
    calI.np <- tmp3.TF[ 9]
    calD.np <- tmp3.TF[10]

    Denom1.a <- Denom1.i <- Denom1.d <-
                Denom2.i <- Denom2.d <- 0  # Denom2.a is unneeded

    if (!is.matrix(eta)) eta <- as.matrix(eta)
    la.mix <- length((a.mix <- as.vector( .a.mix )))
    li.mix <- length((i.mix <- as.vector( .i.mix )))
    ld.mix <- length((d.mix <- as.vector( .d.mix )))
    la.mlm <- length((a.mlm <- as.vector( .a.mlm )))
    li.mlm <- length((i.mlm <- as.vector( .i.mlm )))
    ld.mlm <- length((d.mlm <- as.vector( .d.mlm )))
    truncate <- as.vector( .truncate )
    max.support <- as.vector( .max.support )



    
    lall.len <- la.mix + li.mix + ld.mix + la.mlm + li.mlm + ld.mlm
    pobs.mix <- pstr.mix <- pdip.mix <- 0  # 4 rowSums()
    pobs.mlm <- pstr.mlm <- pdip.mlm <- 0  # matrix(0, NROW(eta), 1)
    shape.p <- cbind(eta2theta(eta[, 1], .lshape.p , .eshape.p ))
    ind.shape.z <- 1  # Points to shape.p only.
    shape.a <- shape.i <-
    shape.d <- shape.p  # Needed; and answer not corrupted

    if (any(tmp3.TF[c(3, 5, 7)])) {  # At least one shape.[aid]
      ind.shape.z <- extra$indeta[c(1, 3, 5, 7), 'launch']  # Vectors
      ind.shape.z <- c(na.omit(ind.shape.z))  # At least one value

      shape.a <- if (!tmp3.TF[ 3]) shape.p else
        eta2theta(eta[, extra$indeta[3, 1]], .lshape.a , .eshape.a )
      shape.i <- if (!tmp3.TF[ 5]) shape.p else
        eta2theta(eta[, extra$indeta[5, 1]], .lshape.i , .eshape.i )
      shape.d <- if (!tmp3.TF[ 7]) shape.p else
        eta2theta(eta[, extra$indeta[7, 1]], .lshape.d , .eshape.d )
    }  # la.mix + li.mix + ld.mix > 0



    if (lall.len) {  # A MLM was fitted.
      allprobs <-
        multilogitlink(eta[, -ind.shape.z, drop = FALSE],
                       refLevel = "(Last)",  # Make sure
                       inverse = TRUE)  # rowSums == 1
      minprob.baseline <- min(allprobs[, ncol(allprobs)], na.rm = TRUE)
      if (anyNA(allprobs))
        warning("there are NAs here in slot linkinv")
      if (min(allprobs) == 0 || max(allprobs) == 1) {
        warning("fitted probabilities numerically 0 or 1 occurred")
      } else
      if (minprob.baseline < 0.10)
        warning("Minimum baseline (reserve) probability close to 0")
      if (control$trace)
        cat("Minimum baseline (reserve) probability = ",
            format(minprob.baseline, digits = 3), "\n")


      Nextone <- 0  # Might not be used actually; 0, not 1
      if (tmp3.TF[ 2])
        pobs.mix <- allprobs[, (Nextone <- Nextone + 1)]
      if (tmp3.TF[ 4])
        pstr.mix <- allprobs[, (Nextone <- Nextone + 1)]
      if (tmp3.TF[ 6])
        pdip.mix <- allprobs[, (Nextone <- Nextone + 1)]

      if (tmp3.TF[ 8]) {
        ind8 <- (Nextone + 1):(Nextone + la.mlm)
        pobs.mlm <- allprobs[, ind8, drop = FALSE]
        dimnames(pobs.mlm) <- list(rownames(eta), as.character(a.mlm))
        Nextone <- Nextone + la.mlm
      }
      if (tmp3.TF[ 9]) {
        ind9 <- (Nextone + 1):(Nextone + li.mlm)
        pstr.mlm <- allprobs[, ind9, drop = FALSE]
        dimnames(pstr.mlm) <- list(rownames(eta), as.character(i.mlm))
        Nextone <- Nextone + li.mlm
      }
      if (tmp3.TF[10]) {
        ind10 <- (Nextone + 1):(Nextone + ld.mlm)
        pdip.mlm <- allprobs[, ind10, drop = FALSE]
        dimnames(pdip.mlm) <- list(rownames(eta), as.character(d.mlm))
        Nextone <- Nextone + ld.mlm  # Not needed
      }
    }  # lall.len


    ltruncat <- length(truncate)
    M1 <- max(extra$indeta, na.rm = TRUE)
    NOS <- ncol(eta) / M1  # extra$NOS
    if (NOS != 1) stop("can only handle 1 response")



    is.a.mixed <- if (tmp3.TF[ 2])
      rowSums(extra$skip.mix.a) > 0 else rep(FALSE, n)
    is.i.mixed <- if (tmp3.TF[ 4])
      rowSums(extra$skip.mix.i) > 0 else rep(FALSE, n)
    is.d.mixed <- if (tmp3.TF[ 6])
      rowSums(extra$skip.mix.d) > 0 else rep(FALSE, n)
    is.a.mlmed <- if (tmp3.TF[ 8])
      rowSums(extra$skip.mlm.a) > 0 else rep(FALSE, n)
    is.i.mlmed <- if (tmp3.TF[ 9])
      rowSums(extra$skip.mlm.i) > 0 else rep(FALSE, n)
    is.d.mlmed <- if (tmp3.TF[10])
      rowSums(extra$skip.mlm.d) > 0 else rep(FALSE, n)

    is.ns <- !is.a.mlmed & !is.i.mlmed  & !is.d.mlmed  &
             !is.a.mixed & !is.i.mixed  & !is.d.mixed  # & !is.truncd


    A8.p <- -1 / log1p(-shape.p)
    A8.a <- -1 / log1p(-shape.a)
    A8.i <- -1 / log1p(-shape.i)
    prob.mlm.a <- if (la.mlm) rowSums(pobs.mlm) else 0  # scalar okay
    prob.mlm.i <- if (li.mlm) rowSums(pstr.mlm) else 0  # scalar okay
    prob.mlm.d <- if (ld.mlm) rowSums(pdip.mlm) else 0  # scalar okay
    pmf.deriv1 <- function(y, shape) {
      A8 <- -1 / log1p(-shape)
      deriv0 <- A8 * (shape^y) / y
      A8 * (shape^(y-1) - deriv0 / (1 - shape))
    }
    pmf.deriv2 <- function(y, shape) {
      A8 <- -1 / log1p(-shape)
      A8prime <- -(A8^2) / (1 - shape)
      deriv0 <- A8 * (shape^y) / y
      deriv1 <- A8 * (shape^(y-1) - deriv0 / (1 - shape))
      A8prime * (shape^(y-1) - deriv0 / (1 - shape)) +
      A8 * ((y - 1) * shape^(y - 2) - deriv0 / (1 - shape)^2 -
            deriv1 / (1 - shape)) 
    }


    sumD.mix.1a.p <- sumD.mix.2a.p <- matrix(0, n, NOS)
    if (la.mix > 0) {  # \calA_p
      DA.mix.0mat.a <-  # Matches naming convention further below
      DA.mix.1mat.a <- matrix(0, n, la.mix)
      for (jay in seq(la.mix)) {
        aval <- a.mix[jay]
        sumD.mix.1a.p <- sumD.mix.1a.p + pmf.deriv1(aval, shape.p)
        sumD.mix.2a.p <- sumD.mix.2a.p + pmf.deriv2(aval, shape.p)
        pmf.a <- dlog(aval, shape.a)
        DA.mix.0mat.a[, jay] <- pmf.a
        DA.mix.1mat.a[, jay] <- pmf.deriv1(aval, shape.a)
      }
      Denom1.a <- rowSums(DA.mix.1mat.a)  # aka sumD.mix.1a.a
    }  # la.mix > 0





    if (li.mix) {
      DI.mix.0mat.i <-  # wrt inflated distribution
      DI.mix.1mat.i <- DI.mix.2mat.i <- matrix(0, n, li.mix)
      DP.mix.0mat.i <-  # wrt parent distribution
      DP.mix.1mat.i <- DP.mix.2mat.i <- matrix(0, n, li.mix)
        for (jay in seq(li.mix)) {
          ival <- i.mix[jay]
          pmf.i <- dlog(ival, shape.i)
          DI.mix.0mat.i[, jay] <- pmf.i
          DI.mix.1mat.i[, jay] <- pmf.deriv1(ival, shape.i)
          DI.mix.2mat.i[, jay] <- pmf.deriv2(ival, shape.i)
          pmf.p <- dlog(ival, shape.p)
          DP.mix.0mat.i[, jay] <- pmf.p
          DP.mix.1mat.i[, jay] <- pmf.deriv1(ival, shape.p)
          DP.mix.2mat.i[, jay] <- pmf.deriv2(ival, shape.p)
        }  # jay
      Denom1.i <- rowSums(DI.mix.1mat.i)
      Denom2.i <- rowSums(DI.mix.2mat.i)
    }  # li.mix


    if (ld.mix) {
      DD.mix.0mat.d <-  # wrt deflated distribution
      DD.mix.1mat.d <- DD.mix.2mat.d <- matrix(0, n, ld.mix)
      DP.mix.0mat.d <-  # wrt parent distribution
      DP.mix.1mat.d <- DP.mix.2mat.d <- matrix(0, n, ld.mix)
        for (jay in seq(ld.mix)) {
          dval <- d.mix[jay]
          pmf.d <- dlog(dval, shape.d)
          DD.mix.0mat.d[, jay] <- pmf.d
          DD.mix.1mat.d[, jay] <- pmf.deriv1(dval, shape.d)
          DD.mix.2mat.d[, jay] <- pmf.deriv2(dval, shape.d)
          pmf.p <- dlog(dval, shape.p)
          DP.mix.0mat.d[, jay] <- pmf.p
          DP.mix.1mat.d[, jay] <- pmf.deriv1(dval, shape.p)
          DP.mix.2mat.d[, jay] <- pmf.deriv2(dval, shape.p)
        }  # jay
      Denom1.d <- rowSums(DD.mix.1mat.d)
      Denom2.d <- rowSums(DD.mix.2mat.d)
    }  # ld.mix

    Bits <- moments.gaitdcombo.log(shape.p,
              pobs.mix = pobs.mix, pobs.mlm = pobs.mlm,
              pstr.mix = pstr.mix, pstr.mlm = pstr.mlm,
              pdip.mix = pdip.mix, pdip.mlm = pdip.mlm,
              a.mix = a.mix, i.mix = i.mix, d.mix = d.mix,
              a.mlm = a.mlm, i.mlm = i.mlm, d.mlm = d.mlm,
              shape.a = shape.a, shape.i = shape.i,
              shape.d = shape.d,
              truncate = truncate, max.support = max.support)


    sumD.mlm.1a.p <- sumD.mlm.2a.p <- matrix(0, n, NOS)
    if (la.mlm)
      for (aval in a.mlm) {
        sumD.mlm.1a.p <- sumD.mlm.1a.p + pmf.deriv1(aval, shape.p)
        sumD.mlm.2a.p <- sumD.mlm.2a.p + pmf.deriv2(aval, shape.p)
      }


    Denom0.p <- c(Bits[["cdf.max.s"]]   - Bits[["SumT0.p"]] -
                  Bits[["SumA0.mix.p"]] - Bits[["SumA0.mlm.p"]])
    Numer <- 1 - pobs.mix - pstr.mix - prob.mlm.a - prob.mlm.i +
                 pdip.mix +            prob.mlm.d
    Denom0.a <- c(Bits[["SumA0.mix.a"]])  # Not .p
    Denom0.i <- c(Bits[["SumI0.mix.i"]])
    Denom0.d <- c(Bits[["SumD0.mix.d"]])


 

    Dp.mlm.0Mat.i <-  # wrt parent distribution
    Dp.mlm.1Mat.i <- Dp.mlm.2Mat.i <- matrix(0, n, NOS)
    if (li.mlm > 0) {
      Dp.mlm.0Mat.i <-  # wrt parent distribution
      Dp.mlm.1Mat.i <- Dp.mlm.2Mat.i <- matrix(0, n, li.mlm)
      for (jay in seq(li.mlm)) {
        ival <- i.mlm[jay]
        pmf.p <- dlog(ival, shape.p)
        Dp.mlm.0Mat.i[, jay] <- pmf.p
        Dp.mlm.1Mat.i[, jay] <- pmf.deriv1(ival, shape.p)
        Dp.mlm.2Mat.i[, jay] <- pmf.deriv2(ival, shape.p)
      }  # jay
    }  # li.mlm



    Dp.mlm.0Mat.d <-  # wrt parent distribution
    Dp.mlm.1Mat.d <- Dp.mlm.2Mat.d <- matrix(0, n, NOS)
    if (ld.mlm > 0) {
      Dp.mlm.0Mat.d <-  # wrt parent distribution
      Dp.mlm.1Mat.d <- Dp.mlm.2Mat.d <- matrix(0, n, ld.mlm)
      for (jay in seq(ld.mlm)) {
        dval <- d.mlm[jay]
        pmf.p <- dlog(dval, shape.p)
        Dp.mlm.0Mat.d[, jay] <- pmf.p
        Dp.mlm.1Mat.d[, jay] <- pmf.deriv1(dval, shape.p)
        Dp.mlm.2Mat.d[, jay] <- pmf.deriv2(dval, shape.p)
      }  # jay
    }  # ld.mlm





    

    sumD.1t.p <- sumD.2t.p <-
    sumD.1t.a <- sumD.2t.a <-
    sumD.1t.i <- sumD.2t.i <-
    sumD.1t.d <- sumD.2t.d <- matrix(0, n, NOS)
    if (ltruncat)
      for (tval in truncate) {
        sumD.1t.p <- sumD.1t.p + pmf.deriv1(tval, shape.p)
        sumD.2t.p <- sumD.2t.p + pmf.deriv2(tval, shape.p)
        sumD.1t.a <- sumD.1t.a + pmf.deriv1(tval, shape.a)
        sumD.2t.a <- sumD.2t.a + pmf.deriv2(tval, shape.a)
        sumD.1t.i <- sumD.1t.i + pmf.deriv1(tval, shape.i)
        sumD.2t.i <- sumD.2t.i + pmf.deriv2(tval, shape.i)
        sumD.1t.d <- sumD.1t.d + pmf.deriv1(tval, shape.d)
        sumD.2t.d <- sumD.2t.d + pmf.deriv2(tval, shape.d)
      }




      
    if (is.finite(max.support)) {
    tmp1.p <- A8.p * (shape.p^max.support -
      (1 - plog(max.support, shape.p))) / (1 - shape.p)
    sumD.1t.p <- sumD.1t.p + tmp1.p
    sumD.2t.p <- sumD.2t.p + (A8.p / (1 - shape.p)) * (
      (shape.p^max.support) / (1 - shape.p) +
      max.support * shape.p^(max.support - 1) -
      (1 - plog(max.support, shape.p)) / (1 - shape.p) - 2 * tmp1.p)


    tmp1.a <- A8.a * (shape.a^max.support -
      (1 - plog(max.support, shape.a))) / (1 - shape.a)
    sumD.1t.a <- sumD.1t.a + tmp1.a
    sumD.2t.a <- sumD.2t.a + (A8.a / (1 - shape.a)) * (
      (shape.a^max.support) / (1 - shape.a) +
      max.support * shape.a^(max.support - 1) -
      (1 - plog(max.support, shape.a)) / (1 - shape.a) - 2 * tmp1.a)
        
    tmp1.i <- A8.i * (shape.i^max.support -
      (1 - plog(max.support, shape.i))) / (1 - shape.i)
    sumD.1t.i <- sumD.1t.i + tmp1.i
    sumD.2t.i <- sumD.2t.i + (A8.i / (1 - shape.i)) * (
      (shape.i^max.support) / (1 - shape.i) +
      max.support * shape.i^(max.support - 1) -
      (1 - plog(max.support, shape.i)) / (1 - shape.i) - 2 * tmp1.i)
    }  # is.finite(max.support)





    Denom1.p <- c(-sumD.1t.p - sumD.mlm.1a.p - sumD.mix.1a.p)
    Denom2.p <- c(-sumD.2t.p - sumD.mlm.2a.p - sumD.mix.2a.p)




    d0B.PI.mlm <- Dp.mlm.0Mat.i / Denom0.p
    d1B.PI.mlm <- Dp.mlm.1Mat.i / Denom0.p -  # This is most general
                  Dp.mlm.0Mat.i * Denom1.p / Denom0.p^2
    d2B.PI.mlm <-     Dp.mlm.2Mat.i / Denom0.p -
                  2 * Dp.mlm.1Mat.i * Denom1.p / Denom0.p^2 -
                      Dp.mlm.0Mat.i * Denom2.p / Denom0.p^2 +
                  2 * Dp.mlm.0Mat.i * (Denom1.p^2) / Denom0.p^3





    d0B.PD.mlm <- Dp.mlm.0Mat.d / Denom0.p
    d1B.PD.mlm <- Dp.mlm.1Mat.d / Denom0.p -  # This is most general
                  Dp.mlm.0Mat.d * Denom1.p / Denom0.p^2
    d2B.PD.mlm <-     Dp.mlm.2Mat.d / Denom0.p -
                  2 * Dp.mlm.1Mat.d * Denom1.p / Denom0.p^2 -
                      Dp.mlm.0Mat.d * Denom2.p / Denom0.p^2 +
                  2 * Dp.mlm.0Mat.d * (Denom1.p^2) / Denom0.p^3


    DELTA.i.mlm <- if (li.mlm > 0) {
      Numer * d0B.PI.mlm + pstr.mlm  # n x li.mlm.
    } else {
      matrix(0, n, 1)  # If li.mlm == 0, for rowSums().
    }

    DELTA.d.mlm <- if (ld.mlm > 0) {
      Numer * d0B.PD.mlm - pdip.mlm  # n x ld.mlm.
    } else {
      matrix(0, n, 1)  # If ld.mlm == 0, for rowSums().
    }


    if (li.mix > 0) {
      d0A.i <- DI.mix.0mat.i / Denom0.i
      d0B.PI.mix <- DP.mix.0mat.i / Denom0.p
      DELTA.i.mix <- Numer * d0B.PI.mix + pstr.mix * d0A.i

      d1A.i <- (DI.mix.1mat.i - DI.mix.0mat.i *
                Denom1.i / Denom0.i) / Denom0.i
      d2A.i <- (DI.mix.2mat.i - (2 * DI.mix.1mat.i * Denom1.i +
                DI.mix.0mat.i * Denom2.i) / Denom0.i +
            2 * DI.mix.0mat.i * (Denom1.i / Denom0.i)^2) / Denom0.i

      d1B.PI.mix <- DP.mix.1mat.i / Denom0.p -
                    DP.mix.0mat.i * Denom1.p / Denom0.p^2
      d2B.PI.mix <-     DP.mix.2mat.i /  Denom0.p -
                    2 * DP.mix.1mat.i *  Denom1.p    / Denom0.p^2 -
                        DP.mix.0mat.i *  Denom2.p    / Denom0.p^2 +
                    2 * DP.mix.0mat.i * (Denom1.p^2) / Denom0.p^3
    }  # li.mix > 0



    if (ld.mix > 0) {
      d0A.d <- DD.mix.0mat.d / Denom0.d
      d0B.PD.mix <- DP.mix.0mat.d / Denom0.p
      DELTA.d.mix <- Numer * d0B.PD.mix - pdip.mix * d0A.d

      d1A.d <- (DD.mix.1mat.d - DD.mix.0mat.d *
                Denom1.d / Denom0.d) / Denom0.d
      d2A.d <- (DD.mix.2mat.d - (2 * DD.mix.1mat.d * Denom1.d +
                DD.mix.0mat.d * Denom2.d) / Denom0.d +
            2 * DD.mix.0mat.d * (Denom1.d / Denom0.d)^2) / Denom0.d

      d1B.PD.mix <- DP.mix.1mat.d / Denom0.p -
                    DP.mix.0mat.d * Denom1.p / Denom0.p^2
      d2B.PD.mix <-     DP.mix.2mat.d /  Denom0.p -
                    2 * DP.mix.1mat.d *  Denom1.p    / Denom0.p^2 -
                        DP.mix.0mat.d *  Denom2.p    / Denom0.p^2 +
                    2 * DP.mix.0mat.d * (Denom1.p^2) / Denom0.p^3
    }  # ld.mix > 0




    if (la.mix) {
      d0A.a <- DA.mix.0mat.a / Denom0.a
      d1A.a <- DA.mix.1mat.a / Denom0.a -
               DA.mix.0mat.a * Denom1.a / Denom0.a^2
    }  # la.mix




    dl.dshape.p <- -A8.p / (1 - shape.p) + y / shape.p
    dl.dshape.p[!is.ns] <- 0  # For is.a.mixed & is.a.mlmed

    dl.dshape.a <- dl.dshape.i <- dl.dshape.d <- numeric(n)
    dl.dpstr.mix <- (-1) / Numer  # \notin A, I, T, D
    dl.dpstr.mix[is.a.mixed] <- 0
    dl.dpstr.mix[is.a.mlmed] <- 0
    dl.dpdip.mix <- (+1) / Numer  # \notin A, I, T, D
    dl.dpdip.mix[is.a.mixed] <- 0
    dl.dpdip.mix[is.a.mlmed] <- 0
    dl.dpobs.mix <- numeric(n)  # 0 for \calA_{np}
    dl.dpobs.mix[is.ns] <- (-1) / Numer[is.ns]
    dl.dpobs.mlm <-
    dl.dpstr.mlm <- matrix(0, n, 1)  # May be unneeded
    dl.dpdip.mlm <- matrix(0, n, max(1, ld.mlm))  # Initzed if used.
    dl.dpdip.mlm[is.ns, ] <- 1 / Numer[is.ns]



    if (tmp3.TF[ 8] && la.mlm) {  # aka \calA_{np}
      dl.dpobs.mlm <- matrix(-1 / Numer, n, la.mlm)  # \notin calS
      dl.dpobs.mlm[!is.ns, ] <- 0  # For a.mix only really
      for (jay in seq(la.mlm)) {
        aval <- a.mlm[jay]
        is.alt.j.mlm <- extra$skip.mlm.a[, jay]  # Logical vector
        tmp7a <- 1 / pobs.mlm[is.alt.j.mlm, jay]
        dl.dpobs.mlm[is.alt.j.mlm, jay] <- tmp7a
      }  # jay
    }  # la.mlm



    dl.dshape.p[is.ns] <- dl.dshape.p[is.ns] -
                           (Denom1.p / Denom0.p)[is.ns]


    
    if (tmp3.TF[ 9] && li.mlm > 0) {  # aka \calI_{np}
      dl.dpstr.mlm <- matrix(-1 / Numer, n, li.mlm)
      dl.dpstr.mlm[!is.ns, ] <- 0  # For a.mlm and a.mix

      for (jay in seq(li.mlm)) {
        is.inf.j.mlm <- extra$skip.mlm.i[, jay]  # Logical vector
        tmp7i <- Numer * d1B.PI.mlm[, jay] / DELTA.i.mlm[, jay]
        dl.dshape.p[is.inf.j.mlm] <- tmp7i[is.inf.j.mlm]


        tmp9i <- d0B.PI.mlm[, jay] / DELTA.i.mlm[, jay]
        n.tmp <- -tmp9i[is.inf.j.mlm]
        p.tmp <- +tmp9i[is.inf.j.mlm]
        if (tmp3.TF[ 8] && la.mlm) dl.dpobs.mlm[is.inf.j.mlm, ] <- n.tmp
        if (tmp3.TF[ 2] && la.mix) dl.dpobs.mix[is.inf.j.mlm  ] <- n.tmp
        if (tmp3.TF[ 4] && li.mix) dl.dpstr.mix[is.inf.j.mlm  ] <- n.tmp
        if (tmp3.TF[10] && ld.mlm) dl.dpdip.mlm[is.inf.j.mlm, ] <- p.tmp
        if (tmp3.TF[ 6] && ld.mix) dl.dpdip.mix[is.inf.j.mlm  ] <- p.tmp


        tmp8 <- (1 - d0B.PI.mlm[, jay]) / DELTA.i.mlm[, jay]
        dl.dpstr.mlm[is.inf.j.mlm, ] <- n.tmp  # tmp9[is.inf.j.mlm]
        dl.dpstr.mlm[is.inf.j.mlm, jay] <- tmp8[is.inf.j.mlm]
      }  # jay
    }  # li.mlm > 0






    if (tmp3.TF[10] && ld.mlm > 0) {  # aka \calD_{np}

      for (jay in seq(ld.mlm)) {
        is.def.j.mlm <- extra$skip.mlm.d[, jay]  # Logical vector
        tmp7d <- Numer * d1B.PD.mlm[, jay] / DELTA.d.mlm[, jay]
        dl.dshape.p[is.def.j.mlm] <- tmp7d[is.def.j.mlm]  # 20211020

 
        tmp9d <- d0B.PD.mlm[, jay] / DELTA.d.mlm[, jay]
        p.tmp <- +tmp9d[is.def.j.mlm]
        n.tmp <- -tmp9d[is.def.j.mlm]
        if (tmp3.TF[ 9] && li.mlm) dl.dpstr.mlm[is.def.j.mlm, ] <- n.tmp
        if (tmp3.TF[ 4] && li.mix) dl.dpstr.mix[is.def.j.mlm  ] <- n.tmp
        if (tmp3.TF[ 8] && la.mlm) dl.dpobs.mlm[is.def.j.mlm, ] <- n.tmp
        if (tmp3.TF[ 2] && la.mix) dl.dpobs.mix[is.def.j.mlm  ] <- n.tmp
        if (tmp3.TF[ 6] && ld.mix) dl.dpdip.mix[is.def.j.mlm  ] <- p.tmp
                                   dl.dpdip.mlm[is.def.j.mlm, ] <- p.tmp
        dl.dpdip.mlm[is.def.j.mlm, jay] <-
        dl.dpdip.mlm[is.def.j.mlm, jay] -
        1 / DELTA.d.mlm[is.def.j.mlm, jay]

      }  # jay
    }  # ld.mlm > 0




    



    if (tmp3.TF[ 2] && la.mix) {  # aka \calA_{p}
      dl.dpobs.mix[is.a.mixed] <- 1 / pobs.mix[is.a.mixed]

      if (tmp3.TF[ 3] && la.mix > 1)
        for (jay in seq(la.mix)) {
          is.alt.j.mix <- extra$skip.mix.a[, jay]  # Logical vector
          tmp2 <- d1A.a[, jay] / d0A.a[, jay]
          dl.dshape.a[is.alt.j.mix] <- tmp2[is.alt.j.mix]  # ccc.
        }  # jay
    }  # la.mix




    if (tmp3.TF[ 4] && li.mix > 0) {  # aka \calI_{p}
      for (jay in seq(li.mix)) {
        ival <- i.mix[jay]
        is.inf.j.mix <- extra$skip.mix.i[, jay]  # Logical vector
        tmp7b <- Numer * d1B.PI.mix[, jay] / DELTA.i.mix[, jay]
        dl.dshape.p[is.inf.j.mix] <- tmp7b[is.inf.j.mix]
        tmp8 <- (d0A.i[, jay] - d0B.PI.mix[, jay]) / DELTA.i.mix[, jay]
        dl.dpstr.mix[is.inf.j.mix] <- tmp8[is.inf.j.mix]
        if (li.mix > 1) {
          tmp2 <- pstr.mix * d1A.i[, jay] / DELTA.i.mix[, jay]
          dl.dshape.i[is.inf.j.mix] <- tmp2[is.inf.j.mix]
        }





        tmp9i <- d0B.PI.mix[, jay] / DELTA.i.mix[, jay]
        n.tmp <- -tmp9i[is.inf.j.mix]
        p.tmp <- +tmp9i[is.inf.j.mix]
        if (tmp3.TF[ 2] && la.mix) dl.dpobs.mix[is.inf.j.mix  ] <- n.tmp
        if (tmp3.TF[ 8] && la.mlm) dl.dpobs.mlm[is.inf.j.mix, ] <- n.tmp
        if (tmp3.TF[ 9] && li.mlm) dl.dpstr.mlm[is.inf.j.mix, ] <- n.tmp
        if (tmp3.TF[10] && ld.mlm) dl.dpdip.mlm[is.inf.j.mix, ] <- p.tmp
        if (tmp3.TF[ 6] && ld.mix) dl.dpdip.mix[is.inf.j.mix  ] <- p.tmp

      }  # jay
    }  # li.mix > 0




    if (tmp3.TF[ 6] && ld.mix > 0) {  # aka \calD_{p}
      for (jay in seq(ld.mix)) {
        dval <- d.mix[jay]
        is.def.j.mix <- extra$skip.mix.d[, jay]  # Logical vector
        tmp7b <- Numer * d1B.PD.mix[, jay] / DELTA.d.mix[, jay]
        dl.dshape.p[is.def.j.mix] <- tmp7b[is.def.j.mix]
        tmp8 <- (d0B.PD.mix[, jay] - d0A.d[, jay]) / DELTA.d.mix[, jay]
        dl.dpdip.mix[is.def.j.mix] <- tmp8[is.def.j.mix]

        if (ld.mix > 1) {
          tmp2 <- (-pdip.mix) * d1A.d[, jay] / DELTA.d.mix[, jay]
          dl.dshape.d[is.def.j.mix] <- tmp2[is.def.j.mix]
        }


        tmp9d <- d0B.PD.mix[, jay] / DELTA.d.mix[, jay]
        n.tmp <- -tmp9d[is.def.j.mix]
        p.tmp <- +tmp9d[is.def.j.mix]
        if (tmp3.TF[ 9] && li.mlm) dl.dpstr.mlm[is.def.j.mix, ] <- n.tmp
        if (tmp3.TF[ 4] && li.mix) dl.dpstr.mix[is.def.j.mix  ] <- n.tmp
        if (tmp3.TF[ 8] && la.mlm) dl.dpobs.mlm[is.def.j.mix, ] <- n.tmp
        if (tmp3.TF[ 2] && la.mix) dl.dpobs.mix[is.def.j.mix  ] <- n.tmp
        if (tmp3.TF[10] && ld.mlm) dl.dpdip.mlm[is.def.j.mix, ] <- p.tmp

      }  # jay
   }  # ld.mix > 0







    new.ansd <- matrix(0, n, M)  # Same dimension as eta
    tmp3.TF <- !is.na(rowSums(extra$indeta))



    if (lall.len) {  # An MLM fitted
      all6.dldp <- cbind(if (tmp3.TF[ 2]) dl.dpobs.mix else NULL,
                         if (tmp3.TF[ 4]) dl.dpstr.mix else NULL,
                         if (tmp3.TF[ 6]) dl.dpdip.mix else NULL,
                         if (tmp3.TF[ 8]) dl.dpobs.mlm else NULL,
                         if (tmp3.TF[ 9]) dl.dpstr.mlm else NULL,
                         if (tmp3.TF[10]) dl.dpdip.mlm else NULL)




      rSs.tmp <- rowSums(allprobs[, -ncol(allprobs), drop = FALSE] *
                         all6.dldp)
      new.ansd[, -ind.shape.z] <- allprobs[, -ncol(allprobs)] *
                                   (all6.dldp - rSs.tmp)
    }  # lall.len


    

    dshape.p.deta <- dtheta.deta(shape.p, .lshape.p , .eshape.p )
    if (tmp3.TF[ 3])
      dshape.a.deta <- dtheta.deta(shape.a, .lshape.a , .eshape.a )
    if (tmp3.TF[ 5])
      dshape.i.deta <- dtheta.deta(shape.i, .lshape.i , .eshape.i )
    if (tmp3.TF[ 7])
      dshape.d.deta <- dtheta.deta(shape.d, .lshape.d , .eshape.d )
    new.ansd[, 1] <- dl.dshape.p * dshape.p.deta
    if (tmp3.TF[ 3])
      new.ansd[, extra$indeta[3, 1]] <- dl.dshape.a * dshape.a.deta
    if (tmp3.TF[ 5])
      new.ansd[, extra$indeta[5, 1]] <- dl.dshape.i * dshape.i.deta
    if (tmp3.TF[ 7])
      new.ansd[, extra$indeta[7, 1]] <- dl.dshape.d * dshape.d.deta
    onecoln.indeta <- extra$indeta[1:7, ]  # One coln params only
    onecoln.indeta <- na.omit(onecoln.indeta)  # Only those present
    allcnames <- c(rownames(onecoln.indeta),
                   as.character(c(a.mlm, i.mlm, d.mlm)))
    colnames(new.ansd) <- allcnames

 


    c(w) * new.ansd
  }), list(
    .lshape.p = lshape.p, .eshape.p = eshape.p,
    .lshape.a = lshape.a, .eshape.a = eshape.a,
    .lshape.i = lshape.i, .eshape.i = eshape.i,
    .lshape.d = lshape.d, .eshape.d = eshape.d,
    .lpstr.mix = lpstr.mix, .lpobs.mix = lpobs.mix,
    .lpdip.mix = lpdip.mix,
    .epstr.mix = epstr.mix, .epobs.mix = epobs.mix,
    .epdip.mix = epdip.mix,
    .a.mix = a.mix, .i.mix = i.mix, .d.mix = d.mix,
    .a.mlm = a.mlm, .i.mlm = i.mlm, .d.mlm = d.mlm,
    .tmp3.TF = tmp3.TF,  # .tmp3 = tmp3,
    .truncate = truncate, .max.support = max.support ))),

  weight = eval(substitute(expression({  # gaitdlog



    wz <- matrix(0, n, M * (M + 1) / 2)  # The complete size
    mean.true.p <- A8.p * shape.p / (1 - shape.p)
    cond.EY.p <- c(mean.true.p - Bits[["SumT1.p"]] -
        Bits[["SumI1.mlm.p"]] - Bits[["SumI1.mix.p"]] -
        Bits[["SumD1.mlm.p"]] - Bits[["SumD1.mix.p"]] -  # 20211109
        Bits[["SumA1.mlm.p"]] - Bits[["SumA1.mix.p"]]) / c(
        Denom0.p -
        Bits[["SumD0.mix.p"]] - Bits[["SumD0.mlm.p"]] -  # 20211109
        Bits[["SumI0.mix.p"]] - Bits[["SumI0.mlm.p"]])
    probns <- Numer * (1 -
      (c(Bits[["SumI0.mix.p"]] + Bits[["SumI0.mlm.p"]]) +
       c(Bits[["SumD0.mix.p"]] + Bits[["SumD0.mlm.p"]])) / Denom0.p)



    if (min(probns) < 0 || 1 < max(probns))
      stop("variable 'probns' for P(nonspecial) is out of range")







    zero0n <- numeric(n)
    ned2l.dpobs.mix.shape.p <- zero0n  # mB overwritten below [4279]
    ned2l.dpobs.mix.shape.a <- zero0n  # Fini; (2, 3) element
    ned2l.dpobs.mix.shape.i <- zero0n  # mB overwritten below
    ned2l.dpobs.mix.shape.d <- zero0n  # mB overwritten below
    ned2l.dpstr.mix.shape.p <- zero0n  # Optional (1, 4) element
    ned2l.dpstr.mix.shape.a <- zero0n  # Final; nothing to do
    ned2l.dpstr.mix.shape.i <- zero0n  # mB overwritten below
    ned2l.dpstr.mix.shape.d <- zero0n  # mB overwritten below
    ned2l.dpdip.mix.shape.p <- zero0n  # Optional (1, 6) element




      posn.pobs.mix <- as.vector(extra$indeta[ 2, 'launch'])
      posn.shape.a <- as.vector(extra$indeta[ 3, 'launch'])
      posn.pstr.mix <- as.vector(extra$indeta[ 4, 'launch'])
      posn.shape.i <- as.vector(extra$indeta[ 5, 'launch'])
      posn.pdip.mix <- as.vector(extra$indeta[ 6, 'launch'])
      posn.shape.d <- as.vector(extra$indeta[ 7, 'launch'])
      posn.pobs.mlm <- as.vector(extra$indeta[ 8, 'launch'])
      posn.pstr.mlm <- as.vector(extra$indeta[ 9, 'launch'])
      posn.pdip.mlm <- as.vector(extra$indeta[10, 'launch'])





    ned2l.dpdip.mix2         <-  # Elt (6, 6)
    ned2l.dpstr.mix2         <-  # Elt (4, 4). Unchanged by deflation.
    ned2l.dpobs.mlm.pstr.mix <-  # Elts (4, >=8). (((09)))
    ned2l.dpobs.mix.pstr.mix <- +probns / Numer^2  # ccc Elt (2, 4)
    if (all(c(la.mix, li.mlm) > 0))  # (((08)))
    ned2l.dpobs.mix.pstr.mlm <- matrix( probns / Numer^2, n, li.mlm)
    if (all(c(li.mix, li.mlm) > 0))  # (((10)))
    ned2l.dpstr.mix.pstr.mlm <- matrix( probns / Numer^2, n, li.mlm)
    if (all(c(ld.mix, ld.mlm) > 0))  # (((21)))
    ned2l.dpdip.mix.pdip.mlm <- matrix( probns / Numer^2, n, ld.mlm)


    ned2l.dpobs.mlm.pdip.mix <-  # Elts (6, >=8). (((19)))
    ned2l.dpstr.mix.pdip.mix <-  # Elt (4, 6)
    ned2l.dpobs.mix.pdip.mix <- -probns / Numer^2  # ccc Elt (2, 6)
    if (all(c(la.mix, ld.mlm) > 0))  # (((17)))
    ned2l.dpobs.mix.pdip.mlm <- matrix(-probns / Numer^2, n, ld.mlm)
    if (all(c(li.mix, ld.mlm) > 0))  # (((18)))
    ned2l.dpstr.mix.pdip.mlm <- matrix(-probns / Numer^2, n, ld.mlm)
    if (all(c(ld.mix, li.mlm) > 0))  # (((20)))
    ned2l.dpdip.mix.pstr.mlm <- matrix(-probns / Numer^2, n, li.mlm)


    


    ned2l.dshape.p2 <- probns * (cond.EY.p / shape.p^2 +  # ccc
                                 A8.p * (1 - A8.p) / (1 - shape.p)^2 +
    Denom2.p / Denom0.p - (Denom1.p / Denom0.p)^2) + 
    (if (tmp3.TF[ 4] && li.mix) Numer *
    rowSums(Numer * (d1B.PI.mix^2) / DELTA.i.mix - d2B.PI.mix) else 0) +
    (if (tmp3.TF[ 9] && li.mlm) Numer *
    rowSums(Numer * (d1B.PI.mlm^2) / DELTA.i.mlm - d2B.PI.mlm) else 0) +
    (if (tmp3.TF[ 6] && ld.mix) Numer *
    rowSums(Numer * (d1B.PD.mix^2) / DELTA.d.mix - d2B.PD.mix) else 0) +
    (if (tmp3.TF[10] && ld.mlm) Numer *  # nnn.
    rowSums(Numer * (d1B.PD.mlm^2) / DELTA.d.mlm - d2B.PD.mlm) else 0)


    wz[, iam(1, 1, M)] <- ned2l.dshape.p2 * dshape.p.deta^2


    ned2l.dpobs.mix2 <- 1 / pobs.mix + probns / Numer^2
    if (tmp3.TF[ 4] && li.mix > 0) {
      ned2l.dpobs.mix2 <-  # More just below, ccc
      ned2l.dpobs.mix2 + rowSums(d0B.PI.mix^2 / DELTA.i.mix)
    }
    if (tmp3.TF[ 9] && li.mlm > 0) {
      ned2l.dpobs.mix2 <-  # ccc.
      ned2l.dpobs.mix2 + rowSums(d0B.PI.mlm^2 / DELTA.i.mlm)
    }
    if (tmp3.TF[ 6] && ld.mix > 0) {
      ned2l.dpobs.mix2 <-  # nnn
      ned2l.dpobs.mix2 + rowSums(d0B.PD.mix^2 / DELTA.d.mix)
    }
    if (tmp3.TF[10] && ld.mlm > 0) {
      ned2l.dpobs.mix2 <-  # nnn
      ned2l.dpobs.mix2 + rowSums(d0B.PD.mlm^2 / DELTA.d.mlm)
    }
    if (tmp3.TF[ 2] && la.mix > 0)
      wz[, iam(2, 2, M)] <- ned2l.dpobs.mix2  # Link done later


    if (tmp3.TF[ 3] && la.mix > 1) {
      ned2l.dshape.a2 <- pobs.mix * (
        rowSums((DA.mix.1mat.a^2) / DA.mix.0mat.a) / Denom0.a -
        (Denom1.a / Denom0.a)^2)  # ccc.
      wz[, iam(3, 3, M)] <- ned2l.dshape.a2 * dshape.a.deta^2
    }




    if (tmp3.TF[ 4] && li.mix > 0) {
      ned2l.dpstr.mix2 <-
      ned2l.dpstr.mix2 +
        rowSums((d0A.i - d0B.PI.mix)^2 / DELTA.i.mix)

      if (tmp3.TF[ 2] && la.mix > 0)
        ned2l.dpobs.mix.shape.p <-
        ned2l.dpobs.mix.shape.p +
          rowSums(d1B.PI.mix * (1 - Numer * d0B.PI.mix / DELTA.i.mix))

      ned2l.dpstr.mix.shape.p <-
      ned2l.dpstr.mix.shape.p + rowSums(
        d1B.PI.mix * (1 + Numer * (d0A.i - d0B.PI.mix) / DELTA.i.mix))

      if (tmp3.TF[ 6])
        ned2l.dpdip.mix.shape.p <-
        ned2l.dpdip.mix.shape.p - rowSums(
          d1B.PI.mix * (1 - Numer * d0B.PI.mix / DELTA.i.mix))

      if (all(tmp3.TF[c(2, 4)]))
        ned2l.dpobs.mix.pstr.mix <-  # ccc
        ned2l.dpobs.mix.pstr.mix +
          rowSums(-d0B.PI.mix * (d0A.i - d0B.PI.mix) / DELTA.i.mix)

      if (all(tmp3.TF[c(4, 6)]))
        ned2l.dpstr.mix.pdip.mix <-
        ned2l.dpstr.mix.pdip.mix + rowSums(
          d0B.PI.mix * (d0A.i - d0B.PI.mix) / DELTA.i.mix)

      if (!is.na(posn.pdip.mix)) {
        ned2l.dpdip.mix2 <-
        ned2l.dpdip.mix2 + rowSums(d0B.PI.mix^2 / DELTA.i.mix)
      }
    }  # (tmp3.TF[ 4] && li.mix > 0)







    if (all(tmp3.TF[c(2, 4,  9)])) {  # was la.mix > 0 & DELTA.i.mix
      ned2l.dpobs.mix.pstr.mix <-  # ccc
      ned2l.dpobs.mix.pstr.mix + rowSums(d0B.PI.mlm^2 / DELTA.i.mlm)
    }
    if (all(tmp3.TF[c(2, 4,  6)])) {  # ==  ld.mix > 0 & DELTA.d.mix
      ned2l.dpobs.mix.pstr.mix <-  # nnn
      ned2l.dpobs.mix.pstr.mix + rowSums(d0B.PD.mix^2 / DELTA.d.mix)
    }
    if (all(tmp3.TF[c(2, 4, 10)])) {  # ==  ld.mlm > 0 & DELTA.d.mlm
      ned2l.dpobs.mix.pstr.mix <-  # nnn.
      ned2l.dpobs.mix.pstr.mix + rowSums(d0B.PD.mlm^2 / DELTA.d.mlm)
    }
    if (!is.na(posn.pobs.mix) && !is.na(posn.pstr.mix))
      wz[, iam(posn.pobs.mix, posn.pstr.mix, M)] <-
        ned2l.dpobs.mix.pstr.mix  # Link done later



    if (all(tmp3.TF[c(2, 6)]))
      ned2l.dpobs.mix.pdip.mix <-  # nnn
      ned2l.dpobs.mix.pdip.mix +
        rowSums( d0B.PD.mix * (d0A.d - d0B.PD.mix) / DELTA.d.mix)

    if (all(tmp3.TF[c(2, 6,  9)])) {  # ==  li.mlm > 0 & DELTA.i.mix
      ned2l.dpobs.mix.pdip.mix <-  # nnn
      ned2l.dpobs.mix.pdip.mix - rowSums(d0B.PI.mlm^2 / DELTA.i.mlm)
    }
    if (all(tmp3.TF[c(2, 6,  4)])) {  # ==  li.mix > 0 & DELTA.i.mix
      ned2l.dpobs.mix.pdip.mix <-  # nnn
      ned2l.dpobs.mix.pdip.mix - rowSums(d0B.PI.mix^2 / DELTA.i.mix)
    }
    if (all(tmp3.TF[c(2, 6, 10)])) {  # ==  ld.mlm > 0 & DELTA.d.mlm
      ned2l.dpobs.mix.pdip.mix <-  # nnn.
      ned2l.dpobs.mix.pdip.mix - rowSums(d0B.PD.mlm^2 / DELTA.d.mlm)
    }
    if (!is.na(posn.pobs.mix) && !is.na(posn.pdip.mix))
      wz[, iam(posn.pobs.mix, posn.pdip.mix, M)] <-
        ned2l.dpobs.mix.pdip.mix  # Link done later





    if (tmp3.TF[ 5] && li.mix > 1) {  # \calI_{p}, includes \theta_i.

      ned2l.dshape.p.shape.i <- pstr.mix * Numer *
        rowSums(d1A.i * d1B.PI.mix / DELTA.i.mix)  # ccc.
      wz[, iam(1, posn.shape.i, M)] <- ned2l.dshape.p.shape.i *
        dshape.p.deta * dshape.i.deta  # All links done here

      ned2l.dshape.i2 <- pstr.mix *
        rowSums(pstr.mix * (d1A.i^2) / DELTA.i.mix - d2A.i)  # ccc.
      wz[, iam(posn.shape.i, posn.shape.i, M)] <-
        ned2l.dshape.i2 * dshape.i.deta^2

      if (tmp3.TF[ 2]) {  # tmp3.TF[ 4] is TRUE, given tmp3.TF[ 5]
        ned2l.dpobs.mix.shape.i <-
          rowSums(-pstr.mix * d1A.i * d0B.PI.mix / DELTA.i.mix)  # ccc.
        wz[, iam(posn.pobs.mix, posn.shape.i, M)] <-
          ned2l.dpobs.mix.shape.i  # * dshape.i.deta done later
      }

      if (tmp3.TF[ 4]) {
        ned2l.dpstr.mix.shape.i <- rowSums(  # ccc.
          d1A.i * (pstr.mix * (d0A.i - d0B.PI.mix) / DELTA.i.mix - 1))
        wz[, iam(posn.pstr.mix, posn.shape.i, M)] <-
          ned2l.dpstr.mix.shape.i  # * dshape.i.deta done later
      }

      if (all(tmp3.TF[c(5, 6)])) {
        ned2l.dpdip.mix.shape.i <- rowSums(
          (-pstr.mix) * d0B.PI.mix * d1A.i / DELTA.i.mix)
        wz[, iam(posn.pdip.mix, posn.shape.i, M)] <-
          ned2l.dpdip.mix.shape.i  # link done later
      }

      if (tmp3.TF[ 8]) {
        ned2l.dpobs.mlm.shape.i <- rowSums(
          -pstr.mix * d0B.PI.mix * d1A.i / DELTA.i.mix)  # ccc.
        for (uuu in seq(la.mlm))
          wz[, iam(posn.pobs.mlm - 1 + uuu, posn.shape.i, M)] <-
            ned2l.dpobs.mlm.shape.i  # * dshape.i.deta done later
      }


    }  # (tmp3.TF[ 5] && li.mix > 1)





    if (tmp3.TF[ 6] && ld.mix > 0) {  # \calD_{p}, maybe w. \theta_d

      if (tmp3.TF[ 2] && la.mix > 0)
        ned2l.dpobs.mix.shape.p <-
        ned2l.dpobs.mix.shape.p +
          rowSums(d1B.PD.mix * (1 - Numer * d0B.PD.mix / DELTA.d.mix))

      ned2l.dpstr.mix.shape.p <-
      ned2l.dpstr.mix.shape.p + rowSums(
        d1B.PD.mix * (1 - Numer * d0B.PD.mix / DELTA.d.mix))

      ned2l.dpdip.mix.shape.p <-
      ned2l.dpdip.mix.shape.p - rowSums(
        d1B.PD.mix * (1 + Numer * (d0A.d - d0B.PD.mix) / DELTA.d.mix))

      if (!is.na(posn.pstr.mix)) {
        ned2l.dpstr.mix2 <-
        ned2l.dpstr.mix2 + rowSums(d0B.PD.mix^2 / DELTA.d.mix)
      }

      if (all(tmp3.TF[c(4, 6)]))
        ned2l.dpstr.mix.pdip.mix <-
        ned2l.dpstr.mix.pdip.mix + rowSums(
          d0B.PD.mix * (d0A.d - d0B.PD.mix) / DELTA.d.mix)

      ned2l.dpdip.mix2 <-
      ned2l.dpdip.mix2 +
        rowSums((d0A.d - d0B.PD.mix)^2 / DELTA.d.mix)

    }  # (tmp3.TF[ 6] && ld.mix > 0)





    if (tmp3.TF[ 7] && ld.mix > 1) {  # \calD_{p}, includes \theta_d


      ned2l.dshape.p.shape.d <- (-pdip.mix) * Numer *
        rowSums(d1A.d * d1B.PD.mix / DELTA.d.mix)  # nnn.
      wz[, iam(1, posn.shape.d, M)] <- ned2l.dshape.p.shape.d *
        dshape.p.deta * dshape.d.deta  # All links done here

      if (tmp3.TF[ 2]) {  # tmp3.TF[ 6] is TRUE, given tmp3.TF[ 7]
        ned2l.dpobs.mix.shape.d <-
          rowSums(pdip.mix * d1A.d * d0B.PD.mix / DELTA.d.mix)  # nnn.
        wz[, iam(posn.pobs.mix, posn.shape.d, M)] <-
          ned2l.dpobs.mix.shape.d  # link done later
      }

      if (tmp3.TF[ 4]) {
        ned2l.dpstr.mix.shape.d <- rowSums(
          pdip.mix * d1A.d * d0B.PD.mix / DELTA.d.mix)
        wz[, iam(posn.pstr.mix, posn.shape.d, M)] <-
          ned2l.dpstr.mix.shape.d  # * dshape.i.deta done later
      }

        ned2l.dpdip.mix.shape.d <- rowSums(
          d1A.d * (1 + pdip.mix * (d0A.d - d0B.PD.mix) / DELTA.d.mix))
        wz[, iam(posn.pdip.mix, posn.shape.d, M)] <-
          ned2l.dpdip.mix.shape.d  # * dshape.d.deta done later

      ned2l.dshape.d2 <- pdip.mix *
        rowSums(pdip.mix * (d1A.d^2) / DELTA.d.mix + d2A.d)  # nnn.
      wz[, iam(posn.shape.d, posn.shape.d, M)] <-
        ned2l.dshape.d2 * dshape.d.deta^2

      if (tmp3.TF[ 8]) {
        ned2l.dpobs.mlm.shape.d <- rowSums(
           pdip.mix * d0B.PD.mix * d1A.d / DELTA.d.mix)  # nnn.
        for (uuu in seq(la.mlm))
          wz[, iam(posn.pobs.mlm - 1 + uuu, posn.shape.d, M)] <-
            ned2l.dpobs.mlm.shape.d  # * dshape.d.deta done later
      }

    }  # (tmp3.TF[ 7] && ld.mix > 1)




        
    if (tmp3.TF[ 9] && li.mlm > 0) {  # \calI_{np}, includes \phi_s.

      if (la.mix && tmp3.TF[ 2])
        ned2l.dpobs.mix.shape.p <-  # ccc
        ned2l.dpobs.mix.shape.p +
          rowSums(d1B.PI.mlm * (1 - Numer * d0B.PI.mlm / DELTA.i.mlm))

      ned2l.dpstr.mix.shape.p <-  # ccc.
      ned2l.dpstr.mix.shape.p + rowSums(
        d1B.PI.mlm * (1 - Numer * d0B.PI.mlm / DELTA.i.mlm))

      if (tmp3.TF[ 6])
        ned2l.dpdip.mix.shape.p <-
        ned2l.dpdip.mix.shape.p - rowSums(
          d1B.PI.mlm * (1 - Numer * d0B.PI.mlm / DELTA.i.mlm))

      if (!is.na(posn.pstr.mix)) {
        ned2l.dpstr.mix2 <-
        ned2l.dpstr.mix2 + rowSums(d0B.PI.mlm^2 / DELTA.i.mlm)
      }

      if (all(tmp3.TF[c(4, 6)]))
        ned2l.dpstr.mix.pdip.mix <-
        ned2l.dpstr.mix.pdip.mix - rowSums(d0B.PI.mlm^2 / DELTA.i.mlm)

      if (!is.na(posn.pdip.mix)) {
        ned2l.dpdip.mix2 <-
        ned2l.dpdip.mix2 + rowSums(d0B.PI.mlm^2 / DELTA.i.mlm)
      }

    }  # tmp3.TF[ 9] && li.mlm > 0



        
    if (tmp3.TF[10] && ld.mlm > 0) {  # \calD_{np}, includes \psi_s.


      if (la.mix && tmp3.TF[ 2])
        ned2l.dpobs.mix.shape.p <-  # nnn.
        ned2l.dpobs.mix.shape.p +
          rowSums(d1B.PD.mlm * (1 - Numer * d0B.PD.mlm / DELTA.d.mlm))

      ned2l.dpstr.mix.shape.p <-  # nnn.
      ned2l.dpstr.mix.shape.p + rowSums(
        d1B.PD.mlm * (1 - Numer * d0B.PD.mlm / DELTA.d.mlm))

      if (tmp3.TF[ 6])
        ned2l.dpdip.mix.shape.p <-
        ned2l.dpdip.mix.shape.p - rowSums(
          d1B.PD.mlm * (1 - Numer * d0B.PD.mlm / DELTA.d.mlm))

      if (!is.na(posn.pstr.mix)) {
        ned2l.dpstr.mix2 <-
        ned2l.dpstr.mix2 + rowSums(d0B.PD.mlm^2 / DELTA.d.mlm)
      }

      if (all(tmp3.TF[c(4, 6)]))
        ned2l.dpstr.mix.pdip.mix <-
        ned2l.dpstr.mix.pdip.mix - rowSums(d0B.PD.mlm^2 / DELTA.d.mlm)

      if (!is.na(posn.pdip.mix)) {
        ned2l.dpdip.mix2 <-
        ned2l.dpdip.mix2 + rowSums(d0B.PD.mlm^2 / DELTA.d.mlm)
      }



    }  # tmp3.TF[10] && ld.mlm > 0






      




    if (!is.na(posn.pobs.mix))  # Optional (1, 2) element:
      wz[, iam(1, posn.pobs.mix, M)] <-
        ned2l.dpobs.mix.shape.p  # One link done later
    if (!is.na(posn.pstr.mix))  # Optional (1, 4) element
      wz[, iam(1, posn.pstr.mix, M)] <-
        ned2l.dpstr.mix.shape.p  # One link done later
    if (!is.na(posn.pdip.mix))  # Optional (1, 6) element
      wz[, iam(1, posn.pdip.mix, M)] <-
        ned2l.dpdip.mix.shape.p  # One link done later
    if (!is.na(posn.pstr.mix) &&
        !is.na(posn.pdip.mix))  # Optional (4, 6) element
      wz[, iam(posn.pstr.mix, posn.pdip.mix, M)] <-
        ned2l.dpstr.mix.pdip.mix  # Links done later zz1

    if (!is.na(posn.pstr.mix))  # Optional (4, 4) element
      wz[, iam(posn.pstr.mix,  # Link done later
               posn.pstr.mix, M)] <- ned2l.dpstr.mix2

    if (!is.na(posn.pdip.mix))  # Optional (6, 6) element
      wz[, iam(posn.pdip.mix,  # Link done later
               posn.pdip.mix, M)] <- ned2l.dpdip.mix2







    if (tmp3.TF[ 8] && la.mlm) {  # \calA_{np}, includes \omega_s
      ofset <- posn.pobs.mlm - 1  # 7 for GAITD combo
      for (uuu in seq(la.mlm)) {  # Diagonal elts only
        wz[, iam(ofset + uuu,
                 ofset + uuu, M)] <- 1 / pobs.mlm[, uuu]
      }  # uuu

      tmp8a <- probns / Numer^2
      if (tmp3.TF[ 4] && li.mix)
        tmp8a <- tmp8a + rowSums((d0B.PI.mix^2) / DELTA.i.mix)
      if (tmp3.TF[ 9] && li.mlm)
        tmp8a <- tmp8a + rowSums((d0B.PI.mlm^2) / DELTA.i.mlm)
      if (tmp3.TF[ 6] && ld.mix)
        tmp8a <- tmp8a + rowSums((d0B.PD.mix^2) / DELTA.d.mix)
      if (tmp3.TF[10] && ld.mlm)
        tmp8a <- tmp8a + rowSums((d0B.PD.mlm^2) / DELTA.d.mlm)
      for (uuu in seq(la.mlm))  # All elts
        for (vvv in uuu:la.mlm)
          wz[, iam(ofset + uuu, ofset + vvv, M)] <-
          wz[, iam(ofset + uuu, ofset + vvv, M)] + tmp8a  # All elts
    }  # la.mlm


 

    if (tmp3.TF[ 8] && la.mlm) {

      init0.i.val <- init0.d.val <- 0
      if (tmp3.TF[ 9] && li.mlm) init0.i.val <-
        rowSums(d1B.PI.mlm * (1 - Numer * d0B.PI.mlm / DELTA.i.mlm))
      if (tmp3.TF[10] && ld.mlm) init0.d.val <-
        rowSums(d1B.PD.mlm * (1 - Numer * d0B.PD.mlm / DELTA.d.mlm))
      ned2l.dpobs.mlm.shape.p <- init0.i.val + init0.d.val  # Vector

      if (tmp3.TF[ 4] && li.mix)
        ned2l.dpobs.mlm.shape.p <-
        ned2l.dpobs.mlm.shape.p + rowSums(
          d1B.PI.mix * (1 - Numer * d0B.PI.mix / DELTA.i.mix))
      if (tmp3.TF[ 6] && ld.mix)
        ned2l.dpobs.mlm.shape.p <-
        ned2l.dpobs.mlm.shape.p + rowSums(  # nnn
          d1B.PD.mix * (1 - Numer * d0B.PD.mix / DELTA.d.mix))

      ofset <- posn.pobs.mlm - 1  # 5 for combo
      for (vvv in seq(la.mlm))  # ccc.
        wz[, iam(1, ofset + vvv, M)] <- ned2l.dpobs.mlm.shape.p
    }  # la.mlm > 0





    if (tmp3.TF[ 9] && li.mlm > 0) {  # \calI_{np}, includes \phi_s
      init0.val <- probns / Numer^2
      if (li.mix)
        init0.val <- init0.val + rowSums((d0B.PI.mix^2) / DELTA.i.mix)
      if (ld.mix)  # nnn
        init0.val <- init0.val + rowSums((d0B.PD.mix^2) / DELTA.d.mix)
      if (ld.mlm)  # nnn
        init0.val <- init0.val + rowSums((d0B.PD.mlm^2) / DELTA.d.mlm)
      ned2l.dpstr.mlm2 <-
        matrix(init0.val, n, li.mlm * (li.mlm + 1) / 2)
      for (uuu in seq(li.mlm))
        for (sss in seq(li.mlm))
          ned2l.dpstr.mlm2[, iam(uuu, uuu, li.mlm)] <-
          ned2l.dpstr.mlm2[, iam(uuu, uuu, li.mlm)] +
            ((sss == uuu) - d0B.PI.mlm[, sss])^2 / DELTA.i.mlm[, sss]
      if (li.mlm > 1) {
        for (uuu in seq(li.mlm - 1))
          for (vvv in (uuu + 1):li.mlm)
            for (sss in seq(li.mlm))
              ned2l.dpstr.mlm2[, iam(uuu, vvv, li.mlm)] <-
              ned2l.dpstr.mlm2[, iam(uuu, vvv, li.mlm)] +
              ((sss == uuu) - d0B.PI.mlm[, sss]) *
              ((sss == vvv) - d0B.PI.mlm[, sss]) / DELTA.i.mlm[, sss]
      }  # if (li.mlm > 1)

      ofset <- posn.pstr.mlm - 1
      for (uuu in seq(li.mlm))
        for (vvv in uuu:li.mlm)
          wz[, iam(ofset + uuu, ofset + vvv, M)] <-
            ned2l.dpstr.mlm2[, iam(uuu, vvv, li.mlm)] 
    }  # li.mlm > 0






    if (tmp3.TF[10] && ld.mlm > 0) {  # \calD_{np}, includes \psi_s
      init0.val <- probns / Numer^2
      if (ld.mix)
        init0.val <- init0.val + rowSums((d0B.PD.mix^2) / DELTA.d.mix)
      if (li.mix)
        init0.val <- init0.val + rowSums((d0B.PI.mix^2) / DELTA.i.mix)
      if (li.mlm)
        init0.val <- init0.val + rowSums((d0B.PI.mlm^2) / DELTA.i.mlm)
      ned2l.dpdip.mlm2 <-
        matrix(init0.val, n, ld.mlm * (ld.mlm + 1) / 2)
      for (uuu in seq(ld.mlm))
        for (sss in seq(ld.mlm))
          ned2l.dpdip.mlm2[, iam(uuu, uuu, ld.mlm)] <-
          ned2l.dpdip.mlm2[, iam(uuu, uuu, ld.mlm)] +
            (d0B.PD.mlm[, sss] - (sss == uuu))^2 / DELTA.d.mlm[, sss]
      if (ld.mlm > 1) {
        for (uuu in seq(ld.mlm - 1))
          for (vvv in (uuu + 1):ld.mlm)
            for (sss in seq(ld.mlm))
              ned2l.dpdip.mlm2[, iam(uuu, vvv, ld.mlm)] <-
              ned2l.dpdip.mlm2[, iam(uuu, vvv, ld.mlm)] +
              (d0B.PD.mlm[, sss] - (sss == uuu)) *
              (d0B.PD.mlm[, sss] - (sss == vvv)) / DELTA.d.mlm[, sss]
      }  # if (ld.mlm > 1)

      ofset <- posn.pdip.mlm - 1
      for (uuu in seq(ld.mlm))
        for (vvv in uuu:ld.mlm)
          wz[, iam(ofset + uuu, ofset + vvv, M)] <-
            ned2l.dpdip.mlm2[, iam(uuu, vvv, ld.mlm)] 
    }  # ld.mlm > 0









    if (tmp3.TF[ 9] && li.mlm > 0) {
      ned2l.dpstr.mlm.theta.p <- matrix(0, n, li.mlm)
      for (vvv in seq(li.mlm))
        for (sss in seq(li.mlm))
          ned2l.dpstr.mlm.theta.p[, vvv] <-
          ned2l.dpstr.mlm.theta.p[, vvv] +
          d1B.PI.mlm[, sss] * (1 + Numer *
          (max(0, sss == vvv) - d0B.PI.mlm[, sss]) / (
          DELTA.i.mlm[, sss]))
      if (li.mix && tmp3.TF[ 4])
        ned2l.dpstr.mlm.theta.p <-
        ned2l.dpstr.mlm.theta.p +
        rowSums(d1B.PI.mix * (1 - Numer * d0B.PI.mix / DELTA.i.mix))
      if (ld.mix && tmp3.TF[ 6])
        ned2l.dpstr.mlm.theta.p <-  # nnn
        ned2l.dpstr.mlm.theta.p +
        rowSums(d1B.PD.mix * (1 - Numer * d0B.PD.mix / DELTA.d.mix))
      if (ld.mlm && tmp3.TF[10])
        ned2l.dpstr.mlm.theta.p <-  # nnn.
        ned2l.dpstr.mlm.theta.p +
        rowSums(d1B.PD.mlm * (1 - Numer * d0B.PD.mlm / DELTA.d.mlm))
      ofset <- posn.pstr.mlm - 1
      for (vvv in seq(li.mlm))  # ccc.
        wz[, iam(1, ofset + vvv, M)] <- ned2l.dpstr.mlm.theta.p[, vvv]
    }  # li.mlm > 0




    if (tmp3.TF[10] && ld.mlm > 0) {
      ned2l.dpdip.mlm.theta.p <- matrix(0, n, ld.mlm)
      for (vvv in seq(ld.mlm))
        for (sss in seq(ld.mlm))
          ned2l.dpdip.mlm.theta.p[, vvv] <-
          ned2l.dpdip.mlm.theta.p[, vvv] -  # Minus
          d1B.PD.mlm[, sss] * (1 + Numer *
          (max(0, sss == vvv) - d0B.PD.mlm[, sss]) / (
          DELTA.d.mlm[, sss]))
      if (ld.mix && tmp3.TF[ 6])
        ned2l.dpdip.mlm.theta.p <-
        ned2l.dpdip.mlm.theta.p -  # Minus
        rowSums(d1B.PD.mix * (1 - Numer * d0B.PD.mix / DELTA.d.mix))
      if (li.mix && tmp3.TF[ 4])
        ned2l.dpdip.mlm.theta.p <-
        ned2l.dpdip.mlm.theta.p -  # Minus
        rowSums(d1B.PI.mix * (1 - Numer * d0B.PI.mix / DELTA.i.mix))
      if (li.mlm && tmp3.TF[ 9])
        ned2l.dpdip.mlm.theta.p <-  # nnn.
        ned2l.dpdip.mlm.theta.p -  # Minus
        rowSums(d1B.PI.mlm * (1 - Numer * d0B.PI.mlm / DELTA.i.mlm))
      ofset <- posn.pdip.mlm - 1
      for (vvv in seq(ld.mlm))  # nnn.
        wz[, iam(1, ofset + vvv, M)] <- ned2l.dpdip.mlm.theta.p[, vvv]
    }  # ld.mlm > 0








    if (li.mlm && li.mix > 1) {

      ned2l.dpstr.mlm.theta.i <-  # Not a matrix, just a vector
        rowSums(-pstr.mix * d0B.PI.mix * d1A.i / DELTA.i.mix)

      for (vvv in seq(li.mlm))
        wz[, iam(posn.shape.i, posn.pstr.mlm - 1 + vvv, M)] <-
          ned2l.dpstr.mlm.theta.i  # ccc.
    }  # li.mlm && li.mix > 1





    if (ld.mlm && ld.mix > 1) {

      ned2l.dpdip.mlm.theta.d <-  # Not a matrix, just a vector
        rowSums(pdip.mix * d0B.PD.mix * d1A.d / DELTA.d.mix)

      for (vvv in seq(ld.mlm))
        wz[, iam(posn.shape.d, posn.pdip.mlm - 1 + vvv, M)] <-
          ned2l.dpdip.mlm.theta.d  # nnn.
    }  # ld.mlm && ld.mix > 1



    if (ld.mlm && li.mix > 1) {

      ned2l.dpdip.mlm.theta.i <-  # Not a matrix, just a vector
        rowSums(-pstr.mix * d0B.PI.mix * d1A.i / DELTA.i.mix)

      for (vvv in seq(ld.mlm))
        wz[, iam(posn.shape.i, posn.pdip.mlm - 1 + vvv, M)] <-
          ned2l.dpdip.mlm.theta.i  # nnn.
    }  # ld.mlm && li.mix > 1



   


    if (li.mlm && ld.mix > 1) {

      ned2l.dpstr.mlm.theta.d <-  # Not a matrix, just a vector
        rowSums(pdip.mix * d0B.PD.mix * d1A.d / DELTA.d.mix)

      for (vvv in seq(li.mlm))
        wz[, iam(posn.shape.d, posn.pstr.mlm - 1 + vvv, M)] <-
          ned2l.dpstr.mlm.theta.d  # nnn.
    }  # li.mlm && ld.mix > 1






    if (all(c(la.mlm, li.mlm) > 0)) {
      ned2l.dpobs.mlm.pstr.mlm <-
        array(probns / Numer^2, c(n, la.mlm, li.mlm))
      for (uuu in seq(la.mlm))
        for (vvv in seq(li.mlm))
          for (sss in seq(li.mlm))
            ned2l.dpobs.mlm.pstr.mlm[, uuu, vvv] <- 
            ned2l.dpobs.mlm.pstr.mlm[, uuu, vvv] - d0B.PI.mlm[, sss] *
              ((sss == vvv) - d0B.PI.mlm[, sss]) / DELTA.i.mlm[, sss]
      if (tmp3.TF[ 4] && li.mix)
        ned2l.dpobs.mlm.pstr.mlm <-
        ned2l.dpobs.mlm.pstr.mlm + rowSums(d0B.PI.mix^2 / DELTA.i.mix)
      if (tmp3.TF[ 6] && ld.mix)
        ned2l.dpobs.mlm.pstr.mlm <-  # nnn
        ned2l.dpobs.mlm.pstr.mlm + rowSums(d0B.PD.mix^2 / DELTA.d.mix)
      if (tmp3.TF[10] && ld.mlm)
        ned2l.dpobs.mlm.pstr.mlm <-  # nnn
        ned2l.dpobs.mlm.pstr.mlm + rowSums(d0B.PD.mlm^2 / DELTA.d.mlm)
      ofset.pobs <- posn.pobs.mlm - 1
      ofset.pstr <- posn.pstr.mlm - 1
      for (uuu in seq(la.mlm))
        for (vvv in seq(li.mlm))
          wz[, iam(ofset.pobs + uuu, ofset.pstr + vvv, M)] <-
            ned2l.dpobs.mlm.pstr.mlm[, uuu, vvv] 
    }  # all(c(la.mlm, li.mlm) > 0)







    if (all(c(li.mlm, ld.mlm) > 0)) {
      ned2l.dpstr.mlm.pdip.mlm <-
        array(-probns / Numer^2, c(n, li.mlm, ld.mlm))
      for (uuu in seq(li.mlm))
        for (vvv in seq(ld.mlm))
          for (sss in seq(li.mlm))
            ned2l.dpstr.mlm.pdip.mlm[, uuu, vvv] <- 
            ned2l.dpstr.mlm.pdip.mlm[, uuu, vvv] + d0B.PI.mlm[, sss] *
              ((sss == uuu) - d0B.PI.mlm[, sss]) / DELTA.i.mlm[, sss]
      for (uuu in seq(li.mlm))
        for (vvv in seq(ld.mlm))
          for (sss in seq(ld.mlm))
            ned2l.dpstr.mlm.pdip.mlm[, uuu, vvv] <- 
            ned2l.dpstr.mlm.pdip.mlm[, uuu, vvv] + d0B.PD.mlm[, sss] *
              ((sss == vvv) - d0B.PD.mlm[, sss]) / DELTA.d.mlm[, sss]
      if (tmp3.TF[ 4] && li.mix)
        ned2l.dpstr.mlm.pdip.mlm <-
        ned2l.dpstr.mlm.pdip.mlm - rowSums(d0B.PI.mix^2 / DELTA.i.mix)
      if (tmp3.TF[ 6] && ld.mix)
        ned2l.dpstr.mlm.pdip.mlm <-  # nnn.
        ned2l.dpstr.mlm.pdip.mlm - rowSums(d0B.PD.mix^2 / DELTA.d.mix)
      ofset.pstr <- posn.pstr.mlm - 1
      ofset.pdip <- posn.pdip.mlm - 1
      for (uuu in seq(li.mlm))
        for (vvv in seq(ld.mlm))
          wz[, iam(ofset.pstr + uuu, ofset.pdip + vvv, M)] <-
            ned2l.dpstr.mlm.pdip.mlm[, uuu, vvv] 
    }  # all(c(li.mlm, ld.mlm) > 0)







    if (all(c(la.mlm, ld.mlm) > 0)) {
      ned2l.dpobs.mlm.pdip.mlm <-
        array(-probns / Numer^2, c(n, la.mlm, ld.mlm))
      for (uuu in seq(la.mlm))
        for (vvv in seq(ld.mlm))
          for (sss in seq(ld.mlm))
            ned2l.dpobs.mlm.pdip.mlm[, uuu, vvv] <- 
            ned2l.dpobs.mlm.pdip.mlm[, uuu, vvv] + d0B.PD.mlm[, sss] *
              ((sss == vvv) - d0B.PD.mlm[, sss]) / DELTA.d.mlm[, sss]
      if (tmp3.TF[ 4] && li.mix)
        ned2l.dpobs.mlm.pdip.mlm <-
        ned2l.dpobs.mlm.pdip.mlm - rowSums(d0B.PI.mix^2 / DELTA.i.mix)
      if (tmp3.TF[ 9] && li.mlm)
        ned2l.dpobs.mlm.pdip.mlm <-
        ned2l.dpobs.mlm.pdip.mlm - rowSums(d0B.PI.mlm^2 / DELTA.i.mlm)
      if (tmp3.TF[ 6] && ld.mix)
        ned2l.dpobs.mlm.pdip.mlm <-
        ned2l.dpobs.mlm.pdip.mlm - rowSums(d0B.PD.mix^2 / DELTA.d.mix)
      ofset.pobs <- posn.pobs.mlm - 1
      ofset.pdip <- posn.pdip.mlm - 1
      for (uuu in seq(la.mlm))
        for (vvv in seq(ld.mlm))
          wz[, iam(ofset.pobs + uuu, ofset.pdip + vvv, M)] <-
            ned2l.dpobs.mlm.pdip.mlm[, uuu, vvv] 
    }  # all(c(la.mlm, li.mlm) > 0)






    if (all(c(la.mix, la.mlm) > 0)) {
      ned2l.dpobs.mix.pobs.mlm <- probns / Numer^2  # Initialize
      if (li.mix)  # tmp3.TF[ 4]
        ned2l.dpobs.mix.pobs.mlm <-
        ned2l.dpobs.mix.pobs.mlm + rowSums(d0B.PI.mix^2 / DELTA.i.mix)
      if (li.mlm)  # tmp3.TF[ 7]
        ned2l.dpobs.mix.pobs.mlm <-
        ned2l.dpobs.mix.pobs.mlm + rowSums(d0B.PI.mlm^2 / DELTA.i.mlm)
      if (ld.mix)  # tmp3.TF[ 6]   nnn
        ned2l.dpobs.mix.pobs.mlm <-
        ned2l.dpobs.mix.pobs.mlm + rowSums(d0B.PD.mix^2 / DELTA.d.mix)
      if (ld.mlm)  # tmp3.TF[10]   nnn
        ned2l.dpobs.mix.pobs.mlm <-
        ned2l.dpobs.mix.pobs.mlm + rowSums(d0B.PD.mlm^2 / DELTA.d.mlm)

      for (uuu in seq(la.mlm))  # ccc.
        wz[, iam(posn.pobs.mix, posn.pobs.mlm - 1 + uuu, M)] <-
          ned2l.dpobs.mix.pobs.mlm  # Link done later
    }



    if (all(c(la.mix, li.mlm) > 0)) {  # all(tmp3.TF[c(2, 9)])
      if (li.mix)  # tmp3.TF[ 4]
        ned2l.dpobs.mix.pstr.mlm <-
        ned2l.dpobs.mix.pstr.mlm + rowSums(d0B.PI.mix^2 / DELTA.i.mix)
      if (ld.mix)  # tmp3.TF[ 6]
        ned2l.dpobs.mix.pstr.mlm <-  # nnn
        ned2l.dpobs.mix.pstr.mlm + rowSums(d0B.PD.mix^2 / DELTA.d.mix)
      if (ld.mlm)  # tmp3.TF[10]
        ned2l.dpobs.mix.pstr.mlm <-  # nnn; + is correct, not -
        ned2l.dpobs.mix.pstr.mlm + rowSums(d0B.PD.mlm^2 / DELTA.d.mlm)
      for (uuu in seq(li.mlm))
        for (sss in seq(li.mlm))
          ned2l.dpobs.mix.pstr.mlm[, uuu] <-
          ned2l.dpobs.mix.pstr.mlm[, uuu] -
            ((sss == uuu) - d0B.PI.mlm[, sss]) *
                            d0B.PI.mlm[, sss] / DELTA.i.mlm[, sss]
      for (uuu in seq(li.mlm))  # ccc.
        wz[, iam(posn.pobs.mix,
                 posn.pstr.mlm - 1 + uuu, M)] <-
          ned2l.dpobs.mix.pstr.mlm[, uuu]  # Link done later
    }  # all(c(la.mix, li.mlm) > 0)






    if (all(c(la.mix, ld.mlm) > 0)) {  # all(tmp3.TF[c(2, 10)])
      if (li.mix)  # tmp3.TF[ 4]
        ned2l.dpobs.mix.pdip.mlm <-
        ned2l.dpobs.mix.pdip.mlm - rowSums(d0B.PI.mix^2 / DELTA.i.mix)
      if (li.mlm)  # tmp3.TF[ 9]
        ned2l.dpobs.mix.pdip.mlm <-
        ned2l.dpobs.mix.pdip.mlm - rowSums(d0B.PI.mlm^2 / DELTA.i.mlm)
      if (ld.mix)  # tmp3.TF[ 6]
        ned2l.dpobs.mix.pdip.mlm <-
        ned2l.dpobs.mix.pdip.mlm - rowSums(d0B.PD.mix^2 / DELTA.d.mix)
      for (uuu in seq(ld.mlm))
        for (sss in seq(ld.mlm))
          ned2l.dpobs.mix.pdip.mlm[, uuu] <-
          ned2l.dpobs.mix.pdip.mlm[, uuu] +
            ((sss == uuu) - d0B.PD.mlm[, sss]) *
                            d0B.PD.mlm[, sss] / DELTA.d.mlm[, sss]
      for (uuu in seq(ld.mlm))  # nnn.
        wz[, iam(posn.pobs.mix,
                 posn.pdip.mlm - 1 + uuu, M)] <-
          ned2l.dpobs.mix.pdip.mlm[, uuu]  # Link done later
    }  # all(c(la.mix, ld.mlm) > 0)





    if (all(c(li.mix, la.mlm) > 0)) {  # all(tmp3.TF[c(4, 8)])
      if (li.mlm)  # tmp3.TF[ 9]
        ned2l.dpobs.mlm.pstr.mix <-
        ned2l.dpobs.mlm.pstr.mix + rowSums(d0B.PI.mlm^2 / DELTA.i.mlm)
      if (ld.mix)  # tmp3.TF[ 6]
        ned2l.dpobs.mlm.pstr.mix <-  # nnn
        ned2l.dpobs.mlm.pstr.mix + rowSums(d0B.PD.mix^2 / DELTA.d.mix)
      if (ld.mlm)  # tmp3.TF[10]
        ned2l.dpobs.mlm.pstr.mix <-  # nnn
        ned2l.dpobs.mlm.pstr.mix + rowSums(d0B.PD.mlm^2 / DELTA.d.mlm)
        ned2l.dpobs.mlm.pstr.mix <-  # tmp3.TF[ 4] && li.mix
        ned2l.dpobs.mlm.pstr.mix -
        rowSums((d0A.i - d0B.PI.mix) * d0B.PI.mix / DELTA.i.mix)

      for (uuu in seq(la.mlm))  # ccc.
        wz[, iam(posn.pstr.mix,
                 posn.pobs.mlm - 1 + uuu, M)] <-
          ned2l.dpobs.mlm.pstr.mix  # Link done later
    }  # all(c(li.mix, la.mlm) > 0







    if (all(c(ld.mix, la.mlm) > 0)) {  # all(tmp3.TF[c(6, 8)])
      if (ld.mlm)  # tmp3.TF[10]
        ned2l.dpobs.mlm.pdip.mix <-
        ned2l.dpobs.mlm.pdip.mix - rowSums(d0B.PD.mlm^2 / DELTA.d.mlm)
      if (li.mix)  # tmp3.TF[ 4]
        ned2l.dpobs.mlm.pdip.mix <-
        ned2l.dpobs.mlm.pdip.mix - rowSums(d0B.PI.mix^2 / DELTA.i.mix)
      if (li.mlm)  # tmp3.TF[ 9]
        ned2l.dpobs.mlm.pdip.mix <-
        ned2l.dpobs.mlm.pdip.mix - rowSums(d0B.PI.mlm^2 / DELTA.i.mlm)
        ned2l.dpobs.mlm.pdip.mix <-  # all(tmp3.TF[c(6, 8)]) 
        ned2l.dpobs.mlm.pdip.mix +
        rowSums((d0A.d - d0B.PD.mix) * d0B.PD.mix / DELTA.d.mix)

      for (uuu in seq(la.mlm))  # nnn.
        wz[, iam(posn.pdip.mix,
                 posn.pobs.mlm - 1 + uuu, M)] <-
          ned2l.dpobs.mlm.pdip.mix  # Link done later
    }  # all(c(ld.mix, la.mlm) > 0






    if (all(c(li.mix, li.mlm) > 0)) {  # all(tmp3.TF[c(4, 9)])
      for (uuu in seq(li.mlm))  # tmp3.TF[ 9]
        for (sss in seq(li.mlm))
          ned2l.dpstr.mix.pstr.mlm[, uuu] <-
          ned2l.dpstr.mix.pstr.mlm[, uuu] -
            ((sss == uuu) - d0B.PI.mlm[, sss]) *
              d0B.PI.mlm[, sss] / DELTA.i.mlm[, sss]
      ned2l.dpstr.mix.pstr.mlm <-
      ned2l.dpstr.mix.pstr.mlm -
      rowSums((d0A.i - d0B.PI.mix) * d0B.PI.mix / DELTA.i.mix)
      if (ld.mix)  # tmp3.TF[ 6]
        ned2l.dpstr.mix.pstr.mlm <-  # nnn
        ned2l.dpstr.mix.pstr.mlm + rowSums(d0B.PD.mix^2 / DELTA.d.mix)
      if (ld.mlm)  # tmp3.TF[10]
        ned2l.dpstr.mix.pstr.mlm <-  # nnn
        ned2l.dpstr.mix.pstr.mlm + rowSums(d0B.PD.mlm^2 / DELTA.d.mlm)

      for (uuu in seq(li.mlm))  # Copy it. ccc.
        wz[, iam(posn.pstr.mix,
                 posn.pstr.mlm - 1 + uuu, M)] <-
          ned2l.dpstr.mix.pstr.mlm[, uuu]  # Link done later
    }  # all(c(li.mix, li.mlm) > 0



    if (all(c(ld.mix, ld.mlm) > 0)) {  # all(tmp3.TF[c(6, 10)])
      for (uuu in seq(ld.mlm))  # tmp3.TF[ 9]
        for (sss in seq(ld.mlm))
          ned2l.dpdip.mix.pdip.mlm[, uuu] <-
          ned2l.dpdip.mix.pdip.mlm[, uuu] -
            ((sss == uuu) - d0B.PD.mlm[, sss]) *
              d0B.PD.mlm[, sss] / DELTA.d.mlm[, sss]
      if (ld.mix)  # tmp3.TF[ 6]
        ned2l.dpdip.mix.pdip.mlm <-
        ned2l.dpdip.mix.pdip.mlm -
        rowSums((d0A.d - d0B.PD.mix) * d0B.PD.mix / DELTA.d.mix)
      if (li.mix)  # tmp3.TF[ 4]
        ned2l.dpdip.mix.pdip.mlm <-
        ned2l.dpdip.mix.pdip.mlm + rowSums(d0B.PI.mix^2 / DELTA.i.mix)
      if (li.mlm)  # tmp3.TF[ 9]
        ned2l.dpdip.mix.pdip.mlm <-
        ned2l.dpdip.mix.pdip.mlm + rowSums(d0B.PI.mlm^2 / DELTA.i.mlm)

      for (uuu in seq(ld.mlm))  # Copy it. ccc.
        wz[, iam(posn.pdip.mix,
                 posn.pdip.mlm - 1 + uuu, M)] <-
          ned2l.dpdip.mix.pdip.mlm[, uuu]  # Link done later
    }  # all(c(ld.mix, ld.mlm) > 0









    if (all(c(ld.mix, li.mlm) > 0)) {  # all(tmp3.TF[c(4, 9)])
      for (uuu in seq(li.mlm))  # tmp3.TF[ 9]
        for (sss in seq(li.mlm))
          ned2l.dpdip.mix.pstr.mlm[, uuu] <-
          ned2l.dpdip.mix.pstr.mlm[, uuu] +
            ((sss == uuu) - d0B.PI.mlm[, sss]) *
              d0B.PI.mlm[, sss] / DELTA.i.mlm[, sss]
      if (ld.mix)  # tmp3.TF[ 6]
        ned2l.dpdip.mix.pstr.mlm <-
        ned2l.dpdip.mix.pstr.mlm +
        rowSums((d0A.d - d0B.PD.mix) * d0B.PD.mix / DELTA.d.mix)
      if (li.mix)  # tmp3.TF[ 4]
        ned2l.dpdip.mix.pstr.mlm <-
        ned2l.dpdip.mix.pstr.mlm - rowSums(d0B.PI.mix^2 / DELTA.i.mix)
      if (ld.mlm)  # tmp3.TF[10]
        ned2l.dpdip.mix.pstr.mlm <-
        ned2l.dpdip.mix.pstr.mlm - rowSums(d0B.PD.mlm^2 / DELTA.d.mlm)

      for (uuu in seq(li.mlm))  # Copy it. ccc.
        wz[, iam(posn.pdip.mix,
                 posn.pstr.mlm - 1 + uuu, M)] <-
          ned2l.dpdip.mix.pstr.mlm[, uuu]  # Link done later
    }  # all(c(ld.mix, li.mlm) > 0



    if (all(c(li.mix, ld.mlm) > 0)) {  # all(tmp3.TF[c(4, 10)])
      for (uuu in seq(ld.mlm))  # tmp3.TF[10]
        for (sss in seq(ld.mlm))
          ned2l.dpstr.mix.pdip.mlm[, uuu] <-
          ned2l.dpstr.mix.pdip.mlm[, uuu] +
            ((sss == uuu) - d0B.PD.mlm[, sss]) *
              d0B.PD.mlm[, sss] / DELTA.d.mlm[, sss]
      if (li.mix)  # tmp3.TF[ 4]
        ned2l.dpstr.mix.pdip.mlm <-
        ned2l.dpstr.mix.pdip.mlm +
        rowSums((d0A.i - d0B.PI.mix) * d0B.PI.mix / DELTA.i.mix)
      if (ld.mix)  # tmp3.TF[ 6]
        ned2l.dpstr.mix.pdip.mlm <-
        ned2l.dpstr.mix.pdip.mlm - rowSums(d0B.PD.mix^2 / DELTA.d.mix)
      if (li.mlm)  # tmp3.TF[ 9]
        ned2l.dpstr.mix.pdip.mlm <-  # nnn.
        ned2l.dpstr.mix.pdip.mlm - rowSums(d0B.PI.mlm^2 / DELTA.i.mlm)

      for (uuu in seq(ld.mlm))  # Copy it. ccc.
        wz[, iam(posn.pstr.mix,
                 posn.pdip.mlm - 1 + uuu, M)] <-
          ned2l.dpstr.mix.pdip.mlm[, uuu]  # Link done later
    }  # all(c(li.mix, ld.mlm) > 0)




 
 
 
    if (lall.len) {
      wz.6 <- matrix(0, n, M * (M + 1) / 2)  # Or == 0 * wz
      ind.rc <- setdiff(1:M, ind.shape.z)  # Contiguous rows and
      lind.rc <- length(ind.rc)  # cols of the DAMLM



 # Copy in the thetas values: the looping is overkill.
      for (uuu in ind.shape.z)
        for (sss in seq(M))
          wz.6[, iam(uuu, sss, M)] <- wz[, iam(uuu, sss, M)]


 


 



      speed.up <- intercept.only && (
                  length(offset) == 1 || all(offset[1] == offset))



      IND.mlm <- iam(NA, NA, lind.rc, both = TRUE, diag = TRUE)
      n.use <- if (speed.up) 2 else n  # For sandwich.mlm





      if (!length(extra$ind.wz.match)) {
        Imat <- matrix(NA, lind.rc, lind.rc)
        for (jay in seq(lind.rc)) {
          iptr <- jay
          for (kay in (ind.rc[jay]):M) {
            if (!any(kay %in% ind.shape.z)) {
              Imat[jay, iptr] <-
                which(extra$index.M$row == ind.rc[jay] &
                      extra$index.M$col == kay)
              iptr <- iptr + 1
            }  # if
          }  # kay
        }  # jay
        ind.wz.match <- Imat[cbind(IND.mlm$row.ind,
                                   IND.mlm$col.ind)]
        extra$ind.wz.match <- ind.wz.match  # Assign it once
      }  # !length(extra$ind.wz.match)


      filling <- if (speed.up)
        wz[1:n.use, extra$ind.wz.match, drop = FALSE] else
        wz[, extra$ind.wz.match, drop = FALSE]









      M.mlm <- lind.rc
      if (is.null(extra$iamlist)) {
        extra$iamlist <- iamlist <- 
          iam(NA, NA, M = M.mlm, both = TRUE)
        if (M.mlm > 1) {  # Offdiagonal elts
          extra$iamlist.nod <- iamlist.nod <-
            iam(NA, NA, M.mlm, both = TRUE, diag = FALSE)
        }
      }  # is.null(extra$iamlist)
      iamlist <- extra$iamlist
      iamlist.nod <- extra$iamlist.nod
      MM12.mlm <- M.mlm * (M.mlm + 1) / 2

      Qf3 <- rowSums(filling[, 1:M.mlm, drop = FALSE] *  # Diag elts
                     (allprobs[1:n.use, 1:M.mlm, drop = FALSE])^2)
      if (M.mlm > 1)  # Offdiagonal elts
        Qf3 <- Qf3 + 2 * rowSums(allprobs[1:n.use, iamlist.nod$row] *
                     filling[, -(1:M.mlm), drop = FALSE] *  # n-vector
                                 allprobs[1:n.use, iamlist.nod$col])
      Qf3 <- matrix(Qf3, n.use, MM12.mlm)



      Qf2rowsums <- matrix(0, n.use, M.mlm)  # rowsums stored columnwise
      for (want in seq(M.mlm)) {  # Want the equivalent of rowSums(Qf2a)
        iamvec <- iam(want, 1:M.mlm, M = M.mlm)  # Diagonals included
        Qf2rowsums[, want] <- rowSums(filling[, iamvec, drop = FALSE] *
                                      allprobs[1:n.use, 1:M.mlm])
      }  # want
      Qf2a <- Qf2rowsums[, iamlist$row]
      Qf2b <- Qf2rowsums[, iamlist$col]


      Qform <- filling - Qf2a - Qf2b + Qf3  # n x MM12.mlm
      Qform <- Qform *
               allprobs[1:n.use, iamlist$row, drop = FALSE] *
               allprobs[1:n.use, iamlist$col, drop = FALSE]



      wz.6[, extra$ind.wz.match] <- if (speed.up)
        matrix(Qform[1, ], n, ncol(Qform), byrow = TRUE) else c(Qform)











 


      dstar.deta <- cbind(dshape.p.deta,
                          if (tmp3.TF[ 3]) dshape.a.deta else NULL,
                          if (tmp3.TF[ 5]) dshape.i.deta else NULL,
                          if (tmp3.TF[ 7]) dshape.d.deta else NULL)
      iptr <- 0
      if (length(ind.shape.z))
      for (uuu in ind.shape.z) {  # Could delete 3 for shape.a (orthog)
        iptr <- iptr + 1
        for (ttt in seq(lind.rc)) {
          wz.6[, iam(uuu, ind.rc[ttt], M)] <- 0  # Initialize
          for (sss in seq(lind.rc)) {
            wz.6[, iam(uuu, ind.rc[ttt], M)] <-
            wz.6[, iam(uuu, ind.rc[ttt], M)] +
              allprobs[, sss] * (max(0, sss == ttt) - allprobs[, ttt]) *
              wz[, iam(uuu, ind.rc[sss], M)] * dstar.deta[, iptr]
          }  # sss
        }  # ttt
      }  # uuu

      wz <- wz.6  # Completed
    }  # lall.len



    if (lall.len) {  # A MLM was fitted
      mytiny <- (allprobs <       sqrt(.Machine$double.eps)) |
                (allprobs > 1.0 - sqrt(.Machine$double.eps))
      atiny <- rowSums(mytiny) > 0
      if (any(atiny)) {
        ind.diags <- setdiff(1:M, ind.shape.z)  # Exclude thetas
        wz[atiny, ind.diags] <- .Machine$double.eps +
        wz[atiny, ind.diags] * (1 + .Machine$double.eps^0.5)
      }
    }  # lall.len




    c(w) * wz
  }), list( .truncate = truncate ))))
}  # gaitdlog


















































 moments.gaitdcombo.binom <-
  function(size.p, prob.p,
           a.mix = NULL, a.mlm = NULL,
           i.mix = NULL, i.mlm = NULL,
           d.mix = NULL, d.mlm = NULL,
           truncate = NULL, max.support = Inf,
           pobs.mix = 0, pobs.mlm = 0,  # Vector and matrix resp.
           pstr.mix = 0, pstr.mlm = 0,  # Ditto
           pdip.mix = 0, pdip.mlm = 0,  # Ditto
           byrow.aid = FALSE,  # For pobs.mlm and pstr.mlm
           size.a = size.p,
           size.i = size.p,
           size.d = size.p,
           prob.a = prob.p,
           prob.i = prob.p,
           prob.d = prob.p,
           type.fitted = "All",  # or "mean"
           moments2 = FALSE) {  # Use this for variances.
  if (is.infinite(max.support)) {
    rmlife1 <- rmlife2 <- numeric(length(size.p))  # 0
  } else {
    stop("currently RML unknown for finite 'max.support'")
    x.use <- max.support + 1
    rmlife1 <- NA
    rmlife2 <- NA
  }  # is.infinite(max.support)



  mylist1 <-
    moments.gaitdcombo.2par(
           theta1.p = size.p, theta2.p = prob.p,
           a.mix = a.mix, a.mlm = a.mlm,
           i.mix = i.mix, i.mlm = i.mlm,
           d.mix = d.mix, d.mlm = d.mlm,
           truncate = truncate, max.support = max.support,
           pobs.mix = pobs.mix, pobs.mlm = pobs.mlm,
           pstr.mix = pstr.mix, pstr.mlm = pstr.mlm,
           pdip.mix = pdip.mix, pdip.mlm = pdip.mlm,
           byrow.aid = byrow.aid,  # type.fitted = type.fitted,
           theta1.a = size.a, theta2.a = prob.a,
           theta1.i = size.i, theta2.i = prob.i,
           theta1.d = size.d, theta2.d = prob.d,
           moments2 = moments2,
           rmlife1 = rmlife1, rmlife2 = rmlife2,
           dfun = "dgaitdbinom")  # do.call() called.


  themean <- with(mylist1,
                  aprd1.mix + iprd1.mix + aprd1.mlm + iprd1.mlm -
                              dprd1.mix             - dprd1.mlm +
                  use.this * (munb.p - SumA1.mix.p -
                              SumA1.mlm.p - SumT1.p) / (
                  cdf.max.s - SumA0.mix.p - SumA0.mlm.p - SumT0.p))

  if (type.fitted == "mean") {
    return(themean)
  }

  ans <- c(mylist1,
           list('rmlife1'     = rmlife1,  # Has the right dimension
                'mean'        = themean))
  if (moments2) {  # Add more info
  ans <- c(ans,
           list('rmlife2'     = rmlife2))
  }
  ans
}  # moments.gaitdcombo.binom










 dgaitdbinom <-
  function(x, size.p, prob.p,
           a.mix = NULL,
           a.mlm = NULL,
           i.mix = NULL,
           i.mlm = NULL,
           d.mix = NULL,
           d.mlm = NULL,
           truncate = NULL,
           pobs.mix = 0,  # vector
           pobs.mlm = 0,  # matrix
           pstr.mix = 0,  # vector
           pstr.mlm = 0,  # matrix
           pdip.mix = 0,  # vector
           pdip.mlm = 0,  # matrix
           byrow.aid = FALSE,  # Applies to 'pobs.mlm' & 'pstr.mlm'
           size.a = size.p, size.i = size.p, size.d = size.p,
           prob.a = prob.p, prob.i = prob.p, prob.d = prob.p,
           log = FALSE,
           ...) {  # ... is for max.support (ignored)
  max.support <- Inf


  log.arg <- log;  rm(log)
  if (!length(max.support))  # Manually
    max.support <- max(size.p, size.a, size.i, na.rm = TRUE)
  lowsup <- 0  # Lower support
  gaitd.errorcheck(a.mix, a.mlm, i.mix, i.mlm,
                   d.mix, d.mlm, truncate, max.support)
  la.mix <- length(a.mix <- sort(a.mix))
  li.mix <- length(i.mix <- sort(i.mix))
  ld.mix <- length(d.mix <- sort(d.mix))
  la.mlm <- length(a.mlm)
  li.mlm <- length(i.mlm)
  ld.mlm <- length(d.mlm)
  ltrunc <- length(truncate)
  if (la.mix + la.mlm + li.mix + li.mlm + ld.mix + ld.mlm +
      ltrunc == 0 &&
      max.support >= max(size.p, na.rm = TRUE))
    return(dbinom(x, size.p, prob.p, log = log.arg))

  if (la.mix == 0) pobs.mix <- 0
  if (la.mlm == 0) pobs.mlm <- 0
  if (li.mix == 0) pstr.mix <- 0
  if (li.mlm == 0) pstr.mlm <- 0
  if (ld.mix == 0) pdip.mix <- 0
  if (ld.mlm == 0) pdip.mlm <- 0
 
  if (any(pobs.mix < 0 | 1 <= pobs.mix, na.rm = TRUE))
    stop("bad input for argument 'pobs.mix'")
  if (any(pobs.mlm < 0 | 1 <= pobs.mlm, na.rm = TRUE))
    stop("bad input for argument 'pobs.mlm'")
  if (any(pstr.mix < 0 | 1 <= pstr.mix, na.rm = TRUE))
    stop("bad input for argument 'pstr.mix'")
  if (any(pstr.mlm < 0 | 1 <= pstr.mlm, na.rm = TRUE))
    stop("bad input for argument 'pstr.mlm'")
  if (any(pdip.mix < 0 | 1 <= pdip.mix, na.rm = TRUE))
    stop("bad input for argument 'pdip.mix'")
  if (any(pdip.mlm < 0 | 1 <= pdip.mlm, na.rm = TRUE))
    stop("bad input for argument 'pdip.mlm'")

   
  LLL <- max(length(x),  
             length(pobs.mix), length(pstr.mix), length(pdip.mix),
             length(size.p),   length(size.a),
             length(size.i),   length(size.d),
             length(prob.p),   length(prob.a),
             length(prob.i),   length(prob.d))
  if (length(x)         < LLL) x        <- rep_len(x,         LLL)
  if (length(size.p)    < LLL) size.p   <- rep_len(size.p,    LLL)
  if (length(size.a)    < LLL) size.a   <- rep_len(size.a,    LLL)
  if (length(size.i)    < LLL) size.i   <- rep_len(size.i,    LLL)
  if (length(size.d)    < LLL) size.d   <- rep_len(size.d,    LLL)
  if (length(prob.p)    < LLL) prob.p   <- rep_len(prob.p,    LLL)
  if (length(prob.a)    < LLL) prob.a   <- rep_len(prob.a,    LLL)
  if (length(prob.i)    < LLL) prob.i   <- rep_len(prob.i,    LLL)
  if (length(prob.d)    < LLL) prob.d   <- rep_len(prob.d,    LLL)
  if (length(pobs.mix)  < LLL) pobs.mix <- rep_len(pobs.mix,  LLL)
  if (length(pstr.mix)  < LLL) pstr.mix <- rep_len(pstr.mix,  LLL)
  if (length(pdip.mix)  < LLL) pdip.mix <- rep_len(pdip.mix,  LLL)



  sumt <- 0  # Initialization to 0 important
  if (ltrunc)
    for (tval in truncate)
      sumt <- sumt + dbinom(tval, size.p, prob.p)
  vecTF.t <- is.finite(x) & ((x %in% truncate) | (max.support < x))
  cdf.max.s <- pbinom(max.support, size.p, prob.p)  # Usually 1
  denom.t <- cdf.max.s - sumt  # No sumt on RHS

    pmf0 <- ifelse(vecTF.t, 0, dbinom(x, size.p, prob.p) / denom.t)


  sum.a <- suma <- 0  # numeric(LLL)
  vecTF.a <- rep_len(FALSE, LLL)
  if (la.mlm) {
    pobs.mlm <-  matrix(pobs.mlm, LLL, la.mlm, byrow = byrow.aid)
    sum.a <-   .rowSums(pobs.mlm, LLL, la.mlm)
    if (any(1 < sum.a, na.rm = TRUE))
      stop("bad input for argument 'pobs.mlm'")  # zz

    for (aval in a.mlm)
      suma <- suma + dbinom(aval, size.p, prob.p)  # Part i

    for (jay in seq(la.mlm)) {
      aval <- a.mlm[jay]
      if (any(vecTF <- is.finite(x) & aval == x)) {
        pmf0[vecTF] <- pobs.mlm[vecTF, jay]
      }
      vecTF.a <- vecTF.a | vecTF  # Cumulative
    }  # jay
  }  # la.mlm



  pmf2.a <- pmf2.i <- pmf2.d <- 0
  if (la.mix) {
    allx.a <- lowsup:max(a.mix)
    pmf2.a <- dgaitdbinom(x,  # Outer distribution---mlm type
                          size.a, prob.a,
                          truncate = setdiff(allx.a, a.mix))
    for (aval in a.mix) {  # Part ii added; cumulative
      suma <- suma + dbinom(aval, size.p, prob.p)
      vecTF <- is.finite(x) & aval == x
      pmf0[vecTF] <- 0  # added; the true values are assigned below
      vecTF.a <- vecTF.a | vecTF  # Cumulative; added
    }
  }  # la.mix

  if (li.mix) {
    allx.i <- lowsup:max(i.mix)
    pmf2.i <- dgaitdbinom(x,  # Outer distribution---mlm type
                          size.i, prob.i,
                          truncate = setdiff(allx.i, i.mix))
  }





  sum.d <- 0  # numeric(LLL)
  if (ld.mlm) {
    pdip.mlm <-  matrix(pdip.mlm, LLL, ld.mlm, byrow = byrow.aid)
    sum.d <-   .rowSums(pdip.mlm, LLL, ld.mlm)
    if (any(1 < sum.d, na.rm = TRUE))
      stop("bad input for argument 'pdip.mlm'")
  }  # ld.mlm


  if (ld.mix) {
    allx.d <- lowsup:max(d.mix)
    pmf2.d <- dgaitdbinom(x, size.p = size.d, prob.p = prob.d,
                          truncate = setdiff(allx.d, d.mix))
  }  # ld.mix










  sum.i <- 0
  if (li.mlm) {
    pstr.mlm <- matrix(pstr.mlm, LLL, li.mlm, byrow = byrow.aid)
    sum.i <-  .rowSums(pstr.mlm, LLL, li.mlm)
    if (any(1 < sum.i, na.rm = TRUE))
      stop("bad input for argument 'pstr.mlm'")
  }  # li.mlm

  skip <- vecTF.t | vecTF.a  # Leave these values alone
  tmp6 <- 1 - sum.a - sum.i - pobs.mix - pstr.mix + sum.d + pdip.mix
  if (li.mlm + ld.mlm) {
    if (any(tmp6[!skip] < 0, na.rm = TRUE)) {
      warning("the vector of normalizing constants contains ",
              "some negative values. Replacing them with NAs")
      tmp6[!skip & tmp6 < 0] <- NA
    }
  }  # li.mlm + ld.mlm

  pmf0[!skip] <- (tmp6 *  # added
    dbinom(x, size.p, prob.p) / (cdf.max.s - suma - sumt))[!skip]


  if (li.mlm) {
    for (jay in seq(li.mlm)) {
      ival <- i.mlm[jay]
      if (any(vecTF <- is.finite(x) & ival == x)) {
        pmf0[vecTF] <- pmf0[vecTF] + pstr.mlm[vecTF, jay]
      }
    }  # jay
  }  # li.mlm




  if (ld.mlm) {
    for (jay in seq(ld.mlm)) {
      dval <- d.mlm[jay]
      if (any(vecTF <- is.finite(x) & dval == x)) {
          pmf0[vecTF] <- pmf0[vecTF] - pdip.mlm[vecTF, jay]
      }
    }  # jay
  }  # ld.mlm


  pmf0 <- pmf0 + pobs.mix * pmf2.a + pstr.mix * pmf2.i -
                 pdip.mix * pmf2.d



  if (log.arg) log(pmf0) else pmf0
}  # dgaitdbinom






 pgaitdbinom <-
  function(q, size.p, prob.p,
           a.mix = NULL,
           a.mlm = NULL,
           i.mix = NULL,
           i.mlm = NULL,
           d.mix = NULL,
           d.mlm = NULL,
           truncate = NULL,
           pobs.mix = 0,
           pobs.mlm = 0,
           pstr.mix = 0,
           pstr.mlm = 0,
           pdip.mix = 0,
           pdip.mlm = 0,
           byrow.aid = FALSE,
           size.a = size.p, size.i = size.p, size.d = size.p,
           prob.a = prob.p, prob.i = prob.p, prob.d = prob.p,
           lower.tail = TRUE,
           ...) {  # ... is for max.support (ignored)
  max.support <- Inf

  if (!length(max.support))  # Manually
    max.support <- max(size.p, size.a, size.i, na.rm = TRUE)
  lowsup <- 0
  gaitd.errorcheck(a.mix, a.mlm, i.mix, i.mlm,
                   d.mix, d.mlm, truncate, max.support)
  la.mix <- length(a.mix <- sort(a.mix))
  li.mix <- length(i.mix <- sort(i.mix))
  ld.mix <- length(d.mix <- sort(d.mix))
  la.mlm <- length(a.mlm)
  li.mlm <- length(i.mlm)
  ld.mlm <- length(d.mlm)
  ltrunc <- length(truncate)
  if (la.mix + la.mlm + li.mix + li.mlm + ld.mix + ld.mlm +
      ltrunc == 0 &&
      max.support >= max(size.p, na.rm = TRUE))
    return(pbinom(q, size.p, prob.p, lower.tail = lower.tail))


  if (la.mix == 0) pobs.mix <- 0
  if (la.mlm == 0) pobs.mlm <- 0
  if (li.mix == 0) pstr.mix <- 0
  if (li.mlm == 0) pstr.mlm <- 0
  if (ld.mix == 0) pdip.mix <- 0
  if (ld.mlm == 0) pdip.mlm <- 0

  if (any(pobs.mix < 0 | 1 <= pobs.mix, na.rm = TRUE))
    stop("bad input for argument 'pobs.mix'")
  if (any(pobs.mlm < 0 | 1 <= pobs.mlm, na.rm = TRUE))
    stop("bad input for argument 'pobs.mlm'")
  if (any(pstr.mix < 0 | 1 <= pstr.mix, na.rm = TRUE))
    stop("bad input for argument 'pstr.mix'")
  if (any(pstr.mlm < 0 | 1 <= pstr.mlm, na.rm = TRUE))
    stop("bad input for argument 'pstr.mlm'")
  if (any(pdip.mix < 0 | 1 <= pdip.mix, na.rm = TRUE))
    stop("bad input for argument 'pdip.mix'")
  if (any(pdip.mlm < 0 | 1 <= pdip.mlm, na.rm = TRUE))
    stop("bad input for argument 'pdip.mlm'")

  LLL <- max(length(q),
             length(pobs.mix), length(pstr.mix), length(pdip.mix),
             length(size.p),   length(size.a),   length(size.i),
             length(size.d),
             length(prob.p),   length(prob.a),   length(prob.i),
             length(prob.d))
  offset.a <- offset.i <- offset.d <-
  Offset.a <- Offset.i <- Offset.d <- numeric(LLL)
  if (length(q)        < LLL) q         <- rep_len(q,         LLL)
  if (length(size.p)   < LLL) size.p    <- rep_len(size.p,    LLL)
  if (length(size.a)   < LLL) size.a    <- rep_len(size.a,    LLL)
  if (length(size.i)   < LLL) size.i    <- rep_len(size.i,    LLL)
  if (length(size.d)   < LLL) size.d    <- rep_len(size.d,    LLL)
  if (length(prob.p)   < LLL) prob.p    <- rep_len(prob.p,    LLL)
  if (length(prob.a)   < LLL) prob.a    <- rep_len(prob.a,    LLL)
  if (length(prob.i)   < LLL) prob.i    <- rep_len(prob.i,    LLL)
  if (length(prob.d)   < LLL) prob.d    <- rep_len(prob.d,    LLL)
  if (length(pobs.mix) < LLL) pobs.mix  <- rep_len(pobs.mix,  LLL)
  if (length(pstr.mix) < LLL) pstr.mix  <- rep_len(pstr.mix,  LLL)
  if (length(pdip.mix) < LLL) pdip.mix  <- rep_len(pdip.mix,  LLL)



  sumt <- 0
  fudge.t <- numeric(LLL)
  cdf.max.s <- pbinom(max.support, size.p, prob.p)  # Usually 1
  if (ltrunc) {
    for (tval in truncate) {
      pmf.p <- dbinom(tval, size.p, prob.p)
      sumt <- sumt + pmf.p
      if (any(vecTF <- is.finite(q) & tval <= q))
        fudge.t[vecTF] <- fudge.t[vecTF] + pmf.p[vecTF]
    }
  }  # ltrunc

  sum.a <- suma <- 0  # numeric(LLL)
  fudge.a <- numeric(LLL)
  if (la.mlm) {
    pobs.mlm <- matrix(pobs.mlm, LLL, la.mlm, byrow = byrow.aid)
    sum.a <-  .rowSums(pobs.mlm, LLL, la.mlm)
    if (any(1 < sum.a, na.rm = TRUE))
      stop("bad input for argument 'pobs.mlm'")

    for (jay in seq(la.mlm)) {
      aval <- a.mlm[jay]
      pmf.p <- dbinom(aval, size.p, prob.p)
      suma <- suma + pmf.p  # cumulative; part i
      if (any(vecTF <- (is.finite(q) & aval <= q))) {
        offset.a[vecTF] <- offset.a[vecTF] + pobs.mlm[vecTF, jay]
        fudge.a[vecTF] <- fudge.a[vecTF] + pmf.p[vecTF]  # cumulative
      }
    }  # jay
  }  # la.mlm

  sum.i <- 0
  if (li.mlm) {
    pstr.mlm <- matrix(pstr.mlm, LLL, li.mlm, byrow = byrow.aid)
    sum.i <-  .rowSums(pstr.mlm, LLL, li.mlm)
    if (any(1 < sum.i, na.rm = TRUE))
      stop("bad input for argument 'pstr.mlm'")

    for (jay in seq(li.mlm)) {
      ival <- i.mlm[jay]
      if (any(vecTF <- (is.finite(q) & ival <= q))) {
        offset.i[vecTF] <- offset.i[vecTF] + pstr.mlm[vecTF, jay]
      }
    }  # jay
  }  # li.mlm



  use.pobs.mix <- 0
  if (la.mix) {
    use.pobs.mix <- matrix(0, LLL, la.mix)
    for (jay in seq(la.mix)) {
      aval <- a.mix[jay]
      pmf.a <- dbinom(aval, size.a, prob.a)
      pmf.p <- dbinom(aval, size.p, prob.p)
      use.pobs.mix[, jay] <- pmf.a
      suma <- suma + pmf.p  # cumulative; part ii
    }
    use.pobs.mix <- pobs.mix *
                    use.pobs.mix / rowSums(use.pobs.mix)

    for (jay in seq(la.mix)) {
      aval <- a.mix[jay]
      pmf.p <- dbinom(aval, size.p, prob.p)
      if (any(vecTF <- (is.finite(q) & aval <= q))) {
        Offset.a[vecTF] <- Offset.a[vecTF] + use.pobs.mix[vecTF, jay]
        fudge.a[vecTF] <- fudge.a[vecTF] + pmf.p[vecTF]  # cumulative
      }
    }  # jay
  }  # la.mix

  use.pstr.mix <- 0
  if (li.mix) {
    use.pstr.mix <- matrix(0, LLL, li.mix)
    for (jay in seq(li.mix)) {
      ival <- i.mix[jay]
      use.pstr.mix[, jay] <- dbinom(ival, size.i, prob.i)
    }
    use.pstr.mix <- pstr.mix *
                    use.pstr.mix / rowSums(use.pstr.mix)

    for (jay in seq(li.mix)) {
      ival <- i.mix[jay]
      pmf.p <- dbinom(ival, size.p, prob.p)
      if (any(vecTF <- (is.finite(q) & ival <= q))) {
        Offset.i[vecTF] <- Offset.i[vecTF] + use.pstr.mix[vecTF, jay]
      }
    }  # jay
  }  # li.mix



  sum.d <- 0
  if (ld.mlm) {
    pdip.mlm <- matrix(pdip.mlm, LLL, ld.mlm, byrow = byrow.aid)
    sum.d <-  .rowSums(pdip.mlm, LLL, ld.mlm)
    if (any(1 < sum.d, na.rm = TRUE))
      stop("bad input for argument 'pdip.mlm'")

    for (jay in seq(ld.mlm)) {
      dval <- d.mlm[jay]
      if (any(vecTF <- (is.finite(q) & dval <= q))) {
        offset.d[vecTF] <- offset.d[vecTF] + pdip.mlm[vecTF, jay]
      }
    }  # jay
  }  # ld.mlm



  use.pdip.mix <- 0
  if (ld.mix) {
    use.pdip.mix <- matrix(0, LLL, ld.mix)
    for (jay in seq(ld.mix)) {
      dval <- d.mix[jay]
      use.pdip.mix[, jay] <- dbinom(dval, size.d, prob.d)
    }
    use.pdip.mix <- pdip.mix *
                    use.pdip.mix / rowSums(use.pdip.mix)

    for (jay in seq(ld.mix)) {
      dval <- d.mix[jay]
      pmf.p <- dbinom(dval, size.p, prob.p)
      if (any(vecTF <- (is.finite(q) & dval <= q))) {
        Offset.d[vecTF] <- Offset.d[vecTF] + use.pdip.mix[vecTF, jay]
      }
    }  # jay
  }  # ld.mix



  numer1 <- 1 - sum.i - sum.a - pstr.mix - pobs.mix + sum.d + pdip.mix
  denom1 <- cdf.max.s - sumt - suma
  ans <- numer1 * (pbinom(q, size.p, prob.p) - fudge.t -
                   fudge.a) / denom1 +
         offset.a + offset.i - offset.d +
         Offset.a + Offset.i - Offset.d
  ans[max.support <= q] <- 1
  ans[ans < 0] <- 0  # Occasional roundoff error
  if (lower.tail) ans else 1 - ans
}  # pgaitdbinom






 qgaitdbinom <-
  function(p, size.p, prob.p,
           a.mix = NULL,
           a.mlm = NULL,
           i.mix = NULL,
           i.mlm = NULL,
           d.mix = NULL,
           d.mlm = NULL,
           truncate = NULL,
           pobs.mix = 0,
           pobs.mlm = 0,
           pstr.mix = 0,
           pstr.mlm = 0,
           pdip.mix = 0,
           pdip.mlm = 0,
           byrow.aid = FALSE,
           size.a = size.p, size.i = size.p, size.d = size.p,
           prob.a = prob.p, prob.i = prob.p, prob.d = prob.p,
           ...) {  # ... is for max.support (ignored)
  max.support <- NULL   # Different from Inf


  if (!length(max.support))  # Manually
    max.support <- max(size.p, size.a, size.i, na.rm = TRUE)
  lowsup <- 0
  gaitd.errorcheck(a.mix, a.mlm, i.mix, i.mlm,
                   d.mix, d.mlm, truncate, max.support)
  la.mix <- length(a.mix <- sort(a.mix))
  li.mix <- length(i.mix <- sort(i.mix))
  ld.mix <- length(d.mix <- sort(d.mix))
  la.mlm <- length(a.mlm)
  li.mlm <- length(i.mlm)
  ld.mlm <- length(d.mlm)
  ltrunc <- length(truncate)
  if (la.mix + la.mlm + li.mix + li.mlm + ld.mix + ld.mlm +
      ltrunc == 0 &&
      max.support >= max(size.p, na.rm = TRUE))
    return(qbinom(p, size.p, prob.p ))  # lower.tail, log.p = FALSE

  if (la.mix == 0) pobs.mix <- 0
  if (la.mlm == 0) pobs.mlm <- 0
  if (li.mix == 0) pstr.mix <- 0
  if (li.mlm == 0) pstr.mlm <- 0
  if (ld.mix == 0) pdip.mix <- 0
  if (ld.mlm == 0) pdip.mlm <- 0
 
  if (any(pobs.mix < 0 | 1 <= pobs.mix, na.rm = TRUE))
    stop("bad input for argument 'pobs.mix'")
  if (any(pobs.mlm < 0 | 1 <= pobs.mlm, na.rm = TRUE))
    stop("bad input for argument 'pobs.mlm'")
  if (any(pstr.mix < 0 | 1 <= pstr.mix, na.rm = TRUE))
    stop("bad input for argument 'pstr.mix'")
  if (any(pstr.mlm < 0 | 1 <= pstr.mlm, na.rm = TRUE))
    stop("bad input for argument 'pstr.mlm'")
  if (any(pdip.mix < 0 | 1 <= pdip.mix, na.rm = TRUE))
    stop("bad input for argument 'pdip.mix'")
  if (any(pdip.mlm < 0 | 1 <= pdip.mlm, na.rm = TRUE))
    stop("bad input for argument 'pdip.mlm'")



  LLL <- max(length(p),
             length(pobs.mix), length(pstr.mix), length(pdip.mix),
             length(size.p),   length(size.a),
             length(size.i),   length(size.d),
             length(prob.p),   length(prob.a),
             length(prob.i),   length(prob.d))
  if (length(p)        < LLL) p        <- rep_len(p,        LLL)
  if (length(size.p)   < LLL) size.p   <- rep_len(size.p,   LLL)
  if (length(size.a)   < LLL) size.a   <- rep_len(size.a,   LLL)
  if (length(size.i)   < LLL) size.i   <- rep_len(size.i,   LLL)
  if (length(size.d)   < LLL) size.d   <- rep_len(size.d,   LLL)
  if (length(prob.p)   < LLL) prob.p   <- rep_len(prob.p,   LLL)
  if (length(prob.a)   < LLL) prob.a   <- rep_len(prob.a,   LLL)
  if (length(prob.i)   < LLL) prob.i   <- rep_len(prob.i,   LLL)
  if (length(prob.d)   < LLL) prob.d   <- rep_len(prob.d,   LLL)
  if (length(pobs.mix) < LLL) pobs.mix <- rep_len(pobs.mix, LLL)
  if (length(pstr.mix) < LLL) pstr.mix <- rep_len(pstr.mix, LLL)
  if (length(pdip.mix) < LLL) pdip.mix <- rep_len(pdip.mix, LLL)

  pobs.mlm <- matrix(pobs.mlm, LLL, max(la.mlm, 1),
                     byrow = byrow.aid)
  pstr.mlm <- matrix(pstr.mlm, LLL, max(li.mlm, 1),
                     byrow = byrow.aid)
  pdip.mlm <- matrix(pdip.mlm, LLL, max(ld.mlm, 1),
                     byrow = byrow.aid)

  min.support <- lowsup  # Usual case; same as lowsup
  min.support.use <- if (ltrunc)
    min(setdiff(min.support:(ltrunc+5), truncate)) else
    min.support
  ans <- p + size.p + size.a + size.i + size.d +
             prob.p + prob.a + prob.i + prob.d

  bad0.p <- !is.finite(size.p) | size.p <= 0 |
            !is.finite(prob.p) | prob.p <= 0 | 1 <= prob.p
  bad0.a <- !is.finite(size.a) | size.a <= 0 |
            !is.finite(prob.a) | prob.a <= 0 | 1 <= prob.a
  bad0.i <- !is.finite(size.i) | size.i <= 0 |
            !is.finite(prob.i) | prob.i <= 0 | 1 <= prob.i
  bad0.d <- !is.finite(size.d) | size.d <= 0 |
            !is.finite(prob.d) | prob.d <= 0 | 1 <= prob.d
  bad0 <- bad0.p | bad0.a | bad0.i | bad0.d
  bad <- bad0 | !is.finite(p) | p <= 0 | 1 <= p

  Lo <- rep_len(min.support.use - 0.5, LLL)
  approx.ans <- Lo  # True at lhs
  Hi <- rep_len(max.support + 0.5, LLL)   # Need finite RHS endpoint
  dont.iterate <- bad
  done <- dont.iterate |
    p <= pgaitdbinom(Hi, size.p, prob.p,
                    a.mix = a.mix, a.mlm = a.mlm,
                    i.mix = i.mix, i.mlm = i.mlm,
                    d.mix = d.mix, d.mlm = d.mlm,
                    truncate = truncate,
                    pstr.mix = pstr.mix, pobs.mix = pobs.mix,
                    pdip.mix = pdip.mix,
                    pstr.mlm = pstr.mlm, pobs.mlm = pobs.mlm,
                    pdip.mlm = pdip.mlm,
                    size.a = size.a, size.i = size.i,
                    size.d = size.d,
                    prob.a = prob.a, prob.i = prob.i,
                    prob.d = prob.d,
                    byrow.aid = FALSE)

  iter <- 0
  max.iter <- round(log2(.Machine$double.xmax)) - 3
  while (!all(done) && iter < max.iter) {
    Lo[!done] <- Hi[!done]
    Hi[!done] <- 2 * Hi[!done] + 10.5  # Bug fixed
    Hi <- pmin(max.support + 0.5, Hi)  # 20190924
    done[!done] <-
      (p[!done] <= pgaitdbinom(Hi[!done],
                       size.p[!done], prob.p[!done],
                       a.mix = a.mix, a.mlm = a.mlm,
                       i.mix = i.mix, i.mlm = i.mlm,
                       d.mix = d.mix, d.mlm = d.mlm,
                       truncate = truncate,
                       pobs.mix = pobs.mix[!done],
                       pstr.mix = pstr.mix[!done],
                       pdip.mix = pdip.mix[!done],
                       pobs.mlm = pobs.mlm[!done, , drop = FALSE],
                       pstr.mlm = pstr.mlm[!done, , drop = FALSE],
                       pdip.mlm = pdip.mlm[!done, , drop = FALSE],
                       size.a = size.a[!done],
                       size.i = size.i[!done],
                       size.d = size.d[!done],
                       prob.a = prob.a[!done],
                       prob.i = prob.i[!done],
                       prob.d = prob.d[!done],
                       byrow.aid = FALSE))
    iter <- iter + 1
  }

      foo <- function(q, size.p, prob.p,
                      a.mix = NULL, a.mlm = NULL,
                      i.mix = NULL, i.mlm = NULL,
                      d.mix = NULL, d.mlm = NULL,
                      truncate = NULL,
                      pobs.mix = 0, pstr.mix = 0, pdip.mix = 0,
                      pobs.mlm = 0, pstr.mlm = 0, pdip.mlm = 0,
                      size.a = size.p, size.i = size.p,
                      size.d = size.p,
                      prob.a = prob.p, prob.i = prob.p,
                      prob.d = prob.p,
                      byrow.aid = FALSE, p)
      pgaitdbinom(q, size.p = size.p, prob.p = prob.p,
                       a.mix = a.mix, a.mlm = a.mlm,
                       i.mix = i.mix, i.mlm = i.mlm,
                       d.mix = d.mix, d.mlm = d.mlm,
                       truncate = truncate,
                       pobs.mix = pobs.mix,
                       pstr.mix = pstr.mix,
                       pdip.mix = pdip.mix,
                       pobs.mlm = pobs.mlm,
                       pstr.mlm = pstr.mlm,
                       pdip.mlm = pdip.mlm,
                       size.a = size.a, prob.a = prob.a,
                       size.i = size.i, prob.i = prob.i,
                       size.d = size.d, prob.d = prob.d,
                       byrow.aid = FALSE) - p

      lhs <- dont.iterate |
        p <= dgaitdbinom(min.support.use,
                        size.p = size.p, prob.p = prob.p,
                        a.mix = a.mix, a.mlm = a.mlm,
                        i.mix = i.mix, i.mlm = i.mlm,
                        d.mix = d.mix, d.mlm = d.mlm,
                        truncate = truncate,
                        pobs.mix = pobs.mix,
                        pstr.mix = pstr.mix,
                        pdip.mix = pdip.mix,
                        pobs.mlm = pobs.mlm,
                        pstr.mlm = pstr.mlm,
                        pdip.mlm = pdip.mlm,
                        size.a = size.a, prob.a = prob.a,
                        size.i = size.i, prob.i = prob.i,
                        size.d = size.d, prob.d = prob.d,
                        byrow.aid = FALSE)

  if (any(!lhs)) {
    approx.ans[!lhs] <-
      bisection.basic(foo, Lo[!lhs], Hi[!lhs], tol = 1/16,
                      size.p = size.p[!lhs],
                      prob.p = prob.p[!lhs],
                      a.mix = a.mix, a.mlm = a.mlm,
                      i.mix = i.mix, i.mlm = i.mlm,
                      d.mix = d.mix, d.mlm = d.mlm,
                      truncate = truncate,
                      pstr.mix = pstr.mix[!lhs],
                      pstr.mlm = pstr.mlm[!lhs, , drop = FALSE],
                      pobs.mix = pobs.mix[!lhs],
                      pobs.mlm = pobs.mlm[!lhs, , drop = FALSE],
                      pdip.mix = pdip.mix[!lhs],
                      pdip.mlm = pdip.mlm[!lhs, , drop = FALSE],
                      size.a = size.a[!lhs],
                      prob.a = prob.a[!lhs],
                      size.i = size.i[!lhs],
                      prob.i = prob.i[!lhs],
                      size.d = size.d[!lhs],
                      prob.d = prob.d[!lhs],
                      byrow.aid = FALSE,
                      p = p[!lhs])
    faa <- floor(approx.ans[!lhs])
    tmp <-
      ifelse(pgaitdbinom(faa,
                       size.p = size.p[!lhs],
                       prob.p = prob.p[!lhs],
                       a.mix = a.mix, a.mlm = a.mlm,
                       i.mix = i.mix, i.mlm = i.mlm,
                       d.mix = d.mix, d.mlm = d.mlm,
                       truncate = truncate,
                       pstr.mix = pstr.mix[!lhs],
                       pstr.mlm = pstr.mlm[!lhs, , drop = FALSE],
                       pobs.mix = pobs.mix[!lhs],
                       pobs.mlm = pobs.mlm[!lhs, , drop = FALSE],
                       pdip.mix = pdip.mix[!lhs],
                       pdip.mlm = pdip.mlm[!lhs, , drop = FALSE],
                       size.a = size.a[!lhs],
                       prob.a = prob.a[!lhs],
                       size.i = size.i[!lhs],
                       prob.i = prob.i[!lhs],
                       size.d = size.d[!lhs],
                       prob.d = prob.d[!lhs],
                       byrow.aid = FALSE) < p[!lhs] &
             p[!lhs] <= pgaitdbinom(faa + 1,
                       size.p = size.p[!lhs],
                       prob.p = prob.p[!lhs],
                       a.mix = a.mix, a.mlm = a.mlm,
                       i.mix = i.mix, i.mlm = i.mlm,
                       d.mix = d.mix, d.mlm = d.mlm,
                       truncate = truncate,
                       pstr.mix = pstr.mix[!lhs],
                       pstr.mlm = pstr.mlm[!lhs, , drop = FALSE],
                       pobs.mix = pobs.mix[!lhs],
                       pobs.mlm = pobs.mlm[!lhs, , drop = FALSE],
                       pdip.mix = pdip.mix[!lhs],
                       pdip.mlm = pdip.mlm[!lhs, , drop = FALSE],
                       size.a = size.a[!lhs],
                       prob.a = prob.a[!lhs],
                       size.i = size.i[!lhs],
                       prob.i = prob.i[!lhs],
                       size.d = size.d[!lhs],
                       prob.d = prob.d[!lhs],
                       byrow.aid = FALSE),
             faa + 1, faa)
    ans[!lhs] <- tmp
  }  # any(!lhs)


  if (ltrunc)
    while (any(vecTF <- !bad & ans %in% truncate))
      ans[vecTF] <- 1 + ans[vecTF]



  vecTF <- !bad0 & !is.na(p) &
           p <= dgaitdbinom(min.support.use, size.p, prob.p,
                       a.mix = a.mix, a.mlm = a.mlm,
                       i.mix = i.mix, i.mlm = i.mlm,
                       d.mix = d.mix, d.mlm = d.mlm,
                       truncate = truncate,
                       pobs.mix = pobs.mix,
                       pstr.mix = pstr.mix,
                       pdip.mix = pdip.mix,
                       pobs.mlm = pobs.mlm,
                       pstr.mlm = pstr.mlm,
                       pdip.mlm = pdip.mlm,
                       size.a = size.a, size.i = size.i,
                       size.d = size.d,
                       prob.a = prob.a, prob.i = prob.i,
                       prob.d = prob.d,
                       byrow.aid = FALSE)
  ans[vecTF] <- min.support.use

  ans[!bad0 & !is.na(p) & p == 0] <- min.support.use
  ans[!bad0 & !is.na(p) & p == 1] <- max.support  # Inf
  ans[!bad0 & !is.na(p) & p <  0] <- NaN
  ans[!bad0 & !is.na(p) & p >  1] <- NaN
  ans[ bad0] <- NaN
  ans
}  # qgaitdbinom






 rgaitdbinom <-
  function(n, size.p, prob.p,
           a.mix = NULL,
           a.mlm = NULL,
           i.mix = NULL,
           i.mlm = NULL,
           d.mix = NULL,
           d.mlm = NULL,
           truncate = NULL,
           pobs.mix = 0,  # vector
           pobs.mlm = 0,  # matrix
           pstr.mix = 0,  # vector
           pstr.mlm = 0,  # matrix
           pdip.mix = 0,  # vector
           pdip.mlm = 0,  # matrix
           byrow.aid = FALSE,
           size.a = size.p, size.i = size.p, size.d = size.p,
           prob.a = prob.p, prob.i = prob.p, prob.d = prob.p,
           ...) {  # ... is for max.support (ignored)
    qgaitdbinom(runif(n), size.p, prob.p,
              a.mix = a.mix,
              a.mlm = a.mlm,
              i.mix = i.mix,
              i.mlm = i.mlm,
              d.mix = d.mix,
              d.mlm = d.mlm,
              truncate = truncate,
              pobs.mix = pobs.mix,
              pobs.mlm = pobs.mlm,
              pstr.mix = pstr.mix,
              pstr.mlm = pstr.mlm,
              pdip.mix = pdip.mix,
              pdip.mlm = pdip.mlm,
              size.a = size.a, size.i = size.i, size.d = size.d,
              prob.a = prob.a, prob.i = prob.i, prob.d = prob.d,
              byrow.aid = byrow.aid)
}  # rgaitdbinom









 gaitd.errorcheck <-
  function(a.mix = NULL, a.mlm = NULL,
           i.mix = NULL, i.mlm = NULL,
           d.mix = NULL, d.mlm = NULL,
           truncate = NULL,
           max.support = Inf,
           min.support = 0, nparams = 1) {
  la.mix <- length(a.mix <- sort(a.mix))
  li.mix <- length(i.mix <- sort(i.mix))
  ld.mix <- length(d.mix <- sort(d.mix))
  la.mlm <- length(a.mlm)
  li.mlm <- length(i.mlm)
  ld.mlm <- length(d.mlm)
  ltrunc <- length(truncate)

  if (!is.numeric(max.support) || is.na(max.support) ||
      length(max.support) != 1 || max.support < min.support ||
      round(max.support) != max.support ||
      (length(truncate) && (
          min(truncate, na.rm = TRUE) < min.support ||
          max.support <= max(truncate, na.rm = TRUE))))
    stop("bad input for argument 'max.support' and/or ",
         "'truncate'")

  allargs <- c(a.mix, a.mlm, i.mix, i.mlm, d.mix, d.mlm)
  allargs <- c(allargs, truncate)  # No NA, NaN, -Inf or Inf allowed
  if (la.mix + la.mlm + li.mix + li.mlm + ld.mix + ld.mlm)
    if (!is.Numeric(allargs, integer.valued = TRUE) ||
        any(allargs < min.support) ||
        any(max.support < allargs))
      stop("bad input for arguments 'a.mix', 'a.mlm', ",
           "'i.mix', 'i.mlm', 'd.mix' and/or 'd.mlm'")
  if (length(unique(allargs)) < la.mix + la.mlm + li.mix + li.mlm +
                                ld.mix + ld.mlm + ltrunc)
      stop("duplicate values found in arguments 'a.mix', ",
           "'a.mlm', 'i.mix', 'i.mlm', 'd.mix', 'd.mlm'",
           " and 'truncate'")

  if (nparams == 2) {
    if(la.mix == 2)
      stop("overfitting: trying to fit a ", nparams, "-parameter ",
           "distribution based on length(a.mix) == ", la.mix, " points")
    if (li.mix == 2)
      stop("overfitting: trying to fit a ", nparams, "-parameter ",
           "distribution based on length(i.mix) == ", li.mix, " points")
    if (ld.mix == 2)
      stop("overfitting: trying to fit a ", nparams, "-parameter ",
           "distribution based on length(d.mix) == ", ld.mix, " points")
  }


}  # gaitd.errorcheck










 moments.gaitdcombo.1par <-
  function(theta.p,
           a.mix = NULL, a.mlm = NULL,
           i.mix = NULL, i.mlm = NULL,
           d.mix = NULL, d.mlm = NULL,
           truncate = NULL, max.support = Inf,
           pobs.mix = 0,  # Vector
           pobs.mlm = 0,  # Matrix
           pstr.mix = 0,
           pstr.mlm = 0,  # Ditto
           pdip.mix = 0,
           pdip.mlm = 0,  # Ditto
           byrow.aid = FALSE,  # For pobs.mlm and pstr.mlm
           theta.a = theta.p, theta.i = theta.p, theta.d = theta.p,
           moments2 = FALSE,  # Use this for variances.
           rmlife1 = 0, rmlife2 = 0,
           dfun = "dpois") {


  NOS <- 1
  nnn <- length(theta.p)
  pfun <- dfun
  substring(pfun, 1) <- "p"  # Replace the "d" by a "p"
  cdf.max.s <- do.call(pfun, list(max.support, theta.p))
  LALT.MIX <- length(a.mix)
  LALT.MLM <- length(a.mlm)
  LINF.MIX <- length(i.mix)
  LINF.MLM <- length(i.mlm)
  LDEF.MIX <- length(d.mix)
  LDEF.MLM <- length(d.mlm)
  LTRUNCAT <- length(truncate)

  if (LALT.MLM == 0) {
    if (!all(pobs.mlm == 0))
      warning("ignoring argument 'pobs.mlm'")
    pobs.mlm <- 0
  }
  if (LINF.MLM == 0) {
    if (!all(pstr.mlm == 0))
      warning("ignoring argument 'pstr.mlm'")
    pstr.mlm <- 0
  }
  if (LDEF.MLM == 0) {
    if (!all(pdip.mlm == 0))
      warning("ignoring argument 'pdip.mlm'")
    pdip.mlm <- 0
  }
  if (LALT.MIX == 0) {
    if (!all(pobs.mix == 0))
      warning("ignoring argument 'pobs.mix'")
    pobs.mix <- 0
  }
  if (LINF.MIX == 0) {
    if (!all(pstr.mix == 0))
      warning("ignoring argument 'pstr.mix'")
    pstr.mix <- 0
  }
  if (LDEF.MIX == 0) {
    if (!all(pdip.mix == 0))
      warning("ignoring argument 'pdip.mix'")
    pdip.mix <- 0
  }

  SumT0.p <- matrix(0, nnn, NOS)  # Does not include upper RHS tail
  SumT1.p <- matrix(rmlife1, nnn, NOS)  # Includes RHS tail
  SumT2.p <- matrix(rmlife2, nnn, NOS)  # Includes RHS tail
  if (LTRUNCAT)
    for (tval in truncate) {
      pmf.p <- do.call(dfun, list(tval, theta.p))
      SumT0.p <- SumT0.p + pmf.p  # Need tval<=max.support
      SumT1.p <- SumT1.p + pmf.p * tval
      if (moments2)
        SumT2.p <- SumT2.p + pmf.p * tval^2
    }

  use.pobs.mix <- use.pobs.mlm <-  # So rowSums() works below.
  use.pstr.mix <- use.pstr.mlm <-
  use.pdip.mix <- use.pdip.mlm <- matrix(0, nnn, 1)
  aprd1.mix <- aprd1.mlm <-    # aprd1.m?? is an innerprod
  aprd2.mix <- aprd2.mlm <- 0  # aprd2.m?? is an innerprod
  SumA0.mix.p <- SumA0.mlm.p <-
  SumA0.mix.a <- SumA0.mlm.a <-
  SumA1.mix.p <- SumA1.mlm.p <-
  SumA1.mix.a <- SumA1.mlm.a <-
                 SumA1.mlm.x <-
  SumA2.mix.p <- SumA2.mlm.p <-
  SumA2.mix.a <- SumA2.mlm.a <-
                 SumA2.mlm.x <- matrix(0, nnn, NOS)
  if (LALT.MIX)
    use.pobs.mix <- matrix(pobs.mix, nnn, 1)
  if (LINF.MIX)
    use.pstr.mix <- matrix(pstr.mix, nnn, 1)
  if (LDEF.MIX)
    use.pdip.mix <- matrix(pdip.mix, nnn, 1)
  if (LALT.MLM)
    use.pobs.mlm <- matrix(pobs.mlm, nnn, LALT.MLM, byrow = byrow.aid)
  if (LINF.MLM)
    use.pstr.mlm <- matrix(pstr.mlm, nnn, LINF.MLM, byrow = byrow.aid)
  if (LDEF.MLM)
    use.pdip.mlm <- matrix(pdip.mlm, nnn, LDEF.MLM, byrow = byrow.aid)


  if (LALT.MIX) {
    for (jay in seq_len(LALT.MIX)) {
      aval <- a.mix[jay]
      pmf.p <- do.call(dfun, list(aval, theta.p))
      pmf.a <- do.call(dfun, list(aval, theta.a))
      SumA0.mix.p <- SumA0.mix.p + pmf.p
      SumA0.mix.a <- SumA0.mix.a + pmf.a
      SumA1.mix.p <- SumA1.mix.p + pmf.p * aval
      SumA1.mix.a <- SumA1.mix.a + pmf.a * aval
      if (moments2) {
        SumA2.mix.p <- SumA2.mix.p + pmf.p * aval^2
        SumA2.mix.a <- SumA2.mix.a + pmf.a * aval^2
      }
    }  # for jay
    aprd1.mix <- use.pobs.mix * SumA1.mix.a / SumA0.mix.a
    if (moments2)
      aprd2.mix <- use.pobs.mix * SumA2.mix.a / SumA0.mix.a
  }  # LALT.MIX


  if (LALT.MLM) {
    for (jay in seq_len(LALT.MLM)) {
      aval <- a.mlm[jay]
      pmf.x <- use.pobs.mlm[, jay]
      pmf.p <- do.call(dfun, list(aval, theta.p))
      pmf.a <- do.call(dfun, list(aval, theta.a))
      SumA0.mlm.p <- SumA0.mlm.p + pmf.p
      SumA0.mlm.a <- SumA0.mlm.a + pmf.a
      SumA1.mlm.p <- SumA1.mlm.p + pmf.p * aval
      SumA1.mlm.a <- SumA1.mlm.a + pmf.a * aval
      SumA1.mlm.x <- SumA1.mlm.x + pmf.x * aval
      if (moments2) {
        SumA2.mlm.p <- SumA2.mlm.p + pmf.p * aval^2
        SumA2.mlm.a <- SumA2.mlm.a + pmf.a * aval^2
        SumA2.mlm.x <- SumA2.mlm.x + pmf.x * aval^2
      }
    }  # for jay
    aprd1.mlm <- SumA1.mlm.x
    if (moments2)
      aprd2.mlm <- SumA2.mlm.x
  }  # LALT.MLM


  iprd1.mix <- iprd1.mlm <-    # iprd1.m?? is an innerprod
  iprd2.mix <- iprd2.mlm <- 0  # iprd2.m?? is an innerprod
  SumI0.mix.p <- SumI0.mlm.p <-
  SumI0.mix.i <- SumI0.mlm.i <-
  SumI1.mix.p <- SumI1.mlm.p <-
  SumI1.mix.i <- SumI1.mlm.i <-
                 SumI1.mlm.x <-
  SumI2.mix.p <- SumI2.mlm.p <-
  SumI2.mix.i <- SumI2.mlm.i <-
                 SumI2.mlm.x <- matrix(0, nnn, NOS)


  dprd1.mix <- dprd1.mlm <-
  dprd2.mix <- dprd2.mlm <- 0
  SumD0.mix.p <- SumD0.mlm.p <-
  SumD0.mix.d <- SumD0.mlm.d <-
  SumD1.mix.p <- SumD1.mlm.p <-
  SumD1.mix.d <- SumD1.mlm.d <-
                 SumD1.mlm.x <-
  SumD2.mix.p <- SumD2.mlm.p <-
  SumD2.mix.d <- SumD2.mlm.d <-
                 SumD2.mlm.x <- matrix(0, nnn, NOS)



  if (LINF.MIX) {
    for (jay in seq_len(LINF.MIX)) {
      ival <- i.mix[jay]
      pmf.p <- do.call(dfun, list(ival, theta.p))
      pmf.i <- do.call(dfun, list(ival, theta.i))
      SumI0.mix.p <- SumI0.mix.p + pmf.p
      SumI0.mix.i <- SumI0.mix.i + pmf.i
      SumI1.mix.p <- SumI1.mix.p + pmf.p * ival
      SumI1.mix.i <- SumI1.mix.i + pmf.i * ival
      if (moments2) {
        SumI2.mix.p <- SumI2.mix.p + pmf.p * ival^2
        SumI2.mix.i <- SumI2.mix.i + pmf.i * ival^2
      }
    }  # for jay
    iprd1.mix <- use.pstr.mix * SumI1.mix.i / SumI0.mix.i
    if (moments2)
      iprd2.mix <- use.pstr.mix * SumI2.mix.i / SumI0.mix.i
  }  # LINF.MIX


  if (LINF.MLM) {
    for (jay in seq_len(LINF.MLM)) {
      ival <- i.mlm[jay]
      pmf.x <- use.pstr.mlm[, jay]
      pmf.p <- do.call(dfun, list(ival, theta.p))
      pmf.i <- do.call(dfun, list(ival, theta.i))
      SumI0.mlm.p <- SumI0.mlm.p + pmf.p
      SumI0.mlm.i <- SumI0.mlm.i + pmf.i
      SumI1.mlm.p <- SumI1.mlm.p + pmf.p * ival
      SumI1.mlm.i <- SumI1.mlm.i + pmf.i * ival
      SumI1.mlm.x <- SumI1.mlm.x + pmf.x * ival
      if (moments2) {
        SumI2.mlm.p <- SumI2.mlm.p + pmf.p * ival^2
        SumI2.mlm.i <- SumI2.mlm.i + pmf.i * ival^2
        SumI2.mlm.x <- SumI2.mlm.x + pmf.x * ival^2
      }
    }  # for jay
    iprd1.mlm <- SumI1.mlm.x
    if (moments2)
      iprd2.mlm <- SumI2.mlm.x
  }  # LINF.MLM






  if (LDEF.MIX) {
    for (jay in seq_len(LDEF.MIX)) {
      dval <- d.mix[jay]
      pmf.p <- do.call(dfun, list(dval, theta.p))
      pmf.d <- do.call(dfun, list(dval, theta.d))
      SumD0.mix.p <- SumD0.mix.p + pmf.p
      SumD0.mix.d <- SumD0.mix.d + pmf.d
      SumD1.mix.p <- SumD1.mix.p + pmf.p * dval
      SumD1.mix.d <- SumD1.mix.d + pmf.d * dval
      if (moments2) {
        SumD2.mix.p <- SumD2.mix.p + pmf.p * dval^2
        SumD2.mix.d <- SumD2.mix.d + pmf.d * dval^2
      }
    }  # for jay
    dprd1.mix <- use.pdip.mix * SumD1.mix.d / SumD0.mix.d
    if (moments2)
      dprd2.mix <- use.pdip.mix * SumD2.mix.d / SumD0.mix.d
  }  # LDEF.MIX






  if (LDEF.MLM) {
    for (jay in seq_len(LDEF.MLM)) {
      dval <- d.mlm[jay]
      pmf.x <- use.pdip.mlm[, jay]
      pmf.p <- do.call(dfun, list(dval, theta.p))
      pmf.d <- do.call(dfun, list(dval, theta.d))
      SumD0.mlm.p <- SumD0.mlm.p + pmf.p
      SumD0.mlm.d <- SumD0.mlm.d + pmf.d
      SumD1.mlm.p <- SumD1.mlm.p + pmf.p * dval
      SumD1.mlm.d <- SumD1.mlm.d + pmf.d * dval
      SumD1.mlm.x <- SumD1.mlm.x + pmf.x * dval
      if (moments2) {
        SumD2.mlm.p <- SumD2.mlm.p + pmf.p * dval^2
        SumD2.mlm.d <- SumD2.mlm.d + pmf.d * dval^2
        SumD2.mlm.x <- SumD2.mlm.x + pmf.x * dval^2
      }
    }  # for jay
    dprd1.mlm <- SumD1.mlm.x
    if (moments2)
      dprd2.mlm <- SumD2.mlm.x
  }  # LDEF.MLM




  use.this <- 1 - rowSums(use.pobs.mlm) - rowSums(use.pstr.mlm) +
                  rowSums(use.pdip.mlm) -
              use.pobs.mix - use.pstr.mix + use.pdip.mix
  ans <- list('cdf.max.s'   = cdf.max.s,
              'SumT0.p'     = SumT0.p,
              'SumT1.p'     = SumT1.p,
              'SumA0.mix.a' = SumA0.mix.a,
              'SumA0.mix.p' = SumA0.mix.p,
              'SumA1.mix.a' = SumA1.mix.a,
              'SumA1.mix.p' = SumA1.mix.p,
              'SumA0.mlm.a' = SumA0.mlm.a,
              'SumA0.mlm.p' = SumA0.mlm.p,
              'SumA1.mlm.a' = SumA1.mlm.a,
              'SumA1.mlm.p' = SumA1.mlm.p,
              'SumI0.mix.i' = SumI0.mix.i,
              'SumI0.mix.p' = SumI0.mix.p,
              'SumI1.mix.i' = SumI1.mix.i,
              'SumI1.mix.p' = SumI1.mix.p,
              'SumI0.mlm.i' = SumI0.mlm.i,
              'SumI0.mlm.p' = SumI0.mlm.p,
              'SumI1.mlm.i' = SumI1.mlm.i,
              'SumI1.mlm.p' = SumI1.mlm.p,
              'SumD0.mix.d' = SumD0.mix.d,
              'SumD0.mix.p' = SumD0.mix.p,
              'SumD1.mix.d' = SumD1.mix.d,
              'SumD1.mix.p' = SumD1.mix.p,
              'SumD0.mlm.d' = SumD0.mlm.d,
              'SumD0.mlm.p' = SumD0.mlm.p,
              'SumD1.mlm.d' = SumD1.mlm.d,
              'SumD1.mlm.p' = SumD1.mlm.p,
              'aprd1.mix'   = aprd1.mix,
              'aprd1.mlm'   = aprd1.mlm,
              'iprd1.mix'   = iprd1.mix,
              'iprd1.mlm'   = iprd1.mlm,
              'dprd1.mix'   = dprd1.mix,    #
              'dprd1.mlm'   = dprd1.mlm,    #
              'use.this'    = use.this)

  if (moments2) {  # Add more info
  ans <- c(ans,
         list(  # 'rmlife2'     = rmlife2,   # May be scalar
              'aprd2.mix'   = aprd2.mix,
              'aprd2.mlm'   = aprd2.mlm,
              'iprd2.mix'   = iprd2.mix,
              'iprd2.mlm'   = iprd2.mlm,
              'dprd2.mix'   = dprd2.mix,    #
              'dprd2.mlm'   = dprd2.mlm,    #
              'SumT2.p'     = SumT2.p,
              'SumA2.mix.p' = SumA2.mix.p,
              'SumA2.mix.a' = SumA2.mix.a,
              'SumI2.mix.p' = SumI2.mix.p,
              'SumI2.mix.i' = SumI2.mix.i,
              'SumD2.mix.p' = SumD2.mix.p,    #
              'SumD2.mix.d' = SumD2.mix.d,    #
              'SumA2.mlm.p' = SumA2.mlm.p,
              'SumA2.mlm.a' = SumA2.mlm.a,
              'SumI2.mlm.p' = SumI2.mlm.p,
              'SumI2.mlm.i' = SumI2.mlm.i,
              'SumD2.mlm.p' = SumD2.mlm.p,    #
              'SumD2.mlm.d' = SumD2.mlm.d))   #
  }
  ans
}  # moments.gaitdcombo.1par 






 moments.gaitdcombo.pois <-
  function(lambda.p,
           a.mix = NULL, a.mlm = NULL,
           i.mix = NULL, i.mlm = NULL,
           d.mix = NULL, d.mlm = NULL,
           truncate = NULL, max.support = Inf,
           pobs.mix = 0,  # Vector
           pobs.mlm = 0,  # Matrix
           pstr.mix = 0,
           pstr.mlm = 0,  # Ditto
           pdip.mix = 0,
           pdip.mlm = 0,  # Ditto
           byrow.aid = FALSE,  # For pobs.mlm and pstr.mlm
           lambda.a = lambda.p, lambda.i = lambda.p,
           lambda.d = lambda.p,
           type.fitted = "All",  # or "mean"
           moments2 = FALSE) {  # Use this for variances.
  rmlife1 <- ppois(max.support - 1, lambda.p, lower.tail = FALSE) *
             lambda.p
  rmlife2 <- ppois(max.support - 2, lambda.p, lower.tail = FALSE) * 
             lambda.p^2 + rmlife1

  mylist1 <- moments.gaitdcombo.1par(theta.p = lambda.p,
           a.mix = a.mix, a.mlm = a.mlm,
           i.mix = i.mix, i.mlm = i.mlm,
           d.mix = d.mix, d.mlm = d.mlm,
           truncate = truncate, max.support = max.support,
           pobs.mix = pobs.mix, pobs.mlm = pobs.mlm,
           pstr.mix = pstr.mix, pstr.mlm = pstr.mlm,
           pdip.mix = pdip.mix, pdip.mlm = pdip.mlm,
           byrow.aid = byrow.aid,  # type.fitted = type.fitted,
           theta.a = lambda.a, theta.i = lambda.i, theta.d = lambda.d,
           moments2 = moments2,
           rmlife1 = rmlife1, rmlife2 = rmlife2,
           dfun = "dpois")


  themean <- with(mylist1,
                  aprd1.mix + iprd1.mix + aprd1.mlm + iprd1.mlm -
                              dprd1.mix             - dprd1.mlm +
                  use.this * (lambda.p - SumA1.mix.p -
                              SumA1.mlm.p - SumT1.p) / (
                  cdf.max.s - SumA0.mix.p - SumA0.mlm.p - SumT0.p))

  if (type.fitted == "mean") {
    return(themean)
  }

  ans <- c(mylist1,
           list('rmlife1'     = rmlife1,  # Has the right dimension
                'mean'        = themean))
  if (moments2) {  # Add more info
  ans <- c(ans,
           list('rmlife2'     = rmlife2))
  }
  ans
}  # moments.gaitdcombo.pois






 moments.gaitdcombo.log <-
  function(shape.p,
           a.mix = NULL, a.mlm = NULL,
           i.mix = NULL, i.mlm = NULL,
           d.mix = NULL, d.mlm = NULL,
           truncate = NULL, max.support = Inf,
           pobs.mix = 0, pobs.mlm = 0,  # Vector and matrix resp.
           pstr.mix = 0, pstr.mlm = 0,  # Ditto
           pdip.mix = 0, pdip.mlm = 0,  # Ditto
           byrow.aid = FALSE,  # For pobs.mlm and pstr.mlm
           shape.a = shape.p, shape.i = shape.p, shape.d = shape.p,
           type.fitted = "All",  # or "mean"
           moments2 = FALSE) {  # Use this for variances.
  A8.p <- -1 / log1p(-shape.p)
  rmlife1 <- A8.p * (shape.p^(max.support + 1)) / (1 - shape.p)
  rmlife2 <- A8.p * ((shape.p^(max.support + 1)) *
                     (max.support + 1 / (1 - shape.p))
                     / (1 - shape.p))

  mylist1 <- moments.gaitdcombo.1par(theta.p = shape.p,
           a.mix = a.mix, a.mlm = a.mlm,
           i.mix = i.mix, i.mlm = i.mlm,
           d.mix = d.mix, d.mlm = d.mlm,
           truncate = truncate, max.support = max.support,
           pobs.mix = pobs.mix, pobs.mlm = pobs.mlm,
           pstr.mix = pstr.mix, pstr.mlm = pstr.mlm,
           pdip.mix = pdip.mix, pdip.mlm = pdip.mlm,
           byrow.aid = byrow.aid,  # type.fitted = type.fitted,
           theta.a = shape.a, theta.i = shape.i, theta.d = shape.d,
           moments2 = moments2,
           rmlife1 = rmlife1, rmlife2 = rmlife2,
           dfun = "dlog")


  themean <- with(mylist1,
                  aprd1.mix + iprd1.mix + aprd1.mlm + iprd1.mlm -
                              dprd1.mix             - dprd1.mlm +
                  use.this *
                  (-shape.p / (log1p(-shape.p) * (1 - shape.p)) -
                  SumA1.mix.p - SumA1.mlm.p - SumT1.p) / (
                  cdf.max.s - SumA0.mix.p - SumA0.mlm.p - SumT0.p))

  if (type.fitted == "mean") {
    return(themean)
  }

  ans <- c(mylist1,
           list('rmlife1'     = rmlife1,  # Has the right dimension
                'mean'        = themean))
  if (moments2) {  # Add more info
  ans <- c(ans,
           list('rmlife2'     = rmlife2))
  }
  ans
}  # moments.gaitdcombo.log






 moments.gaitdcombo.zeta <-
  function(shape.p,
           a.mix = NULL, a.mlm = NULL,
           i.mix = NULL, i.mlm = NULL,
           d.mix = NULL, d.mlm = NULL,
           truncate = NULL, max.support = Inf,
           pobs.mix = 0, pobs.mlm = 0,  # Vector and matrix resp.
           pstr.mix = 0, pstr.mlm = 0,  # Ditto
           pdip.mix = 0, pdip.mlm = 0,  # Ditto
           byrow.aid = FALSE,  # For pobs.mlm and pstr.mlm
           shape.a = shape.p, shape.i = shape.p, shape.d = shape.p,
           type.fitted = "All",  # or "mean"
           moments2 = FALSE) {  # Use this for variances.
  rmlife1 <- if (is.finite(max.support)) zeta(shape.p) * (1 -
    pzeta(max.support, shape.p - 1)) / zeta(shape.p + 1) else
    numeric(length(shape.p))
  rmlife1[shape.p <= 1] <- NA  # NA or Inf, not sure
  rmlife2 <- if (is.finite(max.support)) zeta(shape.p - 1) * (1 -
    pzeta(max.support, shape.p - 2)) / zeta(shape.p + 1) else
    numeric(length(shape.p))
  rmlife2[shape.p <= 2] <- NA  # NA or Inf, not sure

  mylist1 <- moments.gaitdcombo.1par(theta.p = shape.p,
           a.mix = a.mix, a.mlm = a.mlm,
           i.mix = i.mix, i.mlm = i.mlm,
           d.mix = d.mix, d.mlm = d.mlm,
           truncate = truncate, max.support = max.support,
           pobs.mix = pobs.mix, pobs.mlm = pobs.mlm,
           pstr.mix = pstr.mix, pstr.mlm = pstr.mlm,
           pdip.mix = pdip.mix, pdip.mlm = pdip.mlm,
           byrow.aid = byrow.aid,  # type.fitted = type.fitted,
           theta.a = shape.a, theta.i = shape.i, theta.d = shape.d,
           moments2 = moments2,
           rmlife1 = rmlife1, rmlife2 = rmlife2,
           dfun = "dzeta")

  themean <-
    with(mylist1,
         aprd1.mix + iprd1.mix + aprd1.mlm + iprd1.mlm -
                     dprd1.mix             - dprd1.mlm +
         use.this *
         (ifelse(shape.p > 1, zeta(shape.p) / zeta(shape.p + 1), NA) -
          SumA1.mix.p - SumA1.mlm.p - SumT1.p) / (
          cdf.max.s - SumA0.mix.p - SumA0.mlm.p - SumT0.p))
  if (type.fitted == "mean") {
    return(themean)
  }

  ans <- c(mylist1,
           list('rmlife1'     = rmlife1,  # Has the right dimension
                'mean'        = themean))
  if (moments2) {  # Add more info
  ans <- c(ans,
           list('rmlife2'     = rmlife2))
  }
  ans
}  # moments.gaitdcombo.zeta














specialsvglm <-
  function(object, ...) {
  infos <- object@family@infos()
  ans <- list(a.mix  = infos$a.mix,
              a.mlm  = infos$a.mlm,
              i.mix  = infos$i.mix,
              i.mlm  = infos$i.mlm,
              d.mix  = infos$d.mix,
              d.mlm  = infos$d.mlm,
              truncate = infos$truncate)
  if (is.numeric(tmp7e <- infos$max.support))
    ans <- c(ans, max.support = tmp7e)
  ans
}  # specialsvglm


if (!isGeneric("specials"))
  setGeneric("specials", function(object, ...)
             standardGeneric("specials"),
             package = "VGAM")
setMethod("specials", "vglm",
          function(object, ...)
          specialsvglm(object, ...))




if (!isGeneric("altered"))
  setGeneric("altered", function(object, ...)
             standardGeneric("altered"),
             package = "VGAM")

setMethod("altered", "vglm",
          function(object, ...) {
          tmp <- specialsvglm(object, ...)
          c(tmp$a.mix, tmp$a.mlm)})




if (!isGeneric("inflated"))
  setGeneric("inflated", function(object, ...)
             standardGeneric("inflated"),
             package = "VGAM")

setMethod("inflated", "vglm",
          function(object, ...) {
          tmp <- specialsvglm(object, ...)
          c(tmp$i.mix, tmp$i.mlm)})


if (!isGeneric("truncated"))
  setGeneric("truncated", function(object, ...)
             standardGeneric("truncated"),
             package = "VGAM")

setMethod("truncated", "vglm",
          function(object, ...) {
          ans <- specialsvglm(object, ...)
          if (any(names(ans) == "max.support"))
            ans[c("truncate", "max.support")] else
            ans[["truncate"]]
          })









  setGeneric("is.altered", function(object, ...)
             standardGeneric("is.altered"),
             package = "VGAM")
setMethod("is.altered", "vglm",
          function(object, ...) {
          tmp <- specialsvglm(object, ...)
          as.logical(length(c(tmp$a.mix, tmp$a.mlm)))})



  setGeneric("is.inflated", function(object, ...)
             standardGeneric("is.inflated"),
             package = "VGAM")
setMethod("is.inflated", "vglm",
          function(object, ...) {
          tmp <- specialsvglm(object, ...)
          as.logical(length(c(tmp$i.mix, tmp$i.mlm)))})



  setGeneric("is.deflated", function(object, ...)
             standardGeneric("is.deflated"),
             package = "VGAM")
setMethod("is.deflated", "vglm",
          function(object, ...) {
          tmp <- specialsvglm(object, ...)
          as.logical(length(c(tmp$d.mix, tmp$d.mlm)))})




  setGeneric("is.truncated", function(object, ...)
             standardGeneric("is.truncated"),
             package = "VGAM")
setMethod("is.truncated", "vglm",
          function(object, ...) {
          tmp <- specialsvglm(object, ...)
          as.logical(length(tmp$truncated)) ||
          (length(tmp$max.support) > 0 && is.finite(tmp$max.support))
          })







 y.gaitcombo.check <-
  function(y, truncate = NULL,
           a.mix = NULL, a.mlm = NULL,
           i.mix = NULL, i.mlm = NULL,
           d.mix = NULL, d.mlm = NULL,
           max.support = Inf, min.support = 0) {
  la.mix <- length(a.mix <- sort(a.mix))
  li.mix <- length(i.mix <- sort(i.mix))
  ld.mix <- length(d.mix <- sort(d.mix))
  la.mlm <- length(a.mlm)
  li.mlm <- length(i.mlm)
  ld.mlm <- length(d.mlm)




  n <- length(y)
   css.mix.a <-  css.mix.i <-  css.mix.d <-
   css.mlm.a <-  css.mlm.i <-  css.mlm.d <- NULL
  skip.mix.a <- skip.mix.i <- skip.mix.d <-  # Default
  skip.mlm.a <- skip.mlm.i <- skip.mlm.d <- NULL

  if (length(truncate) && any(y %in% truncate))
    stop("some response values == values in argument 'truncate'")
  if (max.support < max(y))
    stop("some response values are greater than the ",
         "'max.support' argument")


  y0.mix.a <- y0.mlm.a <-
  y0.mix.i <- y0.mlm.i <-
  y0.mix.d <- y0.mlm.d <- NULL


  if (la.mix > 0) {  # ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
    y0.mix.a <- matrix(0, n, la.mix)
    for (jay in seq(la.mix))
      y0.mix.a[, jay] <- as.numeric(y == a.mix[jay])
    skip.mix.a <- matrix(as.logical(y0.mix.a), n, la.mix)  # dim lost
    if (any((css.mix.a <- colSums(skip.mix.a)) == 0))
      stop("some 'a.mix' argument values have no response values: ",
           paste(a.mix[css.mix.a == 0], collapse = ", "))          
  }  # la.mix

  if (la.mlm > 0) {  # ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
    y0.mlm.a <- matrix(0, n, la.mlm)
    for (jay in seq(la.mlm))
      y0.mlm.a[, jay] <- as.numeric(y == a.mlm[jay])
    skip.mlm.a <- matrix(as.logical(y0.mlm.a), n, la.mlm)  # dim lost
    if (any((css.mlm.a <- colSums(skip.mlm.a)) == 0))
      stop("some 'a.mlm' argument values have no response values: ",
           paste(a.mlm[css.mlm.a == 0], collapse = ", "))          
  }  # la.mlm


  if (li.mix > 0) {  # ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
    y0.mix.i <- matrix(0, n, li.mix)
    for (jay in seq(li.mix))
      y0.mix.i[, jay] <- as.numeric(y == i.mix[jay])
    skip.mix.i <- matrix(as.logical(y0.mix.i), n, li.mix)  # dim lost
      if (any((css.mix.i <- colSums(skip.mix.i)) == 0))
      stop("some 'i.mix' argument values have no response values: ",
           paste(i.mix[css.mix.i == 0], collapse = ", "))          
  }  # li.mix

      
  if (li.mlm > 0) {  # ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
    y0.mlm.i <- matrix(0, n, li.mlm)
    for (jay in seq(li.mlm))
      y0.mlm.i[, jay] <- as.numeric(y == i.mlm[jay])
    skip.mlm.i <- matrix(as.logical(y0.mlm.i), n, li.mlm)  # dim lost
      if (any((css.mlm.i <- colSums(skip.mlm.i)) == 0))
      stop("some 'i.mlm' argument values have no response values: ",
           paste(i.mlm[css.mlm.i == 0], collapse = ", "))          
  }  # li.mlm



  if (ld.mix > 0) {  # ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
    y0.mix.d <- matrix(0, n, ld.mix)
    for (jay in seq(ld.mix))
      y0.mix.d[, jay] <- as.numeric(y == d.mix[jay])
    skip.mix.d <- matrix(as.logical(y0.mix.d), n, ld.mix)  # dim lost
      if (any((css.mix.d <- colSums(skip.mix.d)) == 0))
      stop("some 'd.mix' argument values have no response values: ",
           paste(d.mix[css.mix.d == 0], collapse = ", "))          
  }  # ld.mix

      
  if (ld.mlm > 0) {  # ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
    y0.mlm.d <- matrix(0, n, ld.mlm)
    for (jay in seq(ld.mlm))
      y0.mlm.d[, jay] <- as.numeric(y == d.mlm[jay])
    skip.mlm.d <- matrix(as.logical(y0.mlm.d), n, ld.mlm)  # dim lost
      if (any((css.mlm.d <- colSums(skip.mlm.d)) == 0))
      stop("some 'd.mlm' argument values have no response values: ",
           paste(d.mlm[css.mlm.d == 0], collapse = ", "))          
  }  # ld.mlm




  list(css.mix.a = css.mix.a, skip.mix.a = skip.mix.a,
       css.mix.i = css.mix.i, skip.mix.i = skip.mix.i,
       css.mix.d = css.mix.d, skip.mix.d = skip.mix.d,
       css.mlm.a = css.mlm.a, skip.mlm.a = skip.mlm.a,
       css.mlm.i = css.mlm.i, skip.mlm.i = skip.mlm.i,
       css.mlm.d = css.mlm.d, skip.mlm.d = skip.mlm.d,
        y0.mix.a =  y0.mix.a,   y0.mlm.a =   y0.mlm.a, 
        y0.mix.i =  y0.mix.i,   y0.mlm.i =   y0.mlm.i,  
        y0.mix.d =  y0.mix.d,   y0.mlm.d =   y0.mlm.d) 
}  # y.gaitcombo.check




















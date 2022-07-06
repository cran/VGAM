# These functions are
# Copyright (C) 1998-2022 T.W. Yee, University of Auckland.
# All rights reserved.











gaitdpoisson.control <-
gaitdlog.control <-
gaitdzeta.control <-
gaitdnbinomial.control <-  # Overwrites the summary() default.
  function(summary.HDEtest = FALSE,
           ...) {
  list(summary.HDEtest = summary.HDEtest)
}








 goffset <-
  function(mux, n,
           a.mix = NULL, i.mix = NULL, d.mix = NULL,
           a.mlm = NULL, i.mlm = NULL, d.mlm = NULL,
           par1or2 = 1
           ) {
  if (!is.Numeric(mux, integer.valued = TRUE, positive = TRUE,
                  length.arg = 1))
    stop("bad input for argument 'mux'")
  if (!is.Numeric(n, integer.valued = TRUE, positive = TRUE,
                  length.arg = 1))
    stop("bad input for argument 'n'")
  if (!is.Numeric(par1or2, integer.valued = TRUE, positive = TRUE,
                  length.arg = 1) || par1or2 > 2)
    stop("bad input for argument 'par1or2'")


  gaitd.errorcheck(a.mix, a.mlm, i.mix, i.mlm,
                   d.mix, d.mlm)  # , truncate


  la.mix <- length(a.mix <- sort(a.mix))
  li.mix <- length(i.mix <- sort(i.mix))
  ld.mix <- length(d.mix <- sort(d.mix))
  la.mlm <- length(a.mlm)
  li.mlm <- length(i.mlm)
  ld.mlm <- length(d.mlm)
  tmp3.TF <- c(rep(TRUE, par1or2),
               la.mix > 0, rep(la.mix > 1, par1or2),
               li.mix > 0, rep(li.mix > 1, par1or2),
               ld.mix > 0, rep(ld.mix > 1, par1or2),
               la.mlm > 0, li.mlm > 0, ld.mlm > 0)
  indeta.finish <- cumsum(c(rep(1, par1or2),
                            1, rep(1, par1or2),
                            1, rep(1, par1or2),
                            1, rep(1, par1or2),
                            la.mlm, li.mlm, ld.mlm,
                            ld.mlm + 1) * c(tmp3.TF, 1))
  indeta.launch <- c(1, 1 + head(indeta.finish, -1))

  indeta.launch <- head(indeta.launch, -1)
  indeta.finish <- head(indeta.finish, -1)
  indeta.launch[!tmp3.TF] <- NA  # Not to be accessed
  indeta.finish[!tmp3.TF] <- NA  # Not to be accessed
  indeta <- cbind(launch = indeta.launch,
                  finish = indeta.finish)
  if (FALSE && par1or2 == 1)
  rownames(indeta) <- c("lambda.p",
                        "pobs.mix", "lambda.a",
                        "pstr.mix", "lambda.i",
                        "pdip.mix", "lambda.d",
                        "pobs.mlm", "pstr.mlm", "pdip.mlm")
  if (FALSE && par1or2 == 2)
  rownames(indeta) <- c("munb.p", "size.p",
                        "pobs.mix", "munb.a", "size.a",
                        "pstr.mix", "munb.i", "size.i",
                        "pdip.mix", "munb.d", "size.d",
                        "pobs.mlm", "pstr.mlm", "pdip.mlm")
  M1 <- max(indeta, na.rm = TRUE)
  Mat <- matrix(0, n, M1)


  colptr <- indeta[if (par1or2 == 1) c(1, 3, 5, 7) else
                   c(1, 4, 7, 10), 'launch']
  colptr <- na.omit(colptr)
  Mat[, colptr] <- log(mux)


  Mat
}  # goffset







 cm3gaitd <-
  function(eq.ap = FALSE, eq.ip = FALSE, eq.dp = FALSE, npar = 1) {
    M <- 4 * npar + 3
    use.mat.mix <- diag(M)  # Full model constraint matrices



    if ( (eq.ap) &&  (eq.ip) &&  (eq.dp))
      use.mat.mix <- matrix(c(1, 0, 0, 0,
                              0, 1, 0, 0,
                              1, 0, 0, 0,
                              0, 0, 1, 0,
                              1, 0, 0, 0,
                              0, 0, 0, 1,
                              1, 0, 0, 0), 7, 4, byrow = TRUE)
    if ( (eq.ap) &&  (eq.ip) && !(eq.dp))
      use.mat.mix <- matrix(c(1, 0, 0, 0, 0,
                              0, 1, 0, 0, 0,
                              1, 0, 0, 0, 0,
                              0, 0, 1, 0, 0,
                              1, 0, 0, 0, 0,
                              0, 0, 0, 1, 0,
                              0, 0, 0, 0, 1), 7, 5, byrow = TRUE)
    if ( (eq.ap) && !(eq.ip) &&  (eq.dp))
      use.mat.mix <- matrix(c(1, 0, 0, 0, 0,
                              0, 1, 0, 0, 0,
                              1, 0, 0, 0, 0,
                              0, 0, 1, 0, 0,
                              0, 0, 0, 0, 1,
                              0, 0, 0, 1, 0,
                              1, 0, 0, 0, 0), 7, 5, byrow = TRUE)
    if ( (eq.ap) && !(eq.ip) && !(eq.dp))
      use.mat.mix <- matrix(c(1, 0, 0, 0, 0, 0,
                              0, 1, 0, 0, 0, 0,
                              1, 0, 0, 0, 0, 0,
                              0, 0, 1, 0, 0, 0,
                              0, 0, 0, 0, 1, 0,
                              0, 0, 0, 1, 0, 0,
                              0, 0, 0, 0, 0, 1), 7, 6, byrow = TRUE)
    if (!(eq.ap) &&  (eq.ip) &&  (eq.dp))
      use.mat.mix <- matrix(c(1, 0, 0, 0, 0,
                              0, 1, 0, 0, 0,
                              0, 0, 0, 0, 1,
                              0, 0, 1, 0, 0,
                              1, 0, 0, 0, 0,
                              0, 0, 0, 1, 0,
                              1, 0, 0, 0, 0), 7, 5, byrow = TRUE)
    if (!(eq.ap) &&  (eq.ip) && !(eq.dp))
      use.mat.mix <- matrix(c(1, 0, 0, 0, 0, 0,
                              0, 1, 0, 0, 0, 0,
                              0, 0, 0, 0, 1, 0,
                              0, 0, 1, 0, 0, 0,
                              1, 0, 0, 0, 0, 0,
                              0, 0, 0, 1, 0, 0,
                              0, 0, 0, 0, 0, 1), 7, 6, byrow = TRUE)
    if (!(eq.ap) && !(eq.ip) &&  (eq.dp))
      use.mat.mix <- matrix(c(1, 0, 0, 0, 0, 0,
                              0, 1, 0, 0, 0, 0,
                              0, 0, 0, 0, 1, 0,
                              0, 0, 1, 0, 0, 0,
                              0, 0, 0, 0, 0, 1,
                              0, 0, 0, 1, 0, 0,
                              1, 0, 0, 0, 0, 0), 7, 6, byrow = TRUE)

    if (eq.ap + eq.ip + eq.dp && npar > 1) {
      use.mat.mix <- kronecker(use.mat.mix, diag(npar))
      ind.keep <- seq(npar + 1) 
      ind.keep <- c(ind.keep, 2*npar + seq(npar + 1))
      ind.keep <- c(ind.keep, 4*npar + seq(npar + 1))
      ind.keep <- c(ind.keep, 6*npar + seq(npar))
      use.mat.mix <- use.mat.mix[ind.keep, ]
    }

  use.mat.mix
}  # cm3gaitd













 gaitdpoisson <-
  function(a.mix = NULL, i.mix = NULL, 
           d.mix = NULL,
           a.mlm = NULL, i.mlm = NULL,  # Unstructured probs are
           d.mlm = NULL,                # contiguous
           truncate = NULL, max.support = Inf,
           zero = c("pobs", "pstr", "pdip"),  # Pruned, handles all 6
           eq.ap = TRUE, eq.ip = TRUE, eq.dp = TRUE,
           parallel.a = FALSE, parallel.i = FALSE,
           parallel.d = FALSE,
           llambda.p = "loglink",
           llambda.a = llambda.p,  # "loglink", 20201117
           llambda.i = llambda.p,  # "loglink", 20201117
           llambda.d = llambda.p,  # "loglink", 20211011
           type.fitted = c("mean", "lambdas",
                           "pobs.mlm", "pstr.mlm", "pdip.mlm",
                           "pobs.mix", "pstr.mix", "pdip.mix",
                           "Pobs.mix", "Pstr.mix", "Pdip.mix",
                           "nonspecial", "Numer", "Denom.p",
                           "sum.mlm.i", "sum.mix.i",
                           "sum.mlm.d", "sum.mix.d",
                           "ptrunc.p", "cdf.max.s"),
           gpstr.mix = ppoints(7) / 3,  # ppoints(9) / 2,
           gpstr.mlm = ppoints(7) / (3 + length(i.mlm)),
           imethod = 1,
           mux.init = c(0.75, 0.5, 0.75),   # Order is A, I, D.
           ilambda.p = NULL, ilambda.a = ilambda.p,
           ilambda.i = ilambda.p, ilambda.d = ilambda.p,
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



  lowsup <- 0
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

  llambda.p <- as.list(substitute(llambda.p))
  elambda.p <- link2list(llambda.p)
  llambda.p <- attr(elambda.p, "function.name")
  llambda.p.save <- llambda.p

  lpobs.mix <- "multilogitlink"  # \omega_p
  epobs.mix <- list()  # zz NULL for now 20200907 coz 'multilogitlink'
  elambda.a <- link2list(llambda.a)
  llambda.a <- attr(elambda.a, "function.name")

  lpstr.mix <- "multilogitlink"  # \phi_p
  epstr.mix <- list()  # zz NULL for now 20200907 coz 'multilogitlink'
  lpdip.mix <- "multilogitlink"  # zz unsure 20211002
  epdip.mix <- list()  # zz unsure 20211002
  elambda.i <- link2list(llambda.i)
  llambda.i <- attr(elambda.i, "function.name")
  elambda.d <- link2list(llambda.d)
  llambda.d <- attr(elambda.d, "function.name")


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
           poissonff(link = .llambda.p.save , zero = NULL),
           list( .llambda.p.save = llambda.p.save))))

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


  if (FALSE) {  # Comment this out to allow default eq.ap=TRUE, etc.
  if (la.mix <= 1 && eq.ap)
    stop("<= one unstructured altered value (no 'lambda.a')",
         ", so setting 'eq.ap = TRUE' is meaningless")
  if (li.mix <= 1 && eq.ip)
    stop("<= one unstructured inflated value (no 'lambda.i')",
            ", so setting 'eq.ip = TRUE' is meaningless")
  if (ld.mix <= 1 && eq.dp)
    stop("<= one unstructured deflated value (no 'lambda.d')",
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
            c("mean", "lambdas",
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
  tmp3 <- c(lambda.p = llambda.p,
            pobs.mix = if (la.mix) "multilogitlink" else NULL,
            lambda.a = if (la.mix > 1) llambda.a else NULL,
            pstr.mix = if (li.mix) "multilogitlink" else NULL,
            lambda.i = if (li.mix > 1) llambda.i else NULL,
            pdip.mix = if (ld.mix) "multilogitlink" else NULL,
            lambda.d = if (ld.mix > 1) llambda.d else NULL,
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
  rownames(indeta) <- c("lambda.p",
                        "pobs.mix", "lambda.a",
                        "pstr.mix", "lambda.i",
                        "pdip.mix", "lambda.d",
                        "pobs.mlm", "pstr.mlm", "pdip.mlm")
  M1 <- max(indeta, na.rm = TRUE)
  predictors.names <- tmp3  # Passed into @infos and @initialize.
      

  blurb1 <- ""   # zz1
  if (la.mlm + la.mix) blurb1 <- "Generally-altered "
  if (li.mlm + li.mix) blurb1 <- "Generally-inflated "
  if (ltrunc.use) blurb1 <- "Generally-truncated "
  if ( (la.mlm + la.mix) &&  (li.mlm + li.mix) && !ltrunc.use)
    blurb1 <- "Generally-altered and -inflated "
  if ( (la.mlm + la.mix) && !(li.mlm + li.mix) &&  ltrunc.use)
    blurb1 <- "Generally-altered and -truncated "
  if (!(la.mlm + la.mix) &&  (li.mlm + li.mix) &&  ltrunc.use)
    blurb1 <- "Generally-inflated and -truncated "
  if ( (la.mlm + la.mix) &&  (li.mlm + li.mix) &&  ltrunc.use)
    blurb1 <- "Generally-altered, -inflated and -truncated "

  if (ld.mlm + ld.mix) blurb1 <-
    c(blurb1,
      if (la.mlm + la.mix + li.mlm + li.mix) "and " else "Generally",
      "-deflated ")



      
  new("vglmff",
  blurb = c(blurb1, "Poisson regression\n",
            "(GAITD-Pois(lambda.p)-",
                   "Pois(lambda.a)-MLM-",
                   "Pois(lambda.i)-MLM-\n",
                   "Pois(lambda.d)-MLM generally)\n\n",
            "Links: ",
            namesof("lambda.p", llambda.p, earg = elambda.p,
                    tag = FALSE),
            if (la.mix > 0) c(", ", "multilogit(pobs.mix)"),
            if (la.mix > 1) c(", ",
            namesof("lambda.a",  llambda.a, elambda.a, tag = FALSE)),
            if (la.mix && li.mix) ", \n       ",
            if (li.mix > 0) c(  if (la.mix) "" else ", ",
            "multilogit(pstr.mix)"),
            if (li.mix > 1) c(", ",
            namesof("lambda.i",  llambda.i, elambda.i, tag = FALSE)),
            if (li.mix && ld.mix) ", \n       ",
            if (ld.mix > 0) c(  if (li.mix) "" else ", ",
            "multilogit(pdip.mix)"),
            if (ld.mix > 1) c(", ",
            namesof("lambda.d",  llambda.d, elambda.d, tag = FALSE)),
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
      Use.mat <- use.mat.mlm <- cbind(M)  # lambda.p only
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
              col(Use.mat) > ncol(use.mat.mix)] <-
          use.mat.mlm[-1, -1]
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
         dpqrfun = "gaitdpois",
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
         parent.name = c("poissonff", "pois"),
         type.fitted  = as.vector( .type.fitted ),
         type.fitted.choices = ( .type.fitted.choices ),
         baseparams.argnames  = "lambda",
         MM1 = 1,  # One parameter for 1 response (lambda). Needed.
         zero = .zero )
  }, list( .zero = zero, .lowsup = lowsup,
           .type.fitted = type.fitted,
           .type.fitted.choices = type.fitted.choices,
           .llambda.p = llambda.p, .elambda.p = elambda.p,
           .llambda.a = llambda.a, .elambda.a = elambda.a,
           .llambda.i = llambda.i, .elambda.i = elambda.i,
           .a.mlm = a.mlm, .i.mlm = i.mlm, .d.mlm = d.mlm,
           .a.mix = a.mix, .i.mix = i.mix, .d.mix = d.mix,
           .la.mlm = la.mlm, .li.mlm = li.mlm, .ld.mlm = ld.mlm,
           .la.mix = la.mix, .li.mix = li.mix, .ld.mix = ld.mix,
           .truncate = truncate, .max.support = max.support,
           .predictors.names = predictors.names,
           .M1 = M1, .lall.len = lall.len
         ))),

  rqresslot = eval(substitute(
    function(mu, y, w, eta, extra = NULL) {
    if (!is.matrix(eta)) eta <- as.matrix(eta)
    la.mix <- length((a.mix <- as.vector( .a.mix )))
    li.mix <- length((i.mix <- as.vector( .i.mix )))
    ld.mix <- length((d.mix <- as.vector( .d.mix )))
    la.mlm <- length((a.mlm <- as.vector( .a.mlm )))
    li.mlm <- length((i.mlm <- as.vector( .i.mlm )))
    ld.mlm <- length((d.mlm <- as.vector( .d.mlm )))
    truncate <- as.vector( .truncate )

    tmp3.TF <- ( .tmp3.TF )   # Logical of length 10.

    lall.len <- la.mix + li.mix + ld.mix + la.mlm + li.mlm + ld.mlm
    pobs.mix <- pstr.mix <- pdip.mix <- 0  # 4 rowSums()
    pobs.mlm <- pstr.mlm <- pdip.mlm <- 0  # matrix(0, NROW(eta), 1)
    lambda.p <- cbind(eta2theta(eta[, 1], .llambda.p , .elambda.p ))
    ind.lambda.z <- 1  # Points to lambda.p only.
    lambda.a <- lambda.i <-
    lambda.d <- lambda.p  # Needed and doesnt corrupt the answer

    if (any(tmp3.TF[c(3, 5, 7)])) {  # At least 1 lambda.[aid]
      ind.lambda.z <- extra$indeta[c(1, 3, 5, 7), 'launch']  # Vecs
      ind.lambda.z <- c(na.omit(ind.lambda.z))  # At least 1 value

      lambda.a <- if (!tmp3.TF[ 3]) lambda.p else
        eta2theta(eta[, extra$indeta[3, 1]], .llambda.a , .elambda.a )
      lambda.i <- if (!tmp3.TF[ 5]) lambda.p else
        eta2theta(eta[, extra$indeta[5, 1]], .llambda.i , .elambda.i )
      lambda.d <- if (!tmp3.TF[ 7]) lambda.p else
        eta2theta(eta[, extra$indeta[7, 1]], .llambda.d , .elambda.d )
    }  # la.mix + li.mix + ld.mix > 0

    if (lall.len) {  # An MLM was fitted
      allprobs <-
        multilogitlink(eta[, -ind.lambda.z, drop = FALSE],
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
        dimnames(pobs.mlm) <- list(rownames(eta),
                                   as.character(a.mlm))
        Nextone <- Nextone + la.mlm
      }
      if (tmp3.TF[ 9]) {
        ind9 <- (Nextone + 1):(Nextone + li.mlm)
        pstr.mlm <- allprobs[, ind9, drop = FALSE]
        dimnames(pstr.mlm) <- list(rownames(eta),
                                   as.character(i.mlm))
        Nextone <- Nextone + li.mlm
      }
      if (tmp3.TF[10]) {
        ind10 <- (Nextone + 1):(Nextone + ld.mlm)
        pdip.mlm <- allprobs[, ind10, drop = FALSE]
        dimnames(pdip.mlm) <- list(rownames(eta),
                                   as.character(d.mlm))
        Nextone <- Nextone + ld.mlm  # Not needed
      }
    }  # lall.len


    scrambleseed <- runif(1)  # To scramble the seed
    qnorm(runif(length(y),
        pgaitdpois(y - 1, lambda.p = lambda.p,
                  a.mix = a.mix, i.mix = i.mix, d.mix = d.mix,
                  a.mlm = a.mlm, i.mlm = i.mlm, d.mlm = d.mlm,
                  truncate = truncate,
                  max.support = as.vector( .max.support ),
                  lambda.a = lambda.a, lambda.i = lambda.i, 
                  lambda.d = lambda.d,
                  pobs.mix = pobs.mix, pstr.mix = pstr.mix,
                  pdip.mix = pdip.mix,
                  pobs.mlm = pobs.mlm, pstr.mlm = pstr.mlm,
                  pdip.mlm = pdip.mlm),
        pgaitdpois(y    , lambda.p = lambda.p,
                  a.mix = a.mix, i.mix = i.mix, d.mix = d.mix,
                  a.mlm = a.mlm, i.mlm = i.mlm, d.mlm = d.mlm,
                  truncate = truncate,
                  max.support = as.vector( .max.support ),
                  lambda.a = lambda.a, lambda.i = lambda.i, 
                  lambda.d = lambda.d,
                  pobs.mix = pobs.mix, pstr.mix = pstr.mix,
                  pdip.mix = pdip.mix,
                  pobs.mlm = pobs.mlm, pstr.mlm = pstr.mlm,
                  pdip.mlm = pdip.mlm)))
  }, list(
    .llambda.p = llambda.p, .elambda.p = elambda.p,
    .llambda.a = llambda.a, .elambda.a = elambda.a,
    .llambda.i = llambda.i, .elambda.i = elambda.i,
    .llambda.d = llambda.d, .elambda.d = elambda.d,
    .lpstr.mix = lpstr.mix, .lpobs.mix = lpobs.mix,
    .lpdip.mix = lpdip.mix,
    .epstr.mix = epstr.mix, .epobs.mix = epobs.mix,
    .epdip.mix = epdip.mix,
    .tmp3.TF = tmp3.TF,
    .a.mix = a.mix, .i.mix = i.mix, .d.mix = d.mix,
    .a.mlm = a.mlm, .i.mlm = i.mlm, .d.mlm = d.mlm,
    .truncate = truncate, .max.support = max.support ))),

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
                               max.support = .max.support )
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
      lambda.a.init <-  lambda.i.init <-  # Needed
      lambda.d.init <-
      lambda.p.init <- Init.mu(y = y, w = w, imethod = .imethod ,
                               imu = .ilambda.p ,  # x = x,
                               ishrinkage = .ishrinkage ,
                               probs.y = .probs.y )
      etastart <- matrix(nrow = n, ncol = M,
        theta2eta(lambda.p.init, .llambda.p , earg = .elambda.p ))


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
        lambda.a.init <- if (length( .ilambda.a ))
          rep_len( .ilambda.a , n) else lambda.p.init  # A vector
        etastart[, 3] <-
          theta2eta(lambda.a.init, .llambda.a , earg = .elambda.a )
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
        lambda.i.init <- if (length( .ilambda.i ))
          rep_len( .ilambda.i , n) else lambda.p.init  # A vector
        etastart[, (extra$indeta[5, 'launch'])] <-
          theta2eta(lambda.i.init, .llambda.i , earg = .elambda.i )
      }  # li.mix > 1


      if (tmp3.TF[ 8]) {  #  la.mlm
        init.pobs.mlm <- if (length( .ipobs.mlm )) {
          matrix( .ipobs.mlm , n, la.mlm, byrow = .byrow.aid )
        } else {
          mux.more.a <- extra$mux.init[1]                            
          init.pobs.mlm <- colSums(c(w) *
                                   extra$skip.mlm.a) / colSums(w)
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





      gaitdpois.Loglikfun1.mix <-
        function(pstr.mix.val, y, x, w, extraargs) {
    sum(c(w) *
        dgaitdpois(y, pstr.mix = pstr.mix.val,
                  pstr.mlm    = extraargs$pstr.mlm,  # Differs here
                  lambda.p    = extraargs$lambda.p,
                  lambda.a    = extraargs$lambda.a,
                  lambda.i    = extraargs$lambda.i,
                  lambda.d    = extraargs$lambda.d,
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

 gaitdpois.Loglikfun1.mlm <-
     function(pstr.mlm.val, y, x, w, extraargs) {
    sum(c(w) *
        dgaitdpois(y, pstr.mlm = pstr.mlm.val,
                   pstr.mix    = extraargs$pstr.mix,  # Differs here
                   lambda.p    = extraargs$lambda.p,
                   lambda.a    = extraargs$lambda.a,
                   lambda.i    = extraargs$lambda.i,
                   lambda.d    = extraargs$lambda.d,
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

 gaitdpois.Loglikfun2 <-
     function(pstr.mix.val, pstr.mlm.val, y, x, w, extraargs) {
    sum(c(w) *
        dgaitdpois(y, pstr.mix = pstr.mix.val,
                   pstr.mlm = pstr.mlm.val,
                   lambda.p    = extraargs$lambda.p,
                   lambda.a    = extraargs$lambda.a,
                   lambda.i    = extraargs$lambda.i,
                   lambda.d    = extraargs$lambda.d,
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
              lambda.p    = lambda.p.init,
              lambda.a    = lambda.a.init,
              lambda.i    = lambda.i.init,
              lambda.d    = lambda.d.init,
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
                         objfun = gaitdpois.Loglikfun2,
                         y = y, w = w, extraargs = extraargs,
                         ret.objfun = TRUE)
          } else if (try.gridsearch.pstr.mix) {
            extraargs$pstr.mlm <- init.pstr.mlm
            grid.search ( .gpstr.mix ,
                         objfun = gaitdpois.Loglikfun1.mix,
                         y = y, w = w, extraargs = extraargs,
                         ret.objfun = TRUE)
          } else if (try.gridsearch.pstr.mlm) {
            extraargs$pstr.mix <- init.pstr.mix
            grid.search ( .gpstr.mlm ,
                         objfun = gaitdpois.Loglikfun1.mlm,
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
        etastart[, (extra$indeta[4, 'launch'])] <-
            etastart.z[, nextone]
        nextone <- nextone + 1
      }
      if (tmp3.TF[ 6]) {  # Coln 2 or 4 or 6
        etastart[, (extra$indeta[6, 'launch'])] <-
            etastart.z[, nextone]
        nextone <- nextone + 1
      }
      if (tmp3.TF[ 8]) {
        ind8 <- (extra$indeta[8, 'launch']):(extra$indeta[8,
                                                          'finish'])
        etastart[, ind8] <- etastart.z[, nextone:(nextone+la.mlm - 1)]
        nextone <- nextone + la.mlm
      }
      if (tmp3.TF[ 9]) {
        ind9 <- (extra$indeta[9, 'launch']):(extra$indeta[9,
                                                          'finish'])
        etastart[, ind9] <- etastart.z[, nextone:(nextone+li.mlm - 1)]
        nextone <- nextone + li.mlm
      }
      if (tmp3.TF[10]) {
        ind0 <- (extra$indeta[10, 'launch']):(extra$indeta[10,
                                                           'finish'])
        etastart[, ind0] <- etastart.z[, nextone:(nextone +
                                                  ld.mlm - 1)]
        if (ncol(etastart.z) != nextone + ld.mlm - 1)
          stop("miscalculation")
      }
    }
  }), list(
    .llambda.p = llambda.p, .elambda.p = elambda.p,
    .llambda.a = llambda.a, .elambda.a = elambda.a,
    .llambda.i = llambda.i, .elambda.i = elambda.i,
    .llambda.d = llambda.d, .elambda.d = elambda.d,
    .ilambda.p = ilambda.p,
    .ilambda.a = ilambda.a,
    .ilambda.i = ilambda.i,
    .ilambda.d = ilambda.d,
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
    .tmp3.TF = tmp3.TF,  # .tmp3 = tmp3,
    .predictors.names = predictors.names,
    .mux.init = mux.init,
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
                c("mean", "lambdas",
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
    lambda.p <- cbind(eta2theta(eta[, 1], .llambda.p , .elambda.p ))
    ind.lambda.z <- 1  # Points to lambda.p only.
    lambda.a <- lambda.i <-
    lambda.d <- lambda.p  # Needed; and answer not corrupted
    tmp3.TF <- ( .tmp3.TF )   # Logical of length 10.

    if (any(tmp3.TF[c(3, 5, 7)])) {  # At least one lambda.[aid]
      ind.lambda.z <- extra$indeta[c(1, 3, 5, 7), 'launch']  # Vectors
      ind.lambda.z <- c(na.omit(ind.lambda.z))  # At least one value
      lambda.a <- if (!tmp3.TF[ 3]) lambda.p else
        eta2theta(eta[, extra$indeta[3, 1]], .llambda.a , .elambda.a )
      lambda.i <- if (!tmp3.TF[ 5]) lambda.p else
        eta2theta(eta[, extra$indeta[5, 1]], .llambda.i , .elambda.i )
      lambda.d <- if (!tmp3.TF[ 7]) lambda.p else
        eta2theta(eta[, extra$indeta[7, 1]], .llambda.d , .elambda.d )
    }  # la.mix + li.mix + ld.mix > 0






    if (lall.len) {  # An MLM was fitted
      allprobs <-
        multilogitlink(eta[, -ind.lambda.z, drop = FALSE],
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
        dimnames(pobs.mlm) <- list(rownames(eta),
                                   as.character(a.mlm))
        Nextone <- Nextone + la.mlm
      }
      if (tmp3.TF[ 9]) {
        ind9 <- (Nextone + 1):(Nextone + li.mlm)
        pstr.mlm <- allprobs[, ind9, drop = FALSE]
        dimnames(pstr.mlm) <- list(rownames(eta),
                                   as.character(i.mlm))
        Nextone <- Nextone + li.mlm
      }
      if (tmp3.TF[10]) {
        ind10 <- (Nextone + 1):(Nextone + ld.mlm)
        pdip.mlm <- allprobs[, ind10, drop = FALSE]
        dimnames(pdip.mlm) <- list(rownames(eta),
                                   as.character(d.mlm))
        Nextone <- Nextone + ld.mlm  # Not needed
      }
    }  # lall.len

    M1 <- max(extra$indeta, na.rm = TRUE)
    NOS <- NCOL(eta) / M1

    Bits <- moments.gaitdcombo.pois(lambda.p,
              pobs.mix = pobs.mix, pobs.mlm = pobs.mlm,
              pstr.mix = pstr.mix, pstr.mlm = pstr.mlm,
              pdip.mix = pdip.mix, pdip.mlm = pdip.mlm,
              a.mix = a.mix, i.mix = i.mix, d.mix = d.mix,
              a.mlm = a.mlm, i.mlm = i.mlm, d.mlm = d.mlm,
              lambda.a = lambda.a, lambda.i = lambda.i,
              lambda.d = lambda.d,
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
        dpois(matrix(a.mix,    n.obs, la.mix, byrow = TRUE),
              matrix(lambda.a, n.obs, la.mix)) / (
        c(Bits[["SumA0.mix.a"]]))
      dim(tmp13) <- c(n.obs, la.mix)
      dimnames(tmp13) <- list(rownames(eta),
                              as.character(a.mix))
      propn.mat.a <- tmp13
    }  # la.mix

    if (li.mix && morework) {
      tmp55 <-  # dpois() does not retain the matrix format
        dpois(matrix(i.mix,    n.obs, li.mix, byrow = TRUE),
              matrix(lambda.i, n.obs, li.mix)) / (
        c(Bits[["SumI0.mix.i"]]))
      dim(tmp55) <- c(n.obs, li.mix)
      dimnames(tmp55) <- list(rownames(eta),
                              as.character(i.mix))
      propn.mat.i <- tmp55  # Correct dimension
    }  # li.mix


    if (ld.mix && morework) {
      tmp55 <-  # dpois() does not retain the matrix format
        dpois(matrix(d.mix,    n.obs, ld.mix, byrow = TRUE),
              matrix(lambda.d, n.obs, ld.mix)) / (
        c(Bits[["SumD0.mix.d"]]))
      dim(tmp55) <- c(n.obs, ld.mix)
      dimnames(tmp55) <- list(rownames(eta),
                              as.character(d.mix))
      propn.mat.d <- tmp55  # Correct dimension
    }  # ld.mix

    ans <- switch(type.fitted,
      "mean"       = Bits[["mean"]],  # Unconditional mean
      "lambdas"    = cbind(lambda.p,
                           if (tmp3.TF[ 3]) lambda.a else NULL,
                           if (tmp3.TF[ 5]) lambda.i else NULL,
                           if (tmp3.TF[ 7]) lambda.d else NULL),
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
             dpois(matrix(i.mlm,    n.obs, li.mlm, byrow = TRUE),
                   matrix(lambda.p, n.obs, li.mlm)) / Denom.p,
      "sum.mlm.d"  = -pdip.mlm + Numer *
             dpois(matrix(d.mlm,    n.obs, ld.mlm, byrow = TRUE),
                   matrix(lambda.p, n.obs, ld.mlm)) / Denom.p,
      "sum.mix.i"  =  c(pstr.mix) * propn.mat.i + Numer *
             dpois(matrix(i.mix,    n.obs, li.mix, byrow = TRUE),
                   matrix(lambda.p, n.obs, li.mix)) / Denom.p,
      "sum.mix.d"  = -c(pdip.mix) * propn.mat.d + Numer *
             dpois(matrix(d.mix,    n.obs, ld.mix, byrow = TRUE),
                   matrix(lambda.p, n.obs, ld.mix)) / Denom.p,
      "ptrunc.p"   = Bits[["SumT0.p"]] + 1 - Bits[["cdf.max.s"]],
      "cdf.max.s"  = Bits[["cdf.max.s"]])  # Pr(y <= max.support)



    ynames.pobs.mlm <- as.character(a.mlm)  # Works with NULLs
    ynames.pstr.mlm <- as.character(i.mlm)  # Works with NULLs
    ynames.pdip.mlm <- as.character(d.mlm)  # Works with NULLs
    if (length(ans))
      label.cols.y(ans, NOS = NOS, colnames.y =
      switch(type.fitted,
             "lambdas"   = c("lambda.p", "lambda.a",  # Some colns NA
                 "lambda.i", "lambda.d")[(tmp3.TF[c(1, 3, 5, 7)])],
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
    .llambda.p = llambda.p, .elambda.p = elambda.p,
    .llambda.a = llambda.a, .elambda.a = elambda.a,
    .llambda.i = llambda.i, .elambda.i = elambda.i,
    .llambda.d = llambda.d, .elambda.d = elambda.d,
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
                earg = list()))  # This isnt perfect; info is lost
    misc$predictors.names <- predictors.names  # Useful for coef()
    misc$link <- link.names  # 
    names(misc$link) <- parameter.names  # 


    misc$earg <- vector("list", M1)
    names(misc$earg) <- names(misc$link)
    misc$earg[[1]] <- ( .elambda.p )  # First one always there
    iptr <- 1
    if (tmp3.TF[ 2])
      misc$earg[[(iptr <- iptr + 1)]] <- list()  # multilogitlink
    if (tmp3.TF[ 3])
      misc$earg[[(iptr <- iptr + 1)]] <- ( .elambda.a )
    if (tmp3.TF[ 4])
      misc$earg[[(iptr <- iptr + 1)]] <- list()  # See below
    if (tmp3.TF[ 5])
      misc$earg[[(iptr <- iptr + 1)]] <- ( .elambda.i )
    if (tmp3.TF[ 6])
      misc$earg[[(iptr <- iptr + 1)]] <- list()  # See below
    if (tmp3.TF[ 7])
      misc$earg[[(iptr <- iptr + 1)]] <- ( .elambda.d )
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
    .llambda.p = llambda.p, .elambda.p = elambda.p,
    .llambda.a = llambda.a, .elambda.a = elambda.a,
    .llambda.i = llambda.i, .elambda.i = elambda.i,
    .llambda.d = llambda.d, .elambda.d = elambda.d,
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
    lambda.p <- cbind(eta2theta(eta[, 1], .llambda.p , .elambda.p ))
    ind.lambda.z <- 1  # Points to lambda.p only.
    lambda.a <- lambda.i <-
    lambda.d <- lambda.p  # Needed and doesnt corrupt the answer

    if (any(tmp3.TF[c(3, 5, 7)])) {  # At least one lambda.[aid]
      ind.lambda.z <- extra$indeta[c(1, 3, 5, 7), 'launch']  # Vectors
      ind.lambda.z <- c(na.omit(ind.lambda.z))  # At least one value

      lambda.a <- if (!tmp3.TF[ 3]) lambda.p else
        eta2theta(eta[, extra$indeta[3, 1]], .llambda.a , .elambda.a )
      lambda.i <- if (!tmp3.TF[ 5]) lambda.p else
        eta2theta(eta[, extra$indeta[5, 1]], .llambda.i , .elambda.i )
      lambda.d <- if (!tmp3.TF[ 7]) lambda.p else
        eta2theta(eta[, extra$indeta[7, 1]], .llambda.d , .elambda.d )
    }  # la.mix + li.mix + ld.mix > 0

    if (lall.len) {  # An MLM was fitted
      allprobs <-
        multilogitlink(eta[, -ind.lambda.z, drop = FALSE],
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
        dimnames(pobs.mlm) <- list(rownames(eta),
                                   as.character(a.mlm))
        Nextone <- Nextone + la.mlm
      }
      if (tmp3.TF[ 9]) {
        ind9 <- (Nextone + 1):(Nextone + li.mlm)
        pstr.mlm <- allprobs[, ind9, drop = FALSE]
        dimnames(pstr.mlm) <- list(rownames(eta),
                                   as.character(i.mlm))
        Nextone <- Nextone + li.mlm
      }
      if (tmp3.TF[10]) {
        ind10 <- (Nextone + 1):(Nextone + ld.mlm)
        pdip.mlm <- allprobs[, ind10, drop = FALSE]
        dimnames(pdip.mlm) <- list(rownames(eta),
                                   as.character(d.mlm))
        Nextone <- Nextone + ld.mlm  # Not needed
      }
    }  # lall.len


    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) *
        dgaitdpois(y, lambda.p, log = TRUE,  # byrow.aid = F,
                  a.mix = a.mix, i.mix = i.mix, d.mix = d.mix,
                  a.mlm = a.mlm, i.mlm = i.mlm, d.mlm = d.mlm,
                  truncate = truncate,
                  max.support = as.vector( .max.support ),
                  lambda.a = lambda.a, lambda.i = lambda.i, 
                  lambda.d = lambda.d,
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
    .llambda.p = llambda.p, .elambda.p = elambda.p,
    .llambda.a = llambda.a, .elambda.a = elambda.a,
    .llambda.i = llambda.i, .elambda.i = elambda.i,
    .llambda.d = llambda.d, .elambda.d = elambda.d,
    .lpstr.mix = lpstr.mix, .lpobs.mix = lpobs.mix,
    .lpdip.mix = lpdip.mix,
    .epstr.mix = epstr.mix, .epobs.mix = epobs.mix,
    .epdip.mix = epdip.mix,
    .a.mix = a.mix, .i.mix = i.mix, .d.mix = d.mix,
    .a.mlm = a.mlm, .i.mlm = i.mlm, .d.mlm = d.mlm,
    .truncate = truncate, .max.support = max.support ))),
  vfamily = c("gaitdpoisson"),
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
    lambda.a <- lambda.i <- lambda.d <- 1  # Needed

    if (!is.matrix(eta)) eta <- as.matrix(eta)
    lambda.p <-
      cbind(eta2theta(eta[, 1], .llambda.p , earg = .elambda.p ))
    ind.lambda.z <- 1  # Points to lambda.p only.
    if (any(tmp3.TF[c(3, 5, 7)])) {  # At least one lambda.[aid]
      ind.lambda.z <- extra$indeta[c(1, 3, 5, 7), 1]  # Vectors
      ind.lambda.z <- c(na.omit(ind.lambda.z))  # At least one value

      lambda.a <- if (!tmp3.TF[ 3]) lambda.p else
        eta2theta(eta[, extra$indeta[3, 1]], .llambda.a , .elambda.a )
      lambda.i <- if (!tmp3.TF[ 5]) lambda.p else
        eta2theta(eta[, extra$indeta[5, 1]], .llambda.i , .elambda.i )
      lambda.d <- if (!tmp3.TF[ 7]) lambda.p else
        eta2theta(eta[, extra$indeta[7, 1]], .llambda.d , .elambda.d )
    }  # la.mix + li.mix + ld.mix > 0

    
    if (lall.len) {  # A MLM was fitted
      allprobs <-
        multilogitlink(eta[, -ind.lambda.z, drop = FALSE],
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
        dimnames(pobs.mlm) <- list(rownames(eta),
                                   as.character(a.mlm))
        Nextone <- Nextone + la.mlm
      }
      if (tmp3.TF[ 9]) {
        ind9 <- (Nextone + 1):(Nextone + li.mlm)
        pstr.mlm <- allprobs[, ind9, drop = FALSE]
        dimnames(pstr.mlm) <- list(rownames(eta),
                                   as.character(i.mlm))
        Nextone <- Nextone + li.mlm
      }
      if (tmp3.TF[10]) {
        ind10 <- (Nextone + 1):(Nextone + ld.mlm)
        pdip.mlm <- allprobs[, ind10, drop = FALSE]
        dimnames(pdip.mlm) <- list(rownames(eta),
                                   as.character(d.mlm))
        Nextone <- Nextone + ld.mlm  # Not needed
      }
    }  # lall.len

    okay.mlm <-
      all(is.finite(pobs.mlm)) && all(0 < pobs.mlm) &&
      all(is.finite(pstr.mlm)) && all(0 < pstr.mlm) &&
      all(is.finite(pdip.mlm)) && all(0 < pdip.mlm)
    okay.mix <-
      all(is.finite(lambda.p)) && all(0 < lambda.p) &&
      all(lambda.p < .max.support ) &&
      all(is.finite(lambda.a)) && all(0 < lambda.a) &&
      all(is.finite(lambda.i)) && all(0 < lambda.i) &&
      all(is.finite(lambda.d)) && all(0 < lambda.d) &&
      all(is.finite(pobs.mix)) && all(0 < pobs.mix) &&
      all(is.finite(pstr.mix)) && all(0 < pstr.mix) &&
      all(is.finite(pdip.mix)) && all(0 < pdip.mix) &&
      all(pobs.mix + pstr.mix + pdip.mix +
          rowSums(pobs.mlm) + rowSums(pstr.mlm) +
          rowSums(pdip.mlm) < 1)  # Combined
    okay.mlm && okay.mix
  }, list(
    .llambda.p = llambda.p, .elambda.p = elambda.p,
    .llambda.a = llambda.a, .elambda.a = elambda.a,
    .llambda.i = llambda.i, .elambda.i = elambda.i,
    .llambda.d = llambda.d, .elambda.d = elambda.d,
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
    lambda.p <- cbind(eta2theta(eta[, 1], .llambda.p , .elambda.p ))
    ind.lambda.z <- 1  # Points to lambda.p only.
    lambda.a <- lambda.i <-
    lambda.d <- lambda.p  # Needed; and answer not corrupted
    tmp3.TF <- ( .tmp3.TF )


    if (any(tmp3.TF[c(3, 5, 7)])) {  # At least one lambda.[aid]
      ind.lambda.z <- extra$indeta[c(1, 3, 5, 7), 'launch']  # Vectors
      ind.lambda.z <- c(na.omit(ind.lambda.z))  # At least one value

      lambda.a <- if (!tmp3.TF[ 3]) lambda.p else
        eta2theta(eta[, extra$indeta[3, 1]], .llambda.a , .elambda.a )
      lambda.i <- if (!tmp3.TF[ 5]) lambda.p else
        eta2theta(eta[, extra$indeta[5, 1]], .llambda.i , .elambda.i )
      lambda.d <- if (!tmp3.TF[ 7]) lambda.p else
        eta2theta(eta[, extra$indeta[7, 1]], .llambda.d , .elambda.d )
    }  # la.mix + li.mix + ld.mix > 0


    if (lall.len) {  # A AMLM was fitted
      allprobs <-
        multilogitlink(eta[, -ind.lambda.z, drop = FALSE],
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
        dimnames(pobs.mlm) <- list(rownames(eta),
                                   as.character(a.mlm))
        Nextone <- Nextone + la.mlm
      }
      if (tmp3.TF[ 9]) {
        ind9 <- (Nextone + 1):(Nextone + li.mlm)
        pstr.mlm <- allprobs[, ind9, drop = FALSE]
        dimnames(pstr.mlm) <- list(rownames(eta),
                                   as.character(i.mlm))
        Nextone <- Nextone + li.mlm
      }
      if (tmp3.TF[10]) {
        ind10 <- (Nextone + 1):(Nextone + ld.mlm)
        pdip.mlm <- allprobs[, ind10, drop = FALSE]
        dimnames(pdip.mlm) <- list(rownames(eta),
                                   as.character(d.mlm))
        Nextone <- Nextone + ld.mlm  # Not needed
      }
    }  # lall.len

    rgaitdpois(nsim * length(lambda.p), lambda.p,
               pobs.mlm = pobs.mlm, pstr.mlm = pstr.mlm,
               pobs.mix = pobs.mix, pstr.mix = pstr.mix,
               pdip.mix = pdip.mix, pdip.mlm = pdip.mlm,
               lambda.a = lambda.a, lambda.i = lambda.i, 
               lambda.d = lambda.d,
               a.mix = a.mix, i.mix = i.mix, d.mix = d.mix,
               a.mlm = a.mlm, i.mlm = i.mlm, d.mlm = d.mlm,
               truncate = .truncate , max.support = .max.support )
  }, list(
    .llambda.p = llambda.p, .elambda.p = elambda.p,
    .llambda.a = llambda.a, .elambda.a = elambda.a,
    .llambda.i = llambda.i, .elambda.i = elambda.i,
    .llambda.d = llambda.d, .elambda.d = elambda.d,
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
    lambda.p <- cbind(eta2theta(eta[, 1], .llambda.p , .elambda.p ))
    ind.lambda.z <- 1  # Points to lambda.p only.
    lambda.a <- lambda.i <-
    lambda.d <- lambda.p  # Needed; and answer not corrupted

    if (any(tmp3.TF[c(3, 5, 7)])) {  # At least one lambda.[aid]
      ind.lambda.z <- extra$indeta[c(1, 3, 5, 7), 'launch']  # Vectors
      ind.lambda.z <- c(na.omit(ind.lambda.z))  # At least one value

      lambda.a <- if (!tmp3.TF[ 3]) lambda.p else
        eta2theta(eta[, extra$indeta[3, 1]], .llambda.a , .elambda.a )
      lambda.i <- if (!tmp3.TF[ 5]) lambda.p else
        eta2theta(eta[, extra$indeta[5, 1]], .llambda.i , .elambda.i )
      lambda.d <- if (!tmp3.TF[ 7]) lambda.p else
        eta2theta(eta[, extra$indeta[7, 1]], .llambda.d , .elambda.d )
    }  # la.mix + li.mix + ld.mix > 0



    if (lall.len) {  # A MLM was fitted.
      allprobs <-
        multilogitlink(eta[, -ind.lambda.z, drop = FALSE],
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
        dimnames(pobs.mlm) <- list(rownames(eta),
                                   as.character(a.mlm))
        Nextone <- Nextone + la.mlm
      }
      if (tmp3.TF[ 9]) {
        ind9 <- (Nextone + 1):(Nextone + li.mlm)
        pstr.mlm <- allprobs[, ind9, drop = FALSE]
        dimnames(pstr.mlm) <- list(rownames(eta),
                                   as.character(i.mlm))
        Nextone <- Nextone + li.mlm
      }
      if (tmp3.TF[10]) {
        ind10 <- (Nextone + 1):(Nextone + ld.mlm)
        pdip.mlm <- allprobs[, ind10, drop = FALSE]
        dimnames(pdip.mlm) <- list(rownames(eta),
                                   as.character(d.mlm))
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

    dl.dlambda.p <- y / lambda.p - 1  # == dl.dlambda.p.usual
    dl.dlambda.p[!is.ns] <- 0  # For is.a.mixed & is.a.mlmed
    prob.mlm.a <- if (la.mlm) rowSums(pobs.mlm) else 0  # scalar okay
    prob.mlm.i <- if (li.mlm) rowSums(pstr.mlm) else 0  # scalar okay
    prob.mlm.d <- if (ld.mlm) rowSums(pdip.mlm) else 0  # scalar okay
    pmf.deriv1 <- function(y, lambda)
      dpois(y-1, lambda) - dpois(y, lambda)
    pmf.deriv2 <- function(y, lambda)
      dpois(y-2, lambda) - 2 * dpois(y-1, lambda) + dpois(y, lambda)


    sumD.mix.1a.p <- sumD.mix.2a.p <- matrix(0, n, NOS)
    if (la.mix > 0) {  # \calA_p
      DA.mix.0mat.a <-  # Matches naming convention further below
      DA.mix.1mat.a <- matrix(0, n, la.mix)
      for (jay in seq(la.mix)) {
        aval <- a.mix[jay]
        sumD.mix.1a.p <- sumD.mix.1a.p + pmf.deriv1(aval, lambda.p)
        sumD.mix.2a.p <- sumD.mix.2a.p + pmf.deriv2(aval, lambda.p)
        pmf.a <- dpois(aval, lambda.a)
        DA.mix.0mat.a[, jay] <- pmf.a
        DA.mix.1mat.a[, jay] <- pmf.deriv1(aval, lambda.a)
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
          pmf.i <- dpois(ival, lambda.i)
          DI.mix.0mat.i[, jay] <- pmf.i
          DI.mix.1mat.i[, jay] <- pmf.deriv1(ival, lambda.i)
          DI.mix.2mat.i[, jay] <- pmf.deriv2(ival, lambda.i)
          pmf.p <- dpois(ival, lambda.p)
          DP.mix.0mat.i[, jay] <- pmf.p
          DP.mix.1mat.i[, jay] <- pmf.deriv1(ival, lambda.p)
          DP.mix.2mat.i[, jay] <- pmf.deriv2(ival, lambda.p)
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
          pmf.d <- dpois(dval, lambda.d)
          DD.mix.0mat.d[, jay] <- pmf.d
          DD.mix.1mat.d[, jay] <- pmf.deriv1(dval, lambda.d)
          DD.mix.2mat.d[, jay] <- pmf.deriv2(dval, lambda.d)
          pmf.p <- dpois(dval, lambda.p)
          DP.mix.0mat.d[, jay] <- pmf.p
          DP.mix.1mat.d[, jay] <- pmf.deriv1(dval, lambda.p)
          DP.mix.2mat.d[, jay] <- pmf.deriv2(dval, lambda.p)
        }  # jay
      Denom1.d <- rowSums(DD.mix.1mat.d)
      Denom2.d <- rowSums(DD.mix.2mat.d)
    }  # ld.mix

    Bits <- moments.gaitdcombo.pois(lambda.p,
              pobs.mix = pobs.mix, pobs.mlm = pobs.mlm,
              pstr.mix = pstr.mix, pstr.mlm = pstr.mlm,
              pdip.mix = pdip.mix, pdip.mlm = pdip.mlm,
              a.mix = a.mix, i.mix = i.mix, d.mix = d.mix,
              a.mlm = a.mlm, i.mlm = i.mlm, d.mlm = d.mlm,
              lambda.a = lambda.a, lambda.i = lambda.i,
              lambda.d = lambda.d,
              truncate = truncate, max.support = max.support)


    sumD.mlm.1a.p <- sumD.mlm.2a.p <- matrix(0, n, NOS)
    if (la.mlm)
      for (aval in a.mlm) {
        sumD.mlm.1a.p <- sumD.mlm.1a.p + pmf.deriv1(aval, lambda.p)
        sumD.mlm.2a.p <- sumD.mlm.2a.p + pmf.deriv2(aval, lambda.p)
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
        pmf.p <- dpois(ival, lambda.p)
        Dp.mlm.0Mat.i[, jay] <- pmf.p
        Dp.mlm.1Mat.i[, jay] <- pmf.deriv1(ival, lambda.p)
        Dp.mlm.2Mat.i[, jay] <- pmf.deriv2(ival, lambda.p)
      }  # jay
    }  # li.mlm



    Dp.mlm.0Mat.d <-  # wrt parent distribution
    Dp.mlm.1Mat.d <- Dp.mlm.2Mat.d <- matrix(0, n, NOS)
    if (ld.mlm > 0) {
      Dp.mlm.0Mat.d <-  # wrt parent distribution
      Dp.mlm.1Mat.d <- Dp.mlm.2Mat.d <- matrix(0, n, ld.mlm)
      for (jay in seq(ld.mlm)) {
        dval <- d.mlm[jay]
        pmf.p <- dpois(dval, lambda.p)
        Dp.mlm.0Mat.d[, jay] <- pmf.p
        Dp.mlm.1Mat.d[, jay] <- pmf.deriv1(dval, lambda.p)
        Dp.mlm.2Mat.d[, jay] <- pmf.deriv2(dval, lambda.p)
      }  # jay
    }  # ld.mlm





    

    sumD.1t.p <- sumD.2t.p <-
    sumD.1t.a <- sumD.2t.a <-
    sumD.1t.i <- sumD.2t.i <-
    sumD.1t.d <- sumD.2t.d <- matrix(0, n, NOS)
    if (ltruncat)
      for (tval in truncate) {
        sumD.1t.p <- sumD.1t.p + pmf.deriv1(tval, lambda.p)
        sumD.2t.p <- sumD.2t.p + pmf.deriv2(tval, lambda.p)
        sumD.1t.a <- sumD.1t.a + pmf.deriv1(tval, lambda.a)
        sumD.2t.a <- sumD.2t.a + pmf.deriv2(tval, lambda.a)
        sumD.1t.i <- sumD.1t.i + pmf.deriv1(tval, lambda.i)
        sumD.2t.i <- sumD.2t.i + pmf.deriv2(tval, lambda.i)
        sumD.1t.d <- sumD.1t.d + pmf.deriv1(tval, lambda.d)
        sumD.2t.d <- sumD.2t.d + pmf.deriv2(tval, lambda.d)
      }
    sumD.1t.p <- sumD.1t.p + dpois(max.support  , lambda.p)
    sumD.2t.p <- sumD.2t.p + dpois(max.support-1, lambda.p) -
                             dpois(max.support  , lambda.p)
    sumD.1t.a <- sumD.1t.a + dpois(max.support  , lambda.a)
    sumD.2t.a <- sumD.2t.a + dpois(max.support-1, lambda.a) -
                             dpois(max.support  , lambda.a)
    sumD.1t.i <- sumD.1t.i + dpois(max.support  , lambda.i)
    sumD.2t.i <- sumD.2t.i + dpois(max.support-1, lambda.i) -
                             dpois(max.support  , lambda.i)
    sumD.1t.d <- sumD.1t.d + dpois(max.support  , lambda.d)
    sumD.2t.d <- sumD.2t.d + dpois(max.support-1, lambda.d) -
                             dpois(max.support  , lambda.d)


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




    dl.dlambda.a <- dl.dlambda.i <- dl.dlambda.d <- numeric(n)
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



    dl.dlambda.p[is.ns] <- dl.dlambda.p[is.ns] -
                           (Denom1.p / Denom0.p)[is.ns]


    
    if (tmp3.TF[ 9] && li.mlm > 0) {  # aka \calI_{np}
      dl.dpstr.mlm <- matrix(-1 / Numer, n, li.mlm)
      dl.dpstr.mlm[!is.ns, ] <- 0  # For a.mlm and a.mix

      for (jay in seq(li.mlm)) {
        is.inf.j.mlm <- extra$skip.mlm.i[, jay]  # Logical vector
        tmp7i <- Numer * d1B.PI.mlm[, jay] / DELTA.i.mlm[, jay]
        dl.dlambda.p[is.inf.j.mlm] <- tmp7i[is.inf.j.mlm]


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
        dl.dlambda.p[is.def.j.mlm] <- tmp7d[is.def.j.mlm]  # 20211020

 
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
          dl.dlambda.a[is.alt.j.mix] <- tmp2[is.alt.j.mix]  # ccc.
        }  # jay
    }  # la.mix




    if (tmp3.TF[ 4] && li.mix > 0) {  # aka \calI_{p}
      for (jay in seq(li.mix)) {
        ival <- i.mix[jay]
        is.inf.j.mix <- extra$skip.mix.i[, jay]  # Logical vector
        tmp7b <- Numer * d1B.PI.mix[, jay] / DELTA.i.mix[, jay]
        dl.dlambda.p[is.inf.j.mix] <- tmp7b[is.inf.j.mix]
        tmp8 <- (d0A.i[, jay] - d0B.PI.mix[, jay]) / DELTA.i.mix[, jay]
        dl.dpstr.mix[is.inf.j.mix] <- tmp8[is.inf.j.mix]
        if (li.mix > 1) {
          tmp2 <- pstr.mix * d1A.i[, jay] / DELTA.i.mix[, jay]
          dl.dlambda.i[is.inf.j.mix] <- tmp2[is.inf.j.mix]
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
        dl.dlambda.p[is.def.j.mix] <- tmp7b[is.def.j.mix]
        tmp8 <- (d0B.PD.mix[, jay] - d0A.d[, jay]) / DELTA.d.mix[, jay]
        dl.dpdip.mix[is.def.j.mix] <- tmp8[is.def.j.mix]

        if (ld.mix > 1) {
          tmp2 <- (-pdip.mix) * d1A.d[, jay] / DELTA.d.mix[, jay]
          dl.dlambda.d[is.def.j.mix] <- tmp2[is.def.j.mix]
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
      new.ansd[, -ind.lambda.z] <- allprobs[, -ncol(allprobs)] *
                                   (all6.dldp - rSs.tmp)
    }  # lall.len


    

    dlambda.p.deta <- dtheta.deta(lambda.p, .llambda.p , .elambda.p )
    if (tmp3.TF[ 3])
      dlambda.a.deta <- dtheta.deta(lambda.a, .llambda.a , .elambda.a )
    if (tmp3.TF[ 5])
      dlambda.i.deta <- dtheta.deta(lambda.i, .llambda.i , .elambda.i )
    if (tmp3.TF[ 7])
      dlambda.d.deta <- dtheta.deta(lambda.d, .llambda.d , .elambda.d )
    new.ansd[, 1] <- dl.dlambda.p * dlambda.p.deta
    if (tmp3.TF[ 3])
      new.ansd[, extra$indeta[3, 1]] <- dl.dlambda.a * dlambda.a.deta
    if (tmp3.TF[ 5])
      new.ansd[, extra$indeta[5, 1]] <- dl.dlambda.i * dlambda.i.deta
    if (tmp3.TF[ 7])
      new.ansd[, extra$indeta[7, 1]] <- dl.dlambda.d * dlambda.d.deta
    onecoln.indeta <- extra$indeta[1:7, ]  # One coln params only
    onecoln.indeta <- na.omit(onecoln.indeta)  # Only those present
    allcnames <- c(rownames(onecoln.indeta),
                   as.character(c(a.mlm, i.mlm, d.mlm)))
    colnames(new.ansd) <- allcnames

 


    c(w) * new.ansd
  }), list(
    .llambda.p = llambda.p, .elambda.p = elambda.p,
    .llambda.a = llambda.a, .elambda.a = elambda.a,
    .llambda.i = llambda.i, .elambda.i = elambda.i,
    .llambda.d = llambda.d, .elambda.d = elambda.d,
    .lpstr.mix = lpstr.mix, .lpobs.mix = lpobs.mix,
    .lpdip.mix = lpdip.mix,
    .epstr.mix = epstr.mix, .epobs.mix = epobs.mix,
    .epdip.mix = epdip.mix,
    .a.mix = a.mix, .i.mix = i.mix, .d.mix = d.mix,
    .a.mlm = a.mlm, .i.mlm = i.mlm, .d.mlm = d.mlm,
    .tmp3.TF = tmp3.TF,  # .tmp3 = tmp3,
    .truncate = truncate, .max.support = max.support ))),

  weight = eval(substitute(expression({  # gaitdpoisson



    wz <- matrix(0, n, M * (M + 1) / 2)  # The complete size
    cond.EY.p <-
      c(lambda.p - Bits[["SumT1.p"]] -
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
    ned2l.dpobs.mix.lambda.p <- zero0n  # mB overwritten below [4279]
    ned2l.dpobs.mix.lambda.a <- zero0n  # Fini; (2, 3) element
    ned2l.dpobs.mix.lambda.i <- zero0n  # mB overwritten below
    ned2l.dpobs.mix.lambda.d <- zero0n  # mB overwritten below
    ned2l.dpstr.mix.lambda.p <- zero0n  # Optional (1, 4) element
    ned2l.dpstr.mix.lambda.a <- zero0n  # Final; nothing to do
    ned2l.dpstr.mix.lambda.i <- zero0n  # mB overwritten below
    ned2l.dpstr.mix.lambda.d <- zero0n  # mB overwritten below
    ned2l.dpdip.mix.lambda.p <- zero0n  # Optional (1, 6) element




      posn.pobs.mix <- as.vector(extra$indeta[ 2, 'launch'])
      posn.lambda.a <- as.vector(extra$indeta[ 3, 'launch'])
      posn.pstr.mix <- as.vector(extra$indeta[ 4, 'launch'])
      posn.lambda.i <- as.vector(extra$indeta[ 5, 'launch'])
      posn.pdip.mix <- as.vector(extra$indeta[ 6, 'launch'])
      posn.lambda.d <- as.vector(extra$indeta[ 7, 'launch'])
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


    


    ned2l.dlambda.p2 <- probns * (cond.EY.p / lambda.p^2 +  # ccc
    Denom2.p / Denom0.p - (Denom1.p / Denom0.p)^2) + 
    (if (tmp3.TF[ 4] && li.mix) Numer *
    rowSums(Numer * (d1B.PI.mix^2) / DELTA.i.mix - d2B.PI.mix) else 0) +
    (if (tmp3.TF[ 9] && li.mlm) Numer *
    rowSums(Numer * (d1B.PI.mlm^2) / DELTA.i.mlm - d2B.PI.mlm) else 0) +
    (if (tmp3.TF[ 6] && ld.mix) Numer *
    rowSums(Numer * (d1B.PD.mix^2) / DELTA.d.mix - d2B.PD.mix) else 0) +
    (if (tmp3.TF[10] && ld.mlm) Numer *  # nnn.
    rowSums(Numer * (d1B.PD.mlm^2) / DELTA.d.mlm - d2B.PD.mlm) else 0)


    wz[, iam(1, 1, M)] <- ned2l.dlambda.p2 * dlambda.p.deta^2


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
      ned2l.dlambda.a2 <- pobs.mix * (
        rowSums((DA.mix.1mat.a^2) / DA.mix.0mat.a) / Denom0.a -
        (Denom1.a / Denom0.a)^2)  # ccc.
      wz[, iam(3, 3, M)] <- ned2l.dlambda.a2 * dlambda.a.deta^2
    }




    if (tmp3.TF[ 4] && li.mix > 0) {
      ned2l.dpstr.mix2 <-
      ned2l.dpstr.mix2 +
        rowSums((d0A.i - d0B.PI.mix)^2 / DELTA.i.mix)

      if (tmp3.TF[ 2] && la.mix > 0)
        ned2l.dpobs.mix.lambda.p <-
        ned2l.dpobs.mix.lambda.p +
          rowSums(d1B.PI.mix * (1 - Numer * d0B.PI.mix / DELTA.i.mix))

      ned2l.dpstr.mix.lambda.p <-
      ned2l.dpstr.mix.lambda.p + rowSums(
        d1B.PI.mix * (1 + Numer * (d0A.i - d0B.PI.mix) / DELTA.i.mix))

      if (tmp3.TF[ 6])
        ned2l.dpdip.mix.lambda.p <-
        ned2l.dpdip.mix.lambda.p - rowSums(
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

      ned2l.dlambda.p.lambda.i <- pstr.mix * Numer *
        rowSums(d1A.i * d1B.PI.mix / DELTA.i.mix)  # ccc.
      wz[, iam(1, posn.lambda.i, M)] <- ned2l.dlambda.p.lambda.i *
        dlambda.p.deta * dlambda.i.deta  # All links done here

      ned2l.dlambda.i2 <- pstr.mix *
        rowSums(pstr.mix * (d1A.i^2) / DELTA.i.mix - d2A.i)  # ccc.
      wz[, iam(posn.lambda.i, posn.lambda.i, M)] <-
        ned2l.dlambda.i2 * dlambda.i.deta^2

      if (tmp3.TF[ 2]) {  # tmp3.TF[ 4] is TRUE, given tmp3.TF[ 5]
        ned2l.dpobs.mix.lambda.i <-
          rowSums(-pstr.mix * d1A.i * d0B.PI.mix / DELTA.i.mix)  # ccc.
        wz[, iam(posn.pobs.mix, posn.lambda.i, M)] <-
          ned2l.dpobs.mix.lambda.i  # * dlambda.i.deta done later
      }

      if (tmp3.TF[ 4]) {
        ned2l.dpstr.mix.lambda.i <- rowSums(  # ccc.
          d1A.i * (pstr.mix * (d0A.i - d0B.PI.mix) / DELTA.i.mix - 1))
        wz[, iam(posn.pstr.mix, posn.lambda.i, M)] <-
          ned2l.dpstr.mix.lambda.i  # * dlambda.i.deta done later
      }

      if (all(tmp3.TF[c(5, 6)])) {
        ned2l.dpdip.mix.lambda.i <- rowSums(
          (-pstr.mix) * d0B.PI.mix * d1A.i / DELTA.i.mix)
        wz[, iam(posn.pdip.mix, posn.lambda.i, M)] <-
          ned2l.dpdip.mix.lambda.i  # link done later
      }

      if (tmp3.TF[ 8]) {
        ned2l.dpobs.mlm.lambda.i <- rowSums(
          -pstr.mix * d0B.PI.mix * d1A.i / DELTA.i.mix)  # ccc.
        for (uuu in seq(la.mlm))
          wz[, iam(posn.pobs.mlm - 1 + uuu, posn.lambda.i, M)] <-
            ned2l.dpobs.mlm.lambda.i  # * dlambda.i.deta done later
      }


    }  # (tmp3.TF[ 5] && li.mix > 1)





    if (tmp3.TF[ 6] && ld.mix > 0) {  # \calD_{p}, maybe w. \theta_d

      if (tmp3.TF[ 2] && la.mix > 0)
        ned2l.dpobs.mix.lambda.p <-
        ned2l.dpobs.mix.lambda.p +
          rowSums(d1B.PD.mix * (1 - Numer * d0B.PD.mix / DELTA.d.mix))

      ned2l.dpstr.mix.lambda.p <-
      ned2l.dpstr.mix.lambda.p + rowSums(
        d1B.PD.mix * (1 - Numer * d0B.PD.mix / DELTA.d.mix))

      ned2l.dpdip.mix.lambda.p <-
      ned2l.dpdip.mix.lambda.p - rowSums(
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


      ned2l.dlambda.p.lambda.d <- (-pdip.mix) * Numer *
        rowSums(d1A.d * d1B.PD.mix / DELTA.d.mix)  # nnn.
      wz[, iam(1, posn.lambda.d, M)] <- ned2l.dlambda.p.lambda.d *
        dlambda.p.deta * dlambda.d.deta  # All links done here

      if (tmp3.TF[ 2]) {  # tmp3.TF[ 6] is TRUE, given tmp3.TF[ 7]
        ned2l.dpobs.mix.lambda.d <-
          rowSums(pdip.mix * d1A.d * d0B.PD.mix / DELTA.d.mix)  # nnn.
        wz[, iam(posn.pobs.mix, posn.lambda.d, M)] <-
          ned2l.dpobs.mix.lambda.d  # link done later
      }

      if (tmp3.TF[ 4]) {
        ned2l.dpstr.mix.lambda.d <- rowSums(
          pdip.mix * d1A.d * d0B.PD.mix / DELTA.d.mix)
        wz[, iam(posn.pstr.mix, posn.lambda.d, M)] <-
          ned2l.dpstr.mix.lambda.d  # * dlambda.i.deta done later
      }

        ned2l.dpdip.mix.lambda.d <- rowSums(
          d1A.d * (1 + pdip.mix * (d0A.d - d0B.PD.mix) / DELTA.d.mix))
        wz[, iam(posn.pdip.mix, posn.lambda.d, M)] <-
          ned2l.dpdip.mix.lambda.d  # * dlambda.d.deta done later

      ned2l.dlambda.d2 <- pdip.mix *
        rowSums(pdip.mix * (d1A.d^2) / DELTA.d.mix + d2A.d)  # nnn.
      wz[, iam(posn.lambda.d, posn.lambda.d, M)] <-
        ned2l.dlambda.d2 * dlambda.d.deta^2

      if (tmp3.TF[ 8]) {
        ned2l.dpobs.mlm.lambda.d <- rowSums(
           pdip.mix * d0B.PD.mix * d1A.d / DELTA.d.mix)  # nnn.
        for (uuu in seq(la.mlm))
          wz[, iam(posn.pobs.mlm - 1 + uuu, posn.lambda.d, M)] <-
            ned2l.dpobs.mlm.lambda.d  # * dlambda.d.deta done later
      }

    }  # (tmp3.TF[ 7] && ld.mix > 1)




        
    if (tmp3.TF[ 9] && li.mlm > 0) {  # \calI_{np}, includes \phi_s.

      if (la.mix && tmp3.TF[ 2])
        ned2l.dpobs.mix.lambda.p <-  # ccc
        ned2l.dpobs.mix.lambda.p +
          rowSums(d1B.PI.mlm * (1 - Numer * d0B.PI.mlm / DELTA.i.mlm))

      ned2l.dpstr.mix.lambda.p <-  # ccc.
      ned2l.dpstr.mix.lambda.p + rowSums(
        d1B.PI.mlm * (1 - Numer * d0B.PI.mlm / DELTA.i.mlm))

      if (tmp3.TF[ 6])
        ned2l.dpdip.mix.lambda.p <-
        ned2l.dpdip.mix.lambda.p - rowSums(
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
        ned2l.dpobs.mix.lambda.p <-  # nnn.
        ned2l.dpobs.mix.lambda.p +
          rowSums(d1B.PD.mlm * (1 - Numer * d0B.PD.mlm / DELTA.d.mlm))

      ned2l.dpstr.mix.lambda.p <-  # nnn.
      ned2l.dpstr.mix.lambda.p + rowSums(
        d1B.PD.mlm * (1 - Numer * d0B.PD.mlm / DELTA.d.mlm))

      if (tmp3.TF[ 6])
        ned2l.dpdip.mix.lambda.p <-
        ned2l.dpdip.mix.lambda.p - rowSums(
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
        ned2l.dpobs.mix.lambda.p  # One link done later
    if (!is.na(posn.pstr.mix))  # Optional (1, 4) element
      wz[, iam(1, posn.pstr.mix, M)] <-
        ned2l.dpstr.mix.lambda.p  # One link done later
    if (!is.na(posn.pdip.mix))  # Optional (1, 6) element
      wz[, iam(1, posn.pdip.mix, M)] <-
        ned2l.dpdip.mix.lambda.p  # One link done later
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
      ned2l.dpobs.mlm.lambda.p <- init0.i.val + init0.d.val  # Vector

      if (tmp3.TF[ 4] && li.mix)
        ned2l.dpobs.mlm.lambda.p <-
        ned2l.dpobs.mlm.lambda.p + rowSums(
          d1B.PI.mix * (1 - Numer * d0B.PI.mix / DELTA.i.mix))
      if (tmp3.TF[ 6] && ld.mix)
        ned2l.dpobs.mlm.lambda.p <-
        ned2l.dpobs.mlm.lambda.p + rowSums(  # nnn
          d1B.PD.mix * (1 - Numer * d0B.PD.mix / DELTA.d.mix))

      ofset <- posn.pobs.mlm - 1  # 5 for combo
      for (vvv in seq(la.mlm))  # ccc.
        wz[, iam(1, ofset + vvv, M)] <- ned2l.dpobs.mlm.lambda.p
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
        wz[, iam(posn.lambda.i, posn.pstr.mlm - 1 + vvv, M)] <-
          ned2l.dpstr.mlm.theta.i  # ccc.
    }  # li.mlm && li.mix > 1





    if (ld.mlm && ld.mix > 1) {

      ned2l.dpdip.mlm.theta.d <-  # Not a matrix, just a vector
        rowSums(pdip.mix * d0B.PD.mix * d1A.d / DELTA.d.mix)

      for (vvv in seq(ld.mlm))
        wz[, iam(posn.lambda.d, posn.pdip.mlm - 1 + vvv, M)] <-
          ned2l.dpdip.mlm.theta.d  # nnn.
    }  # ld.mlm && ld.mix > 1



    if (ld.mlm && li.mix > 1) {

      ned2l.dpdip.mlm.theta.i <-  # Not a matrix, just a vector
        rowSums(-pstr.mix * d0B.PI.mix * d1A.i / DELTA.i.mix)

      for (vvv in seq(ld.mlm))
        wz[, iam(posn.lambda.i, posn.pdip.mlm - 1 + vvv, M)] <-
          ned2l.dpdip.mlm.theta.i  # nnn.
    }  # ld.mlm && li.mix > 1



   


    if (li.mlm && ld.mix > 1) {

      ned2l.dpstr.mlm.theta.d <-  # Not a matrix, just a vector
        rowSums(pdip.mix * d0B.PD.mix * d1A.d / DELTA.d.mix)

      for (vvv in seq(li.mlm))
        wz[, iam(posn.lambda.d, posn.pstr.mlm - 1 + vvv, M)] <-
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
      ind.rc <- setdiff(1:M, ind.lambda.z)  # Contiguous rows and
      lind.rc <- length(ind.rc)  # cols of the DAMLM



 # Copy in the thetas values: the looping is overkill.
      for (uuu in ind.lambda.z)
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
            if (!any(kay %in% ind.lambda.z)) {
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











 


      dstar.deta <- cbind(dlambda.p.deta,
                          if (tmp3.TF[ 3]) dlambda.a.deta else NULL,
                          if (tmp3.TF[ 5]) dlambda.i.deta else NULL,
                          if (tmp3.TF[ 7]) dlambda.d.deta else NULL)
      iptr <- 0
      if (length(ind.lambda.z))
      for (uuu in ind.lambda.z) {  # Could delete 3 for lambda.a (orthog)
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
        ind.diags <- setdiff(1:M, ind.lambda.z)  # Exclude thetas
        wz[atiny, ind.diags] <- .Machine$double.eps +
        wz[atiny, ind.diags] * (1 + .Machine$double.eps^0.5)
      }
    }  # lall.len




    c(w) * wz
  }), list( .truncate = truncate ))))
}  # gaitdpoisson











 dgaitdpois <-
  function(x, lambda.p,
           a.mix = NULL,
           a.mlm = NULL,
           i.mix = NULL,
           i.mlm = NULL,
           d.mix = NULL,
           d.mlm = NULL,
           truncate = NULL, max.support = Inf,
           pobs.mix = 0,  # vector
           pobs.mlm = 0,  # matrix
           pstr.mix = 0,  # vector
           pstr.mlm = 0,  # matrix
           pdip.mix = 0,  # vector
           pdip.mlm = 0,  # matrix
           byrow.aid = FALSE,  # Applies to 'pobs.mlm' & 'pstr.mlm'
           lambda.a = lambda.p, lambda.i = lambda.p,
           lambda.d = lambda.p,
           log = FALSE) {
  log.arg <- log;  rm(log)
  lowsup <- 0  # Lower support
  gaitd.errorcheck(a.mix, a.mlm, i.mix, i.mlm,
                   d.mix, d.mlm, truncate, max.support)
  la.mix <- length(a.mix <- sort(a.mix))
  li.mix <- length(i.mix <- sort(i.mix))
  ld.mix <- length(d.mix <- sort(d.mix))
  la.mlm <- length(a.mlm)
  li.mlm <- length(i.mlm)
  ld.mlm <- length(d.mlm)
  ltrunc <- length(truncate <- sort(truncate))
  if (la.mix + la.mlm + li.mix + li.mlm + ld.mix + ld.mlm +
      ltrunc == 0 &&
      is.infinite(max.support))
    return(dpois(x, lambda.p, log = log.arg))


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
    
  if (any(pdip.mix < 0 | 1 <= pdip.mix, na.rm = TRUE)) {
    stop("bad input for argument 'pdip.mix'")
  }
    
  if (any(pdip.mlm < 0 | 1 <= pdip.mlm, na.rm = TRUE)) {
    stop("bad input for argument 'pdip.mlm'")
  }
   
   

  LLL <- max(length(x),  
             length(pobs.mix), length(pstr.mix), length(pdip.mix),
             length(lambda.p), length(lambda.a),
             length(lambda.i), length(lambda.d))
  if (length(x)        < LLL) x        <- rep_len(x,        LLL)
  if (length(lambda.p) < LLL) lambda.p <- rep_len(lambda.p, LLL)
  if (length(lambda.a) < LLL) lambda.a <- rep_len(lambda.a, LLL)
  if (length(lambda.i) < LLL) lambda.i <- rep_len(lambda.i, LLL)
  if (length(lambda.d) < LLL) lambda.d <- rep_len(lambda.d, LLL)
  if (length(pobs.mix) < LLL) pobs.mix <- rep_len(pobs.mix, LLL)
  if (length(pstr.mix) < LLL) pstr.mix <- rep_len(pstr.mix, LLL)
  if (length(pdip.mix) < LLL) pdip.mix <- rep_len(pdip.mix, LLL)



  sumt <- 0  # Initialization to 0 important
  if (ltrunc)
    for (tval in truncate)
      sumt <- sumt + dpois(tval, lambda.p)  # Need tval <= max.support
  vecTF.t <- is.finite(x) & ((x %in% truncate) | (max.support < x))
  cdf.max.s <- ppois(max.support, lambda.p)  # Usually 1
  denom.t <- cdf.max.s - sumt  # No sumt on RHS

    pmf0 <- ifelse(vecTF.t, 0, dpois(x, lambda.p) / denom.t)


  sum.a <- suma <- 0  # numeric(LLL)
  vecTF.a <- rep_len(FALSE, LLL)
  if (la.mlm) {
    pobs.mlm <-  matrix(pobs.mlm, LLL, la.mlm, byrow = byrow.aid)
    sum.a <-   .rowSums(pobs.mlm, LLL, la.mlm)
    if (any(1 < sum.a, na.rm = TRUE))
      stop("bad input for argument 'pobs.mlm'")  # zz

    for (aval in a.mlm)
      suma <- suma + dpois(aval, lambda.p)  # Part i

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
    pmf2.a <- dgaitdpois(x, lambda.a,  # Outer distribution---mlm type
                         truncate = setdiff(allx.a, a.mix),
                         max.support = max(a.mix))
    for (aval in a.mix) {
      suma <- suma + dpois(aval, lambda.p)  # Part ii added; cumulative
      vecTF <- is.finite(x) & aval == x
      pmf0[vecTF] <- 0  # added; the true values are assigned below
      vecTF.a <- vecTF.a | vecTF  # Cumulative; added
    }
  }  # la.mix


  if (li.mix) {
    allx.i <- lowsup:max(i.mix)
    pmf2.i <- dgaitdpois(x, lambda.i,  # Outer distn---mlm type
                         truncate = setdiff(allx.i, i.mix),
                         max.support = max(i.mix))
  }  # li.mix




  sum.d <- 0  # numeric(LLL)
  if (ld.mlm) {
    pdip.mlm <- matrix(pdip.mlm, LLL, ld.mlm, byrow = byrow.aid)
    sum.d <-  .rowSums(pdip.mlm, LLL, ld.mlm)
    if (any(1 < sum.d, na.rm = TRUE))
      stop("bad input for argument 'pdip.mlm'")
  }  # ld.mlm


  if (ld.mix) {
    allx.d <- lowsup:max(d.mix)
    pmf2.d <- dgaitdpois(x, lambda.d,  # Outer distn---mlm type
                         truncate = setdiff(allx.d, d.mix),
                         max.support = max(d.mix))
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
    dpois(x, lambda.p) / (cdf.max.s - suma - sumt))[!skip]


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
}  # dgaitdpois









 pgaitdpois <-
  function(q, lambda.p,
           a.mix = NULL,
           a.mlm = NULL,
           i.mix = NULL,
           i.mlm = NULL,
           d.mix = NULL,
           d.mlm = NULL,
           truncate = NULL, max.support = Inf,
           pobs.mix = 0,
           pobs.mlm = 0,
           pstr.mix = 0,
           pstr.mlm = 0,
           pdip.mix = 0,
           pdip.mlm = 0,
           byrow.aid = FALSE,
           lambda.a = lambda.p, lambda.i = lambda.p,
           lambda.d = lambda.p, lower.tail = TRUE) {
  lowsup <- 0
  gaitd.errorcheck(a.mix, a.mlm, i.mix, i.mlm,
                   d.mix, d.mlm, truncate, max.support)
  la.mix <- length(a.mix <- sort(a.mix))
  li.mix <- length(i.mix <- sort(i.mix))
  ld.mix <- length(d.mix <- sort(d.mix))
  la.mlm <- length(a.mlm)
  li.mlm <- length(i.mlm)
  ld.mlm <- length(d.mlm)
  ltrunc <- length(truncate <- sort(truncate))
  if (la.mix + la.mlm + li.mix + li.mlm + ld.mix + ld.mlm +
      ltrunc == 0 &&
      is.infinite(max.support))
    return(ppois(q, lambda.p, lower.tail = lower.tail))  # log.p

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
             length(lambda.p), length(lambda.a),
             length(lambda.i), length(lambda.d))
  offset.a <- offset.i <- offset.d <-
  Offset.a <- Offset.i <- Offset.d <- numeric(LLL)
  if (length(q)        < LLL) q        <- rep_len(q,        LLL)
  if (length(lambda.p) < LLL) lambda.p <- rep_len(lambda.p, LLL)
  if (length(lambda.a) < LLL) lambda.a <- rep_len(lambda.a, LLL)
  if (length(lambda.i) < LLL) lambda.i <- rep_len(lambda.i, LLL)
  if (length(lambda.d) < LLL) lambda.d <- rep_len(lambda.d, LLL)
  if (length(pobs.mix) < LLL) pobs.mix <- rep_len(pobs.mix, LLL)
  if (length(pstr.mix) < LLL) pstr.mix <- rep_len(pstr.mix, LLL)
  if (length(pdip.mix) < LLL) pdip.mix <- rep_len(pdip.mix, LLL)

  sumt <- 0
  fudge.t <- numeric(LLL)
  cdf.max.s <- ppois(max.support, lambda.p)  # Usually 1
  if (ltrunc) {
    for (tval in truncate) {
      pmf.p <- dpois(tval, lambda.p)
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
      pmf.p <- dpois(aval, lambda.p)
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
      pmf.a <- dpois(aval, lambda.a)
      pmf.p <- dpois(aval, lambda.p)
      use.pobs.mix[, jay] <- pmf.a
      suma <- suma + pmf.p  # cumulative; part ii
    }
    use.pobs.mix <- pobs.mix *
                    use.pobs.mix / rowSums(use.pobs.mix)

    for (jay in seq(la.mix)) {
      aval <- a.mix[jay]
      pmf.p <- dpois(aval, lambda.p)
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
      use.pstr.mix[, jay] <- dpois(ival, lambda.i)
    }
    use.pstr.mix <- pstr.mix *
                    use.pstr.mix / rowSums(use.pstr.mix)

    for (jay in seq(li.mix)) {
      ival <- i.mix[jay]
      pmf.p <- dpois(ival, lambda.p)
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
      use.pdip.mix[, jay] <- dpois(dval, lambda.d)
    }
    use.pdip.mix <- pdip.mix *
                    use.pdip.mix / rowSums(use.pdip.mix)

    for (jay in seq(ld.mix)) {
      dval <- d.mix[jay]
      pmf.p <- dpois(dval, lambda.p)
      if (any(vecTF <- (is.finite(q) & dval <= q))) {
        Offset.d[vecTF] <- Offset.d[vecTF] + use.pdip.mix[vecTF, jay]
      }
    }  # jay
  }  # ld.mix





  numer1 <- 1 - sum.i - sum.a - pstr.mix - pobs.mix + sum.d + pdip.mix
  denom1 <- cdf.max.s - sumt - suma
  ans <- numer1 * (ppois(q, lambda.p) - fudge.t - fudge.a) / denom1 +
         offset.a + offset.i - offset.d +
         Offset.a + Offset.i - Offset.d
  ans[max.support <= q] <- 1
  ans[ans < 0] <- 0  # Occasional roundoff error
  if (lower.tail) ans else 1 - ans
}  # pgaitdpois








 qgaitdpois <-
  function(p, lambda.p,
           a.mix = NULL,
           a.mlm = NULL,
           i.mix = NULL,
           i.mlm = NULL,
           d.mix = NULL,
           d.mlm = NULL,
           truncate = NULL, max.support = Inf,
           pobs.mix = 0,
           pobs.mlm = 0,
           pstr.mix = 0,
           pstr.mlm = 0,
           pdip.mix = 0,
           pdip.mlm = 0,
           byrow.aid = FALSE,
           lambda.a = lambda.p, lambda.i = lambda.p,
           lambda.d = lambda.p) {


  lowsup <- 0
  gaitd.errorcheck(a.mix, a.mlm, i.mix, i.mlm,
                   d.mix, d.mlm, truncate, max.support)
  la.mix <- length(a.mix <- sort(a.mix))
  li.mix <- length(i.mix <- sort(i.mix))
  ld.mix <- length(d.mix <- sort(d.mix))
  la.mlm <- length(a.mlm)
  li.mlm <- length(i.mlm)
  ld.mlm <- length(d.mlm)
  ltrunc <- length(truncate <- sort(truncate))
  if (la.mix + la.mlm + li.mix + li.mlm + ld.mix + ld.mlm +
      ltrunc == 0 &&
      is.infinite(max.support))
    return(qpois(p, lambda.p))  # lower.tail = TRUE, log.p = FALSE

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
             length(lambda.p), length(lambda.a), 
             length(lambda.i), length(lambda.d))
  if (length(p)        < LLL) p        <- rep_len(p,        LLL)
  if (length(lambda.p) < LLL) lambda.p <- rep_len(lambda.p, LLL)
  if (length(lambda.a) < LLL) lambda.a <- rep_len(lambda.a, LLL)
  if (length(lambda.i) < LLL) lambda.i <- rep_len(lambda.i, LLL)
  if (length(lambda.d) < LLL) lambda.d <- rep_len(lambda.d, LLL)
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
    min(setdiff(min.support:(ltrunc+5), truncate)) else min.support
  ans <- p + lambda.p

  bad0 <- !is.finite(lambda.p) | lambda.p <= 0 |
          !is.finite(lambda.a) | lambda.a <= 0 |
          !is.finite(lambda.i) | lambda.i <= 0 |
          !is.finite(lambda.d) | lambda.d <= 0
  bad <- bad0 | !is.finite(p) | p <= 0 | 1 <= p

  Lo <- rep_len(min.support.use - 0.5, LLL)
  approx.ans <- Lo  # True at lhs
  Hi <- if (is.finite(max.support))
    rep(max.support + 0.5, LLL) else 2 * Lo + 10.5
  dont.iterate <- bad
  done <- dont.iterate |
    p <= pgaitdpois(Hi, lambda.p,
                   a.mix = a.mix, a.mlm = a.mlm,
                   i.mix = i.mix, i.mlm = i.mlm,
                   d.mix = d.mix, d.mlm = d.mlm,
                   truncate = truncate, max.support = max.support,
                   pstr.mix = pstr.mix, pobs.mix = pobs.mix,
                   pdip.mix = pdip.mix,
                   pstr.mlm = pstr.mlm, pobs.mlm = pobs.mlm,
                   pdip.mlm = pdip.mlm,
                   lambda.a = lambda.a, lambda.i = lambda.i,
                   lambda.d = lambda.d,
                   byrow.aid = FALSE)

  iter <- 0
  max.iter <- round(log2(.Machine$double.xmax)) - 3
  while (!all(done) && iter < max.iter) {
    Lo[!done] <- Hi[!done]
    Hi[!done] <- 2 * Hi[!done] + 10.5  # Bug fixed
    Hi <- pmin(max.support + 0.5, Hi)  # 20190924
    done[!done] <-
      (p[!done] <= pgaitdpois(Hi[!done], lambda.p[!done],
                       a.mix = a.mix, a.mlm = a.mlm,
                       i.mix = i.mix, i.mlm = i.mlm,
                       d.mix = d.mix, d.mlm = d.mlm,
                       truncate = truncate,
                       max.support = max.support,
                       pobs.mix = pobs.mix[!done],
                       pstr.mix = pstr.mix[!done],
                       pdip.mix = pdip.mix[!done],
                       pobs.mlm = pobs.mlm[!done, , drop = FALSE],
                       pstr.mlm = pstr.mlm[!done, , drop = FALSE],
                       pdip.mlm = pdip.mlm[!done, , drop = FALSE],
                       lambda.a = lambda.a[!done],
                       lambda.i = lambda.i[!done],
                       lambda.d = lambda.d[!done],
                       byrow.aid = FALSE))
    iter <- iter + 1
  }

      foo <- function(q, lambda.p,
                      a.mix = NULL, a.mlm = NULL,
                      i.mix = NULL, i.mlm = NULL,
                      d.mix = NULL, d.mlm = NULL,
                      truncate = NULL, max.support = Inf,
                      pobs.mix = 0, pstr.mix = 0, pdip.mix = 0,
                      pobs.mlm = 0, pstr.mlm = 0, pdip.mlm = 0,
                      lambda.a = lambda.p, lambda.i = lambda.p,
                      lambda.d = lambda.p,
                      byrow.aid = FALSE, p)
      pgaitdpois(q, lambda.p = lambda.p,
                       a.mix = a.mix, a.mlm = a.mlm,
                       i.mix = i.mix, i.mlm = i.mlm,
                       d.mix = d.mix, d.mlm = d.mlm,
                       truncate = truncate,
                       max.support = max.support,
                       pobs.mix = pobs.mix,
                       pstr.mix = pstr.mix,
                       pdip.mix = pdip.mix,
                       pobs.mlm = pobs.mlm,
                       pstr.mlm = pstr.mlm,
                       pdip.mlm = pdip.mlm,
                       lambda.a = lambda.a, lambda.i = lambda.i,
                       lambda.d = lambda.d,
                       byrow.aid = FALSE) - p

      lhs <- dont.iterate |
        p <= dgaitdpois(min.support.use, lambda.p = lambda.p,
                       a.mix = a.mix, a.mlm = a.mlm,
                       i.mix = i.mix, i.mlm = i.mlm,
                       d.mix = d.mix, d.mlm = d.mlm,
                       truncate = truncate,
                       max.support = max.support,
                       pobs.mix = pobs.mix,
                       pstr.mix = pstr.mix,
                       pdip.mix = pdip.mix,
                       pobs.mlm = pobs.mlm,
                       pstr.mlm = pstr.mlm,
                       pdip.mlm = pdip.mlm,
                       lambda.a = lambda.a, lambda.i = lambda.i,
                       lambda.d = lambda.d,
                       byrow.aid = FALSE)

  if (any(!lhs)) {
    approx.ans[!lhs] <-
      bisection.basic(foo, Lo[!lhs], Hi[!lhs], tol = 1/16,
                      lambda.p = lambda.p[!lhs],
                      a.mix = a.mix, a.mlm = a.mlm,
                      i.mix = i.mix, i.mlm = i.mlm,
                      d.mix = d.mix, d.mlm = d.mlm,
                      truncate = truncate,
                      max.support = max.support,
                      pstr.mix = pstr.mix[!lhs],
                      pstr.mlm = pstr.mlm[!lhs, , drop = FALSE],
                      pobs.mix = pobs.mix[!lhs],
                      pobs.mlm = pobs.mlm[!lhs, , drop = FALSE],
                      pdip.mix = pdip.mix[!lhs],
                      pdip.mlm = pdip.mlm[!lhs, , drop = FALSE],
                      lambda.a = lambda.a[!lhs],
                      lambda.i = lambda.i[!lhs],
                      lambda.d = lambda.d[!lhs],
                      byrow.aid = FALSE,
                      p = p[!lhs])
    faa <- floor(approx.ans[!lhs])
    tmp <-
      ifelse(pgaitdpois(faa, lambda.p[!lhs],
                       a.mix = a.mix, a.mlm = a.mlm,
                       i.mix = i.mix, i.mlm = i.mlm,
                       d.mix = d.mix, d.mlm = d.mlm,
                       truncate = truncate,
                       max.support = max.support,
                       pstr.mix = pstr.mix[!lhs],
                       pstr.mlm = pstr.mlm[!lhs, , drop = FALSE],
                       pobs.mix = pobs.mix[!lhs],
                       pobs.mlm = pobs.mlm[!lhs, , drop = FALSE],
                       pdip.mix = pdip.mix[!lhs],
                       pdip.mlm = pdip.mlm[!lhs, , drop = FALSE],
                       lambda.a = lambda.a[!lhs],
                       lambda.i = lambda.i[!lhs],
                       lambda.d = lambda.d[!lhs],
                       byrow.aid = FALSE) < p[!lhs] &
             p[!lhs] <= pgaitdpois(faa + 1, lambda.p[!lhs],
                       a.mix = a.mix, a.mlm = a.mlm,
                       i.mix = i.mix, i.mlm = i.mlm,
                       d.mix = d.mix, d.mlm = d.mlm,
                       truncate = truncate,
                       max.support = max.support,
                       pstr.mix = pstr.mix[!lhs],
                       pstr.mlm = pstr.mlm[!lhs, , drop = FALSE],
                       pobs.mix = pobs.mix[!lhs],
                       pobs.mlm = pobs.mlm[!lhs, , drop = FALSE],
                       pdip.mix = pdip.mix[!lhs],
                       pdip.mlm = pdip.mlm[!lhs, , drop = FALSE],
                       lambda.a = lambda.a[!lhs],
                       lambda.i = lambda.i[!lhs],
                       lambda.d = lambda.d[!lhs],
                       byrow.aid = FALSE),
             faa + 1, faa)
    ans[!lhs] <- tmp
  }  # any(!lhs)


  if (ltrunc)
    while (any(vecTF <- !bad & ans %in% truncate))
      ans[vecTF] <- 1 + ans[vecTF]



  vecTF <- !bad0 & !is.na(p) &
           p <= dgaitdpois(min.support.use, lambda.p,
                       a.mix = a.mix, a.mlm = a.mlm,
                       i.mix = i.mix, i.mlm = i.mlm,
                       d.mix = d.mix, d.mlm = d.mlm,
                       truncate = truncate,
                       max.support = max.support,
                       pobs.mix = pobs.mix,
                       pstr.mix = pstr.mix,
                       pdip.mix = pdip.mix,
                       pobs.mlm = pobs.mlm,
                       pstr.mlm = pstr.mlm,
                       pdip.mlm = pdip.mlm,
                       lambda.a = lambda.a, lambda.i = lambda.i,
                       lambda.d = lambda.d,
                       byrow.aid = FALSE)
  ans[vecTF] <- min.support.use

  ans[!bad0 & !is.na(p) & p == 0] <- min.support.use
  ans[!bad0 & !is.na(p) & p == 1] <- max.support  # Inf
  ans[!bad0 & !is.na(p) & p <  0] <- NaN
  ans[!bad0 & !is.na(p) & p >  1] <- NaN
  ans[ bad0] <- NaN
  ans
}  # qgaitdpois





 rgaitdpois <-
  function(n, lambda.p,
           a.mix = NULL,
           a.mlm = NULL,
           i.mix = NULL,
           i.mlm = NULL,
           d.mix = NULL,
           d.mlm = NULL,
           truncate = NULL, max.support = Inf,
           pobs.mix = 0,  # vector
           pobs.mlm = 0,  # matrix
           pstr.mix = 0,  # vector
           pstr.mlm = 0,  # matrix
           pdip.mix = 0,  # vector
           pdip.mlm = 0,  # matrix
           byrow.aid = FALSE,
           lambda.a = lambda.p, lambda.i = lambda.p,
           lambda.d = lambda.p) {
    qgaitdpois(runif(n), lambda.p,
              a.mix = a.mix,
              a.mlm = a.mlm,
              i.mix = i.mix,
              i.mlm = i.mlm,
              d.mix = d.mix,
              d.mlm = d.mlm,
              truncate = truncate, max.support = max.support,
              pobs.mix = pobs.mix,
              pobs.mlm = pobs.mlm,
              pstr.mix = pstr.mix,
              pstr.mlm = pstr.mlm,
              pdip.mix = pdip.mix,
              pdip.mlm = pdip.mlm,
              lambda.a = lambda.a, lambda.i = lambda.i,
              lambda.d = lambda.p,
              byrow.aid = byrow.aid)
}  # rgaitdpois










 dgaitdplot <-
  function(   # xx, pmf,
           theta.p,  # scalar, else a 2-vector
           fam = "pois",
           a.mix = NULL, i.mix = NULL,  # Unstructured probs are
           d.mix = NULL,
           a.mlm = NULL, i.mlm = NULL,  # contiguous
           d.mlm = NULL,
           truncate = NULL, max.support = Inf,
           pobs.mix = 0,  # scalar
           pobs.mlm = 0,  # vector of length a.mlm
           pstr.mix = 0,  # scalar
           pstr.mlm = 0,  # vector of length a.mlm
           pdip.mix = 0,  # scalar
           pdip.mlm = 0,  # vector of length a.mlm
           byrow.aid = FALSE,  # Applies to 'pobs.mlm' & 'pstr.mlm'
           theta.a = theta.p,  # scalar, else a 2-vector
           theta.i = theta.p,  # scalar, else a 2-vector
           theta.d = theta.p,  # scalar, else a 2-vector
           deflation = FALSE,  # Single logical FALSE, TRUE
           plot.it = TRUE, new.plot = TRUE,
           offset.x = ifelse(new.plot, 0, 0.25),
           type.plot = "h",  # Matches 'type' argument
           xlim = c(0, min(100, max.support + 2)),
           ylim = NULL,
           xlab = "",  # Was "y" prior to using oma
           ylab = "Probability",
           main = "",
           cex.main = 1.2, posn.main = NULL,
           all.col = NULL, all.lty = NULL, all.lwd = NULL, 
           lty.p     = "solid",  # longdash, dashed, twodash, solid
           lty.a.mix = "longdash",
           lty.a.mlm = "longdash",
           lty.i.mix = "dashed",
           lty.i.mlm = "dashed",
           lty.d.mix = "solid",
           lty.d.mlm = "solid",
           lty.d.dip = "dashed",
           col.p = "pink2",  # peach "salmon"  # salmon1,..., salmon4
           col.a.mix = artichoke.col,  # amazon.col,    # azure.col,  
           col.a.mlm = asparagus.col,  # avocado.col,   # "blue",
           col.i.mix = indigo.col,
           col.i.mlm = iris.col,   # "purple",   # 
           col.d.mix = deer.col,
           col.d.mlm = dirt.col,
           col.d.dip = desire.col,   # "red", # "orangered" maybe
           col.t = turquoise.col,   # "tan",
           cex.p = 1,
           lwd.p = NULL, lwd.a = NULL,
           lwd.i = NULL, lwd.d = NULL,  # Default: par()$lwd
           iontop = TRUE, dontop = TRUE,
           las = 0, lend = "round",  # "round", "butt", "square", 0:2
           axes.x = TRUE, axes.y = TRUE,
           Plot.trunc = TRUE, cex.t = 1, pch.t = 1,
           baseparams.argnames = NULL,   #,  # Optional safety
           nparams = 1,
           flip.args = FALSE,  # Set TRUE for "gaitdnbinomial"
           ...
           ) {  # ... ignored currently.



  if (!length(lwd.p)) lwd.p <- 1
  if (!length(lwd.a)) lwd.a <- 1
  if (!length(lwd.i)) lwd.i <- 1
  if (!length(lwd.d)) lwd.d <- 1

  if (length(all.col))
    col.p <- col.a.mix <-  col.a.mlm <- col.i.mix <- col.i.mlm <- 
    col.d.mix <- col.d.mlm <- col.t <- all.col
  if (length(all.lty))
    lty.p <- lty.a.mix <-  lty.a.mlm <- lty.i.mix <- lty.i.mlm <- 
    lty.d.mix <- lty.d.mlm <- all.lty
  if (length(all.lwd))
    lwd.p <- lwd.a <-  lwd.i <- lwd.d <- all.lwd

  gaitd.errorcheck(a.mix, a.mlm, i.mix, i.mlm,
                   d.mix, d.mlm, truncate, max.support,
                   nparams = nparams)  # length(theta.p) might be okay



  la.mix <- length(a.mix <- sort(a.mix))
  li.mix <- length(i.mix <- sort(i.mix))
  ld.mix <- length(d.mix <- sort(d.mix))
  la.mlm <- length(a.mlm)
  li.mlm <- length(i.mlm)
  ld.mlm <- length(d.mlm)
  ltrunc <- length(truncate <- sort(truncate))




  MM <- length(theta.p)
  if (MM != 1 && MM !=  2)
    stop("can only handle 1 or 2 parameters")

  xx <- seq(xlim[1], xlim[2])
  ind.A.mlm <- ind.A.mix <- ind.trunc <-
  ind.I.mlm <- ind.I.mix <-
  ind.D.mlm <- ind.D.mix <- FALSE

  if (length(truncate))
    ind.trunc <- (xx %in% truncate | max.support < xx)
  if (length(a.mix)) {
    ind.A.mix <- xx %in% a.mix
  } else {
    pobs.mix <- 0  # Make sure
  }
  if (length(a.mlm)) {
    ind.A.mlm <- xx %in% a.mlm
  } else {
    pobs.mlm <- 0  # Make sure
  }
  if (length(i.mix)) {
    ind.I.mix <- xx %in% i.mix
  } else {
    pstr.mix <- 0  # Make sure
  }
  if (length(i.mlm)) {
    ind.I.mlm <- xx %in% i.mlm
  } else {
    pstr.mlm <- 0  # Make sure
  }
  if (length(d.mix)) {
    ind.D.mix <- xx %in% d.mix
  } else {
    pdip.mix <- 0  # Make sure
  }
  if (length(d.mlm)) {
    ind.D.mlm <- xx %in% d.mlm
  } else {
    pdip.mlm <- 0  # Make sure
  }

  special.xx <- ind.A.mix | ind.A.mlm | ind.I.mix | ind.I.mlm |
                ind.D.mix | ind.D.mlm | ind.trunc
  if (length(pobs.mix) != 1) stop("bad input for argument 'pobs.mix'")
  if (length(pstr.mix) != 1) stop("bad input for argument 'pstr.mix'")
  if (length(pdip.mix) != 1) stop("bad input for argument 'pdip.mix'")
  if (length(a.mlm) && length(pobs.mlm) > length(a.mlm))
    warning("bad input for argument 'pobs.mlm'?")
  if (length(i.mlm) && length(pstr.mlm) > length(i.mlm))
    warning("bad input for argument 'pstr.mlm'?")
  if (length(d.mlm) && length(pdip.mlm) > length(d.mlm))
    warning("bad input for argument 'pdip.mlm'?")

      
  if (any(ind.A.mlm))
    pobs.mlm <- matrix(pobs.mlm, 1,  # length(xx),
                       length(a.mlm), byrow = byrow.aid)
  if (any(ind.I.mlm))
    pstr.mlm <- matrix(pstr.mlm, 1,  # length(xx),
                       length(i.mlm), byrow = byrow.aid)
  if (any(ind.D.mlm))
    pdip.mlm <- matrix(pdip.mlm, 1,  # length(xx),
                       length(d.mlm), byrow = byrow.aid)

      
  dfun <- paste0("dgaitd", fam)


  pmf.p <-
    if (MM == 1)
      do.call(dfun, list(x =  xx, theta.p)) else
    if (MM == 2) {
      if (flip.args)
        do.call(dfun, list(x =  xx, theta.p[2], theta.p[1])) else
        do.call(dfun, list(x =  xx, theta.p[1], theta.p[2]))
    }



  alist <- list(  # x = xx,  # theta.p,
            a.mix = a.mix, a.mlm = a.mlm,
            i.mix = i.mix, i.mlm = i.mlm,
            d.mix = d.mix, d.mlm = d.mlm,
            truncate = truncate, max.support = max.support,
            pobs.mix = pobs.mix, pobs.mlm = pobs.mlm,
            pstr.mix = pstr.mix, pstr.mlm = pstr.mlm,
            pdip.mix = pdip.mix, pdip.mlm = pdip.mlm,
            byrow.aid = byrow.aid)

  if (length(baseparams.argnames)) {
    alist[[paste0(baseparams.argnames[1], ".p")]] <- theta.p[1]
    alist[[paste0(baseparams.argnames[1], ".a")]] <- theta.a[1]
    alist[[paste0(baseparams.argnames[1], ".i")]] <- theta.i[1]
    alist[[paste0(baseparams.argnames[1], ".d")]] <- theta.d[1]
    if (MM == 2) {
      alist[[paste0(baseparams.argnames[2], ".p")]] <- theta.p[2]
      alist[[paste0(baseparams.argnames[2], ".a")]] <- theta.a[2]
      alist[[paste0(baseparams.argnames[2], ".i")]] <- theta.i[2]
      alist[[paste0(baseparams.argnames[2], ".d")]] <- theta.d[2]
    }  # else
  } else {
    if (MM == 1) {
      alist <- c(alist,  # Unnamed, for lambda.p, etc.:
                 list(theta.p, theta.a, theta.i, theta.d))
    } else {  # MM == 2
      alist <- if (flip.args)
               c(alist,
                 list(theta.p[2], theta.p[1],  # Order is crucial.
                      theta.a[2], theta.i[2], theta.d[2],
                      theta.a[1], theta.i[1], theta.d[1])) else
               c(alist,
                 list(theta.p[1], theta.p[2],  # Order is crucial.
                      theta.a[1], theta.i[1], theta.d[1],
                      theta.a[2], theta.i[2], theta.d[2]))
    }
  }


  dlist <- alist
  dlist$x <- xx
  dlist$log <- FALSE
  pmf.z <- do.call(dfun, dlist)
  if (!all(is.finite(pmf.z)))
    warning("some PMF values are not finite")
  if ((rp.z <- range(pmf.z, na.rm = TRUE))[1] < 0 || rp.z[2] > 1)
    warning("some computed PMF values not in [0, 1]")


  mlist <- alist
  mlist$type.fitted <- "All"
  mlist$moments2 <- TRUE
  mom.fun <- paste0("moments.gaitdcombo.", fam)
  Bits <- do.call(mom.fun, mlist)

  myylim <- if (is.null(ylim))
    c(0, max(0, pmf.z, na.rm = TRUE) * 1.04) else ylim


  if (plot.it) {
    if (new.plot)
      plot(xx[!special.xx] + offset.x,
           pmf.z[!special.xx], type = "n",
           las = las, axes = axes.x && axes.y,  # axes,
           xlim = xlim,
           ylim = myylim,  # c(0, 0.11),  # Fix ylim here

         
           xlab = xlab, ylab = ylab)  # Parent distribution

    lines(xx[!special.xx] + offset.x,
          pmf.z[!special.xx], type = type.plot, col = col.p,
          lwd = lwd.p,  # Ordinary points
          cex = cex.p, lend = lend, lty = lty.p)

    if (length(posn.main)) {
      posn.main <- rep(posn.main, 2)
      text(posn.main[1], posn.main[2], labels = main,
           cex = cex.main)
    } else {
      title(main = main)  # , cex = cex.main
    }
  }  # plot.it

      
  if (plot.it && !(axes.x && axes.y)) {
    box()
    axis(1, labels = FALSE, tick = TRUE)
    axis(2, labels = FALSE, tick = TRUE)
  }
  if (plot.it && !axes.x && axes.y) {  # Avoid clutter
    axis(1, labels = FALSE)
    axis(2, las = 1)
  }
  if (plot.it && axes.x && !axes.y) {  # Avoid clutter
    axis(1)
    axis(2, labels = FALSE)
  }



  if (plot.it && length(a.mlm))  # Altered mlm
    lines(xx[ ind.A.mlm] + offset.x, pmf.z[ind.A.mlm],
          lwd = lwd.a, type = type.plot,
          col = col.a.mlm, lty = lty.a.mlm, lend = lend)


  if (plot.it && length(a.mix))  # Altered mix
    lines(xx[ ind.A.mix] + offset.x, pmf.z[ind.A.mix],
          lwd = lwd.a, type = type.plot,
          col = col.a.mix, lty = lty.a.mix, lend = lend)




    Denom.p <- as.vector(Bits[["cdf.max.s"]] - Bits[["SumT0.p"]] -
                         Bits[["SumA0.mix.p"]] - Bits[["SumA0.mlm.p"]])
    if (any(Denom.p == 0))
      stop("0s found in the denominator (variable 'Denom.p')")
    Numer <- as.vector(1 -
               (if (length(a.mix)) pobs.mix else 0) -
               (if (length(i.mix)) pstr.mix else 0) +
               (if (length(d.mix)) pdip.mix else 0) -
               (if (length(a.mlm)) rowSums(rbind(pobs.mlm)) else 0) -
               (if (length(i.mlm)) rowSums(rbind(pstr.mlm)) else 0) +
               (if (length(d.mlm)) rowSums(rbind(pdip.mlm)) else 0))
    if (any(Numer < 0))
      warning("variable 'Numer' has negative values")
    Numer <- Numer[1]




  if (any(ind.I.mix)) {  # Inflated mix

    spikes.mix <- pmf.z[ind.I.mix]  # Top of the spike
    start.pt.i.mix <- Numer * pmf.p[ind.I.mix] / Denom.p  # Skin

  if (plot.it)
    segments(xx[ind.I.mix] + offset.x,
             if (iontop) start.pt.i.mix else   # Outer distn
             spikes.mix - start.pt.i.mix,
             xx[ind.I.mix] + offset.x,
             spikes.mix,  # This value is unchanged
             lwd = if (iontop) lwd.i else lwd.p,
             lty = if (iontop) lty.i.mix else lty.p,
             col = if (iontop) col.i.mix else col.p, lend = lend)


  if (plot.it)
    lines(xx[ind.I.mix] + offset.x,
          if (iontop) start.pt.i.mix else
          spikes.mix - start.pt.i.mix,
          lend = lend,
          lwd = if (iontop) lwd.p else lwd.i,
          lty = if (iontop) lty.p else lty.i.mix,
          col = if (iontop) col.p else col.i.mix,
          type = type.plot)  # Blend in
  }  # ind.I.mix




  
  if (any(ind.I.mlm)) {  # Inflated mlm
    start.pt.i.mlm <- Numer * pmf.p[ind.I.mlm] / Denom.p
    spikes.mlm <- pmf.z[ind.I.mlm]  # Top of the spike

  if (plot.it)
    segments(xx[ind.I.mlm] + offset.x,
             if (iontop) start.pt.i.mlm else  # Outer distn
             spikes.mlm - start.pt.i.mlm,
             xx[ind.I.mlm] + offset.x,
             spikes.mlm,  # This value is unchanged
             lend = lend,
             lwd = if (iontop) lwd.i else lwd.p,
             col = if (iontop) col.i.mlm else col.p,
             lty = if (iontop) lty.i.mlm else lty.p)


  if (plot.it)
    lines(xx[ind.I.mlm] + offset.x,
          if (iontop) start.pt.i.mlm else
          spikes.mlm - start.pt.i.mlm,
          lend = lend,
          lwd = if (iontop) lwd.p else lwd.i,
          lty = if (iontop) lty.p else lty.i.mlm,
          col = if (iontop) col.p else col.i.mlm,
          type = type.plot)  # Blend in
  }  # ind.I.mlm

      





  if (any(ind.D.mix)) {  # Deflated mix

    start.pt.d.mix <- Numer * pmf.p[ind.D.mix] / Denom.p  # Skin
    dips.mix <- pmf.z[ind.D.mix]  # Bottom of the dip

  if (plot.it && deflation)
    segments(xx[ind.D.mix] + offset.x,
             start.pt.d.mix,  # This value is unchanged; top
             xx[ind.D.mix] + offset.x,
             if (dontop) dips.mix else
             start.pt.d.mix - dips.mix,
             lwd = if (dontop) lwd.d else lwd.d,
             lty = if (dontop) lty.d.dip else lty.d.mix,
             col = if (dontop) col.d.dip else col.d.mix,
             lend = lend)


  if (plot.it)
    lines(xx[ind.D.mix] + offset.x,
          if (deflation) {
            if (dontop) dips.mix else start.pt.d.mix - dips.mix
          } else dips.mix,
          lend = lend,
          lwd = if (deflation && !dontop) lwd.d else lwd.d,
          lty = if (deflation && !dontop) lty.d.dip else lty.d.mix,
          col = if (deflation && !dontop) col.d.dip else col.d.mix,
          type = type.plot)  # Blend in
  }  # ind.D.mix





  
  if (any(ind.D.mlm)) {  # Inflated mlm
    start.pt.d.mlm <- Numer * pmf.p[ind.D.mlm] / Denom.p  # Skin
    dips.mlm <- pmf.z[ind.D.mlm]  # Bottom of the dip; midstream

  if (plot.it && deflation)
    segments(xx[ind.D.mlm] + offset.x,
             start.pt.d.mlm,  # This value is unchanged; top
             xx[ind.D.mlm] + offset.x,
             if (dontop) dips.mlm else start.pt.d.mlm - dips.mlm,
             lwd = if (dontop) lwd.d else lwd.d,
             lty = if (dontop) lty.d.dip else lty.d.mlm,
             col = if (dontop) col.d.dip else col.d.mlm,
             lend = lend)


  if (plot.it)
    lines(xx[ind.D.mlm] + offset.x,
          if (deflation && !dontop) start.pt.d.mlm - dips.mlm else 
          dips.mlm,
          lend = lend,
          lwd = if (deflation && !dontop) lwd.d else lwd.d,
          lty = if (deflation && !dontop) lty.d.dip else lty.d.mlm,
          col = if (deflation && !dontop) col.d.dip else col.d.mlm,
          type = type.plot)  # Blend in

    if (any((if (deflation && !dontop) start.pt.d.mlm - dips.mlm else 
             dips.mlm) < 0))
      stop("negative probabilities; 'pdip.mlm' too large?")
  }  # ind.D.mlm

      




  if (Plot.trunc && plot.it && any(ind.trunc)) {
    lhs.tvec <- xx[ind.trunc] + offset.x * 0
    lhs.tvec <- lhs.tvec[ceiling(par()$usr[1]) + 0.5 < lhs.tvec]
    if (length(lhs.tvec))
      points(lhs.tvec, numeric(length(lhs.tvec)),
             col = col.t, cex = cex.t, pch = pch.t)
  }
  if (Plot.trunc && plot.it) {
    rhs.tvec <- if (floor(par()$usr[2]) > max.support + 2)
                (max.support + 1):(floor(par()$usr[2]) - 1) else NULL
    if (length(rhs.tvec))
      points(rhs.tvec + offset.x * 0,
             numeric(length(rhs.tvec)),
             col = col.t, cex = cex.t, pch = pch.t)
  }






  names(pmf.p) <- as.character(xx)
  invisible(list(x           = xx,
                 pmf.z       = pmf.z,  # Top profile == PMF.
                 sc.parent   = Numer * pmf.p / Denom.p,  # Skin.
                 unsc.parent = pmf.p))   # FYI only.
}  # dgaitdplot












plotdgaitd.vglm <-
  function(object, ...) {

  infos.list <- object@family@infos()
  specvals <- specials(object)


  Inside <- sapply(specvals, is.null)
  if (length(Inside) == 7 && all(Inside))
    stop("'object' has no special values. ",
         "Is it a GAITD regression object?")
  if (length(Inside) == 8 && all(Inside[1:7]) &&
      infos.list$max.support == infos.list$Support[2])
    stop("'object' has no special values. ",
         "Is it really a GAITD regression object?")
  use.max.support <- if (is.numeric(infos.list$max.support))
    infos.list$max.support else Inf

  if (!is.numeric(MM1 <- infos.list$MM1))
    MM1 <- 1  # Default really
  if (MM1 > 2)
    stop("Can only handle 1- or 2-parameter distributions")

  etamat <- predict(object)  # n x M
  eta.p <- etamat[, 1:MM1, drop = FALSE]  # n x MM1
  theta.p1 <- as.vector(eta2theta(eta.p[, 1], linkfun(object)[1]))
  theta.p2 <- if (MM1 == 2)
    as.vector(eta2theta(eta.p[, 2], linkfun(object)[2])) else NULL
  theta.p <- cbind(theta.p1, theta.p2)
  colnames(theta.p) <- paste0(infos.list$baseparams.argnames, ".p")
  if (!is.logical(intercept.only <- object@misc$intercept.only))
    stop("cannot determine whether 'object' is intercept-only")
  if (!intercept.only)
    warning("argument 'object' is not intercept-only")


  pobs.mix <- if (length(specvals$a.mix))
    fitted(object, type.fitted = "pobs.mix") else cbind(0, 0)
  pobs.mlm <- if (length(specvals$a.mlm))
    fitted(object, type.fitted = "pobs.mlm") else cbind(0, 0)
  pstr.mix <- if (length(specvals$i.mix))
    fitted(object, type.fitted = "pstr.mix") else cbind(0, 0)
  pstr.mlm <- if (length(specvals$i.mlm))
    fitted(object, type.fitted = "pstr.mlm") else cbind(0, 0)
  pdip.mix <- if (length(specvals$d.mix))
    fitted(object, type.fitted = "pdip.mix") else cbind(0, 0)
  pdip.mlm <- if (length(specvals$d.mlm))
    fitted(object, type.fitted = "pdip.mlm") else cbind(0, 0)


  indeta <- object@extra$indeta
  if (MM1 == 1) {
    theta.a <- if (any(is.na(indeta[ 3, ]))) theta.p else
                 as.vector(eta2theta(etamat[, (indeta[ 3, 1])],
                           linkfun(object)[(indeta[ 3, 1])]))
    theta.i <- if (any(is.na(indeta[ 5, ]))) theta.p else
                 as.vector(eta2theta(etamat[, (indeta[ 5, 1])],
                           linkfun(object)[(indeta[ 5, 1])]))
    theta.d <- if (any(is.na(indeta[ 7, ]))) theta.p else
                 as.vector(eta2theta(etamat[, (indeta[ 7, 1])],
                           linkfun(object)[(indeta[ 7, 1])]))

    theta.a <- cbind(theta.a)
    theta.i <- cbind(theta.i)
    theta.d <- cbind(theta.d)
  } else {
    theta.a <- if (any(is.na(indeta[ 4, ]))) theta.p else
                 cbind(eta2theta(etamat[, (indeta[ 4, 1])],
                                 linkfun(object)[(indeta[ 4, 1])]),
                       eta2theta(etamat[, (indeta[ 5, 1])],
                                 linkfun(object)[(indeta[ 5, 1])]))
    colnames(theta.a) <- paste0(infos.list$baseparams.argnames, ".a")
    theta.i <- if (any(is.na(indeta[ 7, ]))) theta.p else
                 cbind(eta2theta(etamat[, (indeta[ 7, 1])],
                                 linkfun(object)[(indeta[ 7, 1])]),
                       eta2theta(etamat[, (indeta[ 8, 1])],
                                 linkfun(object)[(indeta[ 8, 1])]))
    colnames(theta.i) <- paste0(infos.list$baseparams.argnames, ".i")
    theta.d <- if (any(is.na(indeta[10, ]))) theta.p else
                 cbind(eta2theta(etamat[, (indeta[10, 1])],
                                 linkfun(object)[(indeta[10, 1])]),
                       eta2theta(etamat[, (indeta[11, 1])],
                                 linkfun(object)[(indeta[11, 1])]))
    colnames(theta.d) <- paste0(infos.list$baseparams.argnames, ".d")
  }


  flip.args <- object@family@infos()$flip.args

  dgaitdplot(theta.p[1, ],  # Reverse ordering may be needed.
       fam = infos.list$parent.name[2],
       a.mix = specvals$a.mix, i.mix = specvals$i.mix, 
       d.mix = specvals$d.mix,
       a.mlm = specvals$a.mlm, i.mlm = specvals$i.mlm,
       d.mlm = specvals$d.mlm,
       truncate = specvals$truncate,
       theta.a = theta.a[1, ],  # Reverse ordering may be needed.
       theta.i = theta.i[1, ],
       theta.d = theta.d[1, ],
       max.support = use.max.support,
       pobs.mix = pobs.mix[1, ],
       pobs.mlm = pobs.mlm[1, ],
       pstr.mix = pstr.mix[1, ],
       pstr.mlm = pstr.mlm[1, ],
       pdip.mix = pdip.mix[1, ],  # 1-coln matrix
       pdip.mlm = pdip.mlm[1, ],
       byrow.aid = TRUE,  # Important really here 20201008
       baseparams.argnames = infos.list$baseparams.argnames,
       nparams = object@family@infos()$MM1,  # Unnecessary?
       flip.args = ifelse(is.logical(flip.args), flip.args, FALSE),
       ...)
}  # plotdgaitd.vglm



 if (!isGeneric("plotdgaitd"))
   setGeneric("plotdgaitd",
              function(object, ...) standardGeneric("plotdgaitd"))

  setMethod("plotdgaitd", signature(object = "vglm"),
            function(object, ...)
            invisible(plotdgaitd.vglm(object, ...)))














 dgaitdnbinom <-
  function(x, size.p,   # prob.p = NULL,
           munb.p,   #  = NULL,
           a.mix = NULL,
           a.mlm = NULL,
           i.mix = NULL,
           i.mlm = NULL,
           d.mix = NULL,
           d.mlm = NULL,
           truncate = NULL, max.support = Inf,
           pobs.mix = 0,  # vector
           pobs.mlm = 0,  # matrix
           pstr.mix = 0,  # vector
           pstr.mlm = 0,  # matrix
           pdip.mix = 0,  # vector
           pdip.mlm = 0,  # matrix
           byrow.aid = FALSE,  # Applies to 'pobs.mlm' & 'pstr.mlm'
           size.a = size.p, size.i = size.p, size.d = size.p,
           munb.a = munb.p, munb.i = munb.p, munb.d = munb.p,
           log = FALSE) {
  prob.p = NULL
  prob.a = prob.p; prob.i = prob.p; prob.d = prob.p

  log.arg <- log;  rm(log)
  lowsup <- 0
  gaitd.errorcheck(a.mix, a.mlm, i.mix, i.mlm,
                   d.mix, d.mlm, truncate, max.support, nparams = 2)

  if ((is.prob <- as.logical(length(prob.p))) &&
      length(munb.p))
    stop("cannot specify both 'prob.p' and 'munb.p' arguments")


  la.mix <- length(a.mix <- sort(a.mix))
  li.mix <- length(i.mix <- sort(i.mix))
  ld.mix <- length(d.mix <- sort(d.mix))
  la.mlm <- length(a.mlm)
  li.mlm <- length(i.mlm)
  ld.mlm <- length(d.mlm)
  ltrunc <- length(truncate <- sort(truncate))
  if (la.mix + la.mlm + li.mix + li.mlm + ld.mix + ld.mlm +
      ltrunc == 0 &&
      is.infinite(max.support))
    return(if (is.prob)
      dnbinom(x, size = size.p, prob = prob.p, log = log.arg) else
      dnbinom(x, size = size.p, mu   = munb.p, log = log.arg))


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

    
  if (any(pdip.mix < 0 | 1 <= pdip.mix, na.rm = TRUE)) {
    stop("bad input for argument 'pdip.mix'")
  }
    
  if (any(pdip.mlm < 0 | 1 <= pdip.mlm, na.rm = TRUE)) {
    stop("bad input for argument 'pdip.mlm'")
  }
   


  LLL <- max(length(x),
             length(pobs.mix), length(pstr.mix), length(pdip.mix),
             length(size.p),   length(size.a),   length(size.i),
             length(size.d),
             length(munb.p),   length(munb.a),   length(munb.i),
             length(munb.d),
             length(prob.p),   length(prob.a),   length(prob.i),
             length(prob.d))
  if (length(x)        < LLL) x        <- rep_len(x,        LLL)
  if (length(size.p)   < LLL) size.p   <- rep_len(size.p,   LLL)
  if (length(size.a)   < LLL) size.a   <- rep_len(size.a,   LLL)
  if (length(size.i)   < LLL) size.i   <- rep_len(size.i,   LLL)
  if (length(size.d)   < LLL) size.d   <- rep_len(size.d,   LLL)
  if (is.prob) {
  if (length(prob.p)   < LLL) prob.p   <- rep_len(prob.p,   LLL)
  if (length(prob.a)   < LLL) prob.a   <- rep_len(prob.a,   LLL)
  if (length(prob.i)   < LLL) prob.i   <- rep_len(prob.i,   LLL)
  if (length(prob.d)   < LLL) prob.d   <- rep_len(prob.d,   LLL)
  } else {
  if (length(munb.p)   < LLL) munb.p   <- rep_len(munb.p,   LLL)
  if (length(munb.a)   < LLL) munb.a   <- rep_len(munb.a,   LLL)
  if (length(munb.i)   < LLL) munb.i   <- rep_len(munb.i,   LLL)
  if (length(munb.d)   < LLL) munb.d   <- rep_len(munb.d,   LLL)
  }
  if (length(pobs.mix) < LLL) pobs.mix <- rep_len(pobs.mix, LLL)
  if (length(pstr.mix) < LLL) pstr.mix <- rep_len(pstr.mix, LLL)
  if (length(pdip.mix) < LLL) pdip.mix <- rep_len(pdip.mix, LLL)



  sumt <- 0  # Initialization to 0 important
  if (ltrunc) {
    if (is.prob) {  # Need tval <= max.support
      for (tval in truncate)
        sumt <- sumt + dnbinom(tval, size.p, prob = prob.p)
    } else {
      for (tval in truncate)
        sumt <- sumt + dnbinom(tval, size.p, mu   = munb.p)
    }
  }
  vecTF.t <- is.finite(x) & ((x %in% truncate) | (max.support < x))
  cdf.max.s <- if (is.prob)
                pnbinom(max.support, size.p, prob = prob.p) else
                pnbinom(max.support, size.p, mu   = munb.p)  # Usually 1
  denom.t <- cdf.max.s - sumt  # No sumt on RHS

    pmf0 <- if (is.prob)  # dgtnbinom
              ifelse(vecTF.t, 0, dnbinom(max.support, size.p,
                                         prob = prob.p) / denom.t) else
              ifelse(vecTF.t, 0, dnbinom(max.support, size.p,
                                         mu   = munb.p) / denom.t)


  sum.a <- suma <- 0  # numeric(LLL)
  vecTF.a <- rep_len(FALSE, LLL)
  if (la.mlm) {
    pobs.mlm <-  matrix(pobs.mlm, LLL, la.mlm, byrow = byrow.aid)
    sum.a <-   .rowSums(pobs.mlm, LLL, la.mlm)
    if (any(1 < sum.a, na.rm = TRUE))
      stop("bad input for argument 'pobs.mlm'")  # zz

    if (is.prob) {  # Part i
      for (aval in a.mlm)
        suma <- suma + dnbinom(aval, size.p, prob = prob.p)
    } else {
      for (aval in a.mlm)
        suma <- suma + dnbinom(aval, size.p, mu   = munb.p)
    }

      
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
    pmf2.a <- dgaitdnbinom(x,  # Outer distribution---mlm type
                          size.p = size.a,
                          munb.p = munb.a,
                          truncate = setdiff(allx.a, a.mix),
                          max.support = max(a.mix))
    for (aval in a.mix) {
      suma <- suma + (if (is.prob)  # Part ii added; cumulative
                      dnbinom(aval, size = size.p, prob = prob.p) else
                      dnbinom(aval, size = size.p, mu   = munb.p))
      vecTF <- is.finite(x) & aval == x
      pmf0[vecTF] <- 0  # added; the true values are assigned below
      vecTF.a <- vecTF.a | vecTF  # Cumulative; added
    }
  }  # la.mix


  if (li.mix) {
    allx.i <- if (length(i.mix)) lowsup:max(i.mix) else NULL
    pmf2.i <- dgaitdnbinom(x,  # Outer distn---mlm type
                           size.p = size.i,
                           munb.p = munb.i,
                           truncate = setdiff(allx.i, i.mix),
                           max.support = max(i.mix))
  }  # li.mix






  sum.d <- 0  # numeric(LLL)
  if (ld.mlm) {
    pdip.mlm <- matrix(pdip.mlm, LLL, ld.mlm, byrow = byrow.aid)
    sum.d <-  .rowSums(pdip.mlm, LLL, ld.mlm)
    if (any(1 < sum.d, na.rm = TRUE))
      stop("bad input for argument 'pdip.mlm'")
  }  # ld.mlm


  if (ld.mix) {
    allx.d <- lowsup:max(d.mix)
    pmf2.d <- dgaitdnbinom(x, size.p = size.d,  # prob.p = prob.d,
                           munb.p = munb.d,
                           truncate = setdiff(allx.d, d.mix),
                           max.support = max(d.mix))
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



  pmf0[!skip] <- (tmp6 * (if (is.prob) 
    dnbinom(x, size.p, prob = prob.p) else
    dnbinom(x, size.p, mu   = munb.p)) / (
    cdf.max.s - suma - sumt))[!skip]  # added


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
}  # dgaitdnbinom







 pgaitdnbinom <-
  function(q, size.p,   # prob.p = NULL,
           munb.p,  # = NULL,
           a.mix = NULL,
           a.mlm = NULL,
           i.mix = NULL,
           i.mlm = NULL,
           d.mix = NULL,
           d.mlm = NULL,
           truncate = NULL, max.support = Inf,
           pobs.mix = 0,
           pobs.mlm = 0,
           pstr.mix = 0,
           pstr.mlm = 0,
           pdip.mix = 0,
           pdip.mlm = 0,
           byrow.aid = FALSE,
           size.a = size.p, size.i = size.p, size.d = size.p,
           munb.a = munb.p, munb.i = munb.p, munb.d = munb.p,
           lower.tail = TRUE) {
  prob.p = NULL
  prob.a = prob.p; prob.i = prob.p; prob.d = prob.p


  lowsup <- 0
  gaitd.errorcheck(a.mix, a.mlm, i.mix, i.mlm,
                   d.mix, d.mlm, truncate, max.support, nparams = 2)

  if ((is.prob <- as.logical(length(prob.p))) &&
      length(munb.p))
    stop("cannot specify both 'prob.p' and 'munb.p' arguments")

  la.mix <- length(a.mix <- sort(a.mix))
  li.mix <- length(i.mix <- sort(i.mix))
  ld.mix <- length(d.mix <- sort(d.mix))
  la.mlm <- length(a.mlm)
  li.mlm <- length(i.mlm)
  ld.mlm <- length(d.mlm)
  ltrunc <- length(truncate <- sort(truncate))
  if (la.mix + la.mlm + li.mix + li.mlm + ld.mix + ld.mlm +
      ltrunc == 0 &&
      is.infinite(max.support))
    return(if (is.prob)
           pnbinom(q, size = size.p, prob = prob.p,
                   lower.tail = lower.tail) else
           pnbinom(q, size = size.p, mu   = munb.p,
                   lower.tail = lower.tail))  # log.p

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
             length(munb.p),   length(munb.a),   length(munb.i),
             length(munb.d),
             length(prob.p),   length(prob.a),   length(prob.i),
             length(prob.d))
  offset.a <- offset.i <- offset.d <-
  Offset.a <- Offset.i <- Offset.d <- numeric(LLL)
  if (length(q)        < LLL) q        <- rep_len(q,        LLL)
  if (length(size.p)   < LLL) size.p   <- rep_len(size.p,   LLL)
  if (length(size.a)   < LLL) size.a   <- rep_len(size.a,   LLL)
  if (length(size.i)   < LLL) size.i   <- rep_len(size.i,   LLL)
  if (length(size.d)   < LLL) size.d   <- rep_len(size.d,   LLL)
  if (is.prob) {
  if (length(prob.p)   < LLL) prob.p   <- rep_len(prob.p,   LLL)
  if (length(prob.a)   < LLL) prob.a   <- rep_len(prob.a,   LLL)
  if (length(prob.i)   < LLL) prob.i   <- rep_len(prob.i,   LLL)
  if (length(prob.d)   < LLL) prob.d   <- rep_len(prob.d,   LLL)
  } else {
  if (length(munb.p)   < LLL) munb.p   <- rep_len(munb.p,   LLL)
  if (length(munb.a)   < LLL) munb.a   <- rep_len(munb.a,   LLL)
  if (length(munb.i)   < LLL) munb.i   <- rep_len(munb.i,   LLL)
  if (length(munb.d)   < LLL) munb.d   <- rep_len(munb.d,   LLL)
  }
  if (length(pobs.mix) < LLL) pobs.mix <- rep_len(pobs.mix, LLL)
  if (length(pstr.mix) < LLL) pstr.mix <- rep_len(pstr.mix, LLL)
  if (length(pdip.mix) < LLL) pdip.mix <- rep_len(pdip.mix, LLL)



  sumt <- 0
  fudge.t <- numeric(LLL)
  cdf.max.s <- if (is.prob)
                pnbinom(max.support, size.p, prob = prob.p) else
                pnbinom(max.support, size.p, mu   = munb.p)  # Usually 1
  if (ltrunc) {
    for (tval in truncate) {
      pmf.p <- if (is.prob) dnbinom(tval, size.p, prob = prob.p) else
                            dnbinom(tval, size.p, mu   = munb.p)
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
      pmf.p <- if (is.prob) dnbinom(aval, size.p, prob = prob.p) else
                            dnbinom(aval, size.p, mu   = munb.p)
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
      pmf.a <- if (is.prob) dnbinom(aval, size.a, prob = prob.a) else
                            dnbinom(aval, size.a, mu   = munb.a)
      pmf.p <- if (is.prob) dnbinom(aval, size.p, prob = prob.p) else
                            dnbinom(aval, size.p, mu   = munb.p)
      use.pobs.mix[, jay] <- pmf.a
      suma <- suma + pmf.p  # cumulative; part ii
    }
    use.pobs.mix <- pobs.mix *
                    use.pobs.mix / rowSums(use.pobs.mix)

    for (jay in seq(la.mix)) {
      aval <- a.mix[jay]
      pmf.p <- if (is.prob) dnbinom(aval, size.p, prob = prob.p) else
                            dnbinom(aval, size.p, mu   = munb.p)
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
      use.pstr.mix[, jay] <- if (is.prob)
        dnbinom(ival, size.i, prob = prob.i) else
        dnbinom(ival, size.i, mu   = munb.i)
    }
    use.pstr.mix <- pstr.mix *
                    use.pstr.mix / rowSums(use.pstr.mix)

    for (jay in seq(li.mix)) {
      ival <- i.mix[jay]
      pmf.p <- if (is.prob) dnbinom(ival, size.p, prob = prob.p) else
                            dnbinom(ival, size.p, mu   = munb.p)
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
      use.pdip.mix[, jay] <- (if (is.prob)
         dnbinom(dval, size.d, prob = prob.d) else
         dnbinom(dval, size.d, mu   = munb.d))
    }
    use.pdip.mix <- pdip.mix *
                    use.pdip.mix / rowSums(use.pdip.mix)

    for (jay in seq(ld.mix)) {
      dval <- d.mix[jay]
      pmf.p <- if (is.prob)
                pnbinom(dval, size.p, prob = prob.p) else
                pnbinom(dval, size.p, mu   = munb.p)
      if (any(vecTF <- (is.finite(q) & dval <= q))) {
        Offset.d[vecTF] <- Offset.d[vecTF] + use.pdip.mix[vecTF, jay]
      }
    }  # jay
  }  # ld.mix







  numer1 <- 1 - sum.i - sum.a - pstr.mix - pobs.mix + sum.d + pdip.mix
  denom1 <- cdf.max.s - sumt - suma
  ans <- numer1 * ((if (is.prob)
         pnbinom(q, size.p, prob = prob.p) else
         pnbinom(q, size.p, mu   = munb.p)) -
         fudge.t - fudge.a) / denom1 +
         offset.a + offset.i - offset.d +
         Offset.a + Offset.i - Offset.d
  ans[max.support <= q] <- 1
  ans[ans < 0] <- 0  # Occasional roundoff error
  if (lower.tail) ans else 1 - ans
}  # pgaitdnbinom








 qgaitdnbinom <-
  function(p, size.p,   # prob.p = NULL,
           munb.p,  # = NULL,
           a.mix = NULL,
           a.mlm = NULL,
           i.mix = NULL,
           i.mlm = NULL,
           d.mix = NULL,
           d.mlm = NULL,
           truncate = NULL, max.support = Inf,
           pobs.mix = 0,
           pobs.mlm = 0,
           pstr.mix = 0,
           pstr.mlm = 0,
           pdip.mix = 0,
           pdip.mlm = 0,
           byrow.aid = FALSE,
           size.a = size.p, size.i = size.p, size.d = size.p,
           munb.a = munb.p, munb.i = munb.p, munb.d = munb.p) {
  prob.p = NULL
  prob.a = prob.p; prob.i = prob.p; prob.d = prob.p

  lowsup <- 0
  gaitd.errorcheck(a.mix, a.mlm, i.mix, i.mlm,
                   d.mix, d.mlm, truncate, max.support, nparams = 2)

  if ((is.prob <- as.logical(length(prob.p))) &&
      length(munb.p))
    stop("cannot specify both 'prob.p' and 'munb.p' arguments")

  la.mix <- length(a.mix <- sort(a.mix))
  li.mix <- length(i.mix <- sort(i.mix))
  ld.mix <- length(d.mix <- sort(d.mix))
  la.mlm <- length(a.mlm)
  li.mlm <- length(i.mlm)
  ld.mlm <- length(d.mlm)
  ltrunc <- length(truncate <- sort(truncate))
  if (la.mix + la.mlm + li.mix + li.mlm + ld.mix + ld.mlm +
      ltrunc == 0 &&
      is.infinite(max.support))
    return(if (is.prob)
      qnbinom(p, size = size.p, prob = prob.p) else
      qnbinom(p, size = size.p, mu   = munb.p))  # lower.t,log.p

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
             length(size.p),   length(size.a),   length(size.i),
             length(size.d),
             length(munb.p),   length(munb.a),   length(munb.i),
             length(munb.d),
             length(prob.p),   length(prob.a),   length(prob.i),
             length(prob.d))
  if (length(p)        < LLL) p        <- rep_len(p,         LLL)
  if (length(size.p)   < LLL) size.p   <- rep_len(size.p,    LLL)
  if (length(size.a)   < LLL) size.a   <- rep_len(size.a,    LLL)
  if (length(size.i)   < LLL) size.i   <- rep_len(size.i,    LLL)
  if (length(size.d)   < LLL) size.d   <- rep_len(size.d,    LLL)
  if (is.prob) {
  if (length(prob.p)   < LLL) prob.p   <- rep_len(prob.p,    LLL)
  if (length(prob.a)   < LLL) prob.a   <- rep_len(prob.a,    LLL)
  if (length(prob.i)   < LLL) prob.i   <- rep_len(prob.i,    LLL)
  if (length(prob.d)   < LLL) prob.d   <- rep_len(prob.d,    LLL)
  } else {
  if (length(munb.p)   < LLL) munb.p   <- rep_len(munb.p,    LLL)
  if (length(munb.a)   < LLL) munb.a   <- rep_len(munb.a,    LLL)
  if (length(munb.i)   < LLL) munb.i   <- rep_len(munb.i,    LLL)
  if (length(munb.d)   < LLL) munb.d   <- rep_len(munb.d,    LLL)
  }
  if (length(pobs.mix) < LLL) pobs.mix <- rep_len(pobs.mix,  LLL)
  if (length(pstr.mix) < LLL) pstr.mix <- rep_len(pstr.mix,  LLL)
  if (length(pdip.mix) < LLL) pdip.mix <- rep_len(pdip.mix,  LLL)


  pobs.mlm <- matrix(pobs.mlm, LLL, max(la.mlm, 1),
                     byrow = byrow.aid)
  pstr.mlm <- matrix(pstr.mlm, LLL, max(li.mlm, 1),
                     byrow = byrow.aid)
  pdip.mlm <- matrix(pdip.mlm, LLL, max(ld.mlm, 1),
                     byrow = byrow.aid)

  min.support <- lowsup  # Usual case; same as lowsup
  min.support.use <- if (ltrunc)
    min(setdiff(min.support:(ltrunc+5), truncate)) else min.support
  ans <- p + size.p + size.a + size.i + size.d + (if (is.prob)
             prob.p + prob.a + prob.i else 
             munb.p + munb.a + munb.i)

  bad0.p <- !is.finite(size.p) | size.p <= 0 |
            (if (is.prob)
              !is.finite(prob.p) | prob.p <= 0 | 1 <= prob.p else
              !is.finite(munb.p) | munb.p <= 0)
  bad0.a <- !is.finite(size.a) | size.a <= 0 |
            (if (is.prob)
              !is.finite(prob.a) | prob.a <= 0 | 1 <= prob.a else
              !is.finite(munb.a) | munb.a <= 0)
  bad0.i <- !is.finite(size.i) | size.i <= 0 |
            (if (is.prob)
              !is.finite(prob.i) | prob.i <= 0 | 1 <= prob.i else
              !is.finite(munb.i) | munb.i <= 0)
  bad0.d <- !is.finite(size.d) | size.d <= 0 |
            (if (is.prob)
              !is.finite(prob.d) | prob.d <= 0 | 1 <= prob.d else
              !is.finite(munb.d) | munb.d <= 0)
  bad0 <- bad0.p | bad0.a | bad0.i | bad0.d
  bad <- bad0 | !is.finite(p) | p <= 0 | 1 <= p

  Lo <- rep_len(min.support.use - 0.5, LLL)
  approx.ans <- Lo  # True at lhs
  Hi <- if (is.finite(max.support))
    rep(max.support + 0.5, LLL) else 2 * Lo + 10.5
  dont.iterate <- bad
  done <- dont.iterate |
  p <= pgaitdnbinom(Hi, size.p,  # prob.p = prob.p,
               munb.p = munb.p,
               a.mix = a.mix, a.mlm = a.mlm,
               i.mix = i.mix, i.mlm = i.mlm,
               d.mix = d.mix, d.mlm = d.mlm,
               truncate = truncate, max.support = max.support,
               pstr.mix = pstr.mix, pobs.mix = pobs.mix,
               pdip.mix = pdip.mix,
               pstr.mlm = pstr.mlm, pobs.mlm = pobs.mlm,
               pdip.mlm = pdip.mlm,
               munb.a = munb.a, munb.i = munb.i, munb.d = munb.d,
               byrow.aid = FALSE)

  iter <- 0
  max.iter <- round(log2(.Machine$double.xmax)) - 3
  while (!all(done) && iter < max.iter) {
    Lo[!done] <- Hi[!done]
    Hi[!done] <- 2 * Hi[!done] + 10.5  # Bug fixed
    Hi <- pmin(max.support + 0.5, Hi)  # 20190924
    done[!done] <-
      (p[!done] <= pgaitdnbinom(Hi[!done], size.p = size.p[!done],
                     munb.p = munb.p[!done],
                     a.mix = a.mix, a.mlm = a.mlm,
                     i.mix = i.mix, i.mlm = i.mlm,
                     d.mix = d.mix, d.mlm = d.mlm,
                     truncate = truncate, max.support = max.support,
                     pobs.mix = pobs.mix[!done],
                     pstr.mix = pstr.mix[!done],
                     pdip.mix = pdip.mix[!done],
                     pobs.mlm = pobs.mlm[!done, , drop = FALSE],
                     pstr.mlm = pstr.mlm[!done, , drop = FALSE],
                     pdip.mlm = pdip.mlm[!done, , drop = FALSE],
                     size.a = size.a[!done], size.i = size.i[!done],
                     size.d = size.d[!done],
                     munb.a = munb.a[!done], munb.i = munb.i[!done],
                     munb.d = munb.d[!done],
                     byrow.aid = FALSE))
    iter <- iter + 1
  }

    foo <- function(q, size.p,   # prob.p = NULL,
                    munb.p,  # = NULL,
                    a.mix = NULL, a.mlm = NULL,
                    i.mix = NULL, i.mlm = NULL,
                    d.mix = NULL, d.mlm = NULL,
                    truncate = NULL, max.support = Inf,
                    pobs.mix = 0, pstr.mix = 0, pdip.mix = 0,
                    pobs.mlm = 0, pstr.mlm = 0, pdip.mlm = 0,
                    size.a = size.p, size.i = size.p, size.d = size.p,
                    munb.a = munb.p, munb.i = munb.p, munb.d = munb.p,
                    byrow.aid = FALSE, p)
          pgaitdnbinom(q, size.p = size.p,
                     munb.p = munb.p,
                     a.mix = a.mix, a.mlm = a.mlm,
                     i.mix = i.mix, i.mlm = i.mlm,
                     d.mix = d.mix, d.mlm = d.mlm,
                     truncate = truncate, max.support = max.support,
                     pobs.mix = pobs.mix,
                     pstr.mix = pstr.mix,
                     pdip.mix = pdip.mix,
                     pobs.mlm = pobs.mlm,
                     pstr.mlm = pstr.mlm,
                     pdip.mlm = pdip.mlm,
                     size.a = size.a, size.i = size.i, size.d = size.d,
                     munb.a = munb.a, munb.i = munb.i, munb.d = munb.d,
                     byrow.aid = FALSE) - p

      lhs <- dont.iterate |
        p <= dgaitdnbinom(min.support.use, size.p = size.p,
                     munb.p = munb.p,
                     a.mix = a.mix, a.mlm = a.mlm,
                     i.mix = i.mix, i.mlm = i.mlm,
                     d.mix = d.mix, d.mlm = d.mlm,
                     truncate = truncate, max.support = max.support,
                     pobs.mix = pobs.mix,
                     pstr.mix = pstr.mix,
                     pdip.mix = pdip.mix,
                     pobs.mlm = pobs.mlm,
                     pstr.mlm = pstr.mlm,
                     pdip.mlm = pdip.mlm,
                     size.a = size.a, size.i = size.i, size.d = size.d,
                     munb.a = munb.a, munb.i = munb.i, munb.d = munb.d,
                     byrow.aid = FALSE)

  if (any(!lhs)) {
    approx.ans[!lhs] <-
      bisection.basic(foo, Lo[!lhs], Hi[!lhs], tol = 1/16,
                      size.p = size.p[!lhs],
                      munb.p = munb.p[!lhs],
                      a.mix = a.mix, a.mlm = a.mlm,
                      i.mix = i.mix, i.mlm = i.mlm,
                      d.mix = d.mix, d.mlm = d.mlm,
                      truncate = truncate, max.support = max.support,
                      pstr.mix = pstr.mix[!lhs],
                      pstr.mlm = pstr.mlm[!lhs, , drop = FALSE],
                      pobs.mix = pobs.mix[!lhs],
                      pobs.mlm = pobs.mlm[!lhs, , drop = FALSE],
                      pdip.mix = pdip.mix[!lhs],
                      pdip.mlm = pdip.mlm[!lhs, , drop = FALSE],
                      size.a = size.a[!lhs], size.i = size.i[!lhs],
                      size.d = size.d[!lhs],
                      munb.a = munb.a[!lhs], munb.i = munb.i[!lhs],
                      munb.d = munb.d[!lhs],
                      byrow.aid = FALSE,
                      p = p[!lhs])
    faa <- floor(approx.ans[!lhs])
    tmp <-
      ifelse(pgaitdnbinom(faa, size.p = size.p[!lhs],
                       munb.p = munb.p[!lhs],
                       a.mix = a.mix, a.mlm = a.mlm,
                       i.mix = i.mix, i.mlm = i.mlm,
                       d.mix = d.mix, d.mlm = d.mlm,
                       truncate = truncate, max.support = max.support,
                       pstr.mix = pstr.mix[!lhs],
                       pstr.mlm = pstr.mlm[!lhs, , drop = FALSE],
                       pobs.mix = pobs.mix[!lhs],
                       pobs.mlm = pobs.mlm[!lhs, , drop = FALSE],
                       pdip.mix = pdip.mix[!lhs],
                       pdip.mlm = pdip.mlm[!lhs, , drop = FALSE],
                       size.a = size.a[!lhs], size.i = size.i[!lhs],
                       size.d = size.d[!lhs],
                       munb.a = munb.a[!lhs], munb.i = munb.i[!lhs],
                       munb.d = munb.d[!lhs],
                       byrow.aid = FALSE) < p[!lhs] &
             p[!lhs] <= pgaitdnbinom(faa + 1, size.p = size.p[!lhs],
                       munb.p = munb.p[!lhs],
                       a.mix = a.mix, a.mlm = a.mlm,
                       i.mix = i.mix, i.mlm = i.mlm,
                       d.mix = d.mix, d.mlm = d.mlm,
                       truncate = truncate, max.support = max.support,
                       pstr.mix = pstr.mix[!lhs],
                       pstr.mlm = pstr.mlm[!lhs, , drop = FALSE],
                       pobs.mix = pobs.mix[!lhs],
                       pobs.mlm = pobs.mlm[!lhs, , drop = FALSE],
                       pdip.mix = pdip.mix[!lhs],
                       pdip.mlm = pdip.mlm[!lhs, , drop = FALSE],
                       size.a = size.a[!lhs], size.i = size.i[!lhs],
                       size.d = size.d[!lhs],
                       munb.a = munb.a[!lhs], munb.i = munb.i[!lhs],
                       munb.d = munb.d[!lhs],
                       byrow.aid = FALSE),
             faa + 1, faa)
    ans[!lhs] <- tmp
  }  # any(!lhs)


  if (ltrunc)
    while (any(vecTF <- !bad & ans %in% truncate))
      ans[vecTF] <- 1 + ans[vecTF]



  vecTF <- !bad0 & !is.na(p) &
           p <= dgaitdnbinom(min.support.use, size.p = size.p,
                     munb.p = munb.p,
                     a.mix = a.mix, a.mlm = a.mlm,
                     i.mix = i.mix, i.mlm = i.mlm,
                     d.mix = d.mix, d.mlm = d.mlm,
                     truncate = truncate, max.support = max.support,
                     pobs.mix = pobs.mix,
                     pstr.mix = pstr.mix,
                     pdip.mix = pdip.mix,
                     pobs.mlm = pobs.mlm,
                     pstr.mlm = pstr.mlm,
                     pdip.mlm = pdip.mlm,
                     size.a = size.a, size.i = size.i, size.d = size.d,
                     munb.a = munb.a, munb.i = munb.i, munb.d = munb.d,
                     byrow.aid = FALSE)
  ans[vecTF] <- min.support.use

  ans[!bad0 & !is.na(p) & p == 0] <- min.support.use
  ans[!bad0 & !is.na(p) & p == 1] <- max.support  # Inf
  ans[!bad0 & !is.na(p) & p <  0] <- NaN
  ans[!bad0 & !is.na(p) & p >  1] <- NaN
  ans[ bad0] <- NaN
  ans
}  # qgaitdnbinom





 rgaitdnbinom <-
  function(n, size.p,  # prob.p = NULL,
           munb.p,  # = NULL,
           a.mix = NULL,
           a.mlm = NULL,
           i.mix = NULL,
           i.mlm = NULL,
           d.mix = NULL,
           d.mlm = NULL,
           truncate = NULL, max.support = Inf,
           pobs.mix = 0,  # vector
           pobs.mlm = 0,  # matrix
           pstr.mix = 0,  # vector
           pstr.mlm = 0,  # matrix
           pdip.mix = 0,  # vector
           pdip.mlm = 0,  # matrix
           byrow.aid = FALSE,
           size.a = size.p, size.i = size.p, size.d = size.p,
           munb.a = munb.p, munb.i = munb.p, munb.d = munb.p) {
    qgaitdnbinom(runif(n),
              size.p = size.p,  # prob.p = prob.p,
              munb.p = munb.p,
              a.mix = a.mix,
              a.mlm = a.mlm,
              i.mix = i.mix,
              i.mlm = i.mlm,
              d.mix = d.mix,
              d.mlm = d.mlm,
              truncate = truncate, max.support = max.support,
              pobs.mix = pobs.mix,
              pobs.mlm = pobs.mlm,
              pstr.mix = pstr.mix,
              pstr.mlm = pstr.mlm,
              pdip.mix = pdip.mix,
              pdip.mlm = pdip.mlm,
              size.a = size.a, munb.a = munb.a,
              size.i = size.i, munb.i = munb.i,
              size.d = size.d, munb.d = munb.d,
              byrow.aid = byrow.aid)
}  # rgaitdnbinom






 moments.gaitdcombo.nbinom <-
  function(size.p, munb.p,
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
           munb.a = munb.p,
           munb.i = munb.p,
           munb.d = munb.p,
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
           theta1.p = size.p, theta2.p = munb.p,
           a.mix = a.mix, a.mlm = a.mlm,
           i.mix = i.mix, i.mlm = i.mlm,
           d.mix = d.mix, d.mlm = d.mlm,
           truncate = truncate, max.support = max.support,
           pobs.mix = pobs.mix, pobs.mlm = pobs.mlm,
           pstr.mix = pstr.mix, pstr.mlm = pstr.mlm,
           pdip.mix = pdip.mix, pdip.mlm = pdip.mlm,
           byrow.aid = byrow.aid,  # type.fitted = type.fitted,
           theta1.a = size.a, theta2.a = munb.a,
           theta1.i = size.i, theta2.i = munb.i,
           theta1.d = size.d, theta2.d = munb.d,
           moments2 = moments2,
           rmlife1 = rmlife1, rmlife2 = rmlife2,
           dfun = "dgaitdnbinom")  # do.call() called.


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
}  # moments.gaitdcombo.nbinom










 moments.gaitdcombo.2par <-
  function(theta1.p, theta2.p,
           a.mix = NULL, a.mlm = NULL,
           i.mix = NULL, i.mlm = NULL,
           d.mix = NULL, d.mlm = NULL,
           truncate = NULL, max.support = Inf,
           pobs.mix = 0, pobs.mlm = 0,  # Vector and matrix resp.
           pstr.mix = 0, pstr.mlm = 0,  # Ditto
           pdip.mix = 0, pdip.mlm = 0,  # Ditto
           byrow.aid = FALSE,  # For pobs.mlm and pstr.mlm
           theta1.a = theta1.p,
           theta1.i = theta1.p,
           theta1.d = theta1.p,
           theta2.a = theta2.p,
           theta2.i = theta2.p,
           theta2.d = theta2.p,
           moments2 = FALSE,  # Use this for variances.
           rmlife1 = 0, rmlife2 = 0,
           dfun = "dgenpois1") {  # "dgaitdnbinom" used

      
  NOS <- 1
  nnn <- length(theta1.p)
  cdf.max.s0 <- pnbinom(max.support, size = theta1.p, mu = theta2.p)
  pfun <- dfun
  substring(pfun, 1) <- "p"  # Replace the "d" by a "p"
  cdf.max.s <- do.call(pfun, list(max.support, theta1.p, theta2.p))
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
      pmf.p <- do.call(dfun, list(tval, theta1.p, theta2.p))
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
  SumA2.mix.p <- SumA2.mlm.p <-
  SumA2.mix.a <- SumA2.mlm.a <-
  SumA0.mix.x <- SumA0.mlm.x <-
  SumA1.mix.x <- SumA1.mlm.x <-
  SumA2.mix.x <- SumA2.mlm.x <- matrix(0, nnn, NOS)
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
      pmf.p <- do.call(dfun, list(aval, theta1.p, theta2.p))
      pmf.a <- do.call(dfun, list(aval, theta1.a, theta2.a))
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
      pmf.p <- do.call(dfun, list(aval, theta1.p, theta2.p))
      pmf.a <- do.call(dfun, list(aval, theta1.a, theta2.a))
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
  SumI2.mix.p <- SumI2.mlm.p <-
  SumI2.mix.i <- SumI2.mlm.i <-
  SumI0.mix.x <- SumI0.mlm.x <-
  SumI1.mix.x <- SumI1.mlm.x <-
  SumI2.mix.x <- SumI2.mlm.x <- matrix(0, nnn, NOS)



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
      pmf.p <- do.call(dfun, list(ival, theta1.p, theta2.p))
      pmf.i <- do.call(dfun, list(ival, theta1.i, theta2.i))
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
      pmf.p <- do.call(dfun, list(ival, theta1.p, theta2.p))
      pmf.i <- do.call(dfun, list(ival, theta1.i, theta2.i))
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
      pmf.p <- do.call(dfun, list(dval, theta1.p, theta2.p))
      pmf.d <- do.call(dfun, list(dval, theta1.d, theta2.d))
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
      pmf.p <- do.call(dfun, list(dval, theta1.p, theta2.p))
      pmf.d <- do.call(dfun, list(dval, theta1.d, theta2.d))
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
  ans <- list('cdf.max.s'    = cdf.max.s,
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
              'SumD0.mix.d' = SumD0.mix.d,  #
              'SumD0.mix.p' = SumD0.mix.p,
              'SumD1.mix.d' = SumD1.mix.d,
              'SumD1.mix.p' = SumD1.mix.p,
              'SumD0.mlm.d' = SumD0.mlm.d,
              'SumD0.mlm.p' = SumD0.mlm.p,
              'SumD1.mlm.d' = SumD1.mlm.d,
              'SumD1.mlm.p' = SumD1.mlm.p,  #
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
              'dprd2.mix'   = dprd2.mix,
              'dprd2.mlm'   = dprd2.mlm,
              'SumT2.p'     = SumT2.p,
              'SumA2.mix.p' = SumA2.mix.p,
              'SumA2.mix.a' = SumA2.mix.a,
              'SumI2.mix.p' = SumI2.mix.p,
              'SumI2.mix.i' = SumI2.mix.i,
              'SumD2.mix.p' = SumD2.mix.p,
              'SumD2.mix.d' = SumD2.mix.d,
              'SumA2.mlm.p' = SumA2.mlm.p,
              'SumA2.mlm.a' = SumA2.mlm.a,
              'SumI2.mlm.p' = SumI2.mlm.p,
              'SumI2.mlm.i' = SumI2.mlm.i,
              'SumD2.mlm.p' = SumD2.mlm.p,    #
              'SumD2.mlm.d' = SumD2.mlm.d))   #
  }
  ans
}  # moments.gaitdcombo.2par






















 gaitdnbinomial <-
  function(a.mix = NULL, i.mix = NULL, d.mix = NULL,
           a.mlm = NULL, i.mlm = NULL, d.mlm = NULL,
           truncate = NULL,    #  max.support = Inf,
           zero = c("size", "pobs", "pstr", "pdip"),
           eq.ap = TRUE, eq.ip = TRUE, eq.dp = TRUE,
           parallel.a = FALSE, parallel.i = FALSE, parallel.d = FALSE,
           lmunb.p = "loglink",
           lmunb.a = lmunb.p, lmunb.i = lmunb.p, lmunb.d = lmunb.p,
           lsize.p = "loglink",
           lsize.a = lsize.p, lsize.i = lsize.p, lsize.d = lsize.p,
           type.fitted = c("mean", "munbs", "sizes",
                           "pobs.mlm", "pstr.mlm", "pdip.mlm",
                           "pobs.mix", "pstr.mix", "pdip.mix",
                           "Pobs.mix", "Pstr.mix", "Pdip.mix",
                           "nonspecial", "Numer", "Denom.p",
                           "sum.mlm.i", "sum.mix.i",
                           "sum.mlm.d", "sum.mix.d",
                           "ptrunc.p", "cdf.max.s"),
           gpstr.mix = ppoints(7) / 3,  # ppoints(9) / 2,
           gpstr.mlm = ppoints(7) / (3 + length(i.mlm)),
           imethod = 1,
           mux.init = c(0.75, 0.5, 0.75, 0.5),  # Order is A, I, D, size
           imunb.p = NULL, imunb.a = imunb.p,
           imunb.i = imunb.p, imunb.d = imunb.p,
           isize.p = NULL,    # NULL, 1 is easy but inflexible
           isize.a = isize.p,   # NULL,    # isize.p,   # 
           isize.i = isize.p, isize.d = isize.p,
           ipobs.mix = NULL, ipstr.mix = NULL,  # 0.25, 
           ipdip.mix = NULL,   # 0.01,    # Easy but inflexible 0.01
           ipobs.mlm = NULL,
           ipstr.mlm = NULL,  # 0.25, 
           ipdip.mlm = NULL,   # 0.01,   # NULL, Easy but inflexible
           byrow.aid = FALSE,
           ishrinkage = 0.95,
           probs.y = 0.35,
           nsimEIM = 500, cutoff.prob = 0.999,  # Maxiter = 5000,
           eps.trig = 1e-7,
           nbd.max.support = 4000,
           max.chunk.MB = 30) {  # max.memory = Inf is allowed


  mux.init <- rep_len(mux.init, 4)
  if (length(a.mix) == 0) a.mix <- NULL
  if (length(i.mix) == 0) i.mix <- NULL
  if (length(d.mix) == 0) d.mix <- NULL
  if (length(a.mlm) == 0) a.mlm <- NULL
  if (length(i.mlm) == 0) i.mlm <- NULL
  if (length(d.mlm) == 0) d.mlm <- NULL
  if (length(truncate) == 0) truncate <- NULL





  max.support <- Inf  # Currently, temporary measure?
  lowsup <- 0
  gaitd.errorcheck(a.mix, a.mlm, i.mix, i.mlm,
                   d.mix, d.mlm, truncate, max.support,
                   min.support = lowsup, nparams = 2)
  la.mix <- length(a.mix <- sort(a.mix))
  li.mix <- length(i.mix <- sort(i.mix))
  ld.mix <- length(d.mix <- sort(d.mix))
  la.mlm <- length(a.mlm)
  li.mlm <- length(i.mlm)
  ld.mlm <- length(d.mlm)
  ltruncat <- length(truncate <- sort(truncate))
  ltrunc.use <- ltruncat > 0 || !is.infinite(max.support) 

  lmunb.p <- as.list(substitute(lmunb.p))
  emunb.p <- link2list(lmunb.p)
  lmunb.p <- attr(emunb.p, "function.name")
  lmunb.p.save <- lmunb.p

  lsize.p <- as.list(substitute(lsize.p))
  esize.p <- link2list(lsize.p)
  lsize.p <- attr(esize.p, "function.name")
  lsize.p.save <- lsize.p

  lpobs.mix <- "multilogitlink"  # \omega_p
  epobs.mix <- list()  # zz NULL for now 20200907 coz 'multilogitlink'
  emunb.a <- link2list(lmunb.a)
  lmunb.a <- attr(emunb.a, "function.name")
  esize.a <- link2list(lsize.a)
  lsize.a <- attr(esize.a, "function.name")

  lpstr.mix <- "multilogitlink"  # \phi_p
  epstr.mix <- list()  # zz NULL for now 20200907 coz 'multilogitlink'
  lpdip.mix <- "multilogitlink"  # zz unsure 20211002
  epdip.mix <- list()  # zz unsure 20211002
  emunb.i <- link2list(lmunb.i)
  lmunb.i <- attr(emunb.i, "function.name")
  esize.i <- link2list(lsize.i)
  lsize.i <- attr(esize.i, "function.name")

  emunb.d <- link2list(lmunb.d)
  lmunb.d <- attr(emunb.d, "function.name")
  esize.d <- link2list(lsize.d)
  lsize.d <- attr(esize.d, "function.name")


  if (is.vector(zero) && is.character(zero) && length(zero) == 4) {
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
           negbinomial(lmu = .lmunb.p.save , lsize = .lsize.p.save ,
                       zero = NULL),
           list( .lmunb.p.save = lmunb.p.save,
                 .lsize.p.save = lsize.p.save))))

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
    stop("<= one unstructured altered value (no 'munb.a')",
         ", so setting 'eq.ap = TRUE' is meaningless")
  if (li.mix <= 1 && eq.ip)
    stop("<= one unstructured inflated value (no 'munb.i')",
            ", so setting 'eq.ip = TRUE' is meaningless")
  if (ld.mix <= 1 && eq.dp)
    stop("<= one unstructured deflated value (no 'munb.d')",
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
            c("mean", "munbs", "sizes",
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
  tmp3 <- c(munb.p   = lmunb.p,
            size.p   = lsize.p,
            pobs.mix = if (la.mix) "multilogitlink" else NULL,
            munb.a   = if (la.mix > 1) lmunb.a else NULL,
            size.a   = if (la.mix > 1) lsize.a else NULL,
            pstr.mix = if (li.mix) "multilogitlink" else NULL,
            munb.i   = if (li.mix > 1) lmunb.i else NULL,
            size.i   = if (li.mix > 1) lsize.i else NULL,
            pdip.mix = if (ld.mix) "multilogitlink" else NULL,
            munb.d   = if (ld.mix > 1) lmunb.d else NULL,
            size.d   = if (ld.mix > 1) lsize.d else NULL,
            if (la.mlm) rep("multilogitlink", la.mlm) else NULL,
            if (li.mlm) rep("multilogitlink", li.mlm) else NULL,
            if (ld.mlm) rep("multilogitlink", ld.mlm) else NULL)
  Ltmp3 <- length(tmp3) 
  if (la.mlm + li.mlm + ld.mlm)
    names(tmp3)[(Ltmp3 - la.mlm - li.mlm - ld.mlm + 1):Ltmp3] <-
      c(tmp7a, tmp7b, tmp7c)
  par1or2 <- 2  # 1
  tmp3.TF <- c(rep(TRUE, par1or2),
               la.mix > 0, rep(la.mix > 1, par1or2),
               li.mix > 0, rep(li.mix > 1, par1or2),
               ld.mix > 0, rep(ld.mix > 1, par1or2),
               la.mlm > 0, li.mlm > 0, ld.mlm > 0)
  indeta.finish <- cumsum(c(rep(1, par1or2),
                            1, rep(1, par1or2),
                            1, rep(1, par1or2),
                            1, rep(1, par1or2),
                            la.mlm, li.mlm, ld.mlm,
                            ld.mlm + 1) * c(tmp3.TF, 1))
  indeta.launch <- c(1, 1 + head(indeta.finish, -1))

  indeta.launch <- head(indeta.launch, -1)
  indeta.finish <- head(indeta.finish, -1)
  indeta.launch[!tmp3.TF] <- NA  # Not to be accessed
  indeta.finish[!tmp3.TF] <- NA  # Not to be accessed
  indeta <- cbind(launch = indeta.launch,
                  finish = indeta.finish)
  rownames(indeta) <- c("munb.p", "size.p",
                        "pobs.mix", "munb.a", "size.a",
                        "pstr.mix", "munb.i", "size.i",
                        "pdip.mix", "munb.d", "size.d",
                        "pobs.mlm", "pstr.mlm", "pdip.mlm")
  M1 <- max(indeta, na.rm = TRUE)
  predictors.names <- tmp3  # Passed into @infos and @initialize.


  blurb1 <- ""
  if (la.mlm + la.mix) blurb1 <- "Generally-altered "
  if (li.mlm + li.mix) blurb1 <- "Generally-inflated "
  if (ltrunc.use) blurb1 <- "Generally-truncated "
  if ( (la.mlm + la.mix) &&  (li.mlm + li.mix) && !ltrunc.use)
    blurb1 <- "Generally-altered & -inflated "
  if ( (la.mlm + la.mix) && !(li.mlm + li.mix) &&  ltrunc.use)
    blurb1 <- "Generally-altered & -truncated "
  if (!(la.mlm + la.mix) &&  (li.mlm + li.mix) &&  ltrunc.use)
    blurb1 <- "Generally-inflated & -truncated "
  if ( (la.mlm + la.mix) &&  (li.mlm + li.mix) &&  ltrunc.use)
    blurb1 <- "Generally-altered, -inflated & -truncated "

  if (ld.mlm + ld.mix) blurb1 <-
    c(blurb1,
      if (la.mlm + la.mix + li.mlm + li.mix) "& " else "Generally",
      "-deflated ")



      
  new("vglmff",
  blurb = c(blurb1, "NB regression\n",
            "(GAITD-NB(munb.p, size.p)-",
                   "NB(munb.a, size.a)-MLM-",
                   "NB(munb.i, size.i)-MLM-\n",
                   "NB(munb.d, size.d)-MLM generally)\n\n",
            "Links: ",
            namesof("munb.p", lmunb.p, earg = emunb.p, tag = FALSE),
            ", ",
            namesof("size.p", lsize.p, earg = esize.p, tag = FALSE),
            if (la.mix) ", \n       ",
            if (la.mix > 0) c("multilogit(pobs.mix)"),
            if (la.mix > 1) c(", ",
            namesof("munb.a",  lmunb.a, emunb.a, tag = FALSE), ", ",
            namesof("size.a",  lsize.a, esize.a, tag = FALSE)),
            if (la.mix && li.mix) ", \n       ",
            if (li.mix > 0) c(  if (la.mix) "" else ", ",
            "multilogit(pstr.mix)"),
            if (li.mix > 1) c(", ",
            namesof("munb.i",  lmunb.i, emunb.i, tag = FALSE), ", ",
            namesof("size.i",  lsize.i, esize.i, tag = FALSE)),
            if (li.mix && ld.mix) ", \n       ",
            if (ld.mix > 0) c(  if (li.mix) "" else ", ",
            "multilogit(pdip.mix)"),
            if (ld.mix > 1) c(", ",
            namesof("munb.d",  lmunb.d, emunb.d, tag = FALSE), ", ",
            namesof("size.d",  lsize.d, esize.d, tag = FALSE)),
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
      warning("20211115; unsure; above vs. below line is right?")
      Use.mat <- use.mat.mlm <- diag(M)  # munb.p only  20211115
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


    use.mat.mix <- cm3gaitd( .eq.ap , .eq.ip , .eq.dp , npar = 2)


 

    tmp3.TF.subset <- tmp3.TF[-(12:14)]  # tmp3.TF[-(8:10)] for Poisson
    use.mat.mix <- use.mat.mix[tmp3.TF.subset, , drop = FALSE]
    notall0 <- function(x) !all(x == 0)
    use.mat.mix <- use.mat.mix[, apply(use.mat.mix, 2, notall0),
                               drop = FALSE]




    if (la.mix + li.mix + ld.mix > 0)
      Use.mat <- use.mat.mix  # Possibly all done







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
         dpqrfun = "gaitdnbinom",
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
  multipleResponses = FALSE,  # negbinomial can be called if TRUE
         parameters.names = names( .predictors.names ),
         parent.name = c("negbinomial", "nbinom"),
         type.fitted  = as.vector( .type.fitted ),
         type.fitted.choices = ( .type.fitted.choices ),
     baseparams.argnames  = c("munb", "size"),  # zz reorder this?
     MM1 = 2,  # 2 parameters for 1 response (munb & size). Needed.
         flip.args = TRUE,  # For dpqr arguments (GAITD plotting).
         zero = .zero )
  }, list( .zero = zero, .lowsup = lowsup,
           .type.fitted = type.fitted,
           .type.fitted.choices = type.fitted.choices,
           .lmunb.p = lmunb.p, .emunb.p = emunb.p,
           .lmunb.a = lmunb.a, .emunb.a = emunb.a,
           .lmunb.i = lmunb.i, .emunb.i = emunb.i,
           .lmunb.d = lmunb.d, .emunb.d = emunb.d,
           .a.mlm = a.mlm, .i.mlm = i.mlm, .d.mlm = d.mlm,
           .a.mix = a.mix, .i.mix = i.mix, .d.mix = d.mix,
           .la.mlm = la.mlm, .li.mlm = li.mlm, .ld.mlm = ld.mlm,
           .la.mix = la.mix, .li.mix = li.mix, .ld.mix = ld.mix,
           .truncate = truncate, .max.support = max.support,
           .predictors.names = predictors.names,
           .M1 = M1, .lall.len = lall.len ))),

  rqresslot = eval(substitute(
    function(mu, y, w, eta, extra = NULL) {
    if (!is.matrix(eta)) eta <- as.matrix(eta)
    la.mix <- length((a.mix <- as.vector( .a.mix )))
    li.mix <- length((i.mix <- as.vector( .i.mix )))
    ld.mix <- length((d.mix <- as.vector( .d.mix )))
    la.mlm <- length((a.mlm <- as.vector( .a.mlm )))
    li.mlm <- length((i.mlm <- as.vector( .i.mlm )))
    ld.mlm <- length((d.mlm <- as.vector( .d.mlm )))
    truncate <- as.vector( .truncate )

    tmp3.TF <- ( .tmp3.TF )   # Logical of length 10.

    lall.len <- la.mix + li.mix + ld.mix + la.mlm + li.mlm + ld.mlm
    pobs.mix <- pstr.mix <- pdip.mix <- 0  # 4 rowSums()
    pobs.mlm <- pstr.mlm <- pdip.mlm <- 0  # matrix(0, NROW(eta), 1)
    munb.p <- cbind(eta2theta(eta[, 1], .lmunb.p , .emunb.p ))
    size.p <- cbind(eta2theta(eta[, 2], .lsize.p , .esize.p ))
    ind.munb.z <- 1:2  # Points to munb.p and size.p only.
    munb.a <- munb.i <- munb.d <- munb.p
    size.a <- size.i <- size.d <- size.p

    if (any(tmp3.TF[c(4, 7, 10)])) {  # At least one munb.[aid]
      ind.munb.z <- extra$indeta[c(1:2, 4:5, 7:8, 10:11), 'launch']
      ind.munb.z <- c(na.omit(ind.munb.z))  # At least one value
      munb.a <- if (!tmp3.TF[ 4]) munb.p else
        eta2theta(eta[, extra$indeta[ 4, 1]], .lmunb.a , .emunb.a )
      munb.i <- if (!tmp3.TF[ 7]) munb.p else
        eta2theta(eta[, extra$indeta[ 7, 1]], .lmunb.i , .emunb.i )
      munb.d <- if (!tmp3.TF[10]) munb.p else
        eta2theta(eta[, extra$indeta[10, 1]], .lmunb.d , .emunb.d )

      size.a <- if (!tmp3.TF[ 5]) size.p else
        eta2theta(eta[, extra$indeta[ 5, 1]], .lsize.a , .esize.a )
      size.i <- if (!tmp3.TF[ 8]) size.p else
        eta2theta(eta[, extra$indeta[ 8, 1]], .lsize.i , .esize.i )
      size.d <- if (!tmp3.TF[11]) size.p else
        eta2theta(eta[, extra$indeta[11, 1]], .lsize.d , .esize.d )
    }  # la.mix + li.mix + ld.mix > 0

    if (lall.len) {  # An MLM was fitted
      allprobs <-
        multilogitlink(eta[, -ind.munb.z, drop = FALSE],
                       refLevel = "(Last)",  # Make sure
                       inverse = TRUE)  # rowSums == 1
      if (anyNA(allprobs))
        warning("there are NAs here in slot linkinv")
      if (min(allprobs) == 0 || max(allprobs) == 1)
        warning("fitted probabilities numerically 0 or 1 occurred")

      Nextone <- 0  # Might not be used actually; 0, not 1
      if (tmp3.TF[ 3])
        pobs.mix <- allprobs[, (Nextone <- Nextone + 1)]
      if (tmp3.TF[ 6])
        pstr.mix <- allprobs[, (Nextone <- Nextone + 1)]
      if (tmp3.TF[ 9])
        pdip.mix <- allprobs[, (Nextone <- Nextone + 1)]

      if (tmp3.TF[12]) {
        ind8 <- (Nextone + 1):(Nextone + la.mlm)
        pobs.mlm <- allprobs[, ind8, drop = FALSE]
        dimnames(pobs.mlm) <- list(rownames(eta),
                                   as.character(a.mlm))
        Nextone <- Nextone + la.mlm
      }
      if (tmp3.TF[13]) {
        ind9 <- (Nextone + 1):(Nextone + li.mlm)
        pstr.mlm <- allprobs[, ind9, drop = FALSE]
        dimnames(pstr.mlm) <- list(rownames(eta),
                                   as.character(i.mlm))
        Nextone <- Nextone + li.mlm
      }
      if (tmp3.TF[14]) {
        ind10 <- (Nextone + 1):(Nextone + ld.mlm)
        pdip.mlm <- allprobs[, ind10, drop = FALSE]
        dimnames(pdip.mlm) <- list(rownames(eta),
                                   as.character(d.mlm))
        Nextone <- Nextone + ld.mlm  # Not needed
      }
    }  # lall.len


    scrambleseed <- runif(1)  # To scramble the seed
    qnorm(runif(length(y),
        pgaitdnbinom(y - 1, munb.p = munb.p, size.p = size.p,
                  a.mix = a.mix, i.mix = i.mix, d.mix = d.mix,
                  a.mlm = a.mlm, i.mlm = i.mlm, d.mlm = d.mlm,
                  truncate = truncate,
                  max.support = as.vector( .max.support ),
                  munb.a = munb.a, munb.i = munb.i,
                  munb.d = munb.d,
                  size.a = size.a, size.i = size.i,
                  size.d = size.d,
                  pobs.mix = pobs.mix, pstr.mix = pstr.mix,
                  pdip.mix = pdip.mix,
                  pobs.mlm = pobs.mlm, pstr.mlm = pstr.mlm,
                  pdip.mlm = pdip.mlm),
        pgaitdnbinom(y    , munb.p = munb.p, size.p = size.p,
                  a.mix = a.mix, i.mix = i.mix, d.mix = d.mix,
                  a.mlm = a.mlm, i.mlm = i.mlm, d.mlm = d.mlm,
                  truncate = truncate,
                  max.support = as.vector( .max.support ),
                  munb.a = munb.a, munb.i = munb.i,
                  munb.d = munb.d,
                  size.a = size.a, size.i = size.i,
                  size.d = size.d,
                  pobs.mix = pobs.mix, pstr.mix = pstr.mix,
                  pdip.mix = pdip.mix,
                  pobs.mlm = pobs.mlm, pstr.mlm = pstr.mlm,
                  pdip.mlm = pdip.mlm)))
  }, list(
    .lmunb.p = lmunb.p, .emunb.p = emunb.p,
    .lmunb.a = lmunb.a, .emunb.a = emunb.a,
    .lmunb.i = lmunb.i, .emunb.i = emunb.i,
    .lmunb.d = lmunb.d, .emunb.d = emunb.d,
    .lsize.p = lsize.p, .esize.p = esize.p,
    .lsize.a = lsize.a, .esize.a = esize.a,
    .lsize.i = lsize.i, .esize.i = esize.i,
    .lsize.d = lsize.d, .esize.d = esize.d,
    .lpstr.mix = lpstr.mix, .lpobs.mix = lpobs.mix,
    .lpdip.mix = lpdip.mix,
    .epstr.mix = epstr.mix, .epobs.mix = epobs.mix,
    .epdip.mix = epdip.mix,
    .tmp3.TF = tmp3.TF,
    .a.mix = a.mix, .i.mix = i.mix, .d.mix = d.mix,
    .a.mlm = a.mlm, .i.mlm = i.mlm, .d.mlm = d.mlm,
    .truncate = truncate, .max.support = max.support ))),

  initialize = eval(substitute(expression({
    extra$indeta <- ( .indeta )  # Avoids recomputing it
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
                               max.support = .max.support )
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
    extra$control.trace <- control$trace  # For summary(fit) postfit.



    predictors.names <- ( .predictors.names )  # Got it, named


    if (!length(etastart)) {
      input.size.p <-
      init.size.a <- init.size.i <- init.size.d <- 1  # Needed
      init.munb.a <- init.munb.i <- init.munb.d <-  # Needed
      init.munb.p <- if (length( .imunb.p )) ( .imunb.p ) else
              Init.mu(y = y, w = w, imethod = .imethod ,
                             imu = .imunb.p ,  # x = x,
                             ishrinkage = .ishrinkage ,
                             probs.y = .probs.y )
      etastart <- matrix(nrow = n, ncol = M,
        theta2eta(init.munb.p, .lmunb.p , earg = .emunb.p ))




      Mom.nb.init <-
        function(y, TFvec = TRUE,
                 w = rep_len(1, length(y)), mux.arg = 10) {
          munb.z.init <- weighted.mean(y[TFvec], w = w[TFvec])
          var.y.wted <- cov.wt(cbind(y[TFvec]), w = w[TFvec])$cov
          size.z.init <- munb.z.init^2 / (var.y.wted - munb.z.init)
          if (!is.Numeric(size.z.init, positive = TRUE))
            size.z.init <- mux.arg * munb.z.init
        c(munb = munb.z.init, size = size.z.init)
      }  # Mom.nb.init



      mux.more.k <- extra$mux.init[4]  # 0.5 and 0.25 seem good.
      if (tmp3.TF[ 2]) {
        init.size.p <-
          if (length( .isize.p )) ( .isize.p ) else {
            keep.Mom.p.mix <- Mom.nb.init(y, w = w)  # TFvec = is.ns
            keep.Mom.p.mix["size"] * mux.more.k
          }
        etastart[, extra$indeta[ 2, 1]] <-
          theta2eta(init.size.p, .lsize.p , earg = .esize.p )
      }  # tmp3.TF[ 2]


      if (tmp3.TF[ 5]) {
        init.size.a <-
          if (length( .isize.a )) ( .isize.a ) else {
            keep.Mom.a.mix <-  # Useful later for the mean.
              Mom.nb.init(y, rowSums(extra$skip.mix.a) > 0, w = w)
            0.5 * (keep.Mom.a.mix["size"] +
                   keep.Mom.p.mix["size"]) * mux.more.k
          }
        etastart[, extra$indeta[ 5, 1]] <-
          theta2eta(init.size.a, .lsize.a , earg = .esize.a )
      }  # tmp3.TF[ 5]



      if (tmp3.TF[ 8]) {
        init.size.i <-
          if (length( .isize.i )) ( .isize.i ) else {
            keep.Mom.i.mix <-  # Useful later for the mean.
              Mom.nb.init(y, rowSums(extra$skip.mix.i) > 0, w = w)
            0.5 * (keep.Mom.i.mix["size"] +
                   keep.Mom.p.mix["size"]) * mux.more.k
          }
        etastart[, extra$indeta[ 8, 1]] <-
          theta2eta(init.size.i, .lsize.i , earg = .esize.i )
      }  # tmp3.TF[ 8]

      if (tmp3.TF[11]) {
        init.size.d <-
          if (length( .isize.d )) ( .isize.d ) else {
            keep.Mom.d.mix <-  # Useful later for the mean.
              Mom.nb.init(y, rowSums(extra$skip.mix.d) > 0, w = w)
            0.5 * (keep.Mom.d.mix["size"] +
                   keep.Mom.p.mix["size"]) * mux.more.k
          }
        etastart[, extra$indeta[11, 1]] <-
          theta2eta(init.size.d, .lsize.d , earg = .esize.d )
      }  # tmp3.TF[11]


 
      mux.more.a <- extra$mux.init[1]   # 0.75 Err to slightly smaller
      init.pobs.mix <- numeric(n)
      if (tmp3.TF[ 3]) {  # la.mix > 0
        init.pobs.mix <- if (length( .ipobs.mix )) {
          rep_len( .ipobs.mix , n)
        } else {
          is.a.mix1 <- rowSums(extra$skip.mix.a) > 0
          rep_len(mux.more.a * sum(w[is.a.mix1]) / sum(w), n) 
        }
      }  # la.mix > 0


      if (tmp3.TF[ 4]) {  # Assign coln 3; la.mix > 1
        init.munb.a <- if (length( .imunb.a ))
          rep_len( .imunb.a , n) else {
          if ( .eq.ap ) init.munb.p else
            rep_len(0.5 * (keep.Mom.a.mix["munb"] +
                           keep.Mom.p.mix["munb"]), n)
        }
        etastart[, 3] <-
          theta2eta(init.munb.a, .lmunb.a , earg = .emunb.a )
      }  # tmp3.TF[ 4]


      init.pstr.mix <- init.pdip.mix <- numeric(n)
      try.gridsearch.pstr.mix <- FALSE
      if (tmp3.TF[ 6]) {  # li.mix > 0
        init.pstr.mix <- if (length( .ipstr.mix )) {
          rep_len( .ipstr.mix , n)
        } else {
          try.gridsearch.pstr.mix <- TRUE
          numeric(n)  # Overwritten by gridsearch
        }
      }  # li.mix > 0

      if (tmp3.TF[ 7]) {  # li.mix > 1
        init.munb.i <- if (length( .imunb.i ))
          rep_len( .imunb.i , n) else {
          if ( .eq.ip ) init.munb.p else
            rep_len( keep.Mom.i.mix["munb"], n)
        }
        etastart[, (extra$indeta[ 7, 'launch'])] <-
          theta2eta(init.munb.i, .lmunb.i , earg = .emunb.i )
      }  # li.mix > 1


      if (tmp3.TF[12]) {  #  la.mlm
        init.pobs.mlm <- if (length( .ipobs.mlm )) {
          matrix( .ipobs.mlm , n, la.mlm, byrow = .byrow.aid )
        } else {
          mux.more.a <- extra$mux.init[1]                            
          init.pobs.mlm <- colSums(c(w) *
                                   extra$skip.mlm.a) / colSums(w)
          init.pobs.mlm <- init.pobs.mlm * as.vector( mux.more.a )
          matrix(init.pobs.mlm, n, la.mlm, byrow = TRUE)
        }
      } else {
        init.pobs.mlm <- matrix(0, n, 1)
      }


      try.gridsearch.pstr.mlm <- FALSE
      if (tmp3.TF[13]) {  #  li.mlm
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





      gaitNBD.Loglikfun1.mix <-
        function(pstr.mix.val, y, x, w, extraargs) {
    sum(c(w) *
        dgaitdnbinom(y, pstr.mix = pstr.mix.val,
                  pstr.mlm    = extraargs$pstr.mlm,  # Differs here
                  munb.p      = extraargs$munb.p,
                  munb.a      = extraargs$munb.a,
                  munb.i      = extraargs$munb.i,
                  munb.d      = extraargs$munb.d,
                  size.p      = extraargs$size.p,
                  size.a      = extraargs$size.a,
                  size.i      = extraargs$size.i,
                  size.d      = extraargs$size.d,
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

 gaitNBD.Loglikfun1.mlm <-
     function(pstr.mlm.val, y, x, w, extraargs) {
    sum(c(w) *
        dgaitdnbinom(y, pstr.mlm = pstr.mlm.val,
                   pstr.mix    = extraargs$pstr.mix,  # Differs here
                   munb.p      = extraargs$munb.p,
                   munb.a      = extraargs$munb.a,
                   munb.i      = extraargs$munb.i,
                   munb.d      = extraargs$munb.d,
                   size.p      = extraargs$size.p,
                   size.a      = extraargs$size.a,
                   size.i      = extraargs$size.i,
                   size.d      = extraargs$size.d,
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

 gaitNBD.Loglikfun2 <-
     function(pstr.mix.val, pstr.mlm.val, y, x, w, extraargs) {
    sum(c(w) *
        dgaitdnbinom(y, pstr.mix = pstr.mix.val,
                   pstr.mlm    = pstr.mlm.val,
                   munb.p      = extraargs$munb.p,
                   munb.a      = extraargs$munb.a,
                   munb.i      = extraargs$munb.i,
                   munb.d      = extraargs$munb.d,
                   size.p      = extraargs$size.p,
                   size.a      = extraargs$size.a,
                   size.i      = extraargs$size.i,
                   size.d      = extraargs$size.d,
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
              munb.p      = init.munb.p,
              munb.a      = init.munb.a,
              munb.i      = init.munb.i,
              munb.d      = init.munb.d,
              size.p      = init.size.p,   # .isize.p ,
              size.a      = init.size.a,   # .isize.a ,
              size.i      = init.size.i,   # .isize.i ,
              size.d      = init.size.d,   # .isize.d ,
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
            grid.search2( .gpstr.mix , .gpstr.mlm ,
                         objfun = gaitNBD.Loglikfun2,
                         y = y, w = w, extraargs = extraargs,
                         ret.objfun = TRUE)
          } else if (try.gridsearch.pstr.mix) {
            extraargs$pstr.mlm <- init.pstr.mlm
            grid.search ( .gpstr.mix ,
                         objfun = gaitNBD.Loglikfun1.mix,
                         y = y, w = w, extraargs = extraargs,
                         ret.objfun = TRUE)
          } else if (try.gridsearch.pstr.mlm) {
            extraargs$pstr.mix <- init.pstr.mix
            grid.search ( .gpstr.mlm ,
                         objfun = gaitNBD.Loglikfun1.mlm,
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
      if (FALSE)
      Numer.init2 <- 1 - rowSums(init.pobs.mlm) -
                         rowSums(init.pstr.mlm) +
                         rowSums(init.pdip.mlm) -
                         init.pobs.mix - init.pstr.mix +
                         init.pdip.mix  # Same as 'Numer'.
        
      etastart.z <- if (lall.len == 0) NULL else {
        tmp.mat <- cbind(if (tmp3.TF[ 3]) init.pobs.mix else NULL,
                         if (tmp3.TF[ 6]) init.pstr.mix else NULL,
                         if (tmp3.TF[ 9]) init.pdip.mix else NULL,
                         if (tmp3.TF[12]) init.pobs.mlm else NULL,
                         if (tmp3.TF[13]) init.pstr.mlm else NULL,
                         if (tmp3.TF[14]) init.pdip.mlm else NULL,
                         Numer.init1)  # Numer.init1  # Numer.init2
        multilogitlink(tmp.mat)
      }  # etastart.z
      if (!is.matrix(etastart.z)) etastart.z <- cbind(etastart.z)



      nextone <- 1  # Might not be used actually
      if (tmp3.TF[ 3]) {
        etastart[, 3] <- etastart.z[, nextone]
        nextone <- nextone + 1
      }
      if (tmp3.TF[ 6]) {  # Coln 3 or 6
        etastart[, (extra$indeta[ 6, 1])] <-
            etastart.z[, nextone]
        nextone <- nextone + 1
      }
      if (tmp3.TF[ 9]) {  # Coln 3 or 6 or 9
        etastart[, (extra$indeta[ 9, 1])] <-
            etastart.z[, nextone]
        nextone <- nextone + 1
      }
      if (tmp3.TF[12]) {
        ind8 <- (extra$indeta[12, 1]):(extra$indeta[12, 2])
        etastart[, ind8] <- etastart.z[, nextone:(nextone +
                                                  la.mlm - 1)]
        nextone <- nextone + la.mlm
      }
      if (tmp3.TF[13]) {
        ind9 <- (extra$indeta[13, 1]):(extra$indeta[13, 2])
        etastart[, ind9] <- etastart.z[, nextone:(nextone +
                                                  li.mlm - 1)]
        nextone <- nextone + li.mlm
      }
      if (tmp3.TF[14]) {
        ind0 <- (extra$indeta[14, 1]):(extra$indeta[14, 2])
        etastart[, ind0] <- etastart.z[, nextone:(nextone +
                                                  ld.mlm - 1)]
        if (ncol(etastart.z) != nextone + ld.mlm - 1)
          stop("miscalculation")
      }
    }
  }), list(
    .lmunb.p = lmunb.p, .emunb.p = emunb.p,
    .lmunb.a = lmunb.a, .emunb.a = emunb.a,
    .lmunb.i = lmunb.i, .emunb.i = emunb.i,
    .lmunb.d = lmunb.d, .emunb.d = emunb.d,
    .lsize.p = lsize.p, .esize.p = esize.p,
    .lsize.a = lsize.a, .esize.a = esize.a,
    .lsize.i = lsize.i, .esize.i = esize.i,
    .lsize.d = lsize.d, .esize.d = esize.d,
    .imunb.p = imunb.p, .isize.p = isize.p,
    .imunb.a = imunb.a, .isize.a = isize.a,
    .imunb.i = imunb.i, .isize.i = isize.i,
    .imunb.d = imunb.d, .isize.d = isize.d,
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
    .tmp3.TF = tmp3.TF,  # .tmp3 = tmp3,
    .predictors.names = predictors.names,
    .mux.init = mux.init,
    .gpstr.mix = gpstr.mix,   # .gpdip.mix = gpdip.mix,
    .gpstr.mlm = gpstr.mlm,   # .gpdip.mlm = gpdip.mlm,
    .ishrinkage = ishrinkage, .probs.y = probs.y,
    .indeta = indeta,
    .eq.ap = eq.ap, .eq.ip = eq.ip, .eq.dp = eq.dp,
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
                c("mean", "munbs", "sizes",
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
    munb.p <- cbind(eta2theta(eta[, 1], .lmunb.p , .emunb.p ))
    size.p <- cbind(eta2theta(eta[, 2], .lsize.p , .esize.p ))
    ind.munb.z <- 1:2  # Points to munb.p and size.p only.
    munb.a <- munb.i <- munb.d <-
              munb.p  # Needed; answer not corrupted
    size.a <- size.i <- size.d <- size.p
    tmp3.TF <- ( .tmp3.TF )   # Logical of length 14.

    if (any(tmp3.TF[c(4, 7, 10)])) {  # At least one munb.[aid]
      ind.munb.z <- extra$indeta[c(1:2, 4:5, 7:8, 10:11), 'launch']
      ind.munb.z <- c(na.omit(ind.munb.z))  # At least one value
      munb.a <- if (!tmp3.TF[ 4]) munb.p else
        eta2theta(eta[, extra$indeta[ 4, 1]], .lmunb.a , .emunb.a )
      munb.i <- if (!tmp3.TF[ 7]) munb.p else
        eta2theta(eta[, extra$indeta[ 7, 1]], .lmunb.i , .emunb.i )
      munb.d <- if (!tmp3.TF[10]) munb.p else
        eta2theta(eta[, extra$indeta[10, 1]], .lmunb.d , .emunb.d )

      size.a <- if (!tmp3.TF[ 5]) size.p else
        eta2theta(eta[, extra$indeta[ 5, 1]], .lsize.a , .esize.a )
      size.i <- if (!tmp3.TF[ 8]) size.p else
        eta2theta(eta[, extra$indeta[ 8, 1]], .lsize.i , .esize.i )
      size.d <- if (!tmp3.TF[11]) size.p else
        eta2theta(eta[, extra$indeta[11, 1]], .lsize.d , .esize.d )
    }  # la.mix + li.mix + ld.mix > 0







    if (lall.len) {  # An MLM was fitted
      allprobs <-
        multilogitlink(eta[, -ind.munb.z, drop = FALSE],
                       inverse = TRUE)  # rowSums == 1
      if (anyNA(allprobs))
        warning("there are NAs here in slot linkinv")
      if (min(allprobs) == 0 || max(allprobs) == 1)
        warning("fitted probabilities numerically 0 or 1 occurred")


      Nextone <- 0  # Might not be used actually
      if (tmp3.TF[ 3])
        pobs.mix <- allprobs[, (Nextone <- Nextone + 1)]
      if (tmp3.TF[ 6])
        pstr.mix <- allprobs[, (Nextone <- Nextone + 1)]
      if (tmp3.TF[ 9])
        pdip.mix <- allprobs[, (Nextone <- Nextone + 1)]

      if (tmp3.TF[12]) {
        ind8 <- (Nextone + 1):(Nextone + la.mlm)
        pobs.mlm <- allprobs[, ind8, drop = FALSE]
        dimnames(pobs.mlm) <- list(rownames(eta),
                                   as.character(a.mlm))
        Nextone <- Nextone + la.mlm
      }
      if (tmp3.TF[13]) {
        ind9 <- (Nextone + 1):(Nextone + li.mlm)
        pstr.mlm <- allprobs[, ind9, drop = FALSE]
        dimnames(pstr.mlm) <- list(rownames(eta),
                                   as.character(i.mlm))
        Nextone <- Nextone + li.mlm
      }
      if (tmp3.TF[14]) {
        ind10 <- (Nextone + 1):(Nextone + ld.mlm)
        pdip.mlm <- allprobs[, ind10, drop = FALSE]
        dimnames(pdip.mlm) <- list(rownames(eta),
                                   as.character(d.mlm))
        Nextone <- Nextone + ld.mlm  # Not needed
      }
    }  # lall.len

    M1 <- max(extra$indeta, na.rm = TRUE)
    NOS <- NCOL(eta) / M1

    Bits <- moments.gaitdcombo.nbinom(
              munb.p = munb.p, size.p = size.p,
              pobs.mix = pobs.mix, pobs.mlm = pobs.mlm,
              pstr.mix = pstr.mix, pstr.mlm = pstr.mlm,
              pdip.mix = pdip.mix, pdip.mlm = pdip.mlm,
              a.mix = a.mix, i.mix = i.mix, d.mix = d.mix,
              a.mlm = a.mlm, i.mlm = i.mlm, d.mlm = d.mlm,
              munb.a = munb.a, munb.i = munb.i, munb.d = munb.d,
              size.a = size.a, size.i = size.i, size.d = size.d,
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
      tmp13 <-  # dnbinom() does not retain the matrix format??
        dnbinom(x    = matrix(a.mix,  n.obs, la.mix, byrow = TRUE),
                size = matrix(size.a, n.obs, la.mix),
                mu   = matrix(munb.a, n.obs, la.mix)) / (
        c(Bits[["SumA0.mix.a"]]))
      dim(tmp13) <- c(n.obs, la.mix)
      dimnames(tmp13) <- list(rownames(eta), as.character(a.mix))
      propn.mat.a <- tmp13
    }  # la.mix

    if (li.mix && morework) {
      tmp55 <-  # dnbinom() does not retain the matrix format??
        dnbinom(x    = matrix(i.mix,  n.obs, li.mix, byrow = TRUE),
                size = matrix(size.i, n.obs, li.mix),
                mu   = matrix(munb.i, n.obs, li.mix)) / (
        c(Bits[["SumI0.mix.i"]]))
      dim(tmp55) <- c(n.obs, li.mix)
      dimnames(tmp55) <- list(rownames(eta), as.character(i.mix))
      propn.mat.i <- tmp55  # Correct dimension
    }  # li.mix


    if (ld.mix && morework) {
      tmp55 <-  # dpois() does not retain the matrix format
        dnbinom(x    = matrix(d.mix,  n.obs, ld.mix, byrow = TRUE),
                size = matrix(size.d, n.obs, ld.mix),
                mu   = matrix(munb.d, n.obs, ld.mix)) / (
        c(Bits[["SumD0.mix.d"]]))
      dim(tmp55) <- c(n.obs, ld.mix)
      dimnames(tmp55) <- list(rownames(eta), as.character(d.mix))
      propn.mat.d <- tmp55  # Correct dimension
    }  # ld.mix

    ans <- switch(type.fitted,
      "mean"       = Bits[["mean"]],  # Unconditional mean
      "munbs"      = cbind(munb.p,
                           if (tmp3.TF[ 4]) munb.a else NULL,
                           if (tmp3.TF[ 7]) munb.i else NULL,
                           if (tmp3.TF[10]) munb.d else NULL),
      "sizes"      = cbind(size.p,
                           if (tmp3.TF[ 4]) size.a else NULL,
                           if (tmp3.TF[ 7]) size.i else NULL,
                           if (tmp3.TF[10]) size.d else NULL),
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
          dnbinom(x    = matrix(i.mlm,  n.obs, li.mlm, byrow = TRUE),
                  size = matrix(size.p, n.obs, li.mlm),
                  mu   = matrix(munb.p, n.obs, li.mlm)) / Denom.p,
      "sum.mlm.d"  = -pdip.mlm + Numer *
          dnbinom(x    = matrix(d.mlm,  n.obs, ld.mlm, byrow = TRUE),
                  size = matrix(size.p, n.obs, ld.mlm),
                  mu   = matrix(munb.p, n.obs, ld.mlm)) / Denom.p,
      "sum.mix.i"  =  c(pstr.mix) * propn.mat.i + Numer *
          dnbinom(x    = matrix(i.mix,  n.obs, li.mix, byrow = TRUE),
                  size = matrix(size.p, n.obs, li.mix),
                  mu   = matrix(munb.p, n.obs, li.mix)) / Denom.p,
      "sum.mix.d"  = -c(pdip.mix) * propn.mat.d + Numer *
          dnbinom(x    = matrix(d.mix,  n.obs, ld.mix, byrow = TRUE),
                  size = matrix(size.p, n.obs, ld.mix),
                  mu   = matrix(munb.p, n.obs, ld.mix)) / Denom.p,
      "ptrunc.p"   = Bits[["SumT0.p"]] + 1 - Bits[["cdf.max.s"]],
      "cdf.max.s"  = Bits[["cdf.max.s"]])  # Pr(y <= max.support)



    ynames.pobs.mlm <- as.character(a.mlm)  # Works with NULLs
    ynames.pstr.mlm <- as.character(i.mlm)  # Works with NULLs
    ynames.pdip.mlm <- as.character(d.mlm)  # Works with NULLs
    if (length(ans))
      label.cols.y(ans, NOS = NOS, colnames.y =
      switch(type.fitted,
             "munbs"   = c("munb.p", "munb.a",  # Some colns NA
                 "munb.i", "munb.d")[(tmp3.TF[c(1, 4, 7, 10)])],
             "sizes"   = c("size.p", "size.a",  # Some colns NA
                 "size.i", "size.d")[(tmp3.TF[c(1, 4, 7, 10)])],
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
    .lmunb.p = lmunb.p, .emunb.p = emunb.p,
    .lmunb.a = lmunb.a, .emunb.a = emunb.a,
    .lmunb.i = lmunb.i, .emunb.i = emunb.i,
    .lmunb.d = lmunb.d, .emunb.d = emunb.d,
    .lsize.p = lsize.p, .esize.p = esize.p,
    .lsize.a = lsize.a, .esize.a = esize.a,
    .lsize.i = lsize.i, .esize.i = esize.i,
    .lsize.d = lsize.d, .esize.d = esize.d,
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
                earg = list()))  # This isnt perfect; info is lost
    misc$predictors.names <- predictors.names  # Useful for coef()
    misc$link <- link.names  # 
    names(misc$link) <- parameter.names  # 


    misc$earg <- vector("list", M1)
    names(misc$earg) <- names(misc$link)
    misc$earg[[1]] <- ( .emunb.p )  # First one always there
    iptr <- 1
    if (tmp3.TF[ 3])
      misc$earg[[(iptr <- iptr + 1)]] <- list()  # multilogitlink
    if (tmp3.TF[ 4])
      misc$earg[[(iptr <- iptr + 1)]] <- ( .emunb.a )
    if (tmp3.TF[ 6])
      misc$earg[[(iptr <- iptr + 1)]] <- list()  # See below
    if (tmp3.TF[ 7])
      misc$earg[[(iptr <- iptr + 1)]] <- ( .emunb.i )
    if (tmp3.TF[ 9])
      misc$earg[[(iptr <- iptr + 1)]] <- list()  # See below
    if (tmp3.TF[10])
      misc$earg[[(iptr <- iptr + 1)]] <- ( .emunb.d )
    if (tmp3.TF[12]) {  # la.mlm
      for (ii in seq(la.mlm)) {
        misc$earg[[(iptr <- iptr + 1)]] <- list()
      }  # ii
    }  # la.mlm
    if (tmp3.TF[13]) {  # li.mlm
      for (ii in seq(li.mlm)) {
        misc$earg[[(iptr <- iptr + 1)]] <- list()
      }  # ii
    }  # li.mlm
    if (tmp3.TF[14]) { # ld.mlm
      for (ii in seq(ld.mlm)) {
        misc$earg[[(iptr <- iptr + 1)]] <- list()
      }  # ii
    }  # ld.mlm
  }), list(
    .lmunb.p = lmunb.p, .emunb.p = emunb.p,
    .lmunb.a = lmunb.a, .emunb.a = emunb.a,
    .lmunb.i = lmunb.i, .emunb.i = emunb.i,
    .lmunb.d = lmunb.d, .emunb.d = emunb.d,
    .lsize.p = lsize.p, .esize.p = esize.p,
    .lsize.a = lsize.a, .esize.a = esize.a,
    .lsize.i = lsize.i, .esize.i = esize.i,
    .lsize.d = lsize.d, .esize.d = esize.d,
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
    munb.p <- cbind(eta2theta(eta[, 1], .lmunb.p , .emunb.p ))
    size.p <- cbind(eta2theta(eta[, 2], .lsize.p , .esize.p ))
    ind.munb.z <- 1:2  # Points to munb.p and size.p only.
    munb.a <- munb.i <- munb.d <- munb.p
    size.a <- size.i <- size.d <- size.p


    if (any(tmp3.TF[c(4, 7, 10)])) {  # At least one munb.[aid]
      ind.munb.z <- extra$indeta[c(1:2, 4:5, 7:8, 10:11), 'launch']
      ind.munb.z <- c(na.omit(ind.munb.z))  # At least one value
      munb.a <- if (!tmp3.TF[ 4]) munb.p else
        eta2theta(eta[, extra$indeta[ 4, 1]], .lmunb.a , .emunb.a )
      munb.i <- if (!tmp3.TF[ 7]) munb.p else
        eta2theta(eta[, extra$indeta[ 7, 1]], .lmunb.i , .emunb.i )
      munb.d <- if (!tmp3.TF[10]) munb.p else
        eta2theta(eta[, extra$indeta[10, 1]], .lmunb.d , .emunb.d )

      size.a <- if (!tmp3.TF[ 5]) size.p else
        eta2theta(eta[, extra$indeta[ 5, 1]], .lsize.a , .esize.a )
      size.i <- if (!tmp3.TF[ 8]) size.p else
        eta2theta(eta[, extra$indeta[ 8, 1]], .lsize.i , .esize.i )
      size.d <- if (!tmp3.TF[11]) size.p else
        eta2theta(eta[, extra$indeta[11, 1]], .lsize.d , .esize.d )
    }  # la.mix + li.mix + ld.mix > 0

    if (lall.len) {  # An MLM was fitted
      allprobs <-
        multilogitlink(eta[, -ind.munb.z, drop = FALSE],
                       refLevel = "(Last)",  # Make sure
                       inverse = TRUE)  # rowSums == 1
      if (anyNA(allprobs))
        warning("there are NAs here in slot linkinv")
      if (min(allprobs) == 0 || max(allprobs) == 1)
        warning("fitted probabilities numerically 0 or 1 occurred")

      Nextone <- 0  # Might not be used actually; 0, not 1
      if (tmp3.TF[ 3])
        pobs.mix <- allprobs[, (Nextone <- Nextone + 1)]
      if (tmp3.TF[ 6])
        pstr.mix <- allprobs[, (Nextone <- Nextone + 1)]
      if (tmp3.TF[ 9])
        pdip.mix <- allprobs[, (Nextone <- Nextone + 1)]

      if (tmp3.TF[12]) {
        ind8 <- (Nextone + 1):(Nextone + la.mlm)
        pobs.mlm <- allprobs[, ind8, drop = FALSE]
        dimnames(pobs.mlm) <- list(rownames(eta),
                                   as.character(a.mlm))
        Nextone <- Nextone + la.mlm
      }
      if (tmp3.TF[13]) {
        ind9 <- (Nextone + 1):(Nextone + li.mlm)
        pstr.mlm <- allprobs[, ind9, drop = FALSE]
        dimnames(pstr.mlm) <- list(rownames(eta),
                                   as.character(i.mlm))
        Nextone <- Nextone + li.mlm
      }
      if (tmp3.TF[14]) {
        ind10 <- (Nextone + 1):(Nextone + ld.mlm)
        pdip.mlm <- allprobs[, ind10, drop = FALSE]
        dimnames(pdip.mlm) <- list(rownames(eta),
                                   as.character(d.mlm))
        Nextone <- Nextone + ld.mlm  # Not needed
      }
    }  # lall.len


    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) *
        dgaitdnbinom(y, munb.p = munb.p, size.p = size.p,
                     log = TRUE,  # byrow.aid = F,
                     a.mix = a.mix, i.mix = i.mix, d.mix = d.mix,
                     a.mlm = a.mlm, i.mlm = i.mlm, d.mlm = d.mlm,
                     truncate = truncate,
                     max.support = as.vector( .max.support ),
                     munb.a = munb.a, munb.i = munb.i,
                     munb.d = munb.d,
                     size.a = size.a, size.i = size.i,
                     size.d = size.d,
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
    .lmunb.p = lmunb.p, .emunb.p = emunb.p,
    .lmunb.a = lmunb.a, .emunb.a = emunb.a,
    .lmunb.i = lmunb.i, .emunb.i = emunb.i,
    .lmunb.d = lmunb.d, .emunb.d = emunb.d,
    .lsize.p = lsize.p, .esize.p = esize.p,
    .lsize.a = lsize.a, .esize.a = esize.a,
    .lsize.i = lsize.i, .esize.i = esize.i,
    .lsize.d = lsize.d, .esize.d = esize.d,
    .lpstr.mix = lpstr.mix, .lpobs.mix = lpobs.mix,
    .lpdip.mix = lpdip.mix,
    .epstr.mix = epstr.mix, .epobs.mix = epobs.mix,
    .epdip.mix = epdip.mix,
    .a.mix = a.mix, .i.mix = i.mix, .d.mix = d.mix,
    .a.mlm = a.mlm, .i.mlm = i.mlm, .d.mlm = d.mlm,
    .truncate = truncate, .max.support = max.support ))),
  vfamily = c("gaitdnbinomial"),
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
    munb.a <- munb.i <- munb.d <- 1  # Needed

    if (!is.matrix(eta)) eta <- as.matrix(eta)
    munb.p <- cbind(eta2theta(eta[, 1], .lmunb.p , .emunb.p ))
    size.p <- cbind(eta2theta(eta[, 2], .lsize.p , .esize.p ))
    ind.munb.z <- 1:2  # Points to munb.p and size.p only.


    if (any(tmp3.TF[c(4, 7, 10)])) {  # At least one munb.[aid]
      ind.munb.z <- extra$indeta[c(1:2, 4:5, 7:8, 10:11), 'launch']
      ind.munb.z <- c(na.omit(ind.munb.z))  # At least one value
      munb.a <- if (!tmp3.TF[ 4]) munb.p else
        eta2theta(eta[, extra$indeta[ 4, 1]], .lmunb.a , .emunb.a )
      munb.i <- if (!tmp3.TF[ 7]) munb.p else
        eta2theta(eta[, extra$indeta[ 7, 1]], .lmunb.i , .emunb.i )
      munb.d <- if (!tmp3.TF[10]) munb.p else
        eta2theta(eta[, extra$indeta[10, 1]], .lmunb.d , .emunb.d )

      size.a <- if (!tmp3.TF[ 5]) size.p else
        eta2theta(eta[, extra$indeta[ 5, 1]], .lsize.a , .esize.a )
      size.i <- if (!tmp3.TF[ 8]) size.p else
        eta2theta(eta[, extra$indeta[ 8, 1]], .lsize.i , .esize.i )
      size.d <- if (!tmp3.TF[11]) size.p else
        eta2theta(eta[, extra$indeta[11, 1]], .lsize.d , .esize.d )
    }  # la.mix + li.mix + ld.mix > 0


    if (lall.len) {  # A MLM was fitted
      allprobs <-
        multilogitlink(eta[, -ind.munb.z, drop = FALSE],
                       inverse = TRUE)  # rowSums == 1

      Nextone <- 0  # Might not be used actually; 0, not 1
      if (tmp3.TF[ 3])
        pobs.mix <- allprobs[, (Nextone <- Nextone + 1)]
      if (tmp3.TF[ 6])
        pstr.mix <- allprobs[, (Nextone <- Nextone + 1)]
      if (tmp3.TF[ 9])
        pdip.mix <- allprobs[, (Nextone <- Nextone + 1)]

      if (tmp3.TF[12]) {
        ind8 <- (Nextone + 1):(Nextone + la.mlm)
        pobs.mlm <- allprobs[, ind8, drop = FALSE]
        dimnames(pobs.mlm) <- list(rownames(eta),
                                   as.character(a.mlm))
        Nextone <- Nextone + la.mlm
      }
      if (tmp3.TF[13]) {
        ind9 <- (Nextone + 1):(Nextone + li.mlm)
        pstr.mlm <- allprobs[, ind9, drop = FALSE]
        dimnames(pstr.mlm) <- list(rownames(eta),
                                   as.character(i.mlm))
        Nextone <- Nextone + li.mlm
      }
      if (tmp3.TF[14]) {
        ind10 <- (Nextone + 1):(Nextone + ld.mlm)
        pdip.mlm <- allprobs[, ind10, drop = FALSE]
        dimnames(pdip.mlm) <- list(rownames(eta),
                                   as.character(d.mlm))
        Nextone <- Nextone + ld.mlm  # Not needed
      }
    }  # lall.len

    okay.mlm <-
      all(is.finite(pobs.mlm)) && all(0 < pobs.mlm) &&
      all(is.finite(pstr.mlm)) && all(0 < pstr.mlm) &&
      all(is.finite(pdip.mlm)) && all(0 < pdip.mlm)
    okay.mix <-
      all(is.finite(munb.p)) && all(0 < munb.p) &&
      all(munb.p < .max.support ) &&
      all(is.finite(munb.a)) && all(0 < munb.a) &&
      all(is.finite(munb.i)) && all(0 < munb.i) &&
      all(is.finite(munb.d)) && all(0 < munb.d) &&
      all(is.finite(pobs.mix)) && all(0 < pobs.mix) &&
      all(is.finite(pstr.mix)) && all(0 < pstr.mix) &&
      all(is.finite(pdip.mix)) && all(0 < pdip.mix) &&
      all(pobs.mix + pstr.mix + pdip.mix +
          rowSums(pobs.mlm) + rowSums(pstr.mlm) +
          rowSums(pdip.mlm) < 1)  # Combined
    okay.mlm && okay.mix
  }, list(
    .lmunb.p = lmunb.p, .emunb.p = emunb.p,
    .lmunb.a = lmunb.a, .emunb.a = emunb.a,
    .lmunb.i = lmunb.i, .emunb.i = emunb.i,
    .lmunb.d = lmunb.d, .emunb.d = emunb.d,
    .lsize.p = lsize.p, .esize.p = esize.p,
    .lsize.a = lsize.a, .esize.a = esize.a,
    .lsize.i = lsize.i, .esize.i = esize.i,
    .lsize.d = lsize.d, .esize.d = esize.d,
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

    munb.p <- cbind(eta2theta(eta[, 1], .lmunb.p , .emunb.p ))
    size.p <- cbind(eta2theta(eta[, 2], .lsize.p , .esize.p ))
    ind.munb.z <- 1:2  # Points to munb.p and size.p only.
    munb.a <- munb.i <- munb.d <- munb.p  # Needed;
    size.a <- size.i <- size.d <- size.p
    tmp3.TF <- ( .tmp3.TF )



    if (any(tmp3.TF[c(4, 7, 10)])) {  # At least one munb.[aid]
      ind.munb.z <- extra$indeta[c(1:2, 4:5, 7:8, 10:11), 'launch']
      ind.munb.z <- c(na.omit(ind.munb.z))  # At least one value
      munb.a <- if (!tmp3.TF[ 4]) munb.p else
        eta2theta(eta[, extra$indeta[ 4, 1]], .lmunb.a , .emunb.a )
      munb.i <- if (!tmp3.TF[ 7]) munb.p else
        eta2theta(eta[, extra$indeta[ 7, 1]], .lmunb.i , .emunb.i )
      munb.d <- if (!tmp3.TF[10]) munb.p else
        eta2theta(eta[, extra$indeta[10, 1]], .lmunb.d , .emunb.d )

      size.a <- if (!tmp3.TF[ 5]) size.p else
        eta2theta(eta[, extra$indeta[ 5, 1]], .lsize.a , .esize.a )
      size.i <- if (!tmp3.TF[ 8]) size.p else
        eta2theta(eta[, extra$indeta[ 8, 1]], .lsize.i , .esize.i )
      size.d <- if (!tmp3.TF[11]) size.p else
        eta2theta(eta[, extra$indeta[11, 1]], .lsize.d , .esize.d )
    }  # la.mix + li.mix + ld.mix > 0



    if (lall.len) {  # A AMLM was fitted
      allprobs <-
        multilogitlink(eta[, -ind.munb.z, drop = FALSE],
                       inverse = TRUE)  # rowSums == 1
      Nextone <- 0  # Might not be used actually; 0, not 1
      if (tmp3.TF[ 3])
        pobs.mix <- allprobs[, (Nextone <- Nextone + 1)]
      if (tmp3.TF[ 6])
        pstr.mix <- allprobs[, (Nextone <- Nextone + 1)]
      if (tmp3.TF[ 9])
        pdip.mix <- allprobs[, (Nextone <- Nextone + 1)]

      if (tmp3.TF[12]) {
        ind8 <- (Nextone + 1):(Nextone + la.mlm)
        pobs.mlm <- allprobs[, ind8, drop = FALSE]
        dimnames(pobs.mlm) <- list(rownames(eta),
                                   as.character(a.mlm))
        Nextone <- Nextone + la.mlm
      }
      if (tmp3.TF[13]) {
        ind9 <- (Nextone + 1):(Nextone + li.mlm)
        pstr.mlm <- allprobs[, ind9, drop = FALSE]
        dimnames(pstr.mlm) <- list(rownames(eta),
                                   as.character(i.mlm))
        Nextone <- Nextone + li.mlm
      }
      if (tmp3.TF[14]) {
        ind10 <- (Nextone + 1):(Nextone + ld.mlm)
        pdip.mlm <- allprobs[, ind10, drop = FALSE]
        dimnames(pdip.mlm) <- list(rownames(eta),
                                   as.character(d.mlm))
        Nextone <- Nextone + ld.mlm  # Not needed
      }
    }  # lall.len

    rgaitdnbinom(nsim * length(munb.p),
                 munb.p = munb.p, size.p = size.p,
                 pobs.mlm = pobs.mlm, pstr.mlm = pstr.mlm,
                 pobs.mix = pobs.mix, pstr.mix = pstr.mix,
                 pdip.mix = pdip.mix, pdip.mlm = pdip.mlm,
                 munb.a = munb.a, munb.i = munb.i, munb.d = munb.d,
                 size.a = size.a, size.i = size.i, size.d = size.d,
                 a.mix = a.mix, i.mix = i.mix, d.mix = d.mix,
                 a.mlm = a.mlm, i.mlm = i.mlm, d.mlm = d.mlm,
                 truncate = .truncate , max.support = .max.support )
  }, list(
    .lmunb.p = lmunb.p, .emunb.p = emunb.p,
    .lmunb.a = lmunb.a, .emunb.a = emunb.a,
    .lmunb.i = lmunb.i, .emunb.i = emunb.i,
    .lmunb.d = lmunb.d, .emunb.d = emunb.d,
    .lsize.p = lsize.p, .esize.p = esize.p,
    .lsize.a = lsize.a, .esize.a = esize.a,
    .lsize.i = lsize.i, .esize.i = esize.i,
    .lsize.d = lsize.d, .esize.d = esize.d,
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
    calA.p  <- tmp3.TF[ 3]
    calI.p  <- tmp3.TF[ 6]
    calD.p  <- tmp3.TF[ 9]
    calA.np <- tmp3.TF[12]
    calI.np <- tmp3.TF[13]
    calD.np <- tmp3.TF[14]

    Denom1.a1 <- Denom1.i1 <- Denom1.d1 <-
                 Denom2.i1 <- Denom2.d1 <- 0  # Denom2.a1 is unneeded

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
    munb.p <- cbind(eta2theta(eta[, 1], .lmunb.p , .emunb.p ))
    size.p <- cbind(eta2theta(eta[, 2], .lsize.p , .esize.p ))
    ind.munb.z <- 1:2  # Points to munb.p and size.p only.
    munb.a <- munb.i <- munb.d <- munb.p  # Needed;
    size.a <- size.i <- size.d <- size.p



    if (any(tmp3.TF[c(4, 7, 10)])) {  # At least one munb.[aid]
      ind.munb.z <- extra$indeta[c(1:2, 4:5, 7:8, 10:11), 'launch']
      ind.munb.z <- c(na.omit(ind.munb.z))  # At least one value
      munb.a <- if (!tmp3.TF[ 4]) munb.p else
        eta2theta(eta[, extra$indeta[ 4, 1]], .lmunb.a , .emunb.a )
      munb.i <- if (!tmp3.TF[ 7]) munb.p else
        eta2theta(eta[, extra$indeta[ 7, 1]], .lmunb.i , .emunb.i )
      munb.d <- if (!tmp3.TF[10]) munb.p else
        eta2theta(eta[, extra$indeta[10, 1]], .lmunb.d , .emunb.d )

      size.a <- if (!tmp3.TF[ 5]) size.p else
        eta2theta(eta[, extra$indeta[ 5, 1]], .lsize.a , .esize.a )
      size.i <- if (!tmp3.TF[ 8]) size.p else
        eta2theta(eta[, extra$indeta[ 8, 1]], .lsize.i , .esize.i )
      size.d <- if (!tmp3.TF[11]) size.p else
        eta2theta(eta[, extra$indeta[11, 1]], .lsize.d , .esize.d )
    }  # la.mix + li.mix + ld.mix > 0


    if (lall.len) {  # A MLM was fitted.
      allprobs <-
        multilogitlink(eta[, -ind.munb.z, drop = FALSE],
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
      if (extra$control.trace)
        cat("Minimum baseline (reserve) probability = ",
            format(minprob.baseline, digits = 3), "\n")


      Nextone <- 0  # Might not be used actually; 0, not 1
      if (tmp3.TF[ 3])
        pobs.mix <- allprobs[, (Nextone <- Nextone + 1)]
      if (tmp3.TF[ 6])
        pstr.mix <- allprobs[, (Nextone <- Nextone + 1)]
      if (tmp3.TF[ 9])
        pdip.mix <- allprobs[, (Nextone <- Nextone + 1)]

      if (tmp3.TF[12]) {
        ind8 <- (Nextone + 1):(Nextone + la.mlm)
        pobs.mlm <- allprobs[, ind8, drop = FALSE]
        dimnames(pobs.mlm) <- list(rownames(eta),
                                   as.character(a.mlm))
        Nextone <- Nextone + la.mlm
      }
      if (tmp3.TF[13]) {
        ind9 <- (Nextone + 1):(Nextone + li.mlm)
        pstr.mlm <- allprobs[, ind9, drop = FALSE]
        dimnames(pstr.mlm) <- list(rownames(eta),
                                   as.character(i.mlm))
        Nextone <- Nextone + li.mlm
      }
      if (tmp3.TF[14]) {
        ind10 <- (Nextone + 1):(Nextone + ld.mlm)
        pdip.mlm <- allprobs[, ind10, drop = FALSE]
        dimnames(pdip.mlm) <- list(rownames(eta),
                                   as.character(d.mlm))
        Nextone <- Nextone + ld.mlm  # Not needed
      }
    }  # lall.len


    ltruncat <- length(truncate)
    M1 <- max(extra$indeta, na.rm = TRUE)
    NOS <- ncol(eta) / M1  # extra$NOS
    if (NOS != 1) stop("can only handle 1 response")



    is.a.mixed <- if (tmp3.TF[ 3])
      rowSums(extra$skip.mix.a) > 0 else rep(FALSE, n)
    is.i.mixed <- if (tmp3.TF[ 6])
      rowSums(extra$skip.mix.i) > 0 else rep(FALSE, n)
    is.d.mixed <- if (tmp3.TF[ 9])
      rowSums(extra$skip.mix.d) > 0 else rep(FALSE, n)
    is.a.mlmed <- if (tmp3.TF[12])
      rowSums(extra$skip.mlm.a) > 0 else rep(FALSE, n)
    is.i.mlmed <- if (tmp3.TF[13])
      rowSums(extra$skip.mlm.i) > 0 else rep(FALSE, n)
    is.d.mlmed <- if (tmp3.TF[14])
      rowSums(extra$skip.mlm.d) > 0 else rep(FALSE, n)

    is.ns <- !is.a.mlmed & !is.i.mlmed  & !is.d.mlmed  &
             !is.a.mixed & !is.i.mixed  & !is.d.mixed  # & !is.truncd

    dl.dmunb.p <- y / munb.p - (1 + y / size.p) / (1 + munb.p / size.p)
    dl.dmunb.p[!is.ns] <- 0  # For is.a.mixed & is.a.mlmed
    dl.dsize.p <- digamma(y + size.p) - digamma(size.p) +
      log1p(-munb.p / (size.p + munb.p)) -
      (y - munb.p) / (size.p + munb.p)
    dl.dsize.p[!is.ns] <- 0


    prob.mlm.a <- if (la.mlm) rowSums(pobs.mlm) else 0  # scalar okay
    prob.mlm.i <- if (li.mlm) rowSums(pstr.mlm) else 0  # scalar okay
    prob.mlm.d <- if (ld.mlm) rowSums(pdip.mlm) else 0  # scalar okay



    nb.munb.der1 <- function(y, munb, size)
      dnbinom(y, size = size, mu = munb) *
      (y / munb - 1) / (1 + munb / size)
    nb.munb.der2 <- function(y, munb, size)
      dnbinom(y, size = size, mu = munb) * (
      (y + size) / (munb + size)^2 - y / munb^2 +
      ((y / munb - 1) / (1 + munb / size))^2)

    nb.size.der1 <- function(y, munb, size)
      dnbinom(y, size = size, mu = munb) *
      (digamma(y + size) - digamma(size) +
      log1p(-munb / (munb + size)) - (y - munb) / (size + munb))
    nb.size.der2 <- function(y, munb, size)
      dnbinom(y, size = size, mu = munb) * (
      trigamma(y + size) - trigamma(size) +
      munb / (size * (size + munb)) +
      (y - munb) / (size + munb)^2 +
      (digamma(y + size) - digamma(size) +
       log1p(-munb / (munb + size)) - (y - munb) / (size + munb))^2)

    nb.musz.der2 <- function(y, munb, size)
      dnbinom(y, size = size, mu = munb) *
      (y - munb) / (size + munb)^2 +
      nb.munb.der1(y, munb, size) *
      nb.size.der1(y, munb, size) / dnbinom(y, size = size, mu = munb)



    sumD.mix.1a.p1 <- sumD.mix.2a.p1 <-
    sumD.mix.1a.p2 <- sumD.mix.2a.p2 <-
                      sumD.mix.2a.p4 <- matrix(0, n, NOS)
    if (la.mix > 0) {  # \calA_p
      DA.mix.0mat.a  <-  # Matches naming convention further below
      DA.mix.1mat.a1 <- DA.mix.1mat.a2 <- matrix(0, n, la.mix)
      for (jay in seq(la.mix)) {
        aval <- a.mix[jay]
        sumD.mix.1a.p1 <- sumD.mix.1a.p1 +
                          nb.munb.der1(aval, munb.p, size.p)
        sumD.mix.2a.p1 <- sumD.mix.2a.p1 +
                          nb.munb.der2(aval, munb.p, size.p)
        sumD.mix.1a.p2 <- sumD.mix.1a.p2 +
                          nb.size.der1(aval, munb.p, size.p)
        sumD.mix.2a.p2 <- sumD.mix.2a.p2 +
                          nb.size.der2(aval, munb.p, size.p)
        sumD.mix.2a.p4 <- sumD.mix.2a.p4 +
                          nb.musz.der2(aval, munb.p, size.p)
        pmf.a <- dnbinom(aval, size = size.a, mu = munb.a)
        DA.mix.0mat.a [, jay] <- pmf.a
        DA.mix.1mat.a1[, jay] <- nb.munb.der1(aval, munb.a, size.a)
        DA.mix.1mat.a2[, jay] <- nb.size.der1(aval, munb.a, size.a)
      }
      Denom1.a1 <- rowSums(DA.mix.1mat.a1)  # aka sumD.mix.1a.a
      Denom1.a2 <- rowSums(DA.mix.1mat.a2)
    }  # la.mix > 0





    if (li.mix) {
      DI.mix.0mat.i  <-  # wrt inflated distribution
      DI.mix.1mat.i1 <- DI.mix.2mat.i1 <-
      DI.mix.1mat.i2 <- DI.mix.2mat.i2 <-
                        DI.mix.2mat.i4 <- matrix(0, n, li.mix)
      DP.mix.0mat.i  <-  # wrt parent distribution
      DP.mix.1mat.i1 <- DP.mix.2mat.i1 <-
      DP.mix.1mat.i2 <- DP.mix.2mat.i2 <-
                        DP.mix.2mat.i4 <- matrix(0, n, li.mix)
      for (jay in seq(li.mix)) {
        ival <- i.mix[jay]
        pmf.i <- dnbinom(ival, size = size.i, mu = munb.i)
        DI.mix.0mat.i [, jay] <- pmf.i
        DI.mix.1mat.i1[, jay] <- nb.munb.der1(ival, munb.i, size.i)
        DI.mix.2mat.i1[, jay] <- nb.munb.der2(ival, munb.i, size.i)
        DI.mix.1mat.i2[, jay] <- nb.size.der1(ival, munb.i, size.i)
        DI.mix.2mat.i2[, jay] <- nb.size.der2(ival, munb.i, size.i)
        DI.mix.2mat.i4[, jay] <- nb.musz.der2(ival, munb.i, size.i)
        pmf.p <- dnbinom(ival, size = size.p, mu = munb.p)
        DP.mix.0mat.i [, jay] <- pmf.p
        DP.mix.1mat.i1[, jay] <- nb.munb.der1(ival, munb.p, size.p)
        DP.mix.2mat.i1[, jay] <- nb.munb.der2(ival, munb.p, size.p)
        DP.mix.1mat.i2[, jay] <- nb.size.der1(ival, munb.p, size.p)
        DP.mix.2mat.i2[, jay] <- nb.size.der2(ival, munb.p, size.p)
        DP.mix.2mat.i4[, jay] <- nb.musz.der2(ival, munb.p, size.p)
      }  # jay
      Denom1.i1 <- rowSums(DI.mix.1mat.i1)
      Denom2.i1 <- rowSums(DI.mix.2mat.i1)
      Denom1.i2 <- rowSums(DI.mix.1mat.i2)
      Denom2.i2 <- rowSums(DI.mix.2mat.i2)
      Denom2.i4 <- rowSums(DI.mix.2mat.i4)
    }  # li.mix


    if (ld.mix) {
      DD.mix.0mat.d  <-  # wrt deflated distribution
      DD.mix.1mat.d1 <- DD.mix.2mat.d1 <-
      DD.mix.1mat.d2 <- DD.mix.2mat.d2 <-
                        DD.mix.2mat.d4 <- matrix(0, n, ld.mix)
      DP.mix.0mat.d  <-  # wrt parent distribution
      DP.mix.1mat.d1 <- DP.mix.2mat.d1 <-
      DP.mix.1mat.d2 <- DP.mix.2mat.d2 <-
                        DP.mix.2mat.d4 <- matrix(0, n, ld.mix)
      for (jay in seq(ld.mix)) {
        dval <- d.mix[jay]
        pmf.d <- dnbinom(dval, size = size.d, mu = munb.d)
        DD.mix.0mat.d [, jay] <- pmf.d
        DD.mix.1mat.d1[, jay] <- nb.munb.der1(dval, munb.d, size.d)
        DD.mix.2mat.d1[, jay] <- nb.munb.der2(dval, munb.d, size.d)
        DD.mix.1mat.d2[, jay] <- nb.size.der1(dval, munb.d, size.d)
        DD.mix.2mat.d2[, jay] <- nb.size.der2(dval, munb.d, size.d)
        DD.mix.2mat.d4[, jay] <- nb.musz.der2(dval, munb.d, size.d)
        pmf.p <- dnbinom(dval, size = size.p, mu = munb.p)
        DP.mix.0mat.d [, jay] <- pmf.p
        DP.mix.1mat.d1[, jay] <- nb.munb.der1(dval, munb.p, size.p)
        DP.mix.2mat.d1[, jay] <- nb.munb.der2(dval, munb.p, size.p)
        DP.mix.1mat.d2[, jay] <- nb.size.der1(dval, munb.p, size.p)
        DP.mix.2mat.d2[, jay] <- nb.size.der2(dval, munb.p, size.p)
        DP.mix.2mat.d4[, jay] <- nb.musz.der2(dval, munb.p, size.p)
      }  # jay
      Denom1.d1 <- rowSums(DD.mix.1mat.d1)
      Denom2.d1 <- rowSums(DD.mix.2mat.d1)
      Denom1.d2 <- rowSums(DD.mix.1mat.d2)
      Denom2.d2 <- rowSums(DD.mix.2mat.d2)
      Denom2.d4 <- rowSums(DD.mix.2mat.d4)
    }  # ld.mix

    Bits <- moments.gaitdcombo.nbinom(
              munb.p = munb.p, size.p = size.p,
              pobs.mix = pobs.mix, pobs.mlm = pobs.mlm,
              pstr.mix = pstr.mix, pstr.mlm = pstr.mlm,
              pdip.mix = pdip.mix, pdip.mlm = pdip.mlm,
              a.mix = a.mix, i.mix = i.mix, d.mix = d.mix,
              a.mlm = a.mlm, i.mlm = i.mlm, d.mlm = d.mlm,
              munb.a = munb.a, munb.i = munb.i, munb.d = munb.d,
              size.a = size.a, size.i = size.i, size.d = size.d,
              truncate = truncate, max.support = max.support)


    sumD.mlm.1a.p1 <- sumD.mlm.2a.p1 <-
    sumD.mlm.1a.p2 <- sumD.mlm.2a.p2 <-
                      sumD.mlm.2a.p4 <- matrix(0, n, NOS)
    if (la.mlm)
      for (aval in a.mlm) {
        sumD.mlm.1a.p1 <- sumD.mlm.1a.p1 +
                          nb.munb.der1(aval, munb.p, size.p)
        sumD.mlm.2a.p1 <- sumD.mlm.2a.p1 +
                          nb.munb.der2(aval, munb.p, size.p)
        sumD.mlm.1a.p2 <- sumD.mlm.1a.p2 +
                          nb.size.der1(aval, munb.p, size.p)
        sumD.mlm.2a.p2 <- sumD.mlm.2a.p2 +
                          nb.size.der2(aval, munb.p, size.p)
        sumD.mlm.2a.p4 <- sumD.mlm.2a.p4 +
                          nb.musz.der2(aval, munb.p, size.p)
      }


    Denom0.p <- c(Bits[["cdf.max.s"]]   - Bits[["SumT0.p"]] -
                  Bits[["SumA0.mix.p"]] - Bits[["SumA0.mlm.p"]])
    Numer <- 1 - pobs.mix - pstr.mix - prob.mlm.a - prob.mlm.i +
                 pdip.mix +            prob.mlm.d
    Denom0.a <- c(Bits[["SumA0.mix.a"]])  # Not .p
    Denom0.i <- c(Bits[["SumI0.mix.i"]])
    Denom0.d <- c(Bits[["SumD0.mix.d"]])



    Dp.mlm.0mat.i  <-  # wrt parent distribution
    Dp.mlm.1mat.i1 <- Dp.mlm.2mat.i1 <-
    Dp.mlm.1mat.i2 <- Dp.mlm.2mat.i2 <-
                      Dp.mlm.2mat.i4 <- matrix(0, n, NOS)
    if (li.mlm > 0) {
      Dp.mlm.0mat.i  <-  # wrt parent distribution
      Dp.mlm.1mat.i1 <- Dp.mlm.2mat.i1 <-
      Dp.mlm.1mat.i2 <- Dp.mlm.2mat.i2 <-
                        Dp.mlm.2mat.i4 <- matrix(0, n, li.mlm)
      for (jay in seq(li.mlm)) {
        ival <- i.mlm[jay]
        pmf.p <- dnbinom(ival, size = size.p, mu = munb.p)
        Dp.mlm.0mat.i [, jay] <- pmf.p
        Dp.mlm.1mat.i1[, jay] <- nb.munb.der1(ival, munb.p, size.p)
        Dp.mlm.2mat.i1[, jay] <- nb.munb.der2(ival, munb.p, size.p)
        Dp.mlm.1mat.i2[, jay] <- nb.size.der1(ival, munb.p, size.p)
        Dp.mlm.2mat.i2[, jay] <- nb.size.der2(ival, munb.p, size.p)
        Dp.mlm.2mat.i4[, jay] <- nb.musz.der2(ival, munb.p, size.p)
      }  # jay
    }  # li.mlm



    Dp.mlm.0mat.d  <-  # wrt parent distribution
    Dp.mlm.1mat.d1 <- Dp.mlm.2mat.d1 <-
    Dp.mlm.1mat.d2 <- Dp.mlm.2mat.d2 <-
                      Dp.mlm.2mat.d4 <- matrix(0, n, NOS)
    if (ld.mlm > 0) {
      Dp.mlm.0mat.d  <-  # wrt parent distribution
      Dp.mlm.1mat.d1 <- Dp.mlm.2mat.d1 <-
      Dp.mlm.1mat.d2 <- Dp.mlm.2mat.d2 <-
                        Dp.mlm.2mat.d4 <- matrix(0, n, ld.mlm)
      for (jay in seq(ld.mlm)) {
        dval <- d.mlm[jay]
        pmf.p <- dnbinom(dval, size = size.p, mu = munb.p)
        Dp.mlm.0mat.d [, jay] <- pmf.p
        Dp.mlm.1mat.d1[, jay] <- nb.munb.der1(dval, munb.p, size.p)
        Dp.mlm.2mat.d1[, jay] <- nb.munb.der2(dval, munb.p, size.p)
        Dp.mlm.1mat.d2[, jay] <- nb.size.der1(dval, munb.p, size.p)
        Dp.mlm.2mat.d2[, jay] <- nb.size.der2(dval, munb.p, size.p)
        Dp.mlm.2mat.d4[, jay] <- nb.musz.der2(dval, munb.p, size.p)
      }  # jay
    }  # ld.mlm




    

    sumD.1t.p1 <- sumD.2t.p1 <-
    sumD.1t.a1 <- sumD.2t.a1 <-
    sumD.1t.i1 <- sumD.2t.i1 <-
    sumD.1t.d1 <- sumD.2t.d1 <- matrix(0, n, NOS)
    sumD.1t.p2 <- sumD.2t.p2 <- sumD.2t.p4 <-
    sumD.1t.a2 <- sumD.2t.a2 <-
    sumD.1t.i2 <- sumD.2t.i2 <-
    sumD.1t.d2 <- sumD.2t.d2 <- matrix(0, n, NOS)
    if (ltruncat)
      for (tval in truncate) {
        sumD.1t.p1 <- sumD.1t.p1 + nb.munb.der1(tval, munb.p, size.p)
        sumD.2t.p1 <- sumD.2t.p1 + nb.munb.der2(tval, munb.p, size.p)
        sumD.1t.a1 <- sumD.1t.a1 + nb.munb.der1(tval, munb.a, size.a)
        sumD.2t.a1 <- sumD.2t.a1 + nb.munb.der2(tval, munb.a, size.a)
        sumD.1t.i1 <- sumD.1t.i1 + nb.munb.der1(tval, munb.i, size.i)
        sumD.2t.i1 <- sumD.2t.i1 + nb.munb.der2(tval, munb.i, size.i)
        sumD.1t.d1 <- sumD.1t.d1 + nb.munb.der1(tval, munb.d, size.d)
        sumD.2t.d1 <- sumD.2t.d1 + nb.munb.der2(tval, munb.d, size.d)

        sumD.1t.p2 <- sumD.1t.p2 + nb.size.der1(tval, munb.p, size.p)
        sumD.2t.p2 <- sumD.2t.p2 + nb.size.der2(tval, munb.p, size.p)
        sumD.1t.a2 <- sumD.1t.a2 + nb.size.der1(tval, munb.a, size.a)
        sumD.2t.a2 <- sumD.2t.a2 + nb.size.der2(tval, munb.a, size.a)
        sumD.1t.i2 <- sumD.1t.i2 + nb.size.der1(tval, munb.i, size.i)
        sumD.2t.i2 <- sumD.2t.i2 + nb.size.der2(tval, munb.i, size.i)
        sumD.1t.d2 <- sumD.1t.d2 + nb.size.der1(tval, munb.d, size.d)
        sumD.2t.d2 <- sumD.2t.d2 + nb.size.der2(tval, munb.d, size.d)

        sumD.2t.p4 <- sumD.2t.p4 + nb.musz.der2(tval, munb.p, size.p)
      }
    if (is.finite(max.support)) {
    stop("upper tail derivs are unavailable for finite 'max.support'")
    }  # is.finite(max.support)


    Denom1.p1  <- c(-sumD.1t.p1 - sumD.mlm.1a.p1 - sumD.mix.1a.p1)
    Denom2.p1  <- c(-sumD.2t.p1 - sumD.mlm.2a.p1 - sumD.mix.2a.p1)
    Denom1.p2  <- c(-sumD.1t.p2 - sumD.mlm.1a.p2 - sumD.mix.1a.p2)
    Denom2.p2  <- c(-sumD.2t.p2 - sumD.mlm.2a.p2 - sumD.mix.2a.p2)
    Denom2.p4  <- c(-sumD.2t.p4 - sumD.mlm.2a.p4 - sumD.mix.2a.p4)




    d0B.PI.mlm  <- Dp.mlm.0mat.i  / Denom0.p
    d1B.PI.mlm1 <- Dp.mlm.1mat.i1 / Denom0.p -  # Most general
                   Dp.mlm.0mat.i  * Denom1.p1 / Denom0.p^2
    d1B.PI.mlm2 <- Dp.mlm.1mat.i2 / Denom0.p -  # Most general
                   Dp.mlm.0mat.i  * Denom1.p2 / Denom0.p^2
    d2B.PI.mlm1 <- Dp.mlm.2mat.i1 / Denom0.p -
               2 * Dp.mlm.1mat.i1 * Denom1.p1 / Denom0.p^2 -
                   Dp.mlm.0mat.i  * Denom2.p1 / Denom0.p^2 +
               2 * Dp.mlm.0mat.i  * (Denom1.p1^2) / Denom0.p^3
    d2B.PI.mlm2 <- Dp.mlm.2mat.i2 / Denom0.p -
               2 * Dp.mlm.1mat.i2 * Denom1.p2 / Denom0.p^2 -
                   Dp.mlm.0mat.i  * Denom2.p2 / Denom0.p^2 +
               2 * Dp.mlm.0mat.i  * (Denom1.p2^2) / Denom0.p^3
    d2B.PI.mlm4 <- Dp.mlm.2mat.i4 / Denom0.p  -
                   Dp.mlm.1mat.i1 * Denom1.p2 / Denom0.p^2 -
                   Dp.mlm.1mat.i2 * Denom1.p1 / Denom0.p^2 -
                   Dp.mlm.0mat.i  * Denom2.p4 / Denom0.p^2 +
               2 * Dp.mlm.0mat.i  * Denom1.p1 * Denom1.p2 / Denom0.p^3




    d0B.PD.mlm  <- Dp.mlm.0mat.d  / Denom0.p
    d1B.PD.mlm1 <- Dp.mlm.1mat.d1 / Denom0.p -  # This is most general
                   Dp.mlm.0mat.d  * Denom1.p1 / Denom0.p^2
    d1B.PD.mlm2 <- Dp.mlm.1mat.d2 / Denom0.p -  # This is most general
                   Dp.mlm.0mat.d  * Denom1.p2 / Denom0.p^2
    d2B.PD.mlm1 <- Dp.mlm.2mat.d1 / Denom0.p -
               2 * Dp.mlm.1mat.d1 * Denom1.p1 / Denom0.p^2 -
                   Dp.mlm.0mat.d  * Denom2.p1 / Denom0.p^2 +
               2 * Dp.mlm.0mat.d  * (Denom1.p1^2) / Denom0.p^3
    d2B.PD.mlm2 <- Dp.mlm.2mat.d2 / Denom0.p -
               2 * Dp.mlm.1mat.d2 * Denom1.p2 / Denom0.p^2 -
                   Dp.mlm.0mat.d  * Denom2.p2 / Denom0.p^2 +
               2 * Dp.mlm.0mat.d  * (Denom1.p2^2) / Denom0.p^3
    d2B.PD.mlm4 <- Dp.mlm.2mat.d4 / Denom0.p -
                   Dp.mlm.1mat.d1 * Denom1.p2 / Denom0.p^2 -
                   Dp.mlm.1mat.d2 * Denom1.p1 / Denom0.p^2 -
                   Dp.mlm.0mat.d  * Denom2.p4 / Denom0.p^2 +
               2 * Dp.mlm.0mat.d  * Denom1.p1 *
                                    Denom1.p2 / Denom0.p^3


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

      d1A.i1 <- (DI.mix.1mat.i1 - DI.mix.0mat.i *
                 Denom1.i1 / Denom0.i) / Denom0.i
      d2A.i1 <- (DI.mix.2mat.i1 - (2 * DI.mix.1mat.i1 * Denom1.i1 +
                 DI.mix.0mat.i  * Denom2.i1) / Denom0.i +
             2 * DI.mix.0mat.i  * (Denom1.i1 / Denom0.i)^2) / Denom0.i

      d1A.i2 <- (DI.mix.1mat.i2 - DI.mix.0mat.i *
                 Denom1.i2 / Denom0.i) / Denom0.i
      d2A.i2 <- (DI.mix.2mat.i2 - (2 * DI.mix.1mat.i2 * Denom1.i2 +
                 DI.mix.0mat.i  * Denom2.i2) / Denom0.i +
             2 * DI.mix.0mat.i  * (Denom1.i2 / Denom0.i)^2) / Denom0.i
      d2A.i4 <-  DI.mix.2mat.i4 /  Denom0.i  -
                 DI.mix.1mat.i1 *  Denom1.i2 / Denom0.i^2 -
                 DI.mix.1mat.i2 *  Denom1.i1 / Denom0.i^2 -
                 DI.mix.0mat.i  *  Denom2.i4 / Denom0.i^2 +
             2 * DI.mix.0mat.i  *  Denom1.i1 *
                                   Denom1.i2 / Denom0.i^3


      d1B.PI.mix1 <- DP.mix.1mat.i1 / Denom0.p -
                     DP.mix.0mat.i  * Denom1.p1 / Denom0.p^2
      d1B.PI.mix2 <- DP.mix.1mat.i2 / Denom0.p -
                     DP.mix.0mat.i  * Denom1.p2 / Denom0.p^2
      d2B.PI.mix1 <-     DP.mix.2mat.i1 /  Denom0.p -
                     2 * DP.mix.1mat.i1 *  Denom1.p1    / Denom0.p^2 -
                         DP.mix.0mat.i  *  Denom2.p1    / Denom0.p^2 +
                     2 * DP.mix.0mat.i  * (Denom1.p1^2) / Denom0.p^3
      d2B.PI.mix2 <-     DP.mix.2mat.i2 /  Denom0.p -
                     2 * DP.mix.1mat.i2 *  Denom1.p2    / Denom0.p^2 -
                         DP.mix.0mat.i  *  Denom2.p2    / Denom0.p^2 +
                     2 * DP.mix.0mat.i  * (Denom1.p2^2) / Denom0.p^3
      d2B.PI.mix4 <-     DP.mix.2mat.i4 /  Denom0.p -
                         DP.mix.1mat.i1 *  Denom1.p2    / Denom0.p^2 -
                         DP.mix.1mat.i2 *  Denom1.p1    / Denom0.p^2 -
                         DP.mix.0mat.i  *  Denom2.p4    / Denom0.p^2 +
                     2 * DP.mix.0mat.i  *  Denom1.p1 *
                                           Denom1.p2    / Denom0.p^3
    }  # li.mix > 0






    if (ld.mix > 0) {
      d0A.d <- DD.mix.0mat.d / Denom0.d
      d0B.PD.mix <- DP.mix.0mat.d / Denom0.p
      DELTA.d.mix <- Numer * d0B.PD.mix - pdip.mix * d0A.d

      d1A.d1 <- (DD.mix.1mat.d1 - DD.mix.0mat.d *
                 Denom1.d1 / Denom0.d) / Denom0.d
      d2A.d1 <- (DD.mix.2mat.d1 - (2 * DD.mix.1mat.d1 * Denom1.d1 +
                 DD.mix.0mat.d * Denom2.d1) / Denom0.d +
             2 * DD.mix.0mat.d * (Denom1.d1 / Denom0.d)^2) / Denom0.d

      d1A.d2 <- (DD.mix.1mat.d2 - DD.mix.0mat.d *
                 Denom1.d2 / Denom0.d) / Denom0.d
      d2A.d2 <- (DD.mix.2mat.d2 - (2 * DD.mix.1mat.d2 * Denom1.d2 +
                 DD.mix.0mat.d * Denom2.d2) / Denom0.d +
             2 * DD.mix.0mat.d * (Denom1.d2 / Denom0.d)^2) / Denom0.d
      d2A.d4 <-  DD.mix.2mat.d4 /  Denom0.d  -
                 DD.mix.1mat.d1 *  Denom1.d2 / Denom0.d^2 -
                 DD.mix.1mat.d2 *  Denom1.d1 / Denom0.d^2 -
                 DD.mix.0mat.d  *  Denom2.d4 / Denom0.d^2 +
             2 * DD.mix.0mat.d  *  Denom1.d1 *
                                   Denom1.d2 / Denom0.d^3

      d1B.PD.mix1 <- DP.mix.1mat.d1 / Denom0.p -
                     DP.mix.0mat.d  * Denom1.p1 / Denom0.p^2
      d2B.PD.mix1 <- DP.mix.2mat.d1 /  Denom0.p -
                 2 * DP.mix.1mat.d1 *  Denom1.p1    / Denom0.p^2 -
                     DP.mix.0mat.d  *  Denom2.p1    / Denom0.p^2 +
                 2 * DP.mix.0mat.d  * (Denom1.p1^2) / Denom0.p^3

      d1B.PD.mix2 <- DP.mix.1mat.d2 / Denom0.p -
                     DP.mix.0mat.d  * Denom1.p2 / Denom0.p^2
      d2B.PD.mix2 <- DP.mix.2mat.d2 /  Denom0.p -
                 2 * DP.mix.1mat.d2 *  Denom1.p2    / Denom0.p^2 -
                     DP.mix.0mat.d  *  Denom2.p2    / Denom0.p^2 +
                 2 * DP.mix.0mat.d  * (Denom1.p2^2) / Denom0.p^3
      d2B.PD.mix4 <- DP.mix.2mat.d4 /  Denom0.p -
                     DP.mix.1mat.d1 *  Denom1.p2    / Denom0.p^2 -
                     DP.mix.1mat.d2 *  Denom1.p1    / Denom0.p^2 -
                     DP.mix.0mat.d  *  Denom2.p4    / Denom0.p^2 +
                 2 * DP.mix.0mat.d  *  Denom1.p1 *
                                       Denom1.p2    / Denom0.p^3
    }  # ld.mix > 0




    if (la.mix) {
      d0A.a  <- DA.mix.0mat.a  / Denom0.a
      d1A.a1 <- DA.mix.1mat.a1 / Denom0.a -
                DA.mix.0mat.a  * Denom1.a1 / Denom0.a^2
      d1A.a2 <- DA.mix.1mat.a2 / Denom0.a -
                DA.mix.0mat.a  * Denom1.a2 / Denom0.a^2
    }  # la.mix




    dl.dmunb.a <- dl.dmunb.i <- dl.dmunb.d <- numeric(n)
    dl.dsize.a <- dl.dsize.i <- dl.dsize.d <- numeric(n)
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



    if (tmp3.TF[12] && la.mlm) {  # aka \calA_{np}
      dl.dpobs.mlm <- matrix(-1 / Numer, n, la.mlm)  # \notin calS
      dl.dpobs.mlm[!is.ns, ] <- 0  # For a.mix only really
      for (jay in seq(la.mlm)) {
        aval <- a.mlm[jay]
        is.alt.j.mlm <- extra$skip.mlm.a[, jay]  # Logical vector
        tmp7a <- 1 / pobs.mlm[is.alt.j.mlm, jay]
        dl.dpobs.mlm[is.alt.j.mlm, jay] <- tmp7a
      }  # jay
    }  # la.mlm



    dl.dmunb.p[is.ns] <- dl.dmunb.p[is.ns] -
                         (Denom1.p1 / Denom0.p)[is.ns]
    dl.dsize.p[is.ns] <- dl.dsize.p[is.ns] -
                         (Denom1.p2 / Denom0.p)[is.ns]


    
    if (tmp3.TF[13] && li.mlm > 0) {  # aka \calI_{np}
      dl.dpstr.mlm <- matrix(-1 / Numer, n, li.mlm)
      dl.dpstr.mlm[!is.ns, ] <- 0  # For a.mlm and a.mix

      for (jay in seq(li.mlm)) {
        is.inf.j.mlm <- extra$skip.mlm.i[, jay]  # Logical vector
        tmp7..m <- Numer * d1B.PI.mlm1[, jay] / DELTA.i.mlm[, jay]
        tmp7..s <- Numer * d1B.PI.mlm2[, jay] / DELTA.i.mlm[, jay]
        dl.dmunb.p[is.inf.j.mlm] <- tmp7..m[is.inf.j.mlm]
        dl.dsize.p[is.inf.j.mlm] <- tmp7..s[is.inf.j.mlm]


        tmp9i <- d0B.PI.mlm[, jay] / DELTA.i.mlm[, jay]
        n.tmp <- -tmp9i[is.inf.j.mlm]
        p.tmp <- +tmp9i[is.inf.j.mlm]
        if (tmp3.TF[12] && la.mlm) dl.dpobs.mlm[is.inf.j.mlm, ] <- n.tmp
        if (tmp3.TF[ 3] && la.mix) dl.dpobs.mix[is.inf.j.mlm  ] <- n.tmp
        if (tmp3.TF[ 6] && li.mix) dl.dpstr.mix[is.inf.j.mlm  ] <- n.tmp
        if (tmp3.TF[14] && ld.mlm) dl.dpdip.mlm[is.inf.j.mlm, ] <- p.tmp
        if (tmp3.TF[ 9] && ld.mix) dl.dpdip.mix[is.inf.j.mlm  ] <- p.tmp


        tmp8 <- (1 - d0B.PI.mlm[, jay]) / DELTA.i.mlm[, jay]
        dl.dpstr.mlm[is.inf.j.mlm, ] <- n.tmp  # tmp9[is.inf.j.mlm]
        dl.dpstr.mlm[is.inf.j.mlm, jay] <- tmp8[is.inf.j.mlm]
      }  # jay
    }  # li.mlm > 0






    if (tmp3.TF[14] && ld.mlm > 0) {  # aka \calD_{np}

      for (jay in seq(ld.mlm)) {
        is.def.j.mlm <- extra$skip.mlm.d[, jay]  # Logical vector
        tmp7..m <- Numer * d1B.PD.mlm1[, jay] / DELTA.d.mlm[, jay]
        tmp7..s <- Numer * d1B.PD.mlm2[, jay] / DELTA.d.mlm[, jay]
        dl.dmunb.p[is.def.j.mlm] <- tmp7..m[is.def.j.mlm]  # 20211020
        dl.dsize.p[is.def.j.mlm] <- tmp7..s[is.def.j.mlm]
 
        tmp9d <- d0B.PD.mlm[, jay] / DELTA.d.mlm[, jay]
        p.tmp <- +tmp9d[is.def.j.mlm]
        n.tmp <- -tmp9d[is.def.j.mlm]
        if (tmp3.TF[13] && li.mlm) dl.dpstr.mlm[is.def.j.mlm, ] <- n.tmp
        if (tmp3.TF[ 6] && li.mix) dl.dpstr.mix[is.def.j.mlm  ] <- n.tmp
        if (tmp3.TF[12] && la.mlm) dl.dpobs.mlm[is.def.j.mlm, ] <- n.tmp
        if (tmp3.TF[ 3] && la.mix) dl.dpobs.mix[is.def.j.mlm  ] <- n.tmp
        if (tmp3.TF[ 9] && ld.mix) dl.dpdip.mix[is.def.j.mlm  ] <- p.tmp
                                   dl.dpdip.mlm[is.def.j.mlm, ] <- p.tmp
        dl.dpdip.mlm[is.def.j.mlm, jay] <-
        dl.dpdip.mlm[is.def.j.mlm, jay] -
        1 / DELTA.d.mlm[is.def.j.mlm, jay]

      }  # jay
    }  # ld.mlm > 0




    



    if (tmp3.TF[ 3] && la.mix) {  # aka \calA_{p}
      dl.dpobs.mix[is.a.mixed] <- 1 / pobs.mix[is.a.mixed]

      if (tmp3.TF[ 4] && la.mix > 1)
        for (jay in seq(la.mix)) {
          is.alt.j.mix <- extra$skip.mix.a[, jay]  # Logical vector
          tmp2..m <- d1A.a1[, jay] / d0A.a[, jay]
          tmp2..s <- d1A.a2[, jay] / d0A.a[, jay]
          dl.dmunb.a[is.alt.j.mix] <- tmp2..m[is.alt.j.mix]  # ccc.
          dl.dsize.a[is.alt.j.mix] <- tmp2..s[is.alt.j.mix]
        }  # jay
    }  # la.mix




    if (tmp3.TF[ 6] && li.mix > 0) {  # aka \calI_{p}
      for (jay in seq(li.mix)) {
        ival <- i.mix[jay]
        is.inf.j.mix <- extra$skip.mix.i[, jay]  # Logical vector
        tmp7..m <- Numer * d1B.PI.mix1[, jay] / DELTA.i.mix[, jay]
        tmp7..s <- Numer * d1B.PI.mix2[, jay] / DELTA.i.mix[, jay]
        dl.dmunb.p[is.inf.j.mix] <- tmp7..m[is.inf.j.mix]
        dl.dsize.p[is.inf.j.mix] <- tmp7..s[is.inf.j.mix]
        tmp8 <- (d0A.i[, jay] - d0B.PI.mix[, jay]) / DELTA.i.mix[, jay]
        dl.dpstr.mix[is.inf.j.mix] <- tmp8[is.inf.j.mix]
        if (li.mix > 1) {
          tmp2..m <- pstr.mix * d1A.i1[, jay] / DELTA.i.mix[, jay]
          tmp2..s <- pstr.mix * d1A.i2[, jay] / DELTA.i.mix[, jay]
          dl.dmunb.i[is.inf.j.mix] <- tmp2..m[is.inf.j.mix]
          dl.dsize.i[is.inf.j.mix] <- tmp2..s[is.inf.j.mix]
        }





        tmp9i <- d0B.PI.mix[, jay] / DELTA.i.mix[, jay]
        n.tmp <- -tmp9i[is.inf.j.mix]
        p.tmp <- +tmp9i[is.inf.j.mix]
        if (tmp3.TF[ 3] && la.mix) dl.dpobs.mix[is.inf.j.mix  ] <- n.tmp
        if (tmp3.TF[12] && la.mlm) dl.dpobs.mlm[is.inf.j.mix, ] <- n.tmp
        if (tmp3.TF[13] && li.mlm) dl.dpstr.mlm[is.inf.j.mix, ] <- n.tmp
        if (tmp3.TF[14] && ld.mlm) dl.dpdip.mlm[is.inf.j.mix, ] <- p.tmp
        if (tmp3.TF[ 9] && ld.mix) dl.dpdip.mix[is.inf.j.mix  ] <- p.tmp

      }  # jay
    }  # li.mix > 0




    if (tmp3.TF[ 9] && ld.mix > 0) {  # aka \calD_{p}
      for (jay in seq(ld.mix)) {
        dval <- d.mix[jay]
        is.def.j.mix <- extra$skip.mix.d[, jay]  # Logical vector
        tmp7..m <- Numer * d1B.PD.mix1[, jay] / DELTA.d.mix[, jay]
        tmp7..s <- Numer * d1B.PD.mix2[, jay] / DELTA.d.mix[, jay]
        dl.dmunb.p[is.def.j.mix] <- tmp7..m[is.def.j.mix]
        dl.dsize.p[is.def.j.mix] <- tmp7..s[is.def.j.mix]
        tmp8 <- (d0B.PD.mix[, jay] - d0A.d[, jay]) / DELTA.d.mix[, jay]
        dl.dpdip.mix[is.def.j.mix] <- tmp8[is.def.j.mix]

        if (ld.mix > 1) {
  if (any(is.na(d1A.d1)))
  stop("NAs found in d1A.d1")
          tmp2..m <- (-pdip.mix) * d1A.d1[, jay] / DELTA.d.mix[, jay]
          tmp2..s <- (-pdip.mix) * d1A.d2[, jay] / DELTA.d.mix[, jay]
          dl.dmunb.d[is.def.j.mix] <- tmp2..m[is.def.j.mix]
          dl.dsize.d[is.def.j.mix] <- tmp2..s[is.def.j.mix]
        }


        tmp9d <- d0B.PD.mix[, jay] / DELTA.d.mix[, jay]
        n.tmp <- -tmp9d[is.def.j.mix]
        p.tmp <- +tmp9d[is.def.j.mix]
        if (tmp3.TF[13] && li.mlm) dl.dpstr.mlm[is.def.j.mix, ] <- n.tmp
        if (tmp3.TF[ 6] && li.mix) dl.dpstr.mix[is.def.j.mix  ] <- n.tmp
        if (tmp3.TF[12] && la.mlm) dl.dpobs.mlm[is.def.j.mix, ] <- n.tmp
        if (tmp3.TF[ 3] && la.mix) dl.dpobs.mix[is.def.j.mix  ] <- n.tmp
        if (tmp3.TF[14] && ld.mlm) dl.dpdip.mlm[is.def.j.mix, ] <- p.tmp

      }  # jay
   }  # ld.mix > 0







    new.ansd <- matrix(0, n, M)  # Same dimension as eta
    tmp3.TF <- !is.na(rowSums(extra$indeta))



    if (lall.len) {  # An MLM fitted
      all6.dldp <- cbind(if (tmp3.TF[ 3]) dl.dpobs.mix else NULL,
                         if (tmp3.TF[ 6]) dl.dpstr.mix else NULL,
                         if (tmp3.TF[ 9]) dl.dpdip.mix else NULL,
                         if (tmp3.TF[12]) dl.dpobs.mlm else NULL,
                         if (tmp3.TF[13]) dl.dpstr.mlm else NULL,
                         if (tmp3.TF[14]) dl.dpdip.mlm else NULL)




      rSs.tmp <- rowSums(allprobs[, -ncol(allprobs), drop = FALSE] *
                         all6.dldp)
      new.ansd[, -ind.munb.z] <- allprobs[, -ncol(allprobs)] *
                                   (all6.dldp - rSs.tmp)
    }  # lall.len


    

    dmunb.p.deta <- dtheta.deta(munb.p, .lmunb.p , .emunb.p )
    dsize.p.deta <- dtheta.deta(size.p, .lsize.p , .esize.p )
    if (tmp3.TF[ 4]) {
      dmunb.a.deta <- dtheta.deta(munb.a, .lmunb.a , .emunb.a )
      dsize.a.deta <- dtheta.deta(size.a, .lsize.a , .esize.a )
    }
    if (tmp3.TF[ 7]) {
      dmunb.i.deta <- dtheta.deta(munb.i, .lmunb.i , .emunb.i )
      dsize.i.deta <- dtheta.deta(size.i, .lsize.i , .esize.i )
    }
    if (tmp3.TF[10]) {
      dmunb.d.deta <- dtheta.deta(munb.d, .lmunb.d , .emunb.d )
      dsize.d.deta <- dtheta.deta(size.d, .lsize.d , .esize.d )
    }
    new.ansd[, 1] <- dl.dmunb.p * dmunb.p.deta
    new.ansd[, 2] <- dl.dsize.p * dsize.p.deta
    if (tmp3.TF[ 4]) {
      new.ansd[, extra$indeta[ 4, 1]] <- dl.dmunb.a * dmunb.a.deta
      new.ansd[, extra$indeta[ 5, 1]] <- dl.dsize.a * dsize.a.deta
    }
    if (tmp3.TF[ 7]) {
      new.ansd[, extra$indeta[ 7, 1]] <- dl.dmunb.i * dmunb.i.deta
      new.ansd[, extra$indeta[ 8, 1]] <- dl.dsize.i * dsize.i.deta
    }
    if (tmp3.TF[10]) {
      new.ansd[, extra$indeta[10, 1]] <- dl.dmunb.d * dmunb.d.deta
      new.ansd[, extra$indeta[11, 1]] <- dl.dsize.d * dsize.d.deta
    }
    onecoln.indeta <- extra$indeta[1:11, ]  # One coln params only
    onecoln.indeta <- na.omit(onecoln.indeta)  # Only those present
    allcnames <- c(rownames(onecoln.indeta),
                   as.character(c(a.mlm, i.mlm, d.mlm)))
    colnames(new.ansd) <- allcnames

 


 if (any(is.na(new.ansd)))
 stop("look here87")
    c(w) * new.ansd
  }), list(
    .lmunb.p = lmunb.p, .emunb.p = emunb.p,
    .lmunb.a = lmunb.a, .emunb.a = emunb.a,
    .lmunb.i = lmunb.i, .emunb.i = emunb.i,
    .lmunb.d = lmunb.d, .emunb.d = emunb.d,
    .lsize.p = lsize.p, .esize.p = esize.p,
    .lsize.a = lsize.a, .esize.a = esize.a,
    .lsize.i = lsize.i, .esize.i = esize.i,
    .lsize.d = lsize.d, .esize.d = esize.d,
    .lpstr.mix = lpstr.mix, .lpobs.mix = lpobs.mix,
    .lpdip.mix = lpdip.mix,
    .epstr.mix = epstr.mix, .epobs.mix = epobs.mix,
    .epdip.mix = epdip.mix,
    .a.mix = a.mix, .i.mix = i.mix, .d.mix = d.mix,
    .a.mlm = a.mlm, .i.mlm = i.mlm, .d.mlm = d.mlm,
    .tmp3.TF = tmp3.TF,  # .tmp3 = tmp3,
    .truncate = truncate, .max.support = max.support ))),

  weight = eval(substitute(expression({  # gaitdnbinomial
      





    wz <- matrix(0, n, M * (M + 1) / 2)  # The complete size
    cond.EY.p <-
      c(munb.p - Bits[["SumT1.p"]] -
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
    ned2l.dpobs.mix.munb.p <- zero0n  # mB overwritten below [4279]
    ned2l.dpobs.mix.munb.a <- zero0n  # Fini; (3, 4) element
    ned2l.dpobs.mix.munb.i <- zero0n  # mB overwritten below
    ned2l.dpobs.mix.munb.d <- zero0n  # mB overwritten below
    ned2l.dpstr.mix.munb.p <- zero0n  # Optional (1, 6) element
    ned2l.dpstr.mix.munb.a <- zero0n  # Final; nothing to do
    ned2l.dpstr.mix.munb.i <- zero0n  # mB overwritten below
    ned2l.dpstr.mix.munb.d <- zero0n  # mB overwritten below
    ned2l.dpdip.mix.munb.p <- zero0n  # Optional (1, 9) element
    ned2l.dpdip.mix.munb.i <- zero0n  # Optional (7, 9) element
    ned2l.dpdip.mix.munb.d <- zero0n  # Optional (9, 10) element


    ned2l.dpobs.mix.size.p <- zero0n  # mB overwritten below [4279]
    ned2l.dpobs.mix.size.a <- zero0n  # Fini; (3, 5) element
    ned2l.dpobs.mix.size.i <- zero0n  # mB overwritten below
    ned2l.dpobs.mix.size.d <- zero0n  # mB overwritten below
    ned2l.dpstr.mix.size.p <- zero0n  # Optional (2, 6) element
    ned2l.dpstr.mix.size.i <- zero0n  # mB overwritten below
    ned2l.dpstr.mix.size.d <- zero0n  # mB overwritten below
    ned2l.dpdip.mix.size.p <- zero0n  # Optional (2, 9) element
    ned2l.dpdip.mix.size.i <- zero0n  # Optional (8, 9) element
    ned2l.dpdip.mix.size.d <- zero0n  # Optional (9, 11) element

    ned2l.dpobs.mlm.size.i <- zero0n  # Optional (8, 12) element
    ned2l.dpobs.mlm.size.d <- zero0n  # Optional (11, 12) element


      posn.pobs.mix <- as.vector(extra$indeta[ 3, 'launch'])
      posn.munb.a <- as.vector(extra$indeta[ 4, 'launch'])
      posn.size.a <- as.vector(extra$indeta[ 5, 'launch'])
      posn.pstr.mix <- as.vector(extra$indeta[ 6, 'launch'])
      posn.munb.i <- as.vector(extra$indeta[ 7, 'launch'])
      posn.size.i <- as.vector(extra$indeta[ 8, 'launch'])
      posn.pdip.mix <- as.vector(extra$indeta[ 9, 'launch'])
      posn.munb.d <- as.vector(extra$indeta[10, 'launch'])
      posn.size.d <- as.vector(extra$indeta[11, 'launch'])
      posn.pobs.mlm <- as.vector(extra$indeta[12, 'launch'])
      posn.pstr.mlm <- as.vector(extra$indeta[13, 'launch'])
      posn.pdip.mlm <- as.vector(extra$indeta[14, 'launch'])





    ned2l.dpdip.mix2         <-  # Elt (9, 9)
    ned2l.dpstr.mix2         <-  # Elt (6, 6). Unchanged by deflation.
    ned2l.dpobs.mlm.pstr.mix <-  # Elts (6, >=12). (((09)))
    ned2l.dpobs.mix.pstr.mix <- +probns / Numer^2  # ccc Elt (3, 6)
    if (all(c(la.mix, li.mlm) > 0))  # (((08)))
    ned2l.dpobs.mix.pstr.mlm <- matrix( probns / Numer^2, n, li.mlm)
    if (all(c(li.mix, li.mlm) > 0))  # (((10)))
    ned2l.dpstr.mix.pstr.mlm <- matrix( probns / Numer^2, n, li.mlm)
    if (all(c(ld.mix, ld.mlm) > 0))  # (((21)))
    ned2l.dpdip.mix.pdip.mlm <- matrix( probns / Numer^2, n, ld.mlm)


    ned2l.dpobs.mlm.pdip.mix <-  # Elts (9, >=12). (((19)))
    ned2l.dpstr.mix.pdip.mix <-  # Elt (6, 9)
    ned2l.dpobs.mix.pdip.mix <- -probns / Numer^2  # ccc Elt (3, 9)
    if (all(c(la.mix, ld.mlm) > 0))  # (((17)))
    ned2l.dpobs.mix.pdip.mlm <- matrix(-probns / Numer^2, n, ld.mlm)
    if (all(c(li.mix, ld.mlm) > 0))  # (((18)))
    ned2l.dpstr.mix.pdip.mlm <- matrix(-probns / Numer^2, n, ld.mlm)
    if (all(c(ld.mix, li.mlm) > 0))  # (((20)))
    ned2l.dpdip.mix.pstr.mlm <- matrix(-probns / Numer^2, n, li.mlm)





    ned2l.dmunb.p2 <- probns * (
       cond.EY.p / munb.p^2 -
      (cond.EY.p + size.p) / (munb.p + size.p)^2 +  # ddd
      Denom2.p1 / Denom0.p - (Denom1.p1 / Denom0.p)^2) + 
      (if (tmp3.TF[ 6] && li.mix) Numer *
      rowSums(Numer *
             (d1B.PI.mix1^2) / DELTA.i.mix - d2B.PI.mix1) else 0) +
      (if (tmp3.TF[13] && li.mlm) Numer *
      rowSums(Numer *
             (d1B.PI.mlm1^2) / DELTA.i.mlm - d2B.PI.mlm1) else 0) +
      (if (tmp3.TF[ 9] && ld.mix) Numer *
      rowSums(Numer *
             (d1B.PD.mix1^2) / DELTA.d.mix - d2B.PD.mix1) else 0) +
      (if (tmp3.TF[14] && ld.mlm) Numer *  # nnn.
      rowSums(Numer *
             (d1B.PD.mlm1^2) / DELTA.d.mlm - d2B.PD.mlm1) else 0)


    wz[, iam(1, 1, M)] <- ned2l.dmunb.p2 * dmunb.p.deta^2



    diff.trig.nbd <- numeric(n)  # Storage
    ind2 <- rep_len(FALSE, n)  # Used for SFS
    max.chunk.MB <- ( .max.chunk.MB )
    eff.p <- sort(c( .cutoff.prob , 1 - .cutoff.prob ))
    Q.mins <- 0
    Q.maxs <- round(qnbinom(p = eff.p[2], mu = munb.p,
                            size = size.p) * 1.1) + 30
    eps.trig <- ( .eps.trig )
    Q.MAXS <- if ( .lsize.p == "loglink")
      pmax(10, ceiling(size.p / sqrt(eps.trig))) else Inf
    Q.maxs <- pmin(Q.maxs, Q.MAXS)
    ind1 <- if (max.chunk.MB > 0)
      (Q.maxs - Q.mins < max.support) else
      stop("argument 'max.chunk.MB' > 0 is needed")
    if ((NN <- sum(ind1)) > 0) {
      Object.Size <- NN * 8 * max(Q.maxs - Q.mins) / (2^20)
      n.chunks <- if (intercept.only) 1 else
                  max(1, ceiling( Object.Size / max.chunk.MB))
      chunk.rows <- ceiling(NN / n.chunks)
      ind2 <- ind1  # Save this
      wind2 <- which(ind1)
      upr.ptr <- 0
      lwr.ptr <- upr.ptr + 1
      while (lwr.ptr <= NN) {
        upr.ptr <- min(upr.ptr + chunk.rows, NN)
        sind2 <- wind2[lwr.ptr:upr.ptr]
        diff.trig.nbd[sind2] <-
          EIM.NB.specialp(mu = munb.p, size = size.p,
                          y.max = max(Q.maxs[sind2]),
                          cutoff.prob = .cutoff.prob ,
                          extra.bit = FALSE,  # extra.bit omitted
                          intercept.only = intercept.only)
        lwr.ptr <- upr.ptr + 1
      }  # while
      cond.diff.trig.nbd <- diff.trig.nbd
      specialvals <- c(truncate, a.mlm, a.mix, i.mlm, i.mix,
                                 d.mlm, d.mix)  # max.support??
      for (sval in specialvals) {
         cond.diff.trig.nbd <- cond.diff.trig.nbd -
           (trigamma(size.p) - trigamma(sval + size.p)) *
           dnbinom(sval, size = size.p, mu = munb.p)
      }  # for sval
      cond.diff.trig.nbd <- cond.diff.trig.nbd / c(
        Denom0.p -
        Bits[["SumD0.mix.p"]] - Bits[["SumD0.mlm.p"]] -  # 20211109
        Bits[["SumI0.mix.p"]] - Bits[["SumI0.mlm.p"]])
    }  # if ((NN <- sum(ind1)) > 0)



    ned2l.dsize.p2 <- probns * (cond.diff.trig.nbd +
      (-munb.p) / (size.p * (size.p + munb.p)) -
      (cond.EY.p - munb.p) / (munb.p + size.p)^2 +
      Denom2.p2 / Denom0.p - (Denom1.p2 / Denom0.p)^2) +
      (if (tmp3.TF[ 6] && li.mix) Numer *
      rowSums(Numer *
             (d1B.PI.mix2^2) / DELTA.i.mix - d2B.PI.mix2) else 0) +
      (if (tmp3.TF[13] && li.mlm) Numer *
      rowSums(Numer *
             (d1B.PI.mlm2^2) / DELTA.i.mlm - d2B.PI.mlm2) else 0) +
      (if (tmp3.TF[ 9] && ld.mix) Numer *
      rowSums(Numer *
             (d1B.PD.mix2^2) / DELTA.d.mix - d2B.PD.mix2) else 0) +
      (if (tmp3.TF[14] && ld.mlm) Numer *  # nnn.
      rowSums(Numer *
             (d1B.PD.mlm2^2) / DELTA.d.mlm - d2B.PD.mlm2) else 0)
    wz[, iam(2, 2, M)] <- ned2l.dsize.p2 * dsize.p.deta^2



    ned2l.dmusz.p2 <- probns * (
      (munb.p - cond.EY.p) / (munb.p + size.p)^2 +
      Denom2.p4 / Denom0.p - Denom1.p1 * Denom1.p2 / Denom0.p^2) +
      (if (tmp3.TF[ 6] && li.mix) Numer *
      rowSums(Numer * d1B.PI.mix1 * d1B.PI.mix2 / DELTA.i.mix -
              d2B.PI.mix4) else 0) +
      (if (tmp3.TF[13] && li.mlm) Numer *
      rowSums(Numer * d1B.PI.mlm1 * d1B.PI.mlm2 / DELTA.i.mlm -
              d2B.PI.mlm4) else 0) +
      (if (tmp3.TF[ 9] && ld.mix) Numer *
      rowSums(Numer * d1B.PD.mix1 * d1B.PD.mix2 / DELTA.d.mix -
              d2B.PD.mix4) else 0) +
      (if (tmp3.TF[14] && ld.mlm) Numer *  # nnn.
      rowSums(Numer * d1B.PD.mlm1 * d1B.PD.mlm2 / DELTA.d.mlm -
              d2B.PD.mlm4) else 0)
    wz[, iam(1, 2, M)] <- ned2l.dmusz.p2 * dmunb.p.deta * dsize.p.deta





    ned2l.dpobs.mix2 <- 1 / pobs.mix + probns / Numer^2
    if (tmp3.TF[ 6] && li.mix > 0) {
      ned2l.dpobs.mix2 <-  # More just below, ccc
      ned2l.dpobs.mix2 + rowSums(d0B.PI.mix^2 / DELTA.i.mix)
    }
    if (tmp3.TF[13] && li.mlm > 0) {
      ned2l.dpobs.mix2 <-  # ccc.
      ned2l.dpobs.mix2 + rowSums(d0B.PI.mlm^2 / DELTA.i.mlm)
    }
    if (tmp3.TF[ 9] && ld.mix > 0) {
      ned2l.dpobs.mix2 <-  # nnn
      ned2l.dpobs.mix2 + rowSums(d0B.PD.mix^2 / DELTA.d.mix)
    }
    if (tmp3.TF[14] && ld.mlm > 0) {
      ned2l.dpobs.mix2 <-  # nnn
      ned2l.dpobs.mix2 + rowSums(d0B.PD.mlm^2 / DELTA.d.mlm)
    }
    if (tmp3.TF[ 3] && la.mix > 0)
      wz[, iam(3, 3, M)] <- ned2l.dpobs.mix2  # Link done later



    if (tmp3.TF[ 4] && la.mix > 1) {
      ned2l.dmunb.a2 <- pobs.mix * (
        rowSums((DA.mix.1mat.a1^2) / DA.mix.0mat.a) / Denom0.a -
        (Denom1.a1 / Denom0.a)^2)  # ccc.
      wz[, iam(4, 4, M)] <- ned2l.dmunb.a2 * dmunb.a.deta^2


      ned2l.dsize.a2 <- pobs.mix * (
        rowSums((DA.mix.1mat.a2^2) / DA.mix.0mat.a) / Denom0.a -
        (Denom1.a2 / Denom0.a)^2)  # ddd.
      wz[, iam(5, 5, M)] <- ned2l.dsize.a2 * dsize.a.deta^2


      ned2l.dmusz.a2 <- pobs.mix * (
        rowSums(DA.mix.1mat.a1 * DA.mix.1mat.a2 / DA.mix.0mat.a) / (
        Denom0.a) -
        Denom1.a1 * Denom1.a2 / Denom0.a^2)  # ddd.
      wz[, iam(4, 5, M)] <- ned2l.dmusz.a2 * dmunb.a.deta * dsize.a.deta
    }  # tmp3.TF[ 4] && la.mix > 1





    if (tmp3.TF[ 6] && li.mix > 0) {
      ned2l.dpstr.mix2 <-
      ned2l.dpstr.mix2 +
        rowSums((d0A.i - d0B.PI.mix)^2 / DELTA.i.mix)

      if (tmp3.TF[ 3] && la.mix > 0)
        ned2l.dpobs.mix.munb.p <-
        ned2l.dpobs.mix.munb.p +
          rowSums(d1B.PI.mix1 * (1 - Numer * d0B.PI.mix / DELTA.i.mix))

      if (tmp3.TF[ 3] && la.mix > 0)
        ned2l.dpobs.mix.size.p <-
        ned2l.dpobs.mix.size.p +
          rowSums(d1B.PI.mix2 * (1 - Numer * d0B.PI.mix / DELTA.i.mix))

      ned2l.dpstr.mix.munb.p <-
      ned2l.dpstr.mix.munb.p + rowSums(
        d1B.PI.mix1 * (1 + Numer * (d0A.i - d0B.PI.mix) / DELTA.i.mix))

      ned2l.dpstr.mix.size.p <-
      ned2l.dpstr.mix.size.p + rowSums(
        d1B.PI.mix2 * (1 + Numer * (d0A.i - d0B.PI.mix) / DELTA.i.mix))

      if (tmp3.TF[ 9])
        ned2l.dpdip.mix.munb.p <-
        ned2l.dpdip.mix.munb.p - rowSums(
          d1B.PI.mix1 * (1 - Numer * d0B.PI.mix / DELTA.i.mix))

      if (tmp3.TF[ 9])
        ned2l.dpdip.mix.size.p <-
        ned2l.dpdip.mix.size.p - rowSums(
          d1B.PI.mix2 * (1 - Numer * d0B.PI.mix / DELTA.i.mix))

      if (all(tmp3.TF[c(3, 6)]))
        ned2l.dpobs.mix.pstr.mix <-  # ccc
        ned2l.dpobs.mix.pstr.mix +
          rowSums(-d0B.PI.mix * (d0A.i - d0B.PI.mix) / DELTA.i.mix)

      if (all(tmp3.TF[c(6, 9)]))
        ned2l.dpstr.mix.pdip.mix <-
        ned2l.dpstr.mix.pdip.mix + rowSums(
          d0B.PI.mix * (d0A.i - d0B.PI.mix) / DELTA.i.mix)

      if (!is.na(posn.pdip.mix)) {
        ned2l.dpdip.mix2 <-
        ned2l.dpdip.mix2 + rowSums(d0B.PI.mix^2 / DELTA.i.mix)
      }
    }  # (tmp3.TF[ 6] && li.mix > 0)





    if (all(tmp3.TF[c(3, 6, 13)])) {  # was la.mix > 0 & DELTA.i.mix
      ned2l.dpobs.mix.pstr.mix <-  # ccc
      ned2l.dpobs.mix.pstr.mix + rowSums(d0B.PI.mlm^2 / DELTA.i.mlm)
    }
    if (all(tmp3.TF[c(3, 6,  9)])) {  # ==  ld.mix > 0 & DELTA.d.mix
      ned2l.dpobs.mix.pstr.mix <-  # nnn
      ned2l.dpobs.mix.pstr.mix + rowSums(d0B.PD.mix^2 / DELTA.d.mix)
    }
    if (all(tmp3.TF[c(3, 6, 14)])) {  # ==  ld.mlm > 0 & DELTA.d.mlm
      ned2l.dpobs.mix.pstr.mix <-  # nnn.
      ned2l.dpobs.mix.pstr.mix + rowSums(d0B.PD.mlm^2 / DELTA.d.mlm)
    }
    if (!is.na(posn.pobs.mix) && !is.na(posn.pstr.mix))
      wz[, iam(posn.pobs.mix, posn.pstr.mix, M)] <-
        ned2l.dpobs.mix.pstr.mix  # Link done later




    if (all(tmp3.TF[c(3, 9)]))
      ned2l.dpobs.mix.pdip.mix <-  # nnn
      ned2l.dpobs.mix.pdip.mix +
        rowSums( d0B.PD.mix * (d0A.d - d0B.PD.mix) / DELTA.d.mix)

    if (all(tmp3.TF[c(3, 9, 13)])) {  # ==  li.mlm > 0 & DELTA.i.mix
      ned2l.dpobs.mix.pdip.mix <-  # nnn
      ned2l.dpobs.mix.pdip.mix - rowSums(d0B.PI.mlm^2 / DELTA.i.mlm)
    }
    if (all(tmp3.TF[c(3, 9,  6)])) {  # ==  li.mix > 0 & DELTA.i.mix
      ned2l.dpobs.mix.pdip.mix <-  # nnn
      ned2l.dpobs.mix.pdip.mix - rowSums(d0B.PI.mix^2 / DELTA.i.mix)
    }
    if (all(tmp3.TF[c(3, 9, 14)])) {  # ==  ld.mlm > 0 & DELTA.d.mlm
      ned2l.dpobs.mix.pdip.mix <-  # nnn.
      ned2l.dpobs.mix.pdip.mix - rowSums(d0B.PD.mlm^2 / DELTA.d.mlm)
    }
    if (!is.na(posn.pobs.mix) && !is.na(posn.pdip.mix))
      wz[, iam(posn.pobs.mix, posn.pdip.mix, M)] <-
        ned2l.dpobs.mix.pdip.mix  # Link done later







    if (tmp3.TF[ 7] && li.mix > 1) {  # \calI_{p}, includes \theta_i.

      ned2l.dmunb.p.munb.i <- pstr.mix * Numer *
        rowSums(d1A.i1 * d1B.PI.mix1 / DELTA.i.mix)  # ccc.
      wz[, iam(1, posn.munb.i, M)] <- ned2l.dmunb.p.munb.i *
        dmunb.p.deta * dmunb.i.deta  # All links done here

      ned2l.dmunb.p.size.i <- pstr.mix * Numer *
        rowSums(d1A.i2 * d1B.PI.mix1 / DELTA.i.mix)
      wz[, iam(1, posn.size.i, M)] <- ned2l.dmunb.p.size.i *
        dmunb.p.deta * dsize.i.deta  # All links done here

      ned2l.dsize.p.munb.i <- pstr.mix * Numer *
        rowSums(d1A.i1 * d1B.PI.mix2 / DELTA.i.mix)
      wz[, iam(2, posn.munb.i, M)] <- ned2l.dsize.p.munb.i *
        dsize.p.deta * dmunb.i.deta  # All links done here

      ned2l.dsize.p.size.i <- pstr.mix * Numer *
        rowSums(d1A.i2 * d1B.PI.mix2 / DELTA.i.mix)
      wz[, iam(2, posn.size.i, M)] <- ned2l.dsize.p.size.i *
        dsize.p.deta * dsize.i.deta  # All links done here

      ned2l.dmunb.i2 <- pstr.mix *
        rowSums(pstr.mix * (d1A.i1^2) / DELTA.i.mix - d2A.i1)  # ccc.
      wz[, iam(posn.munb.i, posn.munb.i, M)] <-
        ned2l.dmunb.i2 * dmunb.i.deta^2

      ned2l.dsize.i2 <- pstr.mix *
        rowSums(pstr.mix * (d1A.i2^2) / DELTA.i.mix - d2A.i2)  # ddd.
      wz[, iam(posn.size.i, posn.size.i, M)] <-
        ned2l.dsize.i2 * dsize.i.deta^2

      ned2l.dmusz.i2 <- pstr.mix *
        rowSums(pstr.mix * d1A.i1 * d1A.i2 / DELTA.i.mix - d2A.i4)
      wz[, iam(posn.munb.i, posn.size.i, M)] <-
        ned2l.dmusz.i2 * dmunb.i.deta * dsize.i.deta

      if (tmp3.TF[ 3]) {  # tmp3.TF[ 6] is TRUE, given tmp3.TF[ 7]
        ned2l.dpobs.mix.munb.i <-
          rowSums(-pstr.mix * d1A.i1 * d0B.PI.mix / DELTA.i.mix)  # ccc.
        wz[, iam(posn.pobs.mix, posn.munb.i, M)] <-
          ned2l.dpobs.mix.munb.i  # * dmunb.i.deta done later
      }

      if (tmp3.TF[ 3]) {  # tmp3.TF[ 6] is TRUE, given tmp3.TF[ 7]
        ned2l.dpobs.mix.size.i <-
          rowSums(-pstr.mix * d1A.i2 * d0B.PI.mix / DELTA.i.mix)  # ddd.
        wz[, iam(posn.pobs.mix, posn.size.i, M)] <-
          ned2l.dpobs.mix.size.i  # * dsize.i.deta done later
      }

      if (tmp3.TF[ 6]) {
        ned2l.dpstr.mix.munb.i <- rowSums(  # ccc.
          d1A.i1 * (pstr.mix * (d0A.i - d0B.PI.mix) / DELTA.i.mix - 1))
        wz[, iam(posn.pstr.mix, posn.munb.i, M)] <-
          ned2l.dpstr.mix.munb.i  # * dmunb.i.deta done later
      }

      if (tmp3.TF[ 6]) {
        ned2l.dpstr.mix.size.i <- rowSums(  # ddd.
          d1A.i2 * (pstr.mix * (d0A.i - d0B.PI.mix) / DELTA.i.mix - 1))
        wz[, iam(posn.pstr.mix, posn.size.i, M)] <-
          ned2l.dpstr.mix.size.i  # * dsize.i.deta done later
      }

      if (all(tmp3.TF[c(7, 9)])) {
        ned2l.dpdip.mix.munb.i <- rowSums(
          (-pstr.mix) * d0B.PI.mix * d1A.i1 / DELTA.i.mix)
        wz[, iam(posn.pdip.mix, posn.munb.i, M)] <-
          ned2l.dpdip.mix.munb.i  # link done later
      }

      if (all(tmp3.TF[c(7, 9)])) {
        ned2l.dpdip.mix.size.i <- rowSums(
          (-pstr.mix) * d0B.PI.mix * d1A.i2 / DELTA.i.mix)
        wz[, iam(posn.pdip.mix, posn.size.i, M)] <-
          ned2l.dpdip.mix.size.i  # link done later
      }

      if (tmp3.TF[12]) {
        ned2l.dpobs.mlm.munb.i <- rowSums(
          -pstr.mix * d0B.PI.mix * d1A.i1 / DELTA.i.mix)  # ccc.
        for (uuu in seq(la.mlm))
          wz[, iam(posn.pobs.mlm - 1 + uuu, posn.munb.i, M)] <-
            ned2l.dpobs.mlm.munb.i  # * dmunb.i.deta done later
      }

      if (tmp3.TF[12]) {
        ned2l.dpobs.mlm.size.i <- rowSums(
          -pstr.mix * d0B.PI.mix * d1A.i2 / DELTA.i.mix)  # ddd.
        for (uuu in seq(la.mlm))
          wz[, iam(posn.pobs.mlm - 1 + uuu, posn.size.i, M)] <-
            ned2l.dpobs.mlm.size.i  # * dsize.i.deta done later
      }
    }  # (tmp3.TF[ 7] && li.mix > 1)






    if (tmp3.TF[ 9] && ld.mix > 0) {  # \calD_{p}, maybe w. \theta_d

      if (tmp3.TF[ 3] && la.mix > 0)
        ned2l.dpobs.mix.munb.p <-
        ned2l.dpobs.mix.munb.p +
          rowSums(d1B.PD.mix1 * (1 - Numer * d0B.PD.mix / DELTA.d.mix))

      if (tmp3.TF[ 3] && la.mix > 0)
        ned2l.dpobs.mix.size.p <-
        ned2l.dpobs.mix.size.p +
          rowSums(d1B.PD.mix2 * (1 - Numer * d0B.PD.mix / DELTA.d.mix))

      ned2l.dpstr.mix.munb.p <-
      ned2l.dpstr.mix.munb.p + rowSums(
        d1B.PD.mix1 * (1 - Numer * d0B.PD.mix / DELTA.d.mix))

      ned2l.dpstr.mix.size.p <-
      ned2l.dpstr.mix.size.p + rowSums(
        d1B.PD.mix2 * (1 - Numer * d0B.PD.mix / DELTA.d.mix))

      ned2l.dpdip.mix.munb.p <-
      ned2l.dpdip.mix.munb.p - rowSums(
        d1B.PD.mix1 * (1 + Numer * (d0A.d - d0B.PD.mix) / DELTA.d.mix))

      ned2l.dpdip.mix.size.p <-
      ned2l.dpdip.mix.size.p - rowSums(
        d1B.PD.mix2 * (1 + Numer * (d0A.d - d0B.PD.mix) / DELTA.d.mix))

      if (!is.na(posn.pstr.mix)) {
        ned2l.dpstr.mix2 <-
        ned2l.dpstr.mix2 + rowSums(d0B.PD.mix^2 / DELTA.d.mix)
      }

      if (all(tmp3.TF[c(6, 9)]))
        ned2l.dpstr.mix.pdip.mix <-
        ned2l.dpstr.mix.pdip.mix + rowSums(
          d0B.PD.mix * (d0A.d - d0B.PD.mix) / DELTA.d.mix)

      ned2l.dpdip.mix2 <-
      ned2l.dpdip.mix2 +
        rowSums((d0A.d - d0B.PD.mix)^2 / DELTA.d.mix)

    }  # (tmp3.TF[ 9] && ld.mix > 0)




    if (tmp3.TF[10] && ld.mix > 1) {  # \calD_{p}, includes \theta_d
      ned2l.dmunb.p.munb.d <- (-pdip.mix) * Numer *
        rowSums(d1A.d1 * d1B.PD.mix1 / DELTA.d.mix)  # nnn.
      wz[, iam(1, posn.munb.d, M)] <- ned2l.dmunb.p.munb.d *
        dmunb.p.deta * dmunb.d.deta  # All links done here

      ned2l.dmunb.p.size.d <- (-pdip.mix) * Numer *
        rowSums(d1A.d2 * d1B.PD.mix1 / DELTA.d.mix)
      wz[, iam(1, posn.size.d, M)] <- ned2l.dmunb.p.size.d *
        dmunb.p.deta * dsize.d.deta  # All links done here

      ned2l.dsize.p.munb.d <- (-pdip.mix) * Numer *
        rowSums(d1A.d1 * d1B.PD.mix2 / DELTA.d.mix)  # ddd.
      wz[, iam(2, posn.munb.d, M)] <- ned2l.dsize.p.munb.d *
        dsize.p.deta * dmunb.d.deta  # All links done here

      ned2l.dsize.p.size.d <- (-pdip.mix) * Numer *
        rowSums(d1A.d2 * d1B.PD.mix2 / DELTA.d.mix)
      wz[, iam(2, posn.size.d, M)] <- ned2l.dsize.p.size.d *
        dsize.p.deta * dsize.d.deta  # All links done here
        
      if (tmp3.TF[ 3]) {  # tmp3.TF[ 9] is TRUE, given tmp3.TF[10]
        ned2l.dpobs.mix.munb.d <-
          rowSums(pdip.mix * d1A.d1 * d0B.PD.mix / DELTA.d.mix)  # nnn.
        wz[, iam(posn.pobs.mix, posn.munb.d, M)] <-
          ned2l.dpobs.mix.munb.d  # link done later
      }
        
      if (tmp3.TF[ 3]) {  # tmp3.TF[ 9] is TRUE, given tmp3.TF[10]
        ned2l.dpobs.mix.size.d <-
          rowSums(pdip.mix * d1A.d2 * d0B.PD.mix / DELTA.d.mix)  # ddd.
        wz[, iam(posn.pobs.mix, posn.size.d, M)] <-
          ned2l.dpobs.mix.size.d  # link done later
      }

      if (tmp3.TF[ 6]) {
        ned2l.dpstr.mix.munb.d <- rowSums(
          pdip.mix * d1A.d1 * d0B.PD.mix / DELTA.d.mix)
        wz[, iam(posn.pstr.mix, posn.munb.d, M)] <-
          ned2l.dpstr.mix.munb.d  # * dmunb.i.deta done later
      }

      if (tmp3.TF[ 6]) {
        ned2l.dpstr.mix.size.d <- rowSums(
          pdip.mix * d1A.d2 * d0B.PD.mix / DELTA.d.mix)
        wz[, iam(posn.pstr.mix, posn.size.d, M)] <-
          ned2l.dpstr.mix.size.d  # * dsize.i.deta done later
      }

        ned2l.dpdip.mix.munb.d <- rowSums(
          d1A.d1 * (1 + pdip.mix * (d0A.d - d0B.PD.mix) / DELTA.d.mix))
        wz[, iam(posn.pdip.mix, posn.munb.d, M)] <-
          ned2l.dpdip.mix.munb.d  # * dmunb.d.deta done later

        ned2l.dpdip.mix.size.d <- rowSums(
          d1A.d2 * (1 + pdip.mix * (d0A.d - d0B.PD.mix) / DELTA.d.mix))
        wz[, iam(posn.pdip.mix, posn.size.d, M)] <-
          ned2l.dpdip.mix.size.d  # * dsize.d.deta done later

      ned2l.dmunb.d2 <- pdip.mix *
        rowSums(pdip.mix * (d1A.d1^2) / DELTA.d.mix + d2A.d1)  # nnn.
      wz[, iam(posn.munb.d, posn.munb.d, M)] <-
        ned2l.dmunb.d2 * dmunb.d.deta^2

      ned2l.dsize.d2 <- pdip.mix *
        rowSums(pdip.mix * (d1A.d2^2) / DELTA.d.mix + d2A.d2)  # ddd.
      wz[, iam(posn.size.d, posn.size.d, M)] <-
        ned2l.dsize.d2 * dsize.d.deta^2

      ned2l.dmusz.d2 <- pdip.mix *
        rowSums(pdip.mix * d1A.d1 * d1A.d2 / DELTA.d.mix + d2A.d4)
      wz[, iam(posn.munb.d, posn.size.d, M)] <-
        ned2l.dmusz.d2 * dmunb.d.deta * dsize.d.deta

      if (tmp3.TF[12]) {
        ned2l.dpobs.mlm.munb.d <- rowSums(
           pdip.mix * d0B.PD.mix * d1A.d1 / DELTA.d.mix)  # ddd.
        for (uuu in seq(la.mlm))
          wz[, iam(posn.pobs.mlm - 1 + uuu, posn.munb.d, M)] <-
            ned2l.dpobs.mlm.munb.d  # * dmunb.d.deta done later
      }

      if (tmp3.TF[12]) {
        ned2l.dpobs.mlm.size.d <- rowSums(
           pdip.mix * d0B.PD.mix * d1A.d2 / DELTA.d.mix)  # ddd.
        for (uuu in seq(la.mlm))
          wz[, iam(posn.pobs.mlm - 1 + uuu, posn.size.d, M)] <-
            ned2l.dpobs.mlm.size.d  # * dsize.d.deta done later
      }

    }  # (tmp3.TF[10] && ld.mix > 1)



        
    if (tmp3.TF[13] && li.mlm > 0) {  # \calI_{np}, includes \phi_s.

      if (la.mix && tmp3.TF[ 3])
        ned2l.dpobs.mix.munb.p <-  # ccc
        ned2l.dpobs.mix.munb.p +
          rowSums(d1B.PI.mlm1 * (1 - Numer * d0B.PI.mlm / DELTA.i.mlm))

      if (la.mix && tmp3.TF[ 3])
        ned2l.dpobs.mix.size.p <-  # ddd
        ned2l.dpobs.mix.size.p +
          rowSums(d1B.PI.mlm2 * (1 - Numer * d0B.PI.mlm / DELTA.i.mlm))

      ned2l.dpstr.mix.munb.p <-  # ccc.
      ned2l.dpstr.mix.munb.p + rowSums(
        d1B.PI.mlm1 * (1 - Numer * d0B.PI.mlm / DELTA.i.mlm))

      ned2l.dpstr.mix.size.p <-  # ddd
      ned2l.dpstr.mix.size.p + rowSums(
        d1B.PI.mlm2 * (1 - Numer * d0B.PI.mlm / DELTA.i.mlm))

      if (tmp3.TF[ 9])
        ned2l.dpdip.mix.munb.p <-
        ned2l.dpdip.mix.munb.p - rowSums(
          d1B.PI.mlm1 * (1 - Numer * d0B.PI.mlm / DELTA.i.mlm))

      if (tmp3.TF[ 9])
        ned2l.dpdip.mix.size.p <-
        ned2l.dpdip.mix.size.p - rowSums(
          d1B.PI.mlm2 * (1 - Numer * d0B.PI.mlm / DELTA.i.mlm))

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

    }  # tmp3.TF[13] && li.mlm > 0




    if (tmp3.TF[14] && ld.mlm > 0) {  # \calD_{np}, includes \psi_s.


      if (la.mix && tmp3.TF[ 3])
        ned2l.dpobs.mix.munb.p <-  # nnn.
        ned2l.dpobs.mix.munb.p +
          rowSums(d1B.PD.mlm1 * (1 - Numer * d0B.PD.mlm / DELTA.d.mlm))

      if (la.mix && tmp3.TF[ 3])
        ned2l.dpobs.mix.size.p <-  # ddd.
        ned2l.dpobs.mix.size.p +
          rowSums(d1B.PD.mlm2 * (1 - Numer * d0B.PD.mlm / DELTA.d.mlm))

      ned2l.dpstr.mix.munb.p <-  # nnn.
      ned2l.dpstr.mix.munb.p + rowSums(
        d1B.PD.mlm1 * (1 - Numer * d0B.PD.mlm / DELTA.d.mlm))

      ned2l.dpstr.mix.size.p <-  # ddd.
      ned2l.dpstr.mix.size.p + rowSums(
        d1B.PD.mlm2 * (1 - Numer * d0B.PD.mlm / DELTA.d.mlm))

      if (tmp3.TF[ 9])
        ned2l.dpdip.mix.munb.p <-
        ned2l.dpdip.mix.munb.p - rowSums(
          d1B.PD.mlm1 * (1 - Numer * d0B.PD.mlm / DELTA.d.mlm))

      if (tmp3.TF[ 9])
        ned2l.dpdip.mix.size.p <-
        ned2l.dpdip.mix.size.p - rowSums(
          d1B.PD.mlm2 * (1 - Numer * d0B.PD.mlm / DELTA.d.mlm))

      if (!is.na(posn.pstr.mix)) {
        ned2l.dpstr.mix2 <-
        ned2l.dpstr.mix2 + rowSums(d0B.PD.mlm^2 / DELTA.d.mlm)
      }

      if (all(tmp3.TF[c(6, 9)]))
        ned2l.dpstr.mix.pdip.mix <-
        ned2l.dpstr.mix.pdip.mix - rowSums(d0B.PD.mlm^2 / DELTA.d.mlm)

      if (!is.na(posn.pdip.mix)) {
        ned2l.dpdip.mix2 <-
        ned2l.dpdip.mix2 + rowSums(d0B.PD.mlm^2 / DELTA.d.mlm)
      }
    }  # tmp3.TF[14] && ld.mlm > 0







      




    if (!is.na(posn.pobs.mix))  # Optional (1, 3) element:
      wz[, iam(1, posn.pobs.mix, M)] <-
        ned2l.dpobs.mix.munb.p  # One link done later

    if (!is.na(posn.pobs.mix))  # Optional (2, 3) element:
      wz[, iam(2, posn.pobs.mix, M)] <-
        ned2l.dpobs.mix.size.p  # One link done later

    if (!is.na(posn.pstr.mix))  # Optional (1, 6) element
      wz[, iam(1, posn.pstr.mix, M)] <-
        ned2l.dpstr.mix.munb.p  # One link done later

    if (!is.na(posn.pstr.mix))  # Optional (2, 6) element
      wz[, iam(2, posn.pstr.mix, M)] <-
        ned2l.dpstr.mix.size.p  # One link done later

    if (!is.na(posn.pdip.mix))  # Optional (1, 9) element
      wz[, iam(1, posn.pdip.mix, M)] <-
        ned2l.dpdip.mix.munb.p  # One link done later

    if (!is.na(posn.pdip.mix))  # Optional (2, 9) element
      wz[, iam(2, posn.pdip.mix, M)] <-
        ned2l.dpdip.mix.size.p  # One link done later

    if (!is.na(posn.pstr.mix) &&
        !is.na(posn.pdip.mix))  # Optional (6, 9) element
      wz[, iam(posn.pstr.mix, posn.pdip.mix, M)] <-
        ned2l.dpstr.mix.pdip.mix  # Links done later

    if (!is.na(posn.pstr.mix))  # Optional (6, 6) element
      wz[, iam(posn.pstr.mix,  # Link done later
               posn.pstr.mix, M)] <- ned2l.dpstr.mix2

    if (!is.na(posn.pdip.mix))  # Optional (9, 9) element
      wz[, iam(posn.pdip.mix,  # Link done later
               posn.pdip.mix, M)] <- ned2l.dpdip.mix2







    if (tmp3.TF[12] && la.mlm) {  # \calA_{np}, includes \omega_s
      ofset <- posn.pobs.mlm - 1  # 11 for GAITD combo
      for (uuu in seq(la.mlm)) {  # Diagonal elts only
        wz[, iam(ofset + uuu,
                 ofset + uuu, M)] <- 1 / pobs.mlm[, uuu]
      }  # uuu

      tmp8a <- probns / Numer^2
      if (tmp3.TF[ 6] && li.mix)
        tmp8a <- tmp8a + rowSums((d0B.PI.mix^2) / DELTA.i.mix)
      if (tmp3.TF[13] && li.mlm)
        tmp8a <- tmp8a + rowSums((d0B.PI.mlm^2) / DELTA.i.mlm)
      if (tmp3.TF[ 9] && ld.mix)
        tmp8a <- tmp8a + rowSums((d0B.PD.mix^2) / DELTA.d.mix)
      if (tmp3.TF[14] && ld.mlm)
        tmp8a <- tmp8a + rowSums((d0B.PD.mlm^2) / DELTA.d.mlm)
      for (uuu in seq(la.mlm))  # All elts
        for (vvv in uuu:la.mlm)
          wz[, iam(ofset + uuu, ofset + vvv, M)] <-
          wz[, iam(ofset + uuu, ofset + vvv, M)] + tmp8a  # All elts
    }  # la.mlm


 

    if (tmp3.TF[12] && la.mlm) {

      init0.i.val <- init0.d.val <- 0
      if (tmp3.TF[13] && li.mlm) init0.i.val <-
        rowSums(d1B.PI.mlm1 * (1 - Numer * d0B.PI.mlm / DELTA.i.mlm))
      if (tmp3.TF[14] && ld.mlm) init0.d.val <-
        rowSums(d1B.PD.mlm1 * (1 - Numer * d0B.PD.mlm / DELTA.d.mlm))
      ned2l.dpobs.mlm.munb.p <- init0.i.val + init0.d.val  # Vector
      if (tmp3.TF[ 6] && li.mix)
        ned2l.dpobs.mlm.munb.p <-
        ned2l.dpobs.mlm.munb.p + rowSums(
          d1B.PI.mix1 * (1 - Numer * d0B.PI.mix / DELTA.i.mix))
      if (tmp3.TF[ 9] && ld.mix)
        ned2l.dpobs.mlm.munb.p <-
        ned2l.dpobs.mlm.munb.p + rowSums(  # nnn
          d1B.PD.mix1 * (1 - Numer * d0B.PD.mix / DELTA.d.mix))
      ofset <- posn.pobs.mlm - 1  # 11 for combo
      for (vvv in seq(la.mlm))  # ccc.
        wz[, iam(1, ofset + vvv, M)] <- ned2l.dpobs.mlm.munb.p



      init0.i.val <- init0.d.val <- 0
      if (tmp3.TF[13] && li.mlm) init0.i.val <-
        rowSums(d1B.PI.mlm2 * (1 - Numer * d0B.PI.mlm / DELTA.i.mlm))
      if (tmp3.TF[14] && ld.mlm) init0.d.val <-
        rowSums(d1B.PD.mlm2 * (1 - Numer * d0B.PD.mlm / DELTA.d.mlm))
      ned2l.dpobs.mlm.size.p <- init0.i.val + init0.d.val  # Vector
      if (tmp3.TF[ 6] && li.mix)
        ned2l.dpobs.mlm.size.p <-
        ned2l.dpobs.mlm.size.p + rowSums(
          d1B.PI.mix2 * (1 - Numer * d0B.PI.mix / DELTA.i.mix))
      if (tmp3.TF[ 9] && ld.mix)
        ned2l.dpobs.mlm.size.p <-
        ned2l.dpobs.mlm.size.p + rowSums(  # ddd
          d1B.PD.mix2 * (1 - Numer * d0B.PD.mix / DELTA.d.mix))
      ofset <- posn.pobs.mlm - 1  # 11 for combo
      for (vvv in seq(la.mlm))  # ccc.
        wz[, iam(2, ofset + vvv, M)] <- ned2l.dpobs.mlm.size.p
    }  # la.mlm > 0






    if (tmp3.TF[13] && li.mlm > 0) {  # \calI_{np}, includes \phi_s
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






    if (tmp3.TF[14] && ld.mlm > 0) {  # \calD_{np}, includes \psi_s
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








    if (tmp3.TF[13] && li.mlm > 0) {
      ned2l.dpstr.mlm.theta1.p <- matrix(0, n, li.mlm)
      for (vvv in seq(li.mlm))
        for (sss in seq(li.mlm))
          ned2l.dpstr.mlm.theta1.p[, vvv] <-
          ned2l.dpstr.mlm.theta1.p[, vvv] +
          d1B.PI.mlm1[, sss] * (1 + Numer *
          (max(0, sss == vvv) - d0B.PI.mlm[, sss]) / (
          DELTA.i.mlm[, sss]))
      if (li.mix && tmp3.TF[ 6])
        ned2l.dpstr.mlm.theta1.p <-
        ned2l.dpstr.mlm.theta1.p +
        rowSums(d1B.PI.mix1 * (1 - Numer * d0B.PI.mix / DELTA.i.mix))
      if (ld.mix && tmp3.TF[ 9])
        ned2l.dpstr.mlm.theta1.p <-  # nnn
        ned2l.dpstr.mlm.theta1.p +
        rowSums(d1B.PD.mix1 * (1 - Numer * d0B.PD.mix / DELTA.d.mix))
      if (ld.mlm && tmp3.TF[14])
        ned2l.dpstr.mlm.theta1.p <-  # nnn.
        ned2l.dpstr.mlm.theta1.p +
        rowSums(d1B.PD.mlm1 * (1 - Numer * d0B.PD.mlm / DELTA.d.mlm))
      ofset <- posn.pstr.mlm - 1
      for (vvv in seq(li.mlm))  # ccc.
        wz[, iam(1, ofset + vvv, M)] <- ned2l.dpstr.mlm.theta1.p[, vvv]


      ned2l.dpstr.mlm.theta2.p <- matrix(0, n, li.mlm)
      for (vvv in seq(li.mlm))
        for (sss in seq(li.mlm))
          ned2l.dpstr.mlm.theta2.p[, vvv] <-
          ned2l.dpstr.mlm.theta2.p[, vvv] +
          d1B.PI.mlm2[, sss] * (1 + Numer *
          (max(0, sss == vvv) - d0B.PI.mlm[, sss]) / (
          DELTA.i.mlm[, sss]))
      if (li.mix && tmp3.TF[ 6])
        ned2l.dpstr.mlm.theta2.p <-
        ned2l.dpstr.mlm.theta2.p +
        rowSums(d1B.PI.mix2 * (1 - Numer * d0B.PI.mix / DELTA.i.mix))
      if (ld.mix && tmp3.TF[ 9])
        ned2l.dpstr.mlm.theta2.p <-  # nnn
        ned2l.dpstr.mlm.theta2.p +
        rowSums(d1B.PD.mix2 * (1 - Numer * d0B.PD.mix / DELTA.d.mix))
      if (ld.mlm && tmp3.TF[14])
        ned2l.dpstr.mlm.theta2.p <-  # nnn.
        ned2l.dpstr.mlm.theta2.p +
        rowSums(d1B.PD.mlm2 * (1 - Numer * d0B.PD.mlm / DELTA.d.mlm))
      ofset <- posn.pstr.mlm - 1
      for (vvv in seq(li.mlm))  # ccc.
        wz[, iam(2, ofset + vvv, M)] <- ned2l.dpstr.mlm.theta2.p[, vvv]
    }  # li.mlm > 0





    if (tmp3.TF[14] && ld.mlm > 0) {
      ned2l.dpdip.mlm.theta1.p <- matrix(0, n, ld.mlm)
      for (vvv in seq(ld.mlm))
        for (sss in seq(ld.mlm))
          ned2l.dpdip.mlm.theta1.p[, vvv] <-
          ned2l.dpdip.mlm.theta1.p[, vvv] -  # Minus
          d1B.PD.mlm1[, sss] * (1 + Numer *
          (max(0, sss == vvv) - d0B.PD.mlm[, sss]) / (
          DELTA.d.mlm[, sss]))
      if (ld.mix && tmp3.TF[ 9])
        ned2l.dpdip.mlm.theta1.p <-
        ned2l.dpdip.mlm.theta1.p -  # Minus
        rowSums(d1B.PD.mix1 * (1 - Numer * d0B.PD.mix / DELTA.d.mix))
      if (li.mix && tmp3.TF[ 6])
        ned2l.dpdip.mlm.theta1.p <-
        ned2l.dpdip.mlm.theta1.p -  # Minus
        rowSums(d1B.PI.mix1 * (1 - Numer * d0B.PI.mix / DELTA.i.mix))
      if (li.mlm && tmp3.TF[13])
        ned2l.dpdip.mlm.theta1.p <-  # nnn.
        ned2l.dpdip.mlm.theta1.p -  # Minus
        rowSums(d1B.PI.mlm1 * (1 - Numer * d0B.PI.mlm / DELTA.i.mlm))
      ofset <- posn.pdip.mlm - 1
      for (vvv in seq(ld.mlm))  # nnn.
        wz[, iam(1, ofset + vvv, M)] <- ned2l.dpdip.mlm.theta1.p[, vvv]


      ned2l.dpdip.mlm.theta2.p <- matrix(0, n, ld.mlm)
      for (vvv in seq(ld.mlm))
        for (sss in seq(ld.mlm))
          ned2l.dpdip.mlm.theta2.p[, vvv] <-
          ned2l.dpdip.mlm.theta2.p[, vvv] -  # Minus
          d1B.PD.mlm2[, sss] * (1 + Numer *
          (max(0, sss == vvv) - d0B.PD.mlm[, sss]) / (
          DELTA.d.mlm[, sss]))
      if (ld.mix && tmp3.TF[ 9])
        ned2l.dpdip.mlm.theta2.p <-
        ned2l.dpdip.mlm.theta2.p -  # Minus
        rowSums(d1B.PD.mix2 * (1 - Numer * d0B.PD.mix / DELTA.d.mix))
      if (li.mix && tmp3.TF[ 6])
        ned2l.dpdip.mlm.theta2.p <-
        ned2l.dpdip.mlm.theta2.p -  # Minus
        rowSums(d1B.PI.mix2 * (1 - Numer * d0B.PI.mix / DELTA.i.mix))
      if (li.mlm && tmp3.TF[13])
        ned2l.dpdip.mlm.theta2.p <-  # ddd.
        ned2l.dpdip.mlm.theta2.p -  # Minus
        rowSums(d1B.PI.mlm2 * (1 - Numer * d0B.PI.mlm / DELTA.i.mlm))
      ofset <- posn.pdip.mlm - 1
      for (vvv in seq(ld.mlm))  # ddd.
        wz[, iam(2, ofset + vvv, M)] <- ned2l.dpdip.mlm.theta2.p[, vvv]
    }  # ld.mlm > 0







    if (li.mlm && li.mix > 1) {

      ned2l.dpstr.mlm.theta1.i <-  # Not a matrix, just a vector
        rowSums(-pstr.mix * d0B.PI.mix * d1A.i1 / DELTA.i.mix)
      for (vvv in seq(li.mlm))
        wz[, iam(posn.munb.i, posn.pstr.mlm - 1 + vvv, M)] <-
          ned2l.dpstr.mlm.theta1.i  # ccc.

      ned2l.dpstr.mlm.theta2.i <-  # Not a matrix, just a vector
        rowSums(-pstr.mix * d0B.PI.mix * d1A.i2 / DELTA.i.mix)
      for (vvv in seq(li.mlm))
        wz[, iam(posn.size.i, posn.pstr.mlm - 1 + vvv, M)] <-
          ned2l.dpstr.mlm.theta2.i  # ddd.
    }  # li.mlm && li.mix > 1






    if (ld.mlm && ld.mix > 1) {

      ned2l.dpdip.mlm.theta1.d <-  # Not a matrix, just a vector
        rowSums(pdip.mix * d0B.PD.mix * d1A.d1 / DELTA.d.mix)
      for (vvv in seq(ld.mlm))
        wz[, iam(posn.munb.d, posn.pdip.mlm - 1 + vvv, M)] <-
          ned2l.dpdip.mlm.theta1.d  # nnn.

      ned2l.dpdip.mlm.theta2.d <-  # Not a matrix, just a vector
        rowSums(pdip.mix * d0B.PD.mix * d1A.d2 / DELTA.d.mix)
      for (vvv in seq(ld.mlm))
        wz[, iam(posn.size.d, posn.pdip.mlm - 1 + vvv, M)] <-
          ned2l.dpdip.mlm.theta2.d  # ddd.
    }  # ld.mlm && ld.mix > 1






    if (ld.mlm && li.mix > 1) {

      ned2l.dpdip.mlm.theta1.i <-  # Not a matrix, just a vector
        rowSums(-pstr.mix * d0B.PI.mix * d1A.i1 / DELTA.i.mix)
      for (vvv in seq(ld.mlm))
        wz[, iam(posn.munb.i, posn.pdip.mlm - 1 + vvv, M)] <-
          ned2l.dpdip.mlm.theta1.i  # nnn.

      ned2l.dpdip.mlm.theta2.i <-  # Not a matrix, just a vector
        rowSums(-pstr.mix * d0B.PI.mix * d1A.i2 / DELTA.i.mix)
      for (vvv in seq(ld.mlm))
        wz[, iam(posn.size.i, posn.pdip.mlm - 1 + vvv, M)] <-
          ned2l.dpdip.mlm.theta2.i  # nnn.
    }  # ld.mlm && li.mix > 1






    if (li.mlm && ld.mix > 1) {

      ned2l.dpstr.mlm.theta1.d <-  # Not a matrix, just a vector
        rowSums(pdip.mix * d0B.PD.mix * d1A.d1 / DELTA.d.mix)
      for (vvv in seq(li.mlm))
        wz[, iam(posn.munb.d, posn.pstr.mlm - 1 + vvv, M)] <-
          ned2l.dpstr.mlm.theta1.d  # nnn.

      ned2l.dpstr.mlm.theta2.d <-  # Not a matrix, just a vector
        rowSums(pdip.mix * d0B.PD.mix * d1A.d2 / DELTA.d.mix)
      for (vvv in seq(li.mlm))
        wz[, iam(posn.size.d, posn.pstr.mlm - 1 + vvv, M)] <-
          ned2l.dpstr.mlm.theta2.d  # nnn.
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
      if (tmp3.TF[ 6] && li.mix)
        ned2l.dpobs.mlm.pstr.mlm <-
        ned2l.dpobs.mlm.pstr.mlm + rowSums(d0B.PI.mix^2 / DELTA.i.mix)
      if (tmp3.TF[ 9] && ld.mix)
        ned2l.dpobs.mlm.pstr.mlm <-  # nnn
        ned2l.dpobs.mlm.pstr.mlm + rowSums(d0B.PD.mix^2 / DELTA.d.mix)
      if (tmp3.TF[14] && ld.mlm)
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
      if (tmp3.TF[ 6] && li.mix)
        ned2l.dpstr.mlm.pdip.mlm <-
        ned2l.dpstr.mlm.pdip.mlm - rowSums(d0B.PI.mix^2 / DELTA.i.mix)
      if (tmp3.TF[ 9] && ld.mix)
        ned2l.dpstr.mlm.pdip.mlm <-  # ddd.
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
      if (tmp3.TF[ 6] && li.mix)
        ned2l.dpobs.mlm.pdip.mlm <-
        ned2l.dpobs.mlm.pdip.mlm - rowSums(d0B.PI.mix^2 / DELTA.i.mix)
      if (tmp3.TF[13] && li.mlm)
        ned2l.dpobs.mlm.pdip.mlm <-
        ned2l.dpobs.mlm.pdip.mlm - rowSums(d0B.PI.mlm^2 / DELTA.i.mlm)
      if (tmp3.TF[ 9] && ld.mix)
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
      if (li.mix)  # tmp3.TF[ 6]
        ned2l.dpobs.mix.pobs.mlm <-
        ned2l.dpobs.mix.pobs.mlm + rowSums(d0B.PI.mix^2 / DELTA.i.mix)
      if (li.mlm)  # tmp3.TF[10]
        ned2l.dpobs.mix.pobs.mlm <-
        ned2l.dpobs.mix.pobs.mlm + rowSums(d0B.PI.mlm^2 / DELTA.i.mlm)
      if (ld.mix)  # tmp3.TF[ 9]   nnn
        ned2l.dpobs.mix.pobs.mlm <-
        ned2l.dpobs.mix.pobs.mlm + rowSums(d0B.PD.mix^2 / DELTA.d.mix)
      if (ld.mlm)  # tmp3.TF[14]   nnn
        ned2l.dpobs.mix.pobs.mlm <-
        ned2l.dpobs.mix.pobs.mlm + rowSums(d0B.PD.mlm^2 / DELTA.d.mlm)

      for (uuu in seq(la.mlm))  # ccc.
        wz[, iam(posn.pobs.mix, posn.pobs.mlm - 1 + uuu, M)] <-
          ned2l.dpobs.mix.pobs.mlm  # Link done later
    }



    if (all(c(la.mix, li.mlm) > 0)) {  # all(tmp3.TF[c(3, 13)])
      if (li.mix)  # tmp3.TF[ 6]
        ned2l.dpobs.mix.pstr.mlm <-
        ned2l.dpobs.mix.pstr.mlm + rowSums(d0B.PI.mix^2 / DELTA.i.mix)
      if (ld.mix)  # tmp3.TF[ 9]
        ned2l.dpobs.mix.pstr.mlm <-  # nnn
        ned2l.dpobs.mix.pstr.mlm + rowSums(d0B.PD.mix^2 / DELTA.d.mix)
      if (ld.mlm)  # tmp3.TF[14]
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





    if (all(c(la.mix, ld.mlm) > 0)) {  # all(tmp3.TF[c(3, 14)])
      if (li.mix)  # tmp3.TF[ 6]
        ned2l.dpobs.mix.pdip.mlm <-
        ned2l.dpobs.mix.pdip.mlm - rowSums(d0B.PI.mix^2 / DELTA.i.mix)
      if (li.mlm)  # tmp3.TF[13]
        ned2l.dpobs.mix.pdip.mlm <-
        ned2l.dpobs.mix.pdip.mlm - rowSums(d0B.PI.mlm^2 / DELTA.i.mlm)
      if (ld.mix)  # tmp3.TF[ 9]
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



    if (all(c(li.mix, la.mlm) > 0)) {  # all(tmp3.TF[c(6, 12)])
      if (li.mlm)  # tmp3.TF[13]
        ned2l.dpobs.mlm.pstr.mix <-
        ned2l.dpobs.mlm.pstr.mix + rowSums(d0B.PI.mlm^2 / DELTA.i.mlm)
      if (ld.mix)  # tmp3.TF[ 9]
        ned2l.dpobs.mlm.pstr.mix <-  # nnn
        ned2l.dpobs.mlm.pstr.mix + rowSums(d0B.PD.mix^2 / DELTA.d.mix)
      if (ld.mlm)  # tmp3.TF[14]
        ned2l.dpobs.mlm.pstr.mix <-  # nnn
        ned2l.dpobs.mlm.pstr.mix + rowSums(d0B.PD.mlm^2 / DELTA.d.mlm)
        ned2l.dpobs.mlm.pstr.mix <-  # tmp3.TF[ 6] && li.mix
        ned2l.dpobs.mlm.pstr.mix -
        rowSums((d0A.i - d0B.PI.mix) * d0B.PI.mix / DELTA.i.mix)

      for (uuu in seq(la.mlm))  # ccc.
        wz[, iam(posn.pstr.mix,
                 posn.pobs.mlm - 1 + uuu, M)] <-
          ned2l.dpobs.mlm.pstr.mix  # Link done later
    }  # all(c(li.mix, la.mlm) > 0








    if (all(c(ld.mix, la.mlm) > 0)) {  # all(tmp3.TF[c(9, 12)])
      if (ld.mlm)  # tmp3.TF[14]
        ned2l.dpobs.mlm.pdip.mix <-
        ned2l.dpobs.mlm.pdip.mix - rowSums(d0B.PD.mlm^2 / DELTA.d.mlm)
      if (li.mix)  # tmp3.TF[ 6]
        ned2l.dpobs.mlm.pdip.mix <-
        ned2l.dpobs.mlm.pdip.mix - rowSums(d0B.PI.mix^2 / DELTA.i.mix)
      if (li.mlm)  # tmp3.TF[13]
        ned2l.dpobs.mlm.pdip.mix <-
        ned2l.dpobs.mlm.pdip.mix - rowSums(d0B.PI.mlm^2 / DELTA.i.mlm)
        ned2l.dpobs.mlm.pdip.mix <-  # all(tmp3.TF[c(9, 12)]) 
        ned2l.dpobs.mlm.pdip.mix +
        rowSums((d0A.d - d0B.PD.mix) * d0B.PD.mix / DELTA.d.mix)

      for (uuu in seq(la.mlm))  # nnn.
        wz[, iam(posn.pdip.mix,
                 posn.pobs.mlm - 1 + uuu, M)] <-
          ned2l.dpobs.mlm.pdip.mix  # Link done later
    }  # all(c(ld.mix, la.mlm) > 0






    if (all(c(li.mix, li.mlm) > 0)) {  # all(tmp3.TF[c(6, 13)])
      for (uuu in seq(li.mlm))  # tmp3.TF[13]
        for (sss in seq(li.mlm))
          ned2l.dpstr.mix.pstr.mlm[, uuu] <-
          ned2l.dpstr.mix.pstr.mlm[, uuu] -
            ((sss == uuu) - d0B.PI.mlm[, sss]) *
              d0B.PI.mlm[, sss] / DELTA.i.mlm[, sss]
      ned2l.dpstr.mix.pstr.mlm <-
      ned2l.dpstr.mix.pstr.mlm -
      rowSums((d0A.i - d0B.PI.mix) * d0B.PI.mix / DELTA.i.mix)
      if (ld.mix)  # tmp3.TF[ 9]
        ned2l.dpstr.mix.pstr.mlm <-  # nnn
        ned2l.dpstr.mix.pstr.mlm + rowSums(d0B.PD.mix^2 / DELTA.d.mix)
      if (ld.mlm)  # tmp3.TF[14]
        ned2l.dpstr.mix.pstr.mlm <-  # nnn
        ned2l.dpstr.mix.pstr.mlm + rowSums(d0B.PD.mlm^2 / DELTA.d.mlm)

      for (uuu in seq(li.mlm))  # Copy it. ccc.
        wz[, iam(posn.pstr.mix,
                 posn.pstr.mlm - 1 + uuu, M)] <-
          ned2l.dpstr.mix.pstr.mlm[, uuu]  # Link done later
    }  # all(c(li.mix, li.mlm) > 0



    if (all(c(ld.mix, ld.mlm) > 0)) {  # all(tmp3.TF[c(9, 14)])

      for (uuu in seq(ld.mlm))  # tmp3.TF[13]
        for (sss in seq(ld.mlm))
          ned2l.dpdip.mix.pdip.mlm[, uuu] <-
          ned2l.dpdip.mix.pdip.mlm[, uuu] -
            ((sss == uuu) - d0B.PD.mlm[, sss]) *
              d0B.PD.mlm[, sss] / DELTA.d.mlm[, sss]
      if (ld.mix)  # tmp3.TF[ 9]
        ned2l.dpdip.mix.pdip.mlm <-
        ned2l.dpdip.mix.pdip.mlm -
        rowSums((d0A.d - d0B.PD.mix) * d0B.PD.mix / DELTA.d.mix)
      if (li.mix)  # tmp3.TF[ 6]
        ned2l.dpdip.mix.pdip.mlm <-
        ned2l.dpdip.mix.pdip.mlm + rowSums(d0B.PI.mix^2 / DELTA.i.mix)
      if (li.mlm)  # tmp3.TF[13]
        ned2l.dpdip.mix.pdip.mlm <-
        ned2l.dpdip.mix.pdip.mlm + rowSums(d0B.PI.mlm^2 / DELTA.i.mlm)

      for (uuu in seq(ld.mlm))  # Copy it. ccc.
        wz[, iam(posn.pdip.mix,
                 posn.pdip.mlm - 1 + uuu, M)] <-
          ned2l.dpdip.mix.pdip.mlm[, uuu]  # Link done later
    }  # all(c(ld.mix, ld.mlm) > 0




    if (all(c(ld.mix, li.mlm) > 0)) {  # all(tmp3.TF[c(9, 13)])

      for (uuu in seq(li.mlm))  # tmp3.TF[13]
        for (sss in seq(li.mlm))
          ned2l.dpdip.mix.pstr.mlm[, uuu] <-
          ned2l.dpdip.mix.pstr.mlm[, uuu] +
            ((sss == uuu) - d0B.PI.mlm[, sss]) *
              d0B.PI.mlm[, sss] / DELTA.i.mlm[, sss]
      if (ld.mix)  # tmp3.TF[ 9]
        ned2l.dpdip.mix.pstr.mlm <-
        ned2l.dpdip.mix.pstr.mlm +
        rowSums((d0A.d - d0B.PD.mix) * d0B.PD.mix / DELTA.d.mix)
      if (li.mix)  # tmp3.TF[ 6]
        ned2l.dpdip.mix.pstr.mlm <-
        ned2l.dpdip.mix.pstr.mlm - rowSums(d0B.PI.mix^2 / DELTA.i.mix)
      if (ld.mlm)  # tmp3.TF[14]
        ned2l.dpdip.mix.pstr.mlm <-
        ned2l.dpdip.mix.pstr.mlm - rowSums(d0B.PD.mlm^2 / DELTA.d.mlm)

      for (uuu in seq(li.mlm))  # Copy it. ccc.
        wz[, iam(posn.pdip.mix,
                 posn.pstr.mlm - 1 + uuu, M)] <-
          ned2l.dpdip.mix.pstr.mlm[, uuu]  # Link done later
    }  # all(c(ld.mix, li.mlm) > 0



    if (all(c(li.mix, ld.mlm) > 0)) {  # all(tmp3.TF[c(6, 14)])

      for (uuu in seq(ld.mlm))  # tmp3.TF[14]
        for (sss in seq(ld.mlm))
          ned2l.dpstr.mix.pdip.mlm[, uuu] <-
          ned2l.dpstr.mix.pdip.mlm[, uuu] +
            ((sss == uuu) - d0B.PD.mlm[, sss]) *
              d0B.PD.mlm[, sss] / DELTA.d.mlm[, sss]
      if (li.mix)  # tmp3.TF[ 6]
        ned2l.dpstr.mix.pdip.mlm <-
        ned2l.dpstr.mix.pdip.mlm +
        rowSums((d0A.i - d0B.PI.mix) * d0B.PI.mix / DELTA.i.mix)
      if (ld.mix)  # tmp3.TF[ 9]
        ned2l.dpstr.mix.pdip.mlm <-
        ned2l.dpstr.mix.pdip.mlm - rowSums(d0B.PD.mix^2 / DELTA.d.mix)
      if (li.mlm)  # tmp3.TF[13]
        ned2l.dpstr.mix.pdip.mlm <-  # ddd.
        ned2l.dpstr.mix.pdip.mlm - rowSums(d0B.PI.mlm^2 / DELTA.i.mlm)

      for (uuu in seq(ld.mlm))  # Copy it. ddd.
        wz[, iam(posn.pstr.mix,
                 posn.pdip.mlm - 1 + uuu, M)] <-
          ned2l.dpstr.mix.pdip.mlm[, uuu]  # Link done later
    }  # all(c(li.mix, ld.mlm) > 0)






 
 
 
    if (lall.len) {
      wz.6 <- matrix(0, n, M * (M + 1) / 2)  # Or == 0 * wz
      ind.rc <- setdiff(1:M, ind.munb.z)  # Contiguous rows and
      lind.rc <- length(ind.rc)  # cols of the DAMLM



 # Copy in the thetas values: the looping is overkill.
      for (uuu in ind.munb.z)
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
            if (!any(kay %in% ind.munb.z)) {
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











 


      dstar.deta <- cbind(dmunb.p.deta, dsize.p.deta,
                          if (tmp3.TF[ 4]) dmunb.a.deta else NULL,
                          if (tmp3.TF[ 5]) dsize.a.deta else NULL,
                          if (tmp3.TF[ 7]) dmunb.i.deta else NULL,
                          if (tmp3.TF[ 8]) dsize.i.deta else NULL,
                          if (tmp3.TF[10]) dmunb.d.deta else NULL,
                          if (tmp3.TF[11]) dsize.d.deta else NULL)
      iptr <- 0
      if (length(ind.munb.z))
      for (uuu in ind.munb.z) {  # Could delete 3 for munb.a (orthog)
        iptr <- iptr + 1
        for (ttt in seq(lind.rc)) {
          wz.6[, iam(uuu, ind.rc[ttt], M)] <- 0  # Initialize
          for (sss in seq(lind.rc)) {
            wz.6[, iam(uuu, ind.rc[ttt], M)] <-
            wz.6[, iam(uuu, ind.rc[ttt], M)] +
                allprobs[, sss] * (max(0, sss == ttt) -
                                   allprobs[, ttt]) *
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
        ind.diags <- setdiff(1:M, ind.munb.z)  # Exclude thetas
        wz[atiny, ind.diags] <- .Machine$double.eps +
        wz[atiny, ind.diags] * (1 + .Machine$double.eps^0.5)
      }
    }  # lall.len




    c(w) * wz
  }), list( .truncate = truncate, .lsize.p = lsize.p ,
            .cutoff.prob = cutoff.prob,
            .nbd.max.support = nbd.max.support,
            .max.chunk.MB = max.chunk.MB,
            .eps.trig = eps.trig,
            .nsimEIM = nsimEIM ))))
}  # gaitdnbinomial






























KLDvglm <-
  function(object, ...) {

  infos.list <- object@family@infos()
  min.support.pmf.p <- infos.list$Support[1]
  specvals <- specials(object)


  Inside <- sapply(specvals, is.null)
  if (length(Inside) == 7 && all(Inside))
    stop("'object' has no special values. ",
         "Is it a GAITD regression object?")
  if (length(Inside) == 8 && all(Inside[1:7]) &&
      infos.list$max.support == infos.list$Support[2])
    stop("'object' has no special values. ",
         "Is it really a GAITD regression object?")
  use.max.support <- if (is.numeric(infos.list$max.support))
    infos.list$max.support else Inf

  if (!is.numeric(MM1 <- infos.list$MM1))
    MM1 <- 1  # Default really
  if (MM1 > 2)
    stop("Can only handle 1- or 2-parameter distributions")

  etamat <- predict(object)  # n x M
  eta.p <- etamat[, 1:MM1, drop = FALSE]  # n x MM1
  theta.p1 <- as.vector(eta2theta(eta.p[, 1], linkfun(object)[1]))
  theta.p2 <- if (MM1 == 2)
    as.vector(eta2theta(eta.p[, 2], linkfun(object)[2])) else NULL
  theta.p <- cbind(theta.p1, theta.p2)
  colnames(theta.p) <- paste0(infos.list$baseparams.argnames, ".p")
  if (!is.logical(intercept.only <- object@misc$intercept.only))
    stop("cannot determine whether 'object' is intercept-only")
  if (!intercept.only)
    stop("argument 'object' is not intercept-only")

  pobs.mix <- if (length(specvals$a.mix))
    fitted(object, type.fitted = "pobs.mix") else cbind(0, 0)
  pobs.mlm <- if (length(specvals$a.mlm))
    fitted(object, type.fitted = "pobs.mlm") else cbind(0, 0)
  pstr.mix <- if (length(specvals$i.mix))
    fitted(object, type.fitted = "pstr.mix") else cbind(0, 0)
  pstr.mlm <- if (length(specvals$i.mlm))
    fitted(object, type.fitted = "pstr.mlm") else cbind(0, 0)
  pdip.mix <- if (length(specvals$d.mix))
    fitted(object, type.fitted = "pdip.mix") else cbind(0, 0)
  pdip.mlm <- if (length(specvals$d.mlm))
    fitted(object, type.fitted = "pdip.mlm") else cbind(0, 0)


  indeta <- object@extra$indeta
  if (MM1 == 1) {
    theta.a <- if (any(is.na(indeta[ 3, ]))) theta.p else
                 as.vector(eta2theta(etamat[, (indeta[ 3, 1])],
                           linkfun(object)[(indeta[ 3, 1])]))
    theta.i <- if (any(is.na(indeta[ 5, ]))) theta.p else
                 as.vector(eta2theta(etamat[, (indeta[ 5, 1])],
                           linkfun(object)[(indeta[ 5, 1])]))
    theta.d <- if (any(is.na(indeta[ 7, ]))) theta.p else
                 as.vector(eta2theta(etamat[, (indeta[ 7, 1])],
                           linkfun(object)[(indeta[ 7, 1])]))

    theta.a <- cbind(theta.a)
    theta.i <- cbind(theta.i)
    theta.d <- cbind(theta.d)
  } else {
    theta.a <- if (any(is.na(indeta[ 4, ]))) theta.p else
                 cbind(eta2theta(etamat[, (indeta[ 4, 1])],
                                 linkfun(object)[(indeta[ 4, 1])]),
                       eta2theta(etamat[, (indeta[ 5, 1])],
                                 linkfun(object)[(indeta[ 5, 1])]))
    colnames(theta.a) <- paste0(infos.list$baseparams.argnames, ".a")
    theta.i <- if (any(is.na(indeta[ 7, ]))) theta.p else
                 cbind(eta2theta(etamat[, (indeta[ 7, 1])],
                                 linkfun(object)[(indeta[ 7, 1])]),
                       eta2theta(etamat[, (indeta[ 8, 1])],
                                 linkfun(object)[(indeta[ 8, 1])]))
    colnames(theta.i) <- paste0(infos.list$baseparams.argnames, ".i")
    theta.d <- if (any(is.na(indeta[10, ]))) theta.p else
                 cbind(eta2theta(etamat[, (indeta[10, 1])],
                                 linkfun(object)[(indeta[10, 1])]),
                       eta2theta(etamat[, (indeta[11, 1])],
                                 linkfun(object)[(indeta[11, 1])]))
    colnames(theta.d) <- paste0(infos.list$baseparams.argnames, ".d")
  }


  flip.args <- object@family@infos()$flip.args
  if (is.null(flip.args)) flip.args <- FALSE
 # zz for checking, delete:

  xlim.vec <- c(min.support.pmf.p,
                min(use.max.support, 
                    max(c(depvar(object)))))  # IMPORTANT
  moreinfo <-
  dgaitdplot(theta.p[1, ],  # Reverse ordering may be needed.
       fam = infos.list$parent.name[2],
       a.mix = specvals$a.mix, i.mix = specvals$i.mix, 
       d.mix = specvals$d.mix,
       a.mlm = specvals$a.mlm, i.mlm = specvals$i.mlm,
       d.mlm = specvals$d.mlm,
       truncate = specvals$truncate,
       theta.a = theta.a[1, ],  # Reverse ordering may be needed.
       theta.i = theta.i[1, ],
       theta.d = theta.d[1, ],
       max.support = use.max.support,
       pobs.mix = pobs.mix[1, ],
       pobs.mlm = pobs.mlm[1, ],
       pstr.mix = pstr.mix[1, ],
       pstr.mlm = pstr.mlm[1, ],
       pdip.mix = pdip.mix[1, ],  # 1-coln matrix
       pdip.mlm = pdip.mlm[1, ],
       byrow.aid = TRUE,  # Important really here 20201008
       baseparams.argnames = infos.list$baseparams.argnames,
       nparams = object@family@infos()$MM1,  # Unnecessary?
       flip.args = ifelse(is.logical(flip.args), flip.args, FALSE),
       xlim = xlim.vec, plot.it = FALSE,  # Both IMPORTANT
       new.plot = FALSE,
       ...)

  pmf.p.hat <- moreinfo$unsc.parent  # This has labels (names() works)
  pmf.z.hat <- moreinfo$pmf.z


  probns    <- c(fitted(object, type.fitted = "nonspecial"))[1]
  Numer     <- c(fitted(object, type.fitted = "Numer"))[1]
  Denom.p   <- c(fitted(object, type.fitted = "Denom.p"))[1]
  Delta <- Numer / Denom.p
  klsum <- Delta * log(Delta) * probns



  if (length(specvals$a.mlm)) {
    pmf.p.mlm.a <- pmf.p.hat[as.character(specvals$a.mlm)]
    ind.mlm.a <- match(specvals$a.mlm, moreinfo$x)
    pobs.mlm  <- pmf.z.hat[ind.mlm.a]  # pmf.z.mlm.a <-
    klsum <- klsum + sum(pobs.mlm * log(pobs.mlm / pmf.p.mlm.a))
  }
  if (length(specvals$i.mlm)) {
    pstr.mlm  <-  (fitted(object, type.fitted = "pstr.mlm"))[1, ]
    pmf.p.mlm.i <- pmf.p.hat[as.character(specvals$i.mlm)]

    ind.mlm.i <- match(specvals$i.mlm, moreinfo$x)
    sum.mlm.i <- pmf.z.hat[ind.mlm.i]


    klsum <- klsum + sum(sum.mlm.i * log(sum.mlm.i / pmf.p.mlm.i))
  }
  if (length(specvals$d.mlm)) {  # zz decreases klsum ??
    pmf.p.mlm.d <- pmf.p.hat[as.character(specvals$d.mlm)]

    ind.mlm.d <- match(specvals$d.mlm, moreinfo$x)
    sum.mlm.d <- pmf.z.hat[ind.mlm.d]


    klsum <- klsum + sum(sum.mlm.d * log(sum.mlm.d / pmf.p.mlm.d))
  }


  if (length(specvals$a.mix)) {
    pmf.p.mix.a <- pmf.p.hat[as.character(specvals$a.mix)]
    ind.mix.a <- match(specvals$a.mix, moreinfo$x)
    Pobs.mix <- pmf.z.hat[ind.mix.a]
    klsum <- klsum + sum(Pobs.mix * log(Pobs.mix / pmf.p.mix.a))
  }
  if (length(specvals$i.mix)) {
    pmf.p.mix.i <- pmf.p.hat[as.character(specvals$i.mix)]
    ind.mix.i <- match(specvals$i.mix, moreinfo$x)
    sum.mix.i <- pmf.z.hat[ind.mix.i]


    klsum <- klsum + sum(sum.mix.i * log(sum.mix.i / pmf.p.mix.i))
  }
  if (length(specvals$d.mix)) {  # zz decreases klsum ??
    pmf.p.mix.d <- pmf.p.hat[as.character(specvals$d.mix)]
    ind.mix.d <- match(specvals$d.mix, moreinfo$x)
    sum.mix.d <- pmf.z.hat[ind.mix.d]


    klsum <- klsum + sum(sum.mix.d * log(sum.mix.d / pmf.p.mix.d))
  }

  klsum
}  # KLDvglm






  setGeneric("KLD", function(object, ...)
             standardGeneric("KLD"),
             package = "VGAM")
setMethod("KLD", "vglm",
          function(object, ...) {
            KLDvglm(object, ...)
          })









Pheapseep <-
  function(object, ...) {
  infos.list <- object@family@infos()
  min.support.pmf.p <- infos.list$Support[1]
  specvals <- specials(object)


  Inside <- sapply(specvals, is.null)
  if (length(Inside) == 7 && all(Inside))
    stop("'object' has no special values. ",
         "Is it a GAITD regression object?")
  if (length(Inside) == 8 && all(Inside[1:7]) &&
      infos.list$max.support == infos.list$Support[2])
    stop("'object' has no special values. ",
         "Is it really a GAITD regression object?")
  use.max.support <- if (is.numeric(infos.list$max.support))
    infos.list$max.support else Inf

  if (!is.numeric(MM1 <- infos.list$MM1))
    MM1 <- 1  # Default really
  if (MM1 > 2)
    stop("Can only handle 1- or 2-parameter distributions")

  etamat <- predict(object)  # n x M
  eta.p <- etamat[, 1:MM1, drop = FALSE]  # n x MM1
  theta.p1 <- as.vector(eta2theta(eta.p[, 1], linkfun(object)[1]))
  theta.p2 <- if (MM1 == 2)
    as.vector(eta2theta(eta.p[, 2], linkfun(object)[2])) else NULL
  theta.p <- cbind(theta.p1, theta.p2)
  colnames(theta.p) <- paste0(infos.list$baseparams.argnames, ".p")
  if (!is.logical(intercept.only <- object@misc$intercept.only))
    stop("cannot determine whether 'object' is intercept-only")
  if (!intercept.only)
    stop("argument 'object' is not intercept-only")

  pobs.mix <- if (length(specvals$a.mix))
    fitted(object, type.fitted = "pobs.mix") else cbind(0, 0)
  pobs.mlm <- if (length(specvals$a.mlm))
    fitted(object, type.fitted = "pobs.mlm") else cbind(0, 0)
  pstr.mix <- if (length(specvals$i.mix))
    fitted(object, type.fitted = "pstr.mix") else cbind(0, 0)
  pstr.mlm <- if (length(specvals$i.mlm))
    fitted(object, type.fitted = "pstr.mlm") else cbind(0, 0)
  pdip.mix <- if (length(specvals$d.mix))
    fitted(object, type.fitted = "pdip.mix") else cbind(0, 0)
  pdip.mlm <- if (length(specvals$d.mlm))
    fitted(object, type.fitted = "pdip.mlm") else cbind(0, 0)


  indeta <- object@extra$indeta
  if (MM1 == 1) {
    theta.a <- if (any(is.na(indeta[ 3, ]))) theta.p else
                 as.vector(eta2theta(etamat[, (indeta[ 3, 1])],
                           linkfun(object)[(indeta[ 3, 1])]))
    theta.i <- if (any(is.na(indeta[ 5, ]))) theta.p else
                 as.vector(eta2theta(etamat[, (indeta[ 5, 1])],
                           linkfun(object)[(indeta[ 5, 1])]))
    theta.d <- if (any(is.na(indeta[ 7, ]))) theta.p else
                 as.vector(eta2theta(etamat[, (indeta[ 7, 1])],
                           linkfun(object)[(indeta[ 7, 1])]))

    theta.a <- cbind(theta.a)
    theta.i <- cbind(theta.i)
    theta.d <- cbind(theta.d)
  } else {
    theta.a <- if (any(is.na(indeta[ 4, ]))) theta.p else
                 cbind(eta2theta(etamat[, (indeta[ 4, 1])],
                                 linkfun(object)[(indeta[ 4, 1])]),
                       eta2theta(etamat[, (indeta[ 5, 1])],
                                 linkfun(object)[(indeta[ 5, 1])]))
    colnames(theta.a) <- paste0(infos.list$baseparams.argnames, ".a")
    theta.i <- if (any(is.na(indeta[ 7, ]))) theta.p else
                 cbind(eta2theta(etamat[, (indeta[ 7, 1])],
                                 linkfun(object)[(indeta[ 7, 1])]),
                       eta2theta(etamat[, (indeta[ 8, 1])],
                                 linkfun(object)[(indeta[ 8, 1])]))
    colnames(theta.i) <- paste0(infos.list$baseparams.argnames, ".i")
    theta.d <- if (any(is.na(indeta[10, ]))) theta.p else
                 cbind(eta2theta(etamat[, (indeta[10, 1])],
                                 linkfun(object)[(indeta[10, 1])]),
                       eta2theta(etamat[, (indeta[11, 1])],
                                 linkfun(object)[(indeta[11, 1])]))
    colnames(theta.d) <- paste0(infos.list$baseparams.argnames, ".d")
  }


  flip.args <- object@family@infos()$flip.args
  if (is.null(flip.args)) flip.args <- FALSE
 # zz for checking, delete:

  xlim.vec <- c(min.support.pmf.p,
                min(use.max.support, 
                    max(c(depvar(object)))))  # IMPORTANT

  moreinfo <-
  dgaitdplot(theta.p[1, ],  # Reverse ordering may be needed.
       fam = infos.list$parent.name[2],
       a.mix = specvals$a.mix, i.mix = specvals$i.mix, 
       d.mix = specvals$d.mix,
       a.mlm = specvals$a.mlm, i.mlm = specvals$i.mlm,
       d.mlm = specvals$d.mlm,
       truncate = specvals$truncate,
       theta.a = theta.a[1, ],  # Reverse ordering may be needed.
       theta.i = theta.i[1, ],
       theta.d = theta.d[1, ],
       max.support = use.max.support,
       pobs.mix = pobs.mix[1, ],
       pobs.mlm = pobs.mlm[1, ],
       pstr.mix = pstr.mix[1, ],
       pstr.mlm = pstr.mlm[1, ],
       pdip.mix = pdip.mix[1, ],  # 1-coln matrix
       pdip.mlm = pdip.mlm[1, ],
       byrow.aid = TRUE,  # Important really here 20201008
       baseparams.argnames = infos.list$baseparams.argnames,
       nparams = object@family@infos()$MM1,  # Unnecessary?
       flip.args = ifelse(is.logical(flip.args), flip.args, FALSE),
       xlim = xlim.vec, plot.it = FALSE,  # Both IMPORTANT
       new.plot = FALSE,
       ...)
  pmf.p.hat <- moreinfo$unsc.parent  # This has labels (names() works)
  pmf.z.hat <- moreinfo$pmf.z


  probns    <- c(fitted(object, type.fitted = "nonspecial"))[1]
  Numer     <- c(fitted(object, type.fitted = "Numer"))[1]
  Denom.p   <- c(fitted(object, type.fitted = "Denom.p"))[1]
  Delta <- Numer / Denom.p
  pheapseepsum.a <- pheapseepsum.i <- pheapseepsum.d <- 0


  
 

  if (length(specvals$a.mlm)) {
    pmf.p.mlm.a <- pmf.p.hat[as.character(specvals$a.mlm)]
    ind.mlm.a <- match(specvals$a.mlm, moreinfo$x)
    pobs.mlm  <- pmf.z.hat[ind.mlm.a]  # pmf.z.mlm.a <-
    pheapseepsum.a <- pheapseepsum.a +
      sum(abs(pobs.mlm - Delta * pmf.p.mlm.a))
  }
  if (length(specvals$i.mlm)) {
    pstr.mlm <-(fitted(object, type.fitted = "pstr.mlm"))[1, ]
    pheapseepsum.i <- pheapseepsum.i + sum(pstr.mlm)
  }
  if (length(specvals$d.mlm)) {  # zz decreases pheapseepsum ??
    pdip.mlm <- c(fitted(object, type.fitted = "pdip.mlm"))[1]
    pheapseepsum.d <- pheapseepsum.d + sum(pdip.mlm)
  }



  if (length(specvals$a.mix)) {
    Pobs.mix2 <-(fitted(object, type.fitted = "Pobs.mix"))[1, ]
    pmf.p.mix.a <- pmf.p.hat[as.character(specvals$a.mix)]
    ind.mix.a <- match(specvals$a.mix, moreinfo$x)
    Pobs.mix <- pmf.z.hat[ind.mix.a]
      pheapseepsum.a <- pheapseepsum.a +
          sum(abs(Pobs.mix - Delta * pmf.p.mix.a))
  }
  if (length(specvals$i.mix)) {
    pstr.mix <- c(fitted(object, type.fitted = "pstr.mix"))[1]
    pheapseepsum.i <- pheapseepsum.i + pstr.mix
  }
  if (length(specvals$d.mix)) {  # zz decreases pheapseepsum ??
    pdip.mix  <- c(fitted(object, type.fitted = "pdip.mix"))[1]
    pheapseepsum.d <- pheapseepsum.d + pdip.mix
  }

  pheapseepsum.a + max(pheapseepsum.i, pheapseepsum.d)
}  # Pheapseep




























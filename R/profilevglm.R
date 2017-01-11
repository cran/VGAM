















profilevglm <-
  function(object, which = 1:p.vlm, alpha = 0.01,
           maxsteps = 10, del = zmax/5, trace = NULL, ...) {




  Pnames <- names(B0 <- coef(object))
  nonA <- !is.na(B0)
  if (any(is.na(B0)))
    stop("currently cannot handle NA-valued regression coefficients")
  pv0 <- t(as.matrix(B0))  # 1 x p.vlm



  p.vlm <- length(Pnames)
  if (is.character(which))
    which <- match(which, Pnames)
  summ <- summary(object)
  std.err <- coef(summ)[, "Std. Error", drop = FALSE]



  M <- npred(object)
  Xm2 <- model.matrix(object, type = "lm2")  # Could be a 0 x 0 matrix
  if (!length(Xm2))
     Xm2 <- NULL  # Make sure. This is safer
  clist <- constraints(object, type = "lm")  # type = c("lm", "term")



  mf <- model.frame(object)

  Y <- model.response(mf)
  if (!is.factor(Y))
    Y <- as.matrix(Y)


  n.lm <- nobs(object, type = "lm")
  OOO <- object@offset
  if (!length(OOO) || all(OOO == 0))
    OOO <- matrix(0, n.lm, M)



  mt <- attr(mf, "terms")




  Wts <- model.weights(mf)
  if (length(Wts) == 0L)
    Wts <- rep(1, n.lm)  # Safest (uses recycling and is a vector)
  Original.de <- deviance(object)  # Could be NULL
  if (!(use.de <- is.Numeric(Original.de)))
    Original.ll <- logLik(object)
  DispersionParameter <- summ@dispersion
  if (!all(DispersionParameter == 1))
    stop("Currently can only handle dispersion parameters ",
         "that are equal to 1")
  X.lm  <- model.matrix(object, type =  "lm")
  X.vlm <- model.matrix(object, type = "vlm")
  fam <- object@family





  quasi.type <- if (length(tmp3 <- fam@infos()$quasi.type))
    tmp3 else FALSE
  if (quasi.type)
    stop("currently this function cannot handle quasi-type models",
         " or models with an estimated dispersion parameter")


  zmax <- sqrt(qchisq(1 - alpha, 1))
  profName <- "z"


  prof <- vector("list", length = length(which))
  names(prof) <- Pnames[which]


  for (i in which) {
    zi <- 0
    pvi <- pv0
    aa <- nonA
    aa[i] <- FALSE
    X.vlm.i <- X.vlm[, aa, drop = FALSE]
    X.lm.i  <-  X.lm  # Try this


 # This is needed by vglm.fit():
    attr(X.vlm.i, "assign") <- attr(X.vlm, "assign")  # zz; this is wrong!
    attr( X.lm.i, "assign") <- attr( X.lm, "assign")


    if (is.logical(trace))
      object@control$trace <- trace


    pnamesi <- Pnames[i]
    for (sgn in c(-1, 1)) {
      if (is.logical(trace) && trace)
        message("\nParameter: ", pnamesi, " ",
                c("down", "up")[(sgn + 1)/2 + 1])
      step <- 0
      zedd <- 0
      LPmat <- matrix(c(X.vlm[, nonA, drop = FALSE] %*% B0[nonA]),
                      n.lm, M, byrow = TRUE) + OOO


      while ((step <- step + 1) < maxsteps &&
             abs(zedd) < zmax) {
        betai <- B0[i] + sgn * step * del * std.err[Pnames[i], 1]
        ooo <- OOO + matrix(X.vlm[, i] * betai, n.lm, M, byrow = TRUE)




        fm <- vglm.fit(x = X.lm.i,  # Possibly use X.lm.i or else X.lm
                       y = Y, w = Wts,
                       X.vlm.arg = X.vlm.i,  # X.vlm,
                       Xm2 = Xm2, Terms = mt,
                       constraints = clist, extra = object@extra,
                       etastart = LPmat,
                       offset = ooo, family = fam,
                       control = object@control)



        fmc <- fm$coefficients
        LPmat <- matrix(X.vlm.i %*% fmc, n.lm, M, byrow = TRUE) + ooo
        ri <- pv0
        ri[, names(fmc)] <- fmc  # coef(fm)
        ri[, pnamesi] <- betai
        pvi <- rbind(pvi, ri)
        zee <- if (use.de) {
          fm$crit.list[["deviance"]] - Original.de
        } else {
          2 * (Original.ll - fm$crit.list[["loglikelihood"]])
        }
        if (zee > -1e-3) {
          zee <- max(zee, 0)
        } else {
          stop("profiling has found a better solution, ",
               "so original fit had not converged")
        }
        zedd <- sgn * sqrt(zee)
        zi <- c(zi, zedd)
      }  # while
    }  # for sgn
    si. <- order(zi)
    prof[[pnamesi]] <- structure(data.frame(zi[si.]), names = profName)
    prof[[pnamesi]]$par.vals <- pvi[si., ,drop = FALSE]
  }  # for i


  val <- structure(prof, original.fit = object, summary = summ)
  class(val) <- c("profile.glm", "profile")
  val
}






if (!isGeneric("profile"))
    setGeneric("profile",
               function(fitted, ...)
               standardGeneric("profile"),
           package = "VGAM")


setMethod("profile", "vglm",
          function(fitted, ...)
          profilevglm(object = fitted, ...))







vplot.profile <-
  function(x, ...) {
  nulls <- sapply(x, is.null)
  if (all(nulls)) return(NULL)
  x <- x[!nulls]
  nm <- names(x)
  nr <- ceiling(sqrt(length(nm)))
  oldpar <- par(mfrow = c(nr, nr))
  on.exit(par(oldpar))
  for (nm in names(x)) {
    tau <- x[[nm]][[1L]]
    parval <- x[[nm]][[2L]][, nm]
    dev.hold()
    plot(parval, tau, xlab = nm, ylab = "tau", type = "n")
    if (sum(tau == 0) == 1) points(parval[tau == 0], 0, pch = 3)
    splineVals <- spline(parval, tau)
    lines(splineVals$x, splineVals$y)
    dev.flush()
  }
}



vpairs.profile <-
function(x, colours = 2:3, ...) {
  parvals <- lapply(x, "[[", "par.vals")


  rng <- apply(do.call("rbind", parvals), 2L, range, na.rm = TRUE)
  Pnames <- colnames(rng)
  npar <- length(Pnames)
  coefs <- coef(attr(x, "original.fit"))
  form <- paste(as.character(formula(attr(x, "original.fit")))[c(2, 1, 3)],
                collapse = "")
  oldpar <- par(mar = c(0, 0, 0, 0), mfrow = c(1, 1),
                oma = c(3, 3, 6, 3), las = 1)
  on.exit(par(oldpar))
  fin <- par("fin")
  dif <- (fin[2L] - fin[1L])/2
  adj <- if (dif > 0) c(dif, 0, dif, 0) else c(0, -dif, 0, -dif)
  par(omi = par("omi") + adj)
  cex <- 1 + 1/npar
  frame()
  mtext(form, side = 3, line = 3, cex = 1.5, outer = TRUE)
  del <- 1/npar
  for (i in 1L:npar) {
    ci <- npar - i
    pi <- Pnames[i]
    for (j in 1L:npar) {
      dev.hold()
      pj <- Pnames[j]
      par(fig = del * c(j - 1, j, ci, ci + 1))
      if (i == j) {
        par(new = TRUE)
        plot(rng[, pj], rng[, pi], axes = FALSE,
             xlab = "", ylab = "", type = "n")
        op <- par(usr = c(-1, 1, -1, 1))
        text(0, 0, pi, cex = cex, adj = 0.5)
        par(op)
      } else {
        col <- colours
        if (i < j) col <- col[2:1]
        if (!is.null(parvals[[pj]])) {
          par(new = TRUE)
          plot(spline(x <- parvals[[pj]][, pj],
                      y <- parvals[[pj]][, pi]),
               type = "l", xlim = rng[, pj],
               ylim = rng[, pi], axes = FALSE,
               xlab = "", ylab = "", col = col[2L])
          pu <- par("usr")
          smidge <- 2/100 * (pu[4L] - pu[3L])
          segments(x, pmax(pu[3L], y - smidge),
                   x, pmin(pu[4L], y + smidge))
        } else
          plot(rng[, pj], rng[, pi], axes = FALSE,
               xlab = "", ylab = "", type = "n")
        if (!is.null(parvals[[pi]])) {
          lines(x <- parvals[[pi]][, pj],
                y <- parvals[[pi]][, pi],
                type = "l", col = col[1L])
          pu <- par("usr")
          smidge <- 2/100 * (pu[2L] - pu[1L])
          segments(pmax(pu[1L], x - smidge), y,
                   pmin(pu[2L], x + smidge), y)
        }
        points(coefs[pj], coefs[pi], pch = 3, cex = 3)
      }
      if (i == npar) axis(1)
      if (j == 1) axis(2)
      if (i == 1) axis(3)
      if (j == npar) axis(4)
      dev.flush()
    }
  }
  par(fig = c(0, 1, 0, 1))
  invisible(x)
}






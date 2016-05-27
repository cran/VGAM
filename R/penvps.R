# These functions are
# Copyright (C) 1998-2016 T.W. Yee, University of Auckland.
# All rights reserved.



Pen.psv <-
  function(constraints = constraints, ps.list = ps.list) {
  assignx <- ps.list$assignx
  nassignx <- names(assignx)
  indexterms <- ps.list$indexterms


  which.X.ps <- ps.list$which.X.ps

  S.arg <- ps.list$S.arg

  lambdalist <- ps.list$lambdalist

  ridge.adj <- ps.list$ridge.adj

  term.labels <- ps.list$term.labels



  index <- numeric()
  lambda.new <- list()
  pen.new.list <- list()
  ncol.X.ps <- sapply(which.X.ps, length)
  ncolHlist.model <- unlist(lapply(constraints, ncol))


  ncolHlist.new <- ncolHlist.model
  if (names(constraints)[[1]] == "(Intercept)") {
    ncolHlist.new <- ncolHlist.new[-1]
    nassignx <- nassignx[-1]
  }



  ncol.H.ps <- ncolHlist.new[indexterms]
  nps <- nassignx[indexterms]



  lambdalen <- sapply(lambdalist, length)





  for (ii in seq_along(ncol.H.ps)) {
    nlambda <- lambdalen[ii]  # lambdalen[[ii]]
    if (nlambda == ncol.H.ps[ii]) {
      lambda.new[[ii]] <- lambdalist[[ii]]
    } else {
      if (nlambda > ncol.H.ps[ii])
        warning("too many lambdas; using the first few")
      lambda.new[[ii]] <- rep_len(lambdalist[[ii]], ncol.H.ps[ii])
    }

    names(lambda.new)[[ii]] <- nps[ii]  # nps[[ii]]






    if (ridge.adj[[ii]] == 0) {
      lambda.diag <- diag(sqrt(lambda.new[[ii]]))
      pen.noridge <- kronecker(lambda.diag, S.arg[[ii]])
      ooo <- matrix(1:(ncol.H.ps[ii] * ncol.X.ps[ii]),
                    ncol = ncol.X.ps[ii], byrow = TRUE)
      pen.new.list[[ii]] <- pen.noridge[, ooo]
      names(pen.new.list)[[ii]] <- nps[ii]  # nps[[ii]]
    } else {
      ioffset <- 0
      joffset <- 0
      Dmat1 <- matrix(0,
                      ncol.H.ps[ii] * (ncol(S.arg[[ii]]) + nrow(S.arg[[ii]])),
                      ncol.H.ps[ii] *  ncol(S.arg[[ii]]))
      for (jay in 1:ncol.H.ps[ii]) {
        pen.set <- sqrt(lambda.new[[ii]][jay]) * S.arg[[ii]]
        pen.ridge <- rbind(pen.set,
                           sqrt(ridge.adj[[ii]]) * diag(ncol(S.arg[[ii]])))
        Dmat1[ioffset + 1:nrow(pen.ridge),
              joffset + 1:ncol(pen.ridge)] <- pen.ridge
        ioffset <- ioffset + nrow(pen.ridge)
        joffset <- joffset + ncol(pen.ridge)
      }  # for jay
      ooo <- matrix(1:(ncol.H.ps[ii] * ncol.X.ps[ii]),
                    ncol = ncol.X.ps[ii], byrow = TRUE)
      pen.new.list[[ii]] <- Dmat1[, ooo]
      names(pen.new.list)[[ii]] <- nps[ii]  # nps[[ii]]
      ioffset <- 0
      joffset <- 0
    }  # if-else ridge.adj
  }  # for




  ncol.allterms <- sapply(assignx, length)

  ncol.model <- if (names(constraints)[[1]] == "(Intercept)")
                  ncol.allterms[-1] else  ncol.allterms
  nrowpen.new.list <- sapply(pen.new.list, nrow)
  nrowPen <- sum(nrowpen.new.list)
  ncolPen <- sum(ncol.allterms * ncolHlist.model)
  iioffset <- 0
  Dmat2 <- matrix(0, nrowPen, ncolPen)
  jay <- 0


  jjoffset <- if (names(constraints)[[1]] == "(Intercept)")
                ncolHlist.model[1] else 0

  for (ii in seq_along(term.labels)) {
    if (indexterms[ii]) {
      jay <- jay + 1
      ind.x <- iioffset + 1:nrow(pen.new.list[[jay]])
      ind.y <- jjoffset + 1:ncol(pen.new.list[[jay]])
      Dmat2[ind.x, ind.y] <- pen.new.list[[jay]]
      iioffset <- iioffset + nrow(pen.new.list[[jay]])
      jjoffset <- jjoffset + ncol(pen.new.list[[jay]])
    } else {
      jjoffset <- jjoffset + ncolHlist.new[ii] * ncol.model[ii]
    }
  }  # ii

  Xvlm.aug <- Dmat2

  attr(Xvlm.aug, "lambda.vlm") <- lambda.new
  Xvlm.aug 
}





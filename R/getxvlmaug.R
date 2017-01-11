# These functions are
# Copyright (C) 1998-2017 T.W. Yee, University of Auckland.
# All rights reserved.





mroot2 <- function(A) {
  if (!isTRUE(all.equal(A, t(A))))
    stop("Supplied matrix not symmetric")

  U <- chol(A, pivot = TRUE, tol = 0)

  opiv <- order(attr(U, "pivot"))
  r <- attr(U, "rank")
  p <- ncol(U)
  if (r < p) U[(r+1):p, (r+1):p] <- 0
  rank <- r
  U <- U[, opiv, drop = FALSE]
  U
}  # mroot2





mroot3 <- function(A, rank = NULL, transpose = FALSE) {
  if (is.null(rank)) rank <- 0 
  if (!isTRUE(all.equal(A, t(A))))
    stop("Supplied matrix not symmetric")
  U <- suppressWarnings(chol(A, pivot = TRUE, tol = 0))
  piv <- order(attr(U, "pivot"))
  r <- attr(U, "rank")
  p <- ncol(U)
  if (r < p) U[(r+1):p, (r+1):p] <- 0
  if (rank < 1) rank <- r
  U <- U[, piv, drop = FALSE]
  if (transpose) t(U[1:rank, , drop = FALSE]) else
                   U[1:rank, , drop = FALSE]
}  # mroot3





get.X.VLM.aug <-
  function(constraints = constraints, sm.osps.list = sm.osps.list) {
  assignx <- sm.osps.list$assignx
  nassignx <- names(assignx)
  indexterms <- sm.osps.list$indexterms


  which.X.sm.osps <- sm.osps.list$which.X.sm.osps

  S.arg <- sm.osps.list$S.arg

  sparlist <- sm.osps.list$sparlist

  ridge.adj <- sm.osps.list$ridge.adj

  term.labels <- sm.osps.list$term.labels



  spar.new <- list()
  pen.new.list <- list()
  ncol.X.sm.osps <- sapply(which.X.sm.osps, length)
  ncolHlist.model <- unlist(lapply(constraints, ncol))


  ncolHlist.new <- ncolHlist.model
  if (names(constraints)[[1]] == "(Intercept)") {
    ncolHlist.new <- ncolHlist.new[-1]
    nassignx <- nassignx[-1]
  }



  ncol.H.sm.osps <- ncolHlist.new[indexterms]
  nsm.osps <- nassignx[indexterms]



  sparlen <- sapply(sparlist, length)






  for (ii in seq_along(ncol.H.sm.osps)) {
    nspar <- sparlen[ii]  # sparlen[[ii]]



    sparlist.use <- sparlist[[ii]]
    sparlist.use[sparlist.use < 0] <- 0



    spar.new[[ii]] <- if (nspar == ncol.H.sm.osps[ii]) {
      sparlist.use
    } else {
      if (ncol.H.sm.osps[ii] < nspar)
        warning("too many 'spar' values; using the first few")
      rep_len(sparlist.use, ncol.H.sm.osps[ii])
    }
   

    names(spar.new)[[ii]] <- nsm.osps[ii]  # nsm.osps[[ii]]







    if (ridge.adj[[ii]] == 0) {
      spar.diag <- diag(sqrt(spar.new[[ii]]))
      pen.noridge <- kronecker(spar.diag, S.arg[[ii]])
      ooo <- matrix(1:(ncol.H.sm.osps[ii] * ncol.X.sm.osps[ii]),
                    ncol = ncol.X.sm.osps[ii], byrow = TRUE)
      pen.new.list[[ii]] <- pen.noridge[, ooo]
      names(pen.new.list)[[ii]] <- nsm.osps[ii]  # nsm.osps[[ii]]
    } else {
      ioffset <- 0
      joffset <- 0
      Dmat1 <-
        matrix(0,
               ncol.H.sm.osps[ii] * (ncol(S.arg[[ii]]) + nrow(S.arg[[ii]])),
               ncol.H.sm.osps[ii] *  ncol(S.arg[[ii]]))
      for (jay in 1:(ncol.H.sm.osps[ii])) {



        pen.set <- mroot2(sqrt(spar.new[[ii]][jay]) * S.arg[[ii]] +
                          sqrt(ridge.adj[[ii]]) * diag(ncol(S.arg[[ii]])))
        pen.ridge <- rbind(pen.set,
                           sqrt(ridge.adj[[ii]]) * diag(ncol(S.arg[[ii]])))
        Dmat1[ioffset + 1:nrow(pen.ridge),
              joffset + 1:ncol(pen.ridge)] <- pen.ridge
        ioffset <- ioffset + nrow(pen.ridge)
        joffset <- joffset + ncol(pen.ridge)
      }  # for jay


      ooo <- matrix(1:(ncol.H.sm.osps[ii] * ncol.X.sm.osps[ii]),
                    nrow = ncol.H.sm.osps[ii],  # Redundant really
                    ncol = ncol.X.sm.osps[ii], byrow = TRUE)
      pen.new.list[[ii]] <- Dmat1[, c(ooo), drop = FALSE]
      names(pen.new.list)[[ii]] <- nsm.osps[ii]  # nsm.osps[[ii]]
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

  attr(Xvlm.aug, "spar.vlm") <- spar.new
  Xvlm.aug
}  # get.X.VLM.aug





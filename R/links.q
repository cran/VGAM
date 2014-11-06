# These functions are
# Copyright (C) 1998-2014 T.W. Yee, University of Auckland.
# All rights reserved.







ToString <- function(x)
  paste(x, collapse = ", ")







 TypicalVGAMfamilyFunction <-
  function(lsigma = "loge",
           isigma = NULL,
           link.list = list("(Default)" = "identitylink",
                            x2          = "loge",
                            x3          = "logoff",
                            x4          = "multilogit",
                            x5          = "multilogit"),
           earg.list = list("(Default)" = list(),
                            x2          = list(),
                            x3          = list(offset = -1),
                            x4          = list(),
                            x5          = list()),
           gsigma = exp(-5:5),
           parallel = TRUE,
           ishrinkage = 0.95,
           nointercept = NULL, imethod = 1,
           type.fitted = c("mean", "pobs0", "pstr0", "onempstr0"),
           probs.x = c(0.15, 0.85),
           probs.y = c(0.25, 0.50, 0.75),
           mv = FALSE, earg.link = FALSE,
           whitespace = FALSE, bred = FALSE, lss = TRUE,
           oim = FALSE, nsimEIM = 100,
           zero = NULL) {
  NULL
}


TypicalVGAMlinkFunction <-
  function(theta,
           someParameter = 0,
           bvalue = NULL,  # .Machine$double.xmin is an alternative
           inverse = FALSE, deriv = 0,
           short = TRUE, tag = FALSE) {
  NULL
}





care.exp <- function(x,
                     thresh = -log( sqrt( .Machine$double.xmin ) )
                     ) {
  x[x >   thresh]  <-  thresh
  x[x < (-thresh)] <- -thresh
  exp(x)
}







 loge <- function(theta,
                  bvalue = NULL,  # .Machine$double.xmin is an alternative
                  inverse = FALSE, deriv = 0,
                  short = TRUE, tag = FALSE) {


  if (is.character(theta)) {
    string <- if (short)
        paste("loge(",  theta, ")", sep = "") else
        paste("loge(",  theta, ")", sep = "")
    if (tag)
      string <- paste("Log:", string)
    return(string)
  }

  if (!inverse && length(bvalue))
    theta[theta <= 0.0] <- bvalue

  if (inverse) {
    if (deriv > 0) {
      1 / Recall(theta = theta,
                 bvalue = bvalue,
                 inverse = FALSE, deriv = deriv)
    } else {
      exp(theta)
    }
  } else {
    switch(deriv + 1, {
       log(theta)},
       theta,
       theta)
  }
}



 logneg <- function(theta,
                    bvalue = NULL,  # .Machine$double.xmin is an alternative
                    inverse = FALSE, deriv = 0,
                    short = TRUE, tag = FALSE) {


  if (is.character(theta)) {
    string <- if (short)
        paste("log(-(",  theta, "))", sep = "") else
        paste("log(-(",  theta, "))", sep = "")
    if (tag)
      string <- paste("Log negative:", string)
    return(string)
  }

  if (!inverse && length(bvalue))
    theta[theta <= 0.0] <- bvalue

  if (inverse) {
    if (deriv > 0) {
      1 / Recall(theta = theta,
                 bvalue = bvalue,
                 inverse = FALSE, deriv = deriv)
    } else {
      -exp(theta)
    }
  } else {
    switch(deriv + 1, {
       log(-theta)},
       theta,
       theta)
  }
}



 logoff <- function(theta,
                    offset = 0,
                    inverse = FALSE, deriv = 0,
                    short = TRUE, tag = FALSE) {

  if (!is.Numeric(offset))
    stop("bad input for argument 'offset'")

  if (is.character(theta)) {
    string <- if (short)
      paste("logoff(", theta,
            ", offset = ", as.character(offset),
            ")", sep = "") else
      paste("log(",
            as.character(offset),
            "+",
            theta,
            ")", sep = "")
    if (tag) 
      string <- paste("Log with offset:", string) 
    return(string)
  }

  if (inverse) {
    if (deriv > 0) {
      1 / Recall(theta = theta,
                 offset = offset,
                 inverse = FALSE, deriv = deriv)
    } else {
      exp(theta) - offset
    }
  } else {
    switch(deriv + 1,
       log(theta + offset),
       theta + offset,
       theta + offset)
  }
}




 identitylink <- function(theta,
                      inverse = FALSE, deriv = 0,
                      short = TRUE, tag = FALSE) {


  if (is.character(theta)) {
    string <- theta
    if (tag)
      string <- paste("Identity:", string)
    return(string)
  }

  if (inverse) {
    if (deriv > 0) {
      1 / Recall(theta = theta,
                 inverse = FALSE, deriv = deriv)
    } else {
      theta
    }
  } else {
    switch(deriv+1,
       theta,
       theta * 0 + 1,
       theta * 0)
  }
}




 negidentity <- function(theta,
                      inverse = FALSE, deriv = 0,
                      short = TRUE, tag = FALSE) {
  if (is.character(theta)) {
    string <- paste("-", theta, sep = "")
    if (tag) 
      string <- paste("Negative-identity:", string) 
    return(string)
  }

  if (inverse) {
    if (deriv > 0) {
      1 / Recall(theta = theta,
                 inverse = FALSE, deriv = deriv)
    } else {
      -theta
    }
  } else {
    switch(deriv+1,
       -theta,
       theta*0 - 1,
       theta*0)
  }
}





 logit <-
  function(theta,
           bvalue = NULL,
           inverse = FALSE, deriv = 0,
           short = TRUE, tag = FALSE) {
  if (is.character(theta)) {
    string <- if (short) 
        paste("logit(", theta, ")", sep = "") else
        paste("log(",   theta, "/(1-", theta, "))", sep = "")
    if (tag) 
      string <- paste("Logit:", string) 
    return(string)
  }

  if (!inverse && length(bvalue)) {
    theta[theta <= 0.0] <- bvalue
    theta[theta >= 1.0] <- 1.0 - bvalue
  }
  if (inverse) {
    if (deriv > 0) {
      1 / Recall(theta = theta, bvalue = bvalue,
                 inverse = FALSE, deriv = deriv)
    } else {
        yy <- theta
        Neg <- (theta <  0) & !is.na(theta)
        yy[ Neg] <- exp(theta[Neg]) / (1 + exp(theta[Neg]))
        Pos <- (theta >= 0) & !is.na(theta)
        yy[Pos] <- 1 / (1 + exp(-theta[Pos]))
        yy
      }
  } else {
    switch(deriv+1, {
       temp2 <- log(theta) - log1p(-theta)
       if (any(near0.5 <- (abs(theta - 0.5) < 0.000125) & !is.na(theta)))
         temp2[near0.5] <- log(theta[near0.5] / (1 - theta[near0.5]))
       temp2
       },
       exp(log(theta) + log1p(-theta)),
       exp(log(theta) + log1p(-theta)) * (1 - 2 * theta))
  }
}






 loglog <- function(theta,
                    bvalue = NULL,  # .Machine$double.eps is an alternative
                    inverse = FALSE, deriv = 0,
                    short = TRUE, tag = FALSE) {
  if (is.character(theta)) {
    string <- if (short) 
        paste("loglog(",  theta, ")",  sep = "") else
        paste("log(log(", theta, "))", sep = "")
    if (tag) 
      string <- paste("Log-Log:", string) 
    return(string)
  }

  if (!inverse && length(bvalue))
    theta[theta <= 1.0] <- bvalue

  if (inverse) {
    if (deriv > 0) {
      1 / Recall(theta = theta,
                 bvalue = bvalue,
                 inverse = FALSE, deriv = deriv)
    } else {
      exp(exp(theta))
    }
  } else {
    switch(deriv+1, {
           log(log(theta))},
           theta * log(theta),
           {  junk <- log(theta)
              -junk^2 / (1 + junk)
           },
           stop("argument 'deriv' unmatched"))
  }
}





 cloglog <- function(theta,
                     bvalue = NULL,  # .Machine$double.eps is an alternative
                     inverse = FALSE, deriv = 0,
                     short = TRUE, tag = FALSE) {
  if (is.character(theta)) {
    string <- if (short) 
        paste("cloglog(",    theta, ")",  sep = "") else
        paste("log(-log(1-", theta, "))", sep = "")
    if (tag) 
      string <- paste("Complementary log-log:", string) 
    return(string)
  }

  if (!inverse && length(bvalue)) {
    theta[theta <= 0.0] <- bvalue
    theta[theta >= 1.0] <- 1.0 - bvalue
  }

  if (inverse) {
    if (deriv > 0) {
      1 / Recall(theta = theta,
                 bvalue = bvalue,
                 inverse = FALSE, deriv = deriv)
    } else {
      junk <- exp(theta)
      -expm1(-junk)
    }
  } else {
    switch(deriv+1, {
           log(-log1p(-theta)) },
           -(1 - theta) * log1p(-theta),
           {  junk <- log1p(-theta)
              -(1 - theta) * (1 + junk) * junk
           },
           stop("argument 'deriv' unmatched"))
  }
}




 probit <- function(theta,
                    bvalue = NULL,  # .Machine$double.eps is an alternative
                    inverse = FALSE, deriv = 0,
                    short = TRUE, tag = FALSE) {
  if (is.character(theta)) {
    string <- if (short) 
        paste("probit(", theta, ")", sep = "") else
        paste("qnorm(",  theta, ")", sep = "")
    if (tag) 
      string <- paste("Probit:", string) 
    return(string)
  }

  if (!inverse && length(bvalue)) {
    theta[theta <= 0.0] <- bvalue
    theta[theta >= 1.0] <- 1 - bvalue
  }

  if (inverse) {
    if (deriv > 0) {
      1 / Recall(theta = theta,
                 bvalue = bvalue,
                 inverse = FALSE, deriv = deriv)
    } else {
      ans <- pnorm(theta)
      if (is.matrix(theta))
        dim(ans) <- dim(theta)
      ans
    }
  } else {
    switch(deriv+1, {
        ans <- qnorm(theta)
        if (is.matrix(theta))
            dim(ans) <- dim(theta)
        ans
     }, {
       if (is.matrix(theta)) {
         ans <- dnorm(qnorm(theta))
         dim(ans) <- dim(theta)
         ans
       } else dnorm(qnorm(as.vector(theta)))
      }, {
        junk <- qnorm(theta)
        ans <- -junk * dnorm(junk)
        if (is.vector(theta)) ans else
        if (is.matrix(theta)) {
          dim(ans) <- dim(theta)
          ans
        } else {
          warning("can only handle vectors and matrices;",
                  " converting to vector")
          ans
        }
      })
  }
}








 explink <- function(theta,
                     bvalue = NULL,  # .Machine$double.eps is an alternative
                     inverse = FALSE, deriv = 0,
                     short = TRUE, tag = FALSE) {
  if (is.character(theta)) {
    string <- if (short) 
        paste("explink(", theta, ")", sep = "") else
        paste("exp(", theta, ")", sep = "")
    if (tag) 
      string <- paste("Exp:", string) 
    return(string)
  }

  if (!inverse && length(bvalue))
    theta[theta <= 0.0] <- bvalue
  if (inverse) {
    if (deriv > 0) {
      1 / Recall(theta = theta,
                 bvalue = bvalue,
                 inverse = FALSE, deriv = deriv)
    } else {
      log(theta)
    }
  } else {
    switch(deriv+1, {
       exp(theta)},
        1 / exp(theta),
       -1 / exp(theta * 2))
  }
}





 reciprocal <- function(theta,
                        bvalue = NULL,  # .Machine$double.eps is an alternative
                        inverse = FALSE, deriv = 0,
                        short = TRUE, tag = FALSE) {
  if (is.character(theta)) {
    string <- paste("1/", theta, sep = "")
    if (tag) 
      string <- paste("Reciprocal:", string) 
    return(string)
  }

  if (!inverse && length(bvalue))
    theta[theta == 0.0] <- bvalue

  if (inverse) {
    if (deriv > 0) {
      1 / Recall(theta = theta,
                 bvalue = bvalue,
                 inverse = FALSE, deriv = deriv)
    } else {
      1/theta
    }
  } else {
    switch(deriv+1, {
       1/theta},
       -theta^2,
       2*theta^3)
  }
}





 negloge <- function(theta,
                   bvalue = NULL,  # .Machine$double.eps is an alternative
                   inverse = FALSE, deriv = 0,
                   short = TRUE, tag = FALSE) {
  if (is.character(theta)) {
      string <- if (short) 
          paste("negloge(", theta, ")", sep = "") else
          paste("-log(",  theta, ")", sep = "")
      if (tag) 
        string <- paste("Negative log:", string) 
      return(string)
  }


  if (!inverse && length(bvalue))
    theta[theta <= 0.0] <- bvalue
  if (inverse) {
    if (deriv > 0) {
      1 / Recall(theta = theta,
                 bvalue = bvalue,
                 inverse = FALSE, deriv = deriv)
    } else {
      exp(-theta)
    }
  } else {
    switch(deriv+1, {
       -log(theta)},
       -theta,
       theta)
  }
}




 negreciprocal <-
  function(theta,
           bvalue = NULL,  # .Machine$double.eps is an alternative
           inverse = FALSE,
           deriv = 0, short = TRUE, tag = FALSE) {
  if (is.character(theta)) {
    string <- paste("-1/", theta, sep = "")
    if (tag) 
      string <- paste("Negative reciprocal:", string) 
    return(string)
  }


  if (!inverse && length(bvalue))
    theta[theta == 0.0] <- bvalue

  if (inverse) {
    if (deriv > 0) {
      1 / Recall(theta,
                 bvalue = bvalue,
                 inverse = FALSE, deriv = deriv)
    } else {
      -1/theta
    }
  } else {
    switch(deriv+1, {
       -1/theta},
       theta^2,
       2*theta^3)
  }
}



 natural.ig <-
  function(theta,
           bvalue = NULL,  # .Machine$double.eps is an alternative
           inverse = FALSE, deriv = 0,
           short = TRUE, tag = FALSE) {

  if (is.character(theta)) {
    string <- paste("-1/", theta, sep = "")
    if (tag) 
      string <- paste("Negative inverse:", string) 
    return(string)
  }

  if (inverse) {
    if (deriv > 0) {
      1 / negreciprocal(theta,
                        bvalue = bvalue,
                        inverse = FALSE, deriv = deriv)
    } else {
      1 / sqrt(-2*theta)
    }
  } else {
    switch(deriv+1,
       -1 / (2 * theta^2),
       theta^3,
       3 * theta^5)
  }
}






 rhobit <- function(theta,
                    bminvalue = NULL,
                    bmaxvalue = NULL,
                    inverse = FALSE, deriv = 0,
                    short = TRUE, tag = FALSE) {
  if (is.character(theta)) {
    string <- if (short) 
        paste("rhobit(", theta, ")", sep = "") else
        paste("log((1+", theta, ")/(1-", theta, "))", sep = "")
    if (tag) 
      string <- paste("Rhobit:", string) 
    return(string)
  }

  if (!inverse) {
   if (length(bminvalue)) theta[theta <= -1.0] <- bminvalue
   if (length(bmaxvalue)) theta[theta >=  1.0] <- bmaxvalue
  }

  if (inverse) {
    if (deriv > 0) {
      1 / Recall(theta = theta,
                 bminvalue = bminvalue,
                 bmaxvalue = bmaxvalue,
                 inverse = FALSE, deriv = deriv)
    } else {
      junk <- exp(theta)
      expm1(theta) / (junk + 1.0)
    }
  } else {
      switch(deriv+1, {
          log1p(theta) - log1p(-theta)},
          (1 - theta^2) / 2,
          (1 - theta^2)^2 / (4*theta))
  }
}




 fisherz <- function(theta,
                     bminvalue = NULL,
                     bmaxvalue = NULL,
                     inverse = FALSE, deriv = 0,
                     short = TRUE, tag = FALSE) {
  if (is.character(theta)) {
    string <- if (short) 
        paste("fisherz(", theta, ")", sep = "") else
        paste("(1/2) * log((1+", theta, ")/(1-", theta, "))", sep = "")
    if (tag) 
      string <- paste("Fisher's Z transformation:", string) 
    return(string)
  }

  if (!inverse) {
     if (length(bminvalue)) theta[theta <= -1.0] <- bminvalue
     if (length(bmaxvalue)) theta[theta >=  1.0] <- bmaxvalue
  }

  if (inverse) {
    if (deriv > 0) {
      1 / Recall(theta = theta,
                 bminvalue = bminvalue,
                 bmaxvalue = bmaxvalue,
                 inverse = FALSE, deriv = deriv)
    } else {
      tanh(theta)
    }
  } else {
      switch(deriv+1,
         atanh(theta),
         1.0 - theta^2,
         (1.0 - theta^2)^2 / (2*theta))
    }
}







 multilogit <-
  function(theta,
           refLevel = "last",
           M = NULL,  # stop("argument 'M' not specified"),
           whitespace = FALSE,
           bvalue = NULL,
           inverse = FALSE, deriv = 0,
           short = TRUE, tag = FALSE) {
 

  fillerChar <- ifelse(whitespace, " ", "")

  if (length(refLevel) != 1)
    stop("the length of 'refLevel' must be one")

  if (is.character(refLevel)) {
    if (refLevel != "last")
      stop('if a character, refLevel must be "last"')
    refLevel <- -1
  } else
  if (is.factor(refLevel)) {
    if (is.ordered(refLevel))
      warning("argument 'refLevel' is from an ordered factor")
    refLevel <- as.character(refLevel) == levels(refLevel)
    refLevel <- (1:length(refLevel))[refLevel]
    if (!is.Numeric(refLevel, length.arg = 1,
                    integer.valued = TRUE, positive = TRUE))
      stop("could not coerce 'refLevel' into a single positive integer")
  } else
  if (!is.Numeric(refLevel, length.arg = 1,
                  integer.valued = TRUE))
    stop("'refLevel' must be a single (positive?) integer")




  if (is.character(theta)) {
    is.M <- is.finite(M) && is.numeric(M)
    string <- if (short)
        paste("multilogit(", theta, ")", sep = "") else {
         if (refLevel < 0) {
           ifelse(whitespace,
             paste("log(", theta, "[,j] / ",
                   theta, "[,",
                   ifelse(is.M, M+1, "M+1"),
                   "]), j = 1:",
                   ifelse(is.M, M, "M"), sep = ""),
             paste("log(", theta, "[,j]/",
                   theta, "[,",
                   ifelse(is.M, M+1, "M+1"),
                   "]), j=1:",
                   ifelse(is.M, M, "M"), sep = ""))
         } else {
             if (refLevel == 1) {
               paste("log(", theta, "[,", "j]",
                   fillerChar, "/", fillerChar,
                     "", theta, "[,", refLevel, "]), j",
                     fillerChar, "=", fillerChar, "2:",
                   ifelse(is.M, (M+1), "(M+1)"),
                     sep = "")
             } else {
               paste("log(", theta, "[,", "j]", fillerChar, "/",
                     "", theta, "[,", refLevel, "]), j",
                     fillerChar, "=", fillerChar,
                     "c(1:", refLevel-1, ",",
                     fillerChar,
                     refLevel+1, ":",
                     ifelse(is.M, (M+1), "(M+1)"),
                     ")", sep = "")
             }
         }
    }
    if (tag)
      string <- paste("Multinomial logit link:", string)
    return(string)
  }



  M.orig <- M
  M <- if (inverse) ncol(cbind(theta)) else
                    ncol(cbind(theta)) - 1
  if (M < 1)
    ifelse(inverse,
           stop("argument 'eta' should have at least one column"),
           stop("argument 'theta' should have at least two columns"))
  if (is.numeric(M.orig) && M != M.orig) {
    warning("argument 'M' does not seem right but using it")
    M <- M.orig
  }



  if (!inverse && length(bvalue))
    theta[theta <= 0.0] <- bvalue
  if (!inverse && length(bvalue))
    theta[theta >= 1.0] <- 1 - bvalue



  foo <- function(eta, refLevel = -1, M) {
    phat <- if ((refLevel < 0) || (refLevel == M+1)) {
      cbind(care.exp(eta), 1.0)
    } else if ( refLevel == 1) {
      cbind(1.0, care.exp(eta))
    } else {
      use.refLevel <- if ( refLevel < 0) M+1 else refLevel
      etamat <- cbind(eta[, 1:( refLevel - 1)],
                      0.0,
                      eta[, ( refLevel ):M])
      care.exp(etamat)
    }
    ans <- phat / rowSums(phat)
    colnames(ans) <- NULL
    ans
  }


  if (inverse) {
    if (deriv > 0) {
      1 / Recall(theta = theta,
                 refLevel = refLevel,
                 bvalue = bvalue,
                 inverse = FALSE, deriv = deriv)
    } else {
       foo(theta, refLevel, M = M)  # log(theta[, -jay] / theta[, jay])
    }
  } else {
    switch(deriv + 1, {
      ans <- if (refLevel < 0) {
        log(theta[, -ncol(theta)] / theta[, ncol(theta)])
      } else {
        use.refLevel <- if (refLevel < 0) ncol(theta) else refLevel
        log(theta[, -( use.refLevel )] / theta[, use.refLevel ])
      }
      colnames(ans) <- NULL
      ans
      },
      care.exp(log(theta) + log1p(-theta)),
      care.exp(log(theta) + log1p(-theta)) * (1 - 2 * theta))
  }
}  # end of multilogit







fsqrt <- function(theta,  #  = NA  , = NULL,
                  min = 0, max = 1, mux = sqrt(2),
                  inverse = FALSE, deriv = 0,
                  short = TRUE, tag = FALSE) {
  if (!is.Numeric(min, length.arg = 1))
    stop("bad input for 'min' component")
  if (!is.Numeric(max, length.arg = 1))
    stop("bad input for 'max' component")
  if (!is.Numeric(mux, length.arg = 1, positive = TRUE))
    stop("bad input for 'mux' component")
  if (min >= max)
    stop("'min' >= 'max' is not allowed")

  if (is.character(theta)) {
    string <- if (short) 
      paste("fsqrt(", theta, ")", sep = "") else {
      if (abs(mux-sqrt(2)) < 1.0e-10)
        paste("sqrt(2*", theta, ") - sqrt(2*(1-", theta, "))",
              sep = "") else
      paste(as.character(mux),
            " * (sqrt(", theta, "-", min, ") - sqrt(",
            max, "-", theta, "))",
            sep = "")
    }
    if (tag) 
      string <- paste("Folded square root:", string) 
    return(string)
  }

  if (inverse) {
    if (deriv > 0) {
      1 / Recall(theta = theta,
                 min = min, max = max, mux = mux,
                 inverse = FALSE, deriv = deriv)
    } else {
      mid <- (min + max) / 2
      boundary <- mux * sqrt(max - min)
      temp <- pmax(0, (theta/mux)^2 * (2*(max-min) - (theta/mux)^2))
      ans <- theta
      if (any(ind5 <- theta <  0))
        ans[ind5] <- mid - 0.5 * sqrt(temp[ind5])
      if (any(ind5 <- theta >= 0))
        ans[ind5] <- mid + 0.5 * sqrt(temp[ind5])
      ans[theta < -boundary] <- NA
      ans[theta >  boundary] <- NA
      ans
    }
  } else {
    switch(deriv+1,
        mux * (sqrt(theta-min) - sqrt(max-theta)),
       (2 / mux) / (1/sqrt(theta-min) + 1/sqrt(max-theta)),
       -(4 / mux) / ((theta-min)^(-3/2) - (max-theta)^(-3/2)))
  }
}




 powerlink <- function(theta,
                       power = 1,
                       inverse = FALSE, deriv = 0,
                       short = TRUE, tag = FALSE) {
    exponent <- power
    if (exponent == 0)
      stop("use the 'loge' link")

  if (is.character(theta)) {
    string <- if (short) 
        paste("powerlink(", theta, ", power = ",
              as.character(exponent), ")",
              sep = "") else
        paste(theta, "^(", as.character(exponent), ")", sep = "")
    if (tag) 
      string <- paste("Power link:", string)
    return(string)
  }

  if (inverse) {
    if (deriv > 0) {
      1 / Recall(theta = theta,
                 power = power,
                 inverse = FALSE, deriv = deriv)
      } else {
          theta^(1/exponent)
      }
  } else {
    switch(deriv+1,
    {
      theta^exponent
    }, {
      (theta^(1-exponent)) / exponent
    }, {
      (theta^(2-exponent)) / (exponent * (exponent-1))
    })
  }
}





 elogit <- function(theta,
                    min = 0, max = 1,
                    bminvalue = NULL,
                    bmaxvalue = NULL,
                    inverse = FALSE, deriv = 0,
                    short = TRUE, tag = FALSE) {

    A <- min
    B <- max
   if (!inverse && length(bminvalue)) theta[theta <= A] <- bminvalue
   if (!inverse && length(bmaxvalue)) theta[theta >= B] <- bmaxvalue



  if (is.character(theta)) {
    string <- if (short) {
      if (A != 0 || B != 1)
        paste("elogit(", theta,
              ", min = ", A,
              ", max = ", B, ")", sep = "") else
        paste("elogit(", theta, ")", sep = "")
    } else {
      paste("log((", theta, "-min)/(max-", theta, "))", sep = "")
    }
    if (tag) 
      string <- paste("Extended logit:", string) 
    return(string)
  }

  if (inverse) {
    if (deriv > 0) {
      1 / Recall(theta = theta,
                 min = min, max = max,
                 bminvalue = bminvalue,
                 bmaxvalue = bmaxvalue,
                 inverse = FALSE, deriv = deriv)
      } else {
        junk <- care.exp(theta)
        (A + B * junk) / (1.0 + junk)
      }
  } else {
    switch(deriv+1, {
           log((theta - A)/(B - theta))},
           (theta - A) * (B - theta) / (B-A),
           (theta - A) * (B - theta) * (B - 2 * theta + A) / (B-A)^2)
  }
}






 logc <- function(theta,
                  bvalue = NULL,  # .Machine$double.xmin is an alternative
                  inverse = FALSE, deriv = 0,
                  short = TRUE, tag = FALSE) {
  if (is.character(theta)) {
    string <- if (short) 
        paste("logc(", theta, ")", sep = "") else
        paste("log(1-", theta, ")", sep = "")
    if (tag) 
      string <- paste("Log Complementary:", string) 
    return(string)
  }


  if (!inverse && length(bvalue)) {
    theta[theta >= 1.0] <- bvalue;
  }
  if (inverse) {
    if (deriv > 0) {
      1 / Recall(theta = theta,
                 bvalue = bvalue,
                 inverse = FALSE, deriv = deriv)
    } else {
        -expm1(theta)
    }
  } else {
    switch(deriv+1, {
           log1p(-theta)},
           -(1.0 - theta),
           -(1.0 - theta)^2)
  }
}








 cauchit <- function(theta,
                     bvalue = .Machine$double.eps,
                     inverse = FALSE, deriv = 0,
                     short = TRUE, tag = FALSE) {
  if (is.character(theta)) {
    string <- if (short) 
        paste("cauchit(", theta, ")", sep = "") else
        paste("tan(pi*(", theta, "-0.5))", sep = "")
    if (tag) 
      string <- paste("Cauchit:", string) 
    return(string)
  }

  if (!inverse && length(bvalue)) {
    theta[theta <= 0.0] <- bvalue
    theta[theta >= 1.0] <- 1.0 - bvalue
  }
  if (inverse) {
    if (deriv > 0) {
      1 / Recall(theta = theta,
                 bvalue = bvalue,
                 inverse = FALSE, deriv = deriv)
      } else {
        0.5 + atan(theta) / pi
      }
  } else {
      switch(deriv+1,
             tan(pi * (theta-0.5)),
             cos(pi * (theta-0.5))^2 / pi,
            -sin(pi * (theta-0.5) * 2)
            )
  }
}





 golf <- function(theta,
                  lambda = 1,
                  cutpoint = NULL,
                  inverse = FALSE, deriv = 0,
                  short = TRUE, tag = FALSE) {



  if (!is.Numeric(lambda, positive = TRUE))
    stop('could not determine lambda or lambda has negative values')
  if (is.Numeric(cutpoint))
    if (any(cutpoint < 0) ||
        !is.Numeric(cutpoint, integer.valued = TRUE))
    warning("argument 'cutpoint' should contain ",
            "non-negative integer values")

  if (is.character(theta)) {
    string <- if (short) {
      lenl <- length(lambda) > 1
      lenc <- length(cutpoint) > 1
      paste("golf(", theta,
            ", lambda = ",
            if (lenl) "c(" else "",
            ToString(lambda),
            if (lenl) ")" else "",
            if (is.Numeric(cutpoint))
      paste(", cutpoint = ",
            if (lenc) "c(" else "",
      ToString(cutpoint),
            if (lenc) ")" else "",
      sep = "") else "",
                  ")", sep = "")
    } else {
      if (is.Numeric(cutpoint)) {
        paste("-3*log(1-qnorm(", theta,
              ")/(3*sqrt(lambda)))",
              " + log(cutpoint)", sep = "")
      } else {
        paste("-3*log(1-qnorm(", theta,
              ")/(3*sqrt(lambda)))", sep = "")
      }
    }
    if (tag) 
      string <- paste("Gamma-ordinal link function:", string) 
    return(string)
  }


  thmat <- cbind(theta)
  lambda <- rep(lambda, len = ncol(thmat))  # Allow recycling for lambda
  if (is.Numeric(cutpoint))
    cutpoint <- rep(cutpoint, len = ncol(thmat))
  if (ncol(thmat) > 1) {
    answer <- thmat
    for (ii in 1:ncol(thmat))
      answer[, ii] <- Recall(theta = thmat[, ii],
                             lambda = lambda[ii],
                             cutpoint = if (is.Numeric(cutpoint))
                                        cutpoint[ii] else NULL,
                            inverse = inverse, deriv = deriv)
    return(answer)
  }


  answer <- if (inverse) {
    if (deriv > 0) {
      1 / Recall(theta = theta,
                 lambda = lambda,
                 cutpoint = cutpoint,
                 inverse = FALSE, deriv = deriv)
    } else {
      if (is.Numeric(cutpoint)) {
        pnorm((1-care.exp(-(theta-log(cutpoint))/3)) * 3 * sqrt(lambda))
      } else {
        pnorm((1-care.exp(-theta/3)) * 3 * sqrt(lambda))
      }
    }
  } else {
    smallno <- 1 * .Machine$double.eps
    Theta <- theta
    Theta <- pmin(Theta, 1 - smallno)  # Since theta == 1 is a possibility
    Theta <- pmax(Theta, smallno)  # Since theta == 0 is a possibility
    Ql <- qnorm(Theta)
    switch(deriv+1, {
        temp <- Ql / (3*sqrt(lambda))
        temp <- pmin(temp, 1.0 - smallno)  # 100 / .Machine$double.eps
        -3*log1p(-temp) +
        if (is.Numeric(cutpoint)) log(cutpoint) else 0},
        (1 - Ql / (3*sqrt(lambda))) * sqrt(lambda) * dnorm(Ql),
        {  stop('cannot handle deriv = 2') },
        stop("argument 'deriv' unmatched"))
  }
  if (!is.Numeric(answer))
    stop("the answer contains some NAs")
  answer
}





 polf <- function(theta,  # = 1,
                  cutpoint = NULL,
                  inverse = FALSE, deriv = 0,
                  short = TRUE, tag = FALSE) {
  if (!is.Numeric(cutpoint))
    stop("could not determine the cutpoint")
  if (any(cutpoint < 0) ||
      !is.Numeric(cutpoint, integer.valued = TRUE))
    warning("argument 'cutpoint' should",
            " contain non-negative integer values")


  if (is.character(theta)) {
    string <- if (short) {
      lenc <- length(cutpoint) > 1
      paste("polf(", theta,
             ", cutpoint = ",
            if (lenc) "c(" else "",
            ToString(cutpoint),
            if (lenc) ")" else "",
            ")", sep = "") 
    } else
      paste("2*log(0.5*qnorm(", theta,
            ") + sqrt(cutpoint+7/8))", sep = "")
    if (tag) 
      string <- paste("Poisson-ordinal link function:", string) 
    return(string)
  }



    thmat <- cbind(theta)
    if (ncol(thmat) > 1) {
        answer <- thmat
        cutpoint <- rep(cutpoint, len = ncol(thmat))  # Reqd for the for loop
        for (ii in 1:ncol(thmat))
            answer[, ii] <- Recall(theta = thmat[, ii],
                                 cutpoint = cutpoint,
                                 inverse = inverse, deriv = deriv)
        return(answer)
    }

  answer <-
  if (inverse) {
      if (deriv > 0) {
          1 / Recall(theta = theta,
                     cutpoint = cutpoint,
                     inverse = FALSE, deriv = deriv)
      } else {
          if (any(cp.index <- cutpoint == 0)) {
              tmp <- theta
              tmp[cp.index] <- 
              cloglog(theta = theta[cp.index],
                      inverse = inverse, deriv = deriv)
              tmp[!cp.index] <-
                pnorm(2 * exp(theta[!cp.index]/2) -
                      2 * sqrt(cutpoint[!cp.index] + 7/8))
              tmp
          } else {
            pnorm(2 * exp(theta/2) - 2 * sqrt(cutpoint + 7/8))
          }
      }
  } else {
    if (any(cp.index <- cutpoint == 0)) {
        cloglog(theta = theta,
                inverse = inverse, deriv = deriv)
    } else {
      smallno <- 1 * .Machine$double.eps
      SMALLNO <- 1 * .Machine$double.xmin
      Theta <- theta
      Theta <- pmin(Theta, 1 - smallno)  # Since theta == 1 is a possibility
      Theta <- pmax(Theta, smallno)  # Since theta == 0 is a possibility
      Ql <- qnorm(Theta)
      switch(deriv+1, {
      temp <- 0.5 * Ql + sqrt(cutpoint + 7/8)
      temp <- pmax(temp, SMALLNO)
      2 * log(temp)},
      (Ql/2 + sqrt(cutpoint + 7/8)) * dnorm(Ql),
      {  stop('cannot handle deriv = 2') },
      stop("argument 'deriv' unmatched"))
    }
  }
  if (!is.Numeric(answer))
    stop("the answer contains some NAs")
  answer
}





 nbolf <- function(theta,
                   cutpoint = NULL,
                   k = NULL,
                   inverse = FALSE, deriv = 0,
                   short = TRUE, tag = FALSE) {

  kay <- k
  if (!is.Numeric(kay, positive = TRUE))
    stop("could not determine 'k' or it is not positive-valued")
  if (!is.Numeric(cutpoint))
    stop("could not determine the cutpoint")
  if (any(cutpoint < 0) ||
      !is.Numeric(cutpoint, integer.valued = TRUE))
    warning("argument 'cutpoint' should",
            " contain non-negative integer values")

  if (is.character(theta)) {
    string <- if (short) {
        lenc <- length(cutpoint) > 1
        lenk <- length(kay) > 1
        paste("nbolf(", theta,
              ", cutpoint = ",
              if (lenc) "c(" else "",
              ToString(cutpoint),
              if (lenc) ")" else "",
              ", k = ",
              if (lenk) "c(" else "",
              ToString(kay),
              if (lenk) ")" else "",
              ")", sep = "")
      } else
        paste("2*log(sqrt(k) * sinh(qnorm(", theta,
              ")/(2*sqrt(k)) + ",
              "asinh(sqrt(cutpoint/k))))", sep = "")
      if (tag) 
        string <- paste("Negative binomial-ordinal link function:",
                        string)
      return(string)
  }


    thmat <- cbind(theta)
    kay <- rep(kay, len = ncol(thmat))  # Allow recycling for kay
    cutpoint <- rep(cutpoint, len = ncol(thmat))  # Allow recycling for cutpoint
    if (ncol(thmat) > 1) {
      answer <- thmat
      for (ii in 1:ncol(thmat))
          answer[, ii] <- Recall(theta = thmat[, ii],
                               cutpoint = cutpoint[ii],
                               k = kay[ii],
                               inverse = inverse, deriv = deriv)
      return(answer)
    }

    answer <-
    if (inverse) {
      if (deriv > 0) {
        1 / Recall(theta = theta,
                   cutpoint = cutpoint,
                   k = kay,
                   inverse = FALSE, deriv = deriv)
      } else {
        if (cutpoint == 0) {
          1.0 - (kay / (kay + care.exp(theta)))^kay
        } else {
            pnorm((asinh(exp(theta/2)/sqrt(kay)) -
                   asinh(sqrt(cutpoint/kay))) * 2 * sqrt(kay))
        }
      }
    } else {
      smallno <- 1 * .Machine$double.eps
      SMALLNO <- 1 * .Machine$double.xmin
      Theta <- theta
      Theta <- pmin(Theta, 1 - smallno)  # Since theta == 1 is a possibility
      Theta <- pmax(Theta, smallno)  # Since theta == 0 is a possibility
      if (cutpoint == 0) {
        switch(deriv+1, {
        temp <- (1 - Theta)^(-1/kay) - 1
        temp <- pmax(temp, SMALLNO)
        log(kay) + log(temp)},
        (kay / (1 - Theta)^(1/kay) - kay) * (1 - Theta)^(kay+1/kay),
        {  stop('cannot handle deriv = 2') },
        stop("argument 'deriv' unmatched"))
      } else {
        Ql <- qnorm(Theta)
        switch(deriv+1, {
              temp <- sqrt(kay) * sinh(Ql/(2*sqrt(kay)) +
                     asinh(sqrt(cutpoint/kay)))
              temp <- pmax(temp, SMALLNO)
              2 * log(temp)}, {
              arg1 <- (Ql/(2*sqrt(kay)) + asinh(sqrt(cutpoint/kay)))
              sqrt(kay) * tanh(arg1) * dnorm(Ql) },
              {  stop('cannot handle deriv = 2') },
              stop("argument 'deriv' unmatched"))
      }
    }
    if (!is.Numeric(answer)) stop("the answer contains some NAs")
    answer
}






 nbolf2 <- function(theta,
                    cutpoint = NULL,
                    k = NULL,
                    inverse = FALSE, deriv = 0,
                    short = TRUE, tag = FALSE) {

  kay <- k
  if (!is.Numeric(kay, positive = TRUE))
    stop("could not determine argument 'k' or ",
         "it is not positive-valued")
  if (!is.Numeric(cutpoint))
    stop("could not determine the cutpoint")
  if (any(cutpoint < 0) ||
      !is.Numeric(cutpoint, integer.valued = TRUE))
    warning("argument 'cutpoint' should ",
            "contain non-negative integer values")

  if (is.character(theta)) {
    string <- if (short) {
      lenc <- length(cutpoint) > 1
      lenk <- length(kay) > 1
      paste("nbolf2(", theta,
            ", earg = list(cutpoint = ",
            if (lenc) "c(" else "",
            ToString(cutpoint),
            if (lenc) ")" else "",
            ", k = ",
            if (lenk) "c(" else "",
            ToString(kay),
            if (lenk) ")" else "",
            "))", sep = "")
  } else {
    paste("3*log(<a complicated expression>)", sep = "")
  }
  if (tag) 
    string <- paste("Negative binomial-ordinal link function 2:",
                   string)
  return(string)
  }


    thmat <- cbind(theta)
    kay <- rep(kay, len = ncol(thmat))  # Allow recycling for kay
    if (ncol(thmat) > 1) {
        answer <- thmat
        for (ii in 1:ncol(thmat))
            answer[, ii] <- Recall(theta = thmat[, ii],
                                 cutpoint = cutpoint[ii],
                                 k = kay[ii],
                                 inverse = inverse, deriv = deriv)
        return(answer)
    }

    answer <-
    if (inverse) {
        if (deriv > 0) {
            1 / Recall(theta = theta,
                       cutpoint = cutpoint,
                       k = kay,
                       inverse = FALSE, deriv = deriv)
        } else {
            if (cutpoint == 0) {
                1.0 - (kay / (kay + care.exp(theta)))^kay
            } else {

            a1 <- -(9*cutpoint+8) / (cutpoint+1)
            a2 <- (9*kay-1) / (kay * (cutpoint+1)^(1/3))
            a3 <- 9 / (kay * (cutpoint+1)^(2/3))
            a4 <- 9 / (cutpoint+1)
            B <- exp(theta/3)
            mymat <- rbind(a1^2*a2^2 + 2*a1*a2^3*B + B^2*a2^4, 0,
                    -2*a1*a2*a3*B - 2*a2^2*a3*B^2 - a1^2*a3 - a2^2*a4, 0,
                    B^2 * a3^2 + a3 * a4)
            ans <- Re(t(apply(mymat, 2, polyroot)))
            theta2 <- invfun <- pnorm(-ans)  # pnorm(-x) = 1-pnorm(x)
            for (ii in 1:4) {
              theta2[, ii] <-
                Recall(theta = theta2[, ii],
                       cutpoint = cutpoint,
                       k = kay,
                       inverse = FALSE, deriv = deriv)
            }
            rankmat <- t(apply(abs(theta2 - theta), 1, rank))
            for (ii in 2:4) {
              if (any(index4 <- (rankmat[, ii] == 1))) {
                invfun[index4, 1] <- invfun[index4, ii]
              }
            }
            invfun[, 1]
            }
        }
    } else {
        smallno <- 1 * .Machine$double.eps
        SMALLNO <- 1 * .Machine$double.xmin
        Theta <- theta
        Theta <- pmin(Theta, 1 - smallno)  # Since theta == 1 is a possibility
        Theta <- pmax(Theta, smallno)  # Since theta == 0 is a possibility
        if (cutpoint == 0) {
            switch(deriv+1, {
            temp <- (1 - Theta)^(-1/kay) - 1
            temp <- pmax(temp, SMALLNO)
            log(kay) + log(temp)},
            (kay / (1 - Theta)^(1/kay) - kay) * (1 - Theta)^(kay+1/kay),
            {  stop("cannot handle 'deriv = 2'") },
            stop("argument 'deriv' unmatched"))
        } else {
            Ql <- qnorm(Theta)
            a1 <- -(9*cutpoint+8) / (cutpoint+1)
            a2 <- (9*kay-1) / (kay * (cutpoint+1)^(1/3))
            a3 <- 9 / (kay * (cutpoint+1)^(2/3))
            a4 <- 9 / (cutpoint+1)
            discrim <- a1^2 * a3 + a2^2 * a4 - Ql^2 * a3 * a4
            denomin <- Ql^2 * a3 - a2^2
            numerat <- (a1*a2 - Ql * sqrt(discrim))
            argmax1 <- numerat / denomin
            switch(deriv+1, {
                argmax2 <- (a1*a2 + Ql * sqrt(discrim)) / denomin
                temp <- ifelse(argmax1 > 0, argmax1, argmax2)
                temp <- pmax(temp, SMALLNO)
                3 * log(temp)}, {
                 BB <- (sqrt(discrim) - Ql^2 * a3 *
                       a4 / sqrt(discrim)) / dnorm(Ql)
                 CC <- 2 * Ql * a3 / dnorm(Ql)
                 dA.dtheta <- (-denomin * BB - numerat * CC) / denomin^2
                 argmax1 / (3 * dA.dtheta)
                },
                {  stop('cannot handle deriv = 2') },
                stop("argument 'deriv' unmatched"))
        }
    }
    if (!is.Numeric(answer)) stop("the answer contains some NAs")
    answer
}



 Cut <- function(y, breaks = c(-Inf, quantile(c(y), prob = (1:4)/4))) {
  y <- as.matrix(y)


  temp <- cut(y, breaks = breaks, labels = FALSE)
  temp <- c(temp)  # integer vector of integers
  if (any(is.na(temp)))
    stop("there are NAs")
  answer <- if (ncol(y) > 1) matrix(temp, nrow(y), ncol(y)) else temp
  if (ncol(y) > 1) {
    ynames <- dimnames(y)[[2]]
    if (!length(ynames))
      ynames <- paste("Y", 1:ncol(y), sep = "")
    xnames <- dimnames(y)[[1]]
    if (!length(xnames)) xnames = as.character(1:nrow(y))
    dimnames(answer) <- list(xnames, ynames)
  }
  attr(answer, "breaks") <- breaks
  answer
}


 checkCut <- function(y) {
  if (!is.Numeric(y, positive = TRUE, integer.valued = TRUE))
    stop("argument 'y' must contain positive integers only")
  uy <- unique(y)
  L <- max(uy)
  oklevels <- 1:L
  if (L == 1)
    stop("only one unique value")
  for (ii in oklevels) {
    if (all(ii != uy))
      stop("there is no ", ii, " value")
  }
  TRUE
}









 nbcanlink <- function(theta,
                       size = NULL,
                       wrt.eta = NULL,
                       bvalue = NULL,
                       inverse = FALSE, deriv = 0,
                       short = TRUE, tag = FALSE) {
  if (is.character(theta)) {
    string <- if (short)
      paste("nbcanlink(", theta, ")", sep = "") else
      paste("log(", theta, " / (", theta, " + size))", sep = "")
    if (tag)
      string <- paste("Nbcanlink:", string)
    return(string)
  }



  kmatrix <- size
  theta <- cbind(theta)
  kmatrix <- cbind(kmatrix)
  if (ncol(kmatrix) != ncol(theta))
    stop("arguments 'theta' and 'size' do not have ",
         "an equal number of cols")
  if (nrow(kmatrix) != nrow(theta))
    stop("arguments 'theta' and 'size' do not have ",
         "an equal number of rows")


  if (deriv > 0) {
    if (!(wrt.eta %in% 1:2))
      stop("argument 'wrt.eta' should be 1 or 2")
  }


  if (!inverse && length(bvalue))
    theta[theta <= 0.0] <- bvalue

  if (inverse) {
    if (deriv > 0) {
      1 / Recall(theta = theta,
                 size = size,
                 wrt.eta = wrt.eta,
                 bvalue = bvalue,
                 inverse = FALSE, deriv = deriv)
    } else {
       ans <- (kmatrix / expm1(-theta))
       if (is.matrix(ans))
         dimnames(ans) <- NULL else
         names(ans) <- NULL
       ans
    }
  } else {
    ans <-
    switch(deriv+1,
        (log(theta / (theta + kmatrix))),
       if (wrt.eta == 1) theta * (theta + kmatrix) / kmatrix else
       -(theta + kmatrix),
       if (wrt.eta == 1)
       -(theta * (theta + kmatrix))^2 / ((2 * theta + kmatrix) *
         kmatrix) else
       (theta + kmatrix)^2)
     if (is.matrix(ans))
       dimnames(ans) <- NULL else
       names(ans) <- NULL
     ans
  }
}





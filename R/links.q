# These functions are
# Copyright (C) 1998-2024 T.W. Yee, University of Auckland.
# All rights reserved.

















ToString <- function(x)
  paste(x, collapse = ", ")












 multilogitlink <-
  function(theta,
           refLevel = "(Last)",
           M = NULL,  # stop("argument 'M' not specified"),
           whitespace = FALSE,
           bvalue = NULL,
           inverse = FALSE, deriv = 0,
           all.derivs = FALSE,
           short = TRUE, tag = FALSE) {




  fillerChar <- ifelse(whitespace, " ", "")

  if (length(refLevel) != 1)
    stop("the length of argument 'refLevel' must be one")

  if (is.character(refLevel)) {
    if (refLevel != "(Last)")
      stop('if a character, refLevel must be "(Last)"')
    refLevel <- -1
  } else
  if (is.factor(refLevel)) {
    if (is.ordered(refLevel))
      warning("argument 'refLevel' is from an ordered factor")
    refLevel <- as.character(refLevel) == levels(refLevel)
    refLevel <- (seq_along(refLevel))[refLevel]
    if (!is.Numeric(refLevel, length.arg = 1,
                    integer.valued = TRUE, positive = TRUE))
      stop("could not coerce 'refLevel' into a single ",
           "positive integer")
  } else
  if (!is.Numeric(refLevel, length.arg = 1,
                  positive = TRUE, integer.valued = TRUE))
    stop("'refLevel' must be a single positive integer")




  if (is.character(theta)) {
    is.M <- is.finite(M) && is.numeric(M)
    string <- if (short) {
        paste0("multilogitlink(", theta, ")")
    } else {
        theta <- as.char.expression(theta)

         if (refLevel < 0) {
           ifelse(whitespace,
             paste0("log(", theta, "[,j] / ",
                   theta, "[,",
                   ifelse(is.M, M+1, "M+1"),
                   "]), j = 1:",
                   ifelse(is.M, M, "M")),
             paste0("log(", theta, "[,j]/",
                   theta, "[,",
                   ifelse(is.M, M+1, "M+1"),
                   "]), j=1:",
                   ifelse(is.M, M, "M")))
         } else {
             if (refLevel == 1) {
               paste0("log(", theta, "[,", "j]",
                   fillerChar, "/", fillerChar,
                     "", theta, "[,", refLevel, "]), j",
                     fillerChar, "=", fillerChar, "2:",
                   ifelse(is.M, (M+1), "(M+1)"))
             } else {
               paste0("log(", theta, "[,", "j]", fillerChar, "/",
                     "", theta, "[,", refLevel, "]), j",
                     fillerChar, "=", fillerChar,
                     "c(1:", refLevel-1, ",",
                     fillerChar,
                     refLevel+1, ":",
                     ifelse(is.M, (M+1), "(M+1)"), ")")
             }
         }
    }
    if (tag)
      string <- paste("Multinomial logit link:", string)
    return(string)
  }



  M.orig <- M
  M <- NCOL(theta) - !(inverse && deriv == 0)
  if (M < 1)
    ifelse(inverse,
       stop("argument 'eta' should have at least one column"),
       stop("argument 'theta' should have at least two columns"))
  if (is.numeric(M.orig) && M != M.orig) {
    warning("argument 'M' does not seem right but using it")
    M <- M.orig
  }
  if (is.numeric(refLevel) && refLevel > M + 1)
    stop("bad input for argument 'refLevel'")



  if (!inverse && length(bvalue))
    theta[theta <= 0.0] <- bvalue
  if (!inverse && length(bvalue))
    theta[theta >= 1.0] <- 1 - bvalue


  foo <- function(eta, refLevel = -1, M) {
    use.refLevel <- if ( refLevel < 0)
                    M+1 else refLevel  # unneeded
    phat <- if ((refLevel < 0) || (refLevel == M+1)) {
      care.exp2(cbind(eta, 0.0))
    } else if ( refLevel == 1) {
      care.exp2(cbind(0.0, eta))
    } else {
      etamat <- cbind(eta[, 1:( refLevel - 1), drop = FALSE],
                      0.0,
                      eta[, ( refLevel ):M, drop = FALSE])
      care.exp2(etamat)
    }

    rSp <- rowSums(phat)
    ans <- phat / rSp
    colnames(ans) <- NULL  # Safest for now
    ans
  }  # foo



  use.refLevel <- if ( refLevel < 0) M+1 else refLevel

  if (inverse) {
    use.refLevel <- if (refLevel < 0) ncol(theta) else refLevel
    switch(as.character(deriv),
      "0" = {
              foo(theta, refLevel,  # refLevel, not use.refLevel
                  M = M)
            },
      "1" = if (all.derivs) {
              index <- iam(NA, NA, M = M, both = TRUE, diag = TRUE)
              theta <- theta[, -use.refLevel, drop = FALSE]  # n x M
              wz <- -theta[, index$row, drop = FALSE] *
                     theta[, index$col, drop = FALSE]
              wz[, 1:M] <- wz[, 1:M] + theta
              wz
            } else {
              theta[, -use.refLevel,  drop = FALSE] *
              theta[,  use.refLevel] / (
              theta[, -use.refLevel,  drop = FALSE] +
              theta[,  use.refLevel])
            },
      "2" = (theta*(1-theta)*(1-2*theta))[, -use.refLevel,
                                          drop = FALSE],
      "3" = {
              temp1 <- theta * (1 - theta)
             (temp1 * (1 - 6 * temp1))[, -use.refLevel,
                                       drop = FALSE]
            },
      stop("argument 'deriv' unmatched"))
  } else {  # Not inverse below here ,,,,,,,,,,,,,,,,,,,

    switch(as.character(deriv),
           "0" = {
      ans <- if (refLevel < 0) {
        log(theta[, -ncol(theta), drop = FALSE] / (
            theta[, ncol(theta)]))
      } else {
        use.refLevel <- if (refLevel < 0)
                          ncol(theta) else refLevel
        log(theta[, -( use.refLevel ), drop = FALSE] / (
            theta[, use.refLevel ]))
      }

      colnames(ans) <- NULL  # Safest for now
      ans
      },
      "1" = care.exp(-log(theta) - log1p(-theta)),
      "2" = (2 * theta - 1) / care.exp(2*log(theta) +
                                       2*log1p(-theta)),
      "3" = {
        temp1 <- care.exp(log(theta) + log1p(-theta))
        2 * (1 - 3 * temp1) / temp1^3
      },
      stop("argument 'deriv' unmatched"))
  }
}  # multilogitlink









 as.char.expression <- function(x) {
  answer <- x
  for (i in length(x)) {
    charvec <- substring(x[i], 1:nchar(x[i]), 1:nchar(x[i]))
    if (!all(is.element(charvec,
                        c(letters,
                          LETTERS,
                          as.character(0:9), ".", "_"))))
      answer[i] <- paste0("(", x[i], ")")
  }
  answer
}



if (FALSE) {
  as.char.expression("a")
  as.char.expression("a+b")
  as.char.expression(c("a", "a+b"))
}






 TypicalVGAMfamilyFunction <-
  function(lsigma = "loglink",
           isigma = NULL,
           zero = NULL,
           gsigma = exp(-5:5),
           eq.mean = FALSE,
           parallel = TRUE,
           imethod = 1,
           vfl = FALSE, Form2 = NULL,
           type.fitted = c("mean", "quantiles", "Qlink",
                           "pobs0", "pstr0", "onempstr0"),
           percentiles = c(25, 50, 75),
           probs.x = c(0.15, 0.85),
           probs.y = c(0.25, 0.50, 0.75),
           multiple.responses = FALSE, earg.link = FALSE,
           ishrinkage = 0.95, nointercept = NULL,
           whitespace = FALSE, bred = FALSE, lss = TRUE,
           oim = FALSE, nsimEIM = 100, byrow.arg = FALSE,
           link.list = list("(Default)" = "identitylink",
                            x2          = "loglink",
                            x3          = "logofflink",
                            x4          = "multilogitlink",
                            x5          = "multilogitlink"),
           earg.list = list("(Default)" = list(),
                            x2          = list(),
                            x3          = list(offset = -1),
                            x4          = list(),
                            x5          = list()),
           Thresh = NULL, nrfs = 1) {
  NULL
}


TypicalVGAMlink <-
  function(theta,
           someParameter = 0,
      bvalue = NULL,  # .Machine$double.xmin
           inverse = FALSE, deriv = 0,
           short = TRUE, tag = FALSE) {
  NULL
}






care.exp <-
  function(x,
           thresh = -log( sqrt( .Machine$double.xmin ) )) {


  x[x >   thresh]  <-  thresh
  x[x < (-thresh)] <- -thresh
  exp(x)
}



care.exp2 <- function(x) {
  if (NCOL(x) == 1)
    x <- cbind(x)
  exp(x - x[cbind(1:NROW(x), max.col(x))])
}










 loglink <-
  function(theta,
           bvalue = NULL,  # .Machine$double.xmin
           inverse = FALSE, deriv = 0,
           short = TRUE, tag = FALSE) {


  if (is.character(theta)) {
    string <- if (short)
        paste0("loglink(",  theta, ")") else
        paste0("log(",  theta, ")")
    if (tag)
      string <- paste("Log:", string)
    return(string)
  }

  if (!inverse && length(bvalue))
    theta[theta <= 0.0] <- bvalue

  if (inverse) {
    switch(as.character(deriv),
           "0" = exp(theta),
           "1" = theta,
           "2" = theta,
           "3" = theta,
           "4" = theta,
           "5" = theta,
           "6" = theta,
           "7" = theta,
           "8" = theta,
           stop("argument 'deriv' unmatched"))
  } else {
    switch(as.character(deriv),
       "0" = log(theta),
       "1" =     1 / theta,
       "2" =    -1 / theta^2,
       "3" =     2 / theta^3,
       "4" =    -6 / theta^4,
       "5" =    24 / theta^5,
       "6" =  -120 / theta^6,
       "7" =   720 / theta^7,
       "8" = -5040 / theta^8,
       stop("argument 'deriv' unmatched"))
  }
}  # loglink






 logneglink <-
  function(theta,
           bvalue = NULL,  # .Machine$double.xmin
           inverse = FALSE, deriv = 0,
           short = TRUE, tag = FALSE) {


  if (is.character(theta)) {
    string <- if (short)
        paste0("logneglink(",  theta, ")") else
        paste0("log(-(",  theta, "))")
    if (tag)
      string <- paste("Log negative:", string)
    return(string)
  }

  if (!inverse && length(bvalue))
    theta[theta <= 0.0] <- bvalue

  if (inverse) {
    switch(as.character(deriv),
           "0" = -exp(theta),
           "1" = theta,
           "2" = theta,
           "3" = theta,
           "4" = theta,
           stop("argument 'deriv' unmatched"))
  } else {
    switch(as.character(deriv),
           "0" = log(-theta),
           "1" =  1 / theta,
           "2" = -1 / theta^2,
           "3" =  2 / theta^3,
           "4" = -6 / theta^4,
           stop("argument 'deriv' unmatched"))
  }
}  # logneglink






 logofflink <-
  function(theta,
           offset = 0,
           inverse = FALSE, deriv = 0,
           short = TRUE, tag = FALSE) {

  if (!is.Numeric(offset))
    stop("bad input for argument 'offset'")

  if (is.character(theta)) {
    string <- if (short)
      paste0("logofflink(", theta,
            ", offset = ", as.character(offset),
            ")") else
      paste0("log(",
            as.character(offset),
            "+",
            as.char.expression(theta),
            ")")
    if (tag)
      string <- paste("Log with offset:", string)
    return(string)
  }

  if (inverse) {
    switch(as.character(deriv),
           "0" = exp(theta) - offset,
           "1" = theta + offset,
           "2" = theta + offset,
           "3" = theta + offset,
           "4" = theta + offset,
           "5" = theta + offset,
           stop("argument 'deriv' unmatched"))
  } else {
    switch(as.character(deriv),
           "0" = log(theta + offset),
           "1" =  1 / (theta + offset),
           "2" = -1 / (theta + offset)^2,
           "3" =  2 / (theta + offset)^3,
           "4" = -6 / (theta + offset)^4,
           "5" = 24 / (theta + offset)^4,
           stop("argument 'deriv' unmatched"))
  }
}  # logofflink








 log1plink <-
  function(theta,
           offset = 0,  # To be left alone.
           inverse = FALSE, deriv = 0,
           short = TRUE, tag = FALSE) {

  if (!is.Numeric(offset))
    stop("bad input for argument 'offset'")
  if (!all(offset == 0))
    stop("'offset' should be left alone because ",
         "it is implicitly unity. ",
         "Use logofflink() instead?")

  if (is.character(theta)) {
    string <- if (short)
      paste0("log1plink(", theta, ")") else
      paste0("log(1+", as.char.expression(theta), ")")
    if (tag)
      string <- paste("Log with unit offset:", string)
    return(string)
  }

  if (inverse) {
    switch(as.character(deriv),
           "0" = expm1(theta),
           "1" = theta + 1,
           "2" = theta + 1,
           "3" = theta + 1,
           "4" = theta + 1,
           "5" = theta + 1,
           stop("argument 'deriv' unmatched"))
  } else {
    switch(as.character(deriv),
           "0" = log1p(theta),
           "1" =  1 / (theta + 1),
           "2" = -1 / (theta + 1)^2,
           "3" =  2 / (theta + 1)^3,
           "4" = -6 / (theta + 1)^4,
           "5" = 24 / (theta + 1)^5,
           stop("argument 'deriv' unmatched"))
  }
}  # log1plink






 identitylink <-
  function(theta,
           inverse = FALSE, deriv = 0,
           short = TRUE, tag = FALSE) {


  if (is.character(theta)) {
    string <- theta
    if (tag)
      string <- paste("Identity:", string)
    return(string)
  }

  switch(as.character(deriv),
         "0" = theta,
         "1" = theta * 0 + 1,
         "2" = theta * 0,  # zz Does not handle Inf and -Inf
         "3" = theta * 0,  # zz Does not handle Inf and -Inf
         "4" = theta * 0,  # zz Does not handle Inf and -Inf
         "5" = theta * 0,  # zz Does not handle Inf and -Inf
         stop("argument 'deriv' unmatched"))
}  # identitylink






 negidentitylink <-
  function(theta,
           inverse = FALSE, deriv = 0,
           short = TRUE, tag = FALSE) {
  if (is.character(theta)) {
    theta <- as.char.expression(theta)
    string <- paste0("-", theta)
    if (tag)
      string <- paste("Negative-identity:", string)
    return(string)
  }

  switch(as.character(deriv),
         "0" = -theta,
         "1" = theta * 0 - 1,
         "2" = theta * 0,  # zz Does not handle Inf and -Inf
         "3" = theta * 0,  # zz Does not handle Inf and -Inf
         "4" = theta * 0,  # zz Does not handle Inf and -Inf
         "5" = theta * 0,  # zz Does not handle Inf and -Inf
         stop("argument 'deriv' unmatched"))
}  # negidentitylink






 logitlink <-
  function(theta,
           bvalue = NULL,
           inverse = FALSE, deriv = 0,
           short = TRUE, tag = FALSE) {
  if (is.character(theta)) {
    string <- if (short)
        paste0("logitlink(",  # "logit(",
               theta,
              ")") else
        paste0("log(",
              as.char.expression(theta),
              "/(1-",
              as.char.expression(theta),
              "))")
    if (tag)
      string <- paste("Logit:", string)
    return(string)
  }

  if (!inverse && length(bvalue)) {
    theta[theta <= 0.0] <- bvalue
    theta[theta >= 1.0] <- 1.0 - bvalue
  }
  if (inverse) {
    switch(as.character(deriv),
   "0" = plogis(theta),
   "1" =    1 / Recall(theta = theta,
                bvalue = bvalue,
                inverse = FALSE,
                deriv = deriv),
   "2" = theta * (1 - theta) * (1 - 2 * theta),
   "3" = (1 - 6 * theta * (1 - theta)) *
          theta * (1 - theta),
        "4" = {
  iD1 <- Recall(theta, deriv = 1, inverse = TRUE)
  iD2 <- Recall(theta, deriv = 2, inverse = TRUE)
  iD3 <- Recall(theta, deriv = 3, inverse = TRUE)
  DD1 <- Recall(theta, deriv = 1, inverse = FALSE)
  DD2 <- Recall(theta, deriv = 2, inverse = FALSE)
  DD3 <- Recall(theta, deriv = 3, inverse = FALSE)
  DD4 <- Recall(theta, deriv = 4, inverse = FALSE)
  (iD1^3) * (15 * iD1 * iD2 * (DD2^2) +
              6 * (iD1^3) * DD2 * DD3 -
              4 * iD2 * DD3 - (iD1^2) * DD4)
        },
           stop("argument 'deriv' unmatched"))
  } else {
    switch(as.character(deriv),
     "0" = qlogis(theta),
     "1" = 1 / (theta * (1 - theta)),
  "2" = (2 * theta - 1) / (theta * (1 - theta))^2,
  "3" = 2 * (1 - 3 * theta *
             (1 - theta)) / (theta * (1 - theta))^3,
     "4" = -6 * (1 - 2 * theta) *
          (1 - 2 * theta *
          (1 - theta)) / (theta * (1 - theta))^4,
     stop("argument 'deriv' unmatched"))
  }
}  # logitlink






 logloglink <-
  function(theta,
           bvalue = NULL,  # .Machine$double.eps
           inverse = FALSE, deriv = 0,
           short = TRUE, tag = FALSE) {
  if (is.character(theta)) {
    string <- if (short)
        paste0("logloglink(",  theta, ")") else
        paste0("log(log(", theta, "))")
    if (tag)
      string <- paste("Log-Log:", string)
    return(string)
  }

  if (!inverse && length(bvalue))
    theta[theta <= 1.0] <- bvalue

  if (inverse) {
    switch(as.character(deriv),
           "0" = exp(exp(theta)),
           "1" = (theta * log(theta)),
           "2" = { junk <- log(theta)
                   theta  * junk * (1 + junk) },
           "3" = { Junk3 <- theta * log(theta)
                   Junk3 * ((1 + log(theta))^2 + Junk3 / theta)
      },
           stop("argument 'deriv' unmatched"))
  } else {
    switch(as.character(deriv),
           "0" = log(log(theta)),
           "1" = 1 / (theta * log(theta)),
           "2" = { junk <- log(theta)
                   -(1 + junk) / (theta * junk)^2
           },
           "3" = { Junk3 <- theta * log(theta)
           (2 * (1 + log(theta))^2 / Junk3 - 1 / theta) / Junk3^2
                 },
           stop("argument 'deriv' unmatched"))
  }
}  # logloglink






 loglogloglink <-
  function(theta,
           bvalue = NULL,  # .Machine$double.eps
           inverse = FALSE, deriv = 0,
           short = TRUE, tag = FALSE) {
  if (is.character(theta)) {
    string <- if (short)
        paste0("loglogloglink(",  theta, ")") else
        paste0("log(log(log(", theta, ")))")
    if (tag)
      string <- paste("Log-Log-Log:", string)
    return(string)
  }

  if (!inverse && length(bvalue))
    theta[theta <= exp(1.0)] <- bvalue

  if (inverse) {
    switch(as.character(deriv),
           "0" = exp(exp(exp(theta))),
           "1" = theta * log(theta) * log(log(theta)),
           "2" = { junk <- log(theta)
                   logjunk <- log(junk)
                   theta * junk * logjunk * (1 + logjunk * (1 + junk))
                 },
           "3" = { junk <- log(theta)
                   logjunk <- log(junk)
                   theta * junk^2 * logjunk^3  * (
                   3 + junk + 1 / junk + 3 / logjunk +
                   3 / (junk * logjunk) +
                   1 / (junk * logjunk^2))
      },
           stop("argument 'deriv' unmatched"))
  } else {
    switch(as.character(deriv),
           "0" = log(log(log(theta))),
           "1" = 1 / (theta * log(theta) * log(log(theta))),
           "2" = { junk <- log(theta)
                   logjunk <- log(junk)
                   (-1 / (theta^2 * junk * logjunk)) *
                   (1 + (1 / junk) * (1 + 1 / logjunk))
           },
           "3" = { junk <- log(theta)
                   logjunk <- log(junk)
                   (3 + 2 * junk + 2 / junk +
                   3 / logjunk +
                   3 / (junk * logjunk) +
                   2 / (junk * logjunk^2)) / (
                   theta^3 * junk^2 * logjunk)
                 },
           stop("argument 'deriv' unmatched"))
  }
}  # loglogloglink









 clogloglink <-
  function(theta,
           bvalue = NULL,  # .Machine$double.eps
           inverse = FALSE, deriv = 0,
           short = TRUE, tag = FALSE) {
  if (is.character(theta)) {
    string <- if (short)
        paste0("clogloglink(", theta, ")") else
        paste0("log(-log(1-",
              as.char.expression(theta), "))")
    if (tag)
      string <- paste("Complementary log-log:", string)
    return(string)
  }

  if (!inverse && length(bvalue)) {
    theta[theta <= 0.0] <- bvalue
    theta[theta >= 1.0] <- 1.0 - bvalue
  }

  if (inverse) {
    switch(as.character(deriv),
    "0" = { -expm1(-exp(theta)) },
    "1" = {
        ans5 <- (-(1 - theta) * log1p(-theta))
        ans5[1 - theta == 0] <- 0  # 20210522; limit
        ans5
    },
    "2" = {
        junk <- log1p(-theta)
        ans6 <- -(1 - theta) * (1 + junk) * junk
        ans6[1 - theta == 0] <- 0  # 20210522; limit
        ans6
     },
     "3" = {
         junk <- log1p(-theta)
         Junk2 <- (1 - theta) * junk
         ans7 <- -Junk2 * (Junk2 / (1 - theta) +
                           (1 + junk)^2)
         ans7[1 - theta == 0] <- 0  # 20210524; limit
         ans7
      },
      stop("argument 'deriv' unmatched"))
  } else {
    switch(as.character(deriv),
    "0" = log(-log1p(-theta)),
    "1" = {
        ans5 <- -1 / ((1 - theta) * log1p(-theta))
        ans5[1 - theta == 0] <- Inf  # 20210522; limit
        ans5
    },
    "2" = {
        junk <- log1p(-theta)
        ans6 <- -(1 + junk) / ((1-theta) * junk)^2
        ans6[1 - theta == 0] <- Inf  # 20210522; limit
        ans6
     },
    "3" = {
        junk <- log1p(-theta)
        Junk3 <- (1 - theta) * junk
        ans7 <- (1 / (1 - theta) -
                2 * (1 + junk)^2 / Junk3) / Junk3^2
        ans7[1 - theta == 0] <- Inf  # 20210524; limit
        ans7
     },
     stop("argument 'deriv' unmatched"))
  }
}  # clogloglink








 cloglink <-
  function(theta,
           bvalue = NULL,  # .Machine$double.eps
           inverse = FALSE, deriv = 0,
           short = TRUE, tag = FALSE) {
  if (is.character(theta)) {
    string <- if (short)
        paste0("cloglink(", theta, ")") else
        paste0("-log(1-",
              as.char.expression(theta), ")")
    if (tag)
      string <- paste("Complementary log:", string)
    return(string)
  }

  if (!inverse && length(bvalue)) {
    theta[theta <= 0.0] <- bvalue
    theta[theta >= 1.0] <- 1.0 - bvalue
  }

  if (inverse) {
    switch(as.character(deriv),
    "0" = -expm1(-theta),
    "1" = 1 - theta,
    "2" = theta - 1,
    "3" = 1 - theta,
    "4" = theta - 1,
    "5" = 1 - theta,
    stop("argument 'deriv' unmatched"))
  } else {
    switch(as.character(deriv),
    "0" = -log1p(-theta),
    "1" = { ans5 <- 1 / (1 - theta)
            ans5[1 - theta == 0] <- Inf  # limit
            ans5 },
    "2" = { ans6 <- 1 / (1 - theta)^2
            ans6[1 - theta == 0] <- Inf  # limit
            ans6 },
    "3" = { ans7 <- 2 / (1 - theta)^3
            ans7[1 - theta == 0] <- Inf  # limit
            ans7 },
    "4" = { ans9 <- 6 / (1 - theta)^4
            ans9[1 - theta == 0] <- Inf  # limit
            ans9 },
    "5" = { ans5 <- 24 / (1 - theta)^5
            ans5[1 - theta == 0] <- Inf  # limit
            ans5 },
     stop("argument 'deriv' unmatched"))
  }
}  # cloglink






 probitlink <-
  function(theta,
           bvalue = NULL,  # .Machine$double.eps
           inverse = FALSE, deriv = 0,
           short = TRUE, tag = FALSE) {
  if (is.character(theta)) {
    string <- if (short)
        paste0("probitlink(", theta, ")") else
        paste0("qnorm(",  theta, ")")
    if (tag)
      string <- paste("Probit:", string)
    return(string)
  }

  if (!inverse && length(bvalue)) {
    theta[theta <= 0.0] <- bvalue
    theta[theta >= 1.0] <- 1 - bvalue
  }

  if (inverse) {
    switch(as.character(deriv),
           "0" = {
                   ans <- pnorm(theta)
                   if (is.matrix(theta))
                     dim(ans) <- dim(theta)
                   ans
                  },
           "1" = {  # 1st deriv
      1 / Recall(theta = theta,
                 bvalue = bvalue,
                 inverse = FALSE, deriv = deriv)
         },
           "2" = {  # 2nd deriv
        Junk2 <- qnorm(theta)
        ans6 <- -Junk2 * dnorm(Junk2)
        ans6[1-theta == 0] <- 0  # 20210525; limit
        if (is.vector(theta)) ans6 else
        if (is.matrix(theta)) {
          dim(ans6) <- dim(theta)
          ans6
        } else {
          warning("can only handle vectors and ",
             "matrices; converting to vector")
          ans6
        }
      },
           "3" = {
             Junk3 <- qnorm(theta)
             junk <- dnorm(Junk3)
             ans7 <- junk * (Junk3^2 - 1)
        ans7[1-theta == 0] <- 0  # 20210525 limit
             ans7
           },
           stop("argument 'deriv' unmatched"))
  } else {
    switch(as.character(deriv),
           "0" = {
        ans <- qnorm(theta)
        if (is.matrix(theta))
          dim(ans) <- dim(theta)
        ans
     },
           "1" = {  # 1st deriv
       if (is.matrix(theta)) {
         ans <- 1 / dnorm(qnorm(theta))
         dim(ans) <- dim(theta)
         ans
       } else {
         1 / dnorm(qnorm(as.vector(theta)))
       }
      },
           "2" = {  # 2nd deriv
        Junk2 <- qnorm(theta)

        ans <- Junk2 / (dnorm(Junk2))^2

        if (is.vector(theta)) ans else
        if (is.matrix(theta)) {
          dim(ans) <- dim(theta)
          ans
        } else {
          warning("can only handle vectors and matrices;",
                  " converting to vector")
          ans
        }
      },
           "3" = {
             Junk3 <- qnorm(theta)
             junk <- dnorm(Junk3)
             (1 + 2 * Junk3^2) / junk^3
           },
           stop("argument 'deriv' unmatched"))
  }
}  # probitlink







 explink <-
  function(theta,
           bvalue = NULL,  # .Machine$double.eps
           inverse = FALSE, deriv = 0,
           short = TRUE, tag = FALSE) {
  if (is.character(theta)) {
    string <- if (short)
        paste0("explink(", theta, ")") else
        paste0("exp(", theta, ")")
    if (tag)
      string <- paste("Exp:", string)
    return(string)
  }

  if (!inverse && length(bvalue))
    theta[theta <= 0.0] <- bvalue
  if (inverse) {
    switch(as.character(deriv),
           "0" = log(theta),
           "1" =      exp(    -theta),
           "2" = -    exp(-2 * theta),  # 20170610 Fixes up a bug
           "3" =  2 * exp(-3 * theta),
           "4" = -6 * exp(-4 * theta),
           stop("argument 'deriv' unmatched"))
  } else {
    switch(as.character(deriv),
           "0" = exp(theta),
           "1" = exp(theta),
           "2" = exp(theta), 
           "3" = exp(theta),
           "4" = exp(theta),
           stop("argument 'deriv' unmatched"))
  }
}  # explink






 reciprocallink <-
  function(theta,
           bvalue = NULL,  # .Machine$double.eps
           inverse = FALSE, deriv = 0,
           short = TRUE, tag = FALSE) {
  if (is.character(theta)) {
    theta <- as.char.expression(theta)
    string <- paste0("1/", theta)
    if (tag)
      string <- paste("Reciprocal:", string)
    return(string)
  }

  if (!inverse && length(bvalue))
    theta[theta == 0.0] <- bvalue

  if (inverse) {
    switch(as.character(deriv),
           "0" = 1 / theta,
           "1" = -   theta^2,
           "2" =  2 * theta^3,
           "3" = -6 * theta^4,
           "4" = 24 * theta^5,
           stop("argument 'deriv' unmatched"))
  } else {
    switch(as.character(deriv),
           "0" =  1 / theta,
           "1" = -1 / theta^2,
           "2" =  2 / theta^3,
           "3" = -6 / theta^4,
           "4" = 24 / theta^5,
           stop("argument 'deriv' unmatched"))
  }
}  # reciprocallink






 negloglink <-
  function(theta,
           bvalue = NULL,  # .Machine$double.eps
           inverse = FALSE, deriv = 0,
           short = TRUE, tag = FALSE) {
  if (is.character(theta)) {
      string <- if (short)
          paste0("negloglink(", theta, ")") else
          paste0("-log(",  theta, ")")
      if (tag)
        string <- paste("Negative log:", string)
      return(string)
  }


  if (!inverse && length(bvalue))
    theta[theta <= 0.0] <- bvalue
  if (inverse) {
    switch(as.character(deriv),
           "0" = exp(-theta),
           "1" = -theta,
           "2" =  theta,
           "3" = -theta,
           "4" =  theta,
           stop("argument 'deriv' unmatched"))
  } else {
    switch(as.character(deriv),
           "0" = -log(theta),
           "1" = -1/theta,
           "2" =  1/theta^2,
           "3" = -2/theta^3,
           "4" =  6/theta^4,
           stop("argument 'deriv' unmatched"))
  }
}  # negloglink






 negreciprocallink <-
  function(theta,
           bvalue = NULL,  # .Machine$double.eps
           inverse = FALSE,
           deriv = 0, short = TRUE, tag = FALSE) {
  if (is.character(theta)) {
    theta <- as.char.expression(theta)
    string <- paste0("-1/", theta)
    if (tag)
      string <- paste("Negative reciprocal:", string)
    return(string)
  }


  if (!inverse && length(bvalue))
    theta[theta == 0.0] <- bvalue

  if (inverse) {
    switch(as.character(deriv),
           "0" = -1 / theta,
           "1" =      theta^2,
           "2" =  2 * theta^3,
           "3" =  6 * theta^4,
           "4" = 24 * theta^5,
           stop("argument 'deriv' unmatched"))
  } else {
    switch(as.character(deriv),
           "0" =  -1 / theta,
           "1" =   1 / theta^2,
           "2" =  -2 / theta^3,
           "3" =   6 / theta^4,
           "4" = -24 / theta^5,
           stop("argument 'deriv' unmatched"))
  }
}  # negreciprocallink






  igcanlink <-
  function(theta,
           bvalue = NULL,  # .Machine$double.eps
           inverse = FALSE, deriv = 0,
           short = TRUE, tag = FALSE) {

  if (is.character(theta)) {
    theta <- as.char.expression(theta)
    string <- paste0("-1/", theta)
    if (tag)
      string <- paste("Negative inverse:", string)
    return(string)
  }

  if (inverse) {
    switch(as.character(deriv),
           "0" =  1 / sqrt(-2*theta),
           "1" =      theta^3,
           "2" =  3 * theta^5, 
           "3" = 15 * theta^7, 
           stop("argument 'deriv' unmatched"))
  } else {
    switch(as.character(deriv),
           "0" =  -1 / (2 * theta^2),
           "1" =   1 / theta^3,
           "2" =  -3 / theta^4,
           "3" =  12 / theta^5,
           "4" = -60 / theta^6,
           stop("argument 'deriv' unmatched"))
  }
}  # igcanlink






 rhobitlink <-
  function(theta,
           bminvalue = NULL,
           bmaxvalue = NULL,
           inverse = FALSE, deriv = 0,
           short = TRUE, tag = FALSE) {
  if (is.character(theta)) {
    string <- if (short)
        paste0("rhobitlink(", theta, ")") else
        paste0("log((1+",
              as.char.expression(theta),
              ")/(1-",
              as.char.expression(theta), "))")
    if (tag)
      string <- paste("Rhobit:", string)
    return(string)
  }

  if (!inverse) {
    if (length(bminvalue))
      theta[theta <= -1.0] <- bminvalue
    if (length(bmaxvalue))
      theta[theta >=  1.0] <- bmaxvalue
  }

  if (inverse) {
    switch(as.character(deriv),
           "0" = {  # junk <- exp(theta)
             expm1(theta) / (exp(theta) + 1.0) },
           "1" = (1 - theta^2) / 2,
       "2" = (-theta / 2) * (1 - theta^2),
       "3" = (3 * theta^2 - 1) * (1 - theta^2) / 4,
             stop("argument 'deriv' unmatched"))
  } else {
    switch(as.character(deriv),
           "0" = {
                 log1p(theta) - log1p(-theta) },
           "1" = 2 / (1 - theta^2),
           "2" = (4*theta) / (1 - theta^2)^2,
           "3" = 4 * (1 + 3 * theta^2) / (1 - theta^2)^3,
           stop("argument 'deriv' unmatched"))
  }
}  # rhobitlink






 fisherzlink <-
  function(theta,
           bminvalue = NULL,
           bmaxvalue = NULL,
           inverse = FALSE, deriv = 0,
           short = TRUE, tag = FALSE) {
  if (is.character(theta)) {
    string <- if (short)
        paste0("fisherzlink(", theta, ")") else
        paste0("(1/2) * log((1+",
              as.char.expression(theta),
              ")/(1-",
              as.char.expression(theta), "))")
    if (tag)
      string <- paste("Fisher's Z transformation:", string)
    return(string)
  }

  if (!inverse) {
     if (length(bminvalue)) theta[theta <= -1.0] <- bminvalue
     if (length(bmaxvalue)) theta[theta >=  1.0] <- bmaxvalue
  }

  if (inverse) {
    switch(as.character(deriv),
           "0" = tanh(theta),
           "1" = 1 - theta^2,
           "2" = 2 * (-theta) * (1 - theta^2),
           "3" = (3 * theta^2 - 1) * (1 - theta^2) * 2,
           stop("argument 'deriv' unmatched"))
  } else {
    switch(as.character(deriv),
           "0" = atanh(theta),
           "1" = 1 / (1 - theta^2),
           "2" = (2*theta) / (1 - theta^2)^2,
           "3" = 2 * (1 + 3 * theta^2) / (1 - theta^2)^3,
           stop("argument 'deriv' unmatched"))
    }
}  # fisherzlink







 foldsqrtlink <-
  function(theta,  #  = NA  , = NULL,
           min = 0, max = 1, mux = sqrt(2),
           inverse = FALSE, deriv = 0,
           short = TRUE, tag = FALSE) {
  if (!is.Numeric(min, length.arg = 1))
    stop("bad input for 'min' component")
  if (!is.Numeric(max, length.arg = 1))
    stop("bad input for 'max' component")
  if (!is.Numeric(mux, length.arg = 1,
                  positive = TRUE))
    stop("bad input for 'mux' component")
  if (min >= max)
    stop("'min' >= 'max' is not allowed")

  if (is.character(theta)) {
    string <- if (short)
      paste0("foldsqrtlink(", theta, ")") else {
    theta <- as.char.expression(theta)
      if (abs(mux-sqrt(2)) < 1.0e-10)
        paste0("sqrt(2*", theta,
               ") - sqrt(2*(1-", theta, "))") else
        paste0(as.character(mux),
               " * (sqrt(", theta, "-", min,
               ") - sqrt(",
            max, "-", theta, "))")
    }
    if (tag)
        string <- paste("Folded square root:",
                        string)
    return(string)
  }

  if (inverse) {
    switch(as.character(deriv),
           "0" = {
      mid <- (min + max) / 2
      boundary <- mux * sqrt(max - min)
      temp <- pmax(0, (theta/mux)^2 *
                  (2*(max-min) - (theta/mux)^2))
      ans <- theta
      if (any(ind5 <- theta <  0))
        ans[ind5] <- mid - 0.5 * sqrt(temp[ind5])
      if (any(ind5 <- theta >= 0))
        ans[ind5] <- mid + 0.5 * sqrt(temp[ind5])
      ans[theta < -boundary] <- NA
      ans[theta >  boundary] <- NA
      ans
        },
      "1" = (2 / mux ) / (1/sqrt(theta-min) +
                          1/sqrt(max-theta)),
      "2" = stop("use the chain rule formula",
                 " to obtain this"),
      "3" = {  # 3rd deriv
        stop("3rd deriv not yet implemented")
      },
      stop("argument 'deriv' unmatched"))
  } else {
    switch(as.character(deriv),
           "0" = mux * (sqrt(theta-min) -
                        sqrt(max-theta)),
           "1" = (1/sqrt(theta-min) +
                  1/sqrt(max-theta)) * mux / 2,
        "2" = -(mux / 4) * ((theta-min)^(-3/2) -
                            (max-theta)^(-3/2)),
        "3" = { #3rd deriv
      stop("3rd deriv not yet implemented")
           },
           stop("argument 'deriv' unmatched"))
  }
}  # foldsqrtlink






 powerlink <-
  function(theta,
           power = 1,
           inverse = FALSE, deriv = 0,
           short = TRUE, tag = FALSE) {
    exponent <- power
    if (exponent == 0)
      stop("use the 'loge' link")

  if (is.character(theta)) {
    string <- if (short)
        paste0("powerlink(", theta, ", power = ",
               as.character(exponent), ")") else
        paste0(as.char.expression(theta),
               "^(", as.character(exponent), ")")
    if (tag)
      string <- paste("Power link:", string)
    return(string)
  }

  if (inverse) {
    switch(as.character(deriv),
           "0" = theta^(1/exponent),
           "1" = (theta^(1-exponent)) / exponent,
           "2" = ((1-exponent) / exponent^2) *
               (theta^(1 - 2*exponent)),
      "3" = { #3rd deriv
      stop("3rd deriv not yet implemented")
      },
           stop("argument 'deriv' unmatched"))
  } else {
    switch(as.character(deriv),
           "0" = theta^exponent,
           "1" = exponent / (theta^(1-exponent)),
      "2" = exponent * (exponent-1) * (theta^(exponent-2)),
      "3" = { #3rd deriv
      stop("3rd deriv not yet implemented")
      },
           stop("argument 'deriv' unmatched"))
  }
}  # powerlink







 extlogitlink <-
  function(theta,
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
        paste0("extlogitlink(", theta,
               ", min = ", A,
               ", max = ", B, ")") else
        paste0("extlogitlink(", theta, ")")
    } else {
      paste0("log((",
             as.char.expression(theta),
             "-min)/(max-",
             as.char.expression(theta), "))")
    }
    if (tag)
      string <- paste("Extended logit:", string)
    return(string)
  }

  if (inverse) {
    switch(as.character(deriv),
           "0" = {
           junk <- care.exp(theta)
           (A + B * junk) / (1.0 + junk) },
           "1" = ((theta - A) * (B - theta)) / (B-A),
           "2" = (A + B - 2 * theta) * (theta - A) *
                 (B - theta) / (B-A)^2,
      "3" = { #3rd deriv
           (theta - A) * (B - theta) * ((2 * theta - A - B)^2 -
       2 * (theta - A) * (B - theta)) / (B - A)^3
      },
           stop("argument 'deriv' unmatched"))
  } else {
    switch(as.character(deriv),
           "0" = {
           log((theta - A)/(B - theta))},
           "1" = (B-A) / ((theta - A) * (B - theta)),
           "2" = ((2 * theta - A - B) * (B-A)) / ((theta - A) *
                 (B - theta))^2,
      "3" = { #3rd deriv
           (B - A) * (2 / ((theta - A) * (B - theta))^2) *
      (1 + (2 * theta - A - B)^2 / ((theta - A) * (B - theta)))
      },
           stop("argument 'deriv' unmatched"))
  }
}  # extlogitlink






 logclink <-
  function(theta,
           bvalue = NULL,  # .Machine$double.xmin
           inverse = FALSE, deriv = 0,
           short = TRUE, tag = FALSE) {
  if (is.character(theta)) {
    string <- if (short)
        paste0("logclink(", theta, ")") else {
        theta <- as.char.expression(theta)
        paste0("log(1-", theta, ")")
    }
    if (tag)
      string <- paste("Log Complementary:", string)
    return(string)
  }


  if (!inverse && length(bvalue)) {
    theta[theta >= 1.0] <- bvalue;
  }
  if (inverse) {
    switch(as.character(deriv),
           "0" = -expm1(theta),
           "1" = theta - 1,
           "2" = theta - 1,
           "3" = theta - 1,
           "4" = theta - 1,
           "5" = theta - 1,
           stop("argument 'deriv' unmatched"))
  } else {
    switch(as.character(deriv),
           "0" = log1p(-theta),
           "1" =  -1 / (1 - theta),
           "2" =  -1 / (1 - theta)^2,
           "3" =  -2 / (1 - theta)^3,
           "4" =  -6 / (1 - theta)^4,
           "5" = -24 / (1 - theta)^5,
           stop("argument 'deriv' unmatched"))
  }
}  # logclink






 cauchitlink <-
  function(theta,
           bvalue = .Machine$double.eps,
           inverse = FALSE, deriv = 0,
           short = TRUE, tag = FALSE) {
  if (is.character(theta)) {
    string <- if (short)
        paste0("cauchitlink(", theta, ")") else {
        theta <- as.char.expression(theta)
        paste0("tan(pi*(", theta, "-0.5))")
    }
    if (tag)
      string <- paste("Cauchit:", string)
    return(string)
  }

  if (!inverse && length(bvalue)) {
    theta[theta <= 0.0] <- bvalue
    theta[theta >= 1.0] <- 1.0 - bvalue
  }
  if (inverse) {
    switch(as.character(deriv),
           "0" = 0.5 + atan(theta) / pi,
           "1" = (cos(pi * (theta-0.5)))^2  / pi,
           "2" = {
             temp2 <- cos(pi * (theta-0.5))
             temp4 <- sin(pi * (theta-0.5))
             -2 * temp4 * temp2^3 / pi
           },
           "3" = {
             temp2 <- cos(pi * (theta-0.5))
             temp5 <- tan(pi * (theta-0.5))
             2 * temp2^6 * (3 * temp5^2 - 1) / pi
           },
           stop("argument 'deriv' unmatched"))
  } else {
    switch(as.character(deriv),
           "0" = tan(pi * (theta-0.5)),
           "1" = pi / (cos(pi * (theta-0.5)))^2,
           "2" =  {
           temp2 <- cos(pi * (theta-0.5))
           temp3 <- tan(pi * (theta-0.5))
           (temp3 * 2 * pi^2) / temp2^2
         },
           "3" =  {
           temp2 <- cos(pi * (theta-0.5))
           temp3 <- tan(pi * (theta-0.5))
           2 * pi^3 * (1 + 3 * temp3^2) / temp2^2
         },
           stop("argument 'deriv' unmatched"))
  }
}  # cauchitlink







 nbcanlink <-
  function(theta,
           size = NULL,
           wrt.param = NULL,
           bvalue = NULL,
           inverse = FALSE, deriv = 0,
           short = TRUE, tag = FALSE) {
  if (is.character(theta)) {
      lastchars1 <-
          substr(theta, nchar(theta),
                 nchar(theta))
      lastchars2 <-
          ifelse(nchar(theta) > 1,
                 substr(theta, nchar(theta) - 1,
                        nchar(theta) - 1),
                rep("", length(theta)))

    size.names <- rep("size", length(theta))
    dig1 <- lastchars1 %in% as.character(0:9)
    dig2 <- lastchars2 %in% as.character(0:9)
    size.names <- ifelse(dig1,
                         paste0("size",
                                lastchars1),
                         size.names)
    size.names <- ifelse(dig2,
                         paste0("size",
                                lastchars2,
                                lastchars1),
                     size.names)

    string <- if (short)
      paste0("nbcanlink(", theta,
        ", ", theta,
        "(", size.names, ")",  # Added 20180803
        ")") else {
      theta <- as.char.expression(theta)
      paste0("log(", theta, " / (", theta, " + ",
             size.names, "))")
    }
    if (tag)
      string <- paste("Nbcanlink:", string)
    return(string)
  }



  kmatrix <- size
  theta <- cbind(theta)
  kmatrix <- cbind(kmatrix)
  if (ncol(kmatrix) != ncol(theta))
    stop("arguments 'theta' & 'size' do not have",
         " an equal number of cols")
  if (nrow(kmatrix) != nrow(theta))
    stop("arguments 'theta' & 'size' do not have",
         " an equal number of rows")


  if (deriv > 0) {
    if (!(wrt.param %in% 1:2))
      stop("argument 'wrt.param' should be 1 or 2")
  }


  if (!inverse && length(bvalue))
    theta[theta <= 0.0] <- bvalue

  if (inverse) {
    switch(as.character(deriv),
           "0" = {
       ans <- (kmatrix / expm1(-theta))
       if (is.matrix(ans))
         dimnames(ans) <- NULL else
         names(ans) <- NULL
       ans
       },

        "1" = if (wrt.param == 1)
     (theta * (theta + kmatrix)) / kmatrix else
          -(theta + kmatrix),

       "2" = if (wrt.param == 1)
       (2 * theta + kmatrix) * theta *
           (theta + kmatrix) / kmatrix^2 else
        theta + kmatrix,
           "3" = { #3rd deriv
      stop("3rd deriv not yet implemented")
      })
  } else {
    ans <-
    switch(as.character(deriv),
           "0" = log(theta / (theta + kmatrix)),

           "1" = if (wrt.param == 1)
                     kmatrix / (theta * (theta +
                     kmatrix)) else
                     -1 / (theta + kmatrix),

       "2" = if (wrt.param == 1)
             (2 * theta + kmatrix) *
             (-kmatrix) / (theta *
             (theta + kmatrix))^2 else
             1 / (theta + kmatrix)^2,
       "3" = { #3rd deriv
      stop("3rd deriv not yet implemented")
      })

     if (is.matrix(ans))
       dimnames(ans) <- NULL else
       names(ans) <- NULL
     ans
  }
}  # nbcanlink








 linkfunvlm <-
    function(object, earg = FALSE, ...) {
  if (!any(slotNames(object) == "extra"))
    stop("no 'extra' slot on the object")
  if (!any(slotNames(object) == "misc"))
    stop("no 'misc' slot on the object")

  M <- npred(object)

  misc <- object@misc
  LINKS1 <- misc$link
  EARGS1 <- misc$earg

  extra <- object@extra
  LINKS2 <- extra$link
  EARGS2 <- extra$earg

  if (length(LINKS1) != M &&
      length(LINKS2) != M) {
    if (LINKS1 != "multilogitlink" &&
        LINKS2 != "multilogitlink")
      warning("the length of the 'links' ",
              "component is not ", M)
  }

  if (length(LINKS1)) {
      if (earg) list(link = LINKS1,
                     earg = EARGS1) else LINKS1
  } else {
      if (earg) list(link = LINKS2,
                     earg = EARGS2) else LINKS2
  }
}  # linkfunvlm



if (!isGeneric("linkfun"))
  setGeneric("linkfun", function(object, ...)
      standardGeneric("linkfun"),
      package = "VGAM")



setMethod("linkfun", "vlm", function(object, ...)
    linkfunvlm(object, ...))







 logitoffsetlink <-
  function(theta,
           offset = 0,
           inverse = FALSE, deriv = 0,
           short = TRUE, tag = FALSE) {
  if (is.character(theta)) {
    string <- if (short)
        paste0("logitoffsetlink(",
               theta,
              ", ", offset[1],
              ")") else
        paste0("log(",
              as.char.expression(theta),
              "/(1-",
              as.char.expression(theta),
              ")",
              " - ", offset[1], ")")
    if (tag)
      string <- paste("Logit-with-offset:", string)
    return(string)
  }




  if (inverse) {
    switch(as.character(deriv),
           "0" = {
      exp.eta <- exp(theta)
      (exp.eta + offset) / (1 + exp.eta + offset)
           },
           "1" = 1 / Recall(theta = theta,
                      offset = offset,
                      inverse = FALSE,
                      deriv = deriv),
   "2" = theta * (1 - theta) * (1 - 2 * theta),
   "3" = { #3rd deriv
      stop("3rd deriv not yet implemented")
      },
           stop("argument 'deriv' unmatched"))
  } else {
    switch(as.character(deriv),
           "0" = {
       temp2 <- log(theta / (1 - theta) - offset)
       temp2
       },
       "1" = 1 / ((1 - theta) *
                  (theta - (1-theta) * offset)),
       "2" = (2 * (theta - offset *
                   (1-theta)) - 1) / (
       (theta - (1-theta)*offset) * (1-theta))^2,
           "3" = { #3rd deriv
      stop("3rd deriv not yet implemented")
      },
       stop("argument 'deriv' unmatched"))
  }
}  # logitoffsetlink






 asinlink <-
  function(theta,
           bvalue = NULL,  # orig.
           inverse = FALSE, deriv = 0,
           short = TRUE, tag = FALSE,
           c10 = c(4, -pi)) {  # === alogitlink()
  T <- TRUE; F <- FALSE
  if (!is.Numeric(c10, length.arg = 2))
    stop("bad input for 'c10'")
  if ((c1 <- c10[1]) <= 0)
    stop("c10[1] must be positive")
  c0 <- c10[2]
  plain <- c1 == 1 && c0 == 0
  fc1 <- format(c1, digits = 4)
  fc0 <- format(c0, digits = 4)
  if (is.character(theta)) {
    string <- if (short) {
      if (plain)
        paste0("asinlink(", theta, ")") else
        paste0(fc1, "*asin(sqrt(", theta, "))",
           ifelse(c0 < 0, fc0,
           ifelse(c0 > 0, paste0("+", fc0), "")))
      } else {
        Theta <- as.char.expression(theta)
        if (plain)
        paste0("sqrt(", Theta, ")") else
        paste0(fc1, "*asin(sqrt(", Theta, "))",
           ifelse(c0 < 0, fc0,
           ifelse(c0 > 0, paste0("+", fc0), "")))
    }
    if (tag)
      string <- paste("Arcsinelink:", string)
    return(string)
  }

  if (TRUE && length(bvalue)) {
    if (inverse && deriv == 0) {
      eps <- 1e-7  # 0
      theta[theta <= c0] <- c0 + eps
      theta[theta >= c0+c1*pi/2] <- c0+c1*pi/2 - eps
    } else {
      eps <- 1e-7  # 0
      theta[theta <= 0] <- 0 + eps
      theta[theta >= 1] <- 1 - eps
    }
  }



  
  if (inverse) {
    switch(as.character(deriv),
      "0" = ifelse(c0 <= theta &
                   theta <= c0+c1*pi/2,
            (sin((theta - c0) / c1))^2,
            NaN),
      "1" = 1 / Recall(theta = theta,
                   bvalue = bvalue,
                   inverse = FALSE,
                   deriv = deriv,
                   c10 = c10),
    "2" = {
  iD1 <- Recall(theta, de = 1, inv = T, c10 = c10)
  DD2 <- Recall(theta, de = 2, inv = F, c10 = c10)
  -(iD1^3) * DD2
          },
    "3" = {
  iD1 <- Recall(theta, de = 1, inv = T, c10 = c10)
  DD2 <- Recall(theta, de = 2, inv = F, c10 = c10)
  DD3 <- Recall(theta, de = 3, inv = F, c10 = c10)
  (iD1^4) * (3 * iD1 * DD2^2 - DD3)
        },
      "4" = {
  iD1 <- Recall(theta, de = 1, inv = T, c10 = c10)
  iD2 <- Recall(theta, de = 2, inv = T, c10 = c10)
  iD3 <- Recall(theta, de = 3, inv = T, c10 = c10)
  DD1 <- Recall(theta, de = 1, inv = F, c10 = c10)
  DD2 <- Recall(theta, de = 2, inv = F, c10 = c10)
  DD3 <- Recall(theta, de = 3, inv = F, c10 = c10)
  DD4 <- Recall(theta, de = 4, inv = F, c10 = c10)
  (iD1^3) * (15 * iD1 * iD2 * (DD2^2) +
              6 * (iD1^3) * DD2 * DD3 -
              4 * iD2 * DD3 - (iD1^2) * DD4)
      },
           stop("argument 'deriv' unmatched"))
  } else {
    switch(as.character(deriv),
     "0" = c1 * asin(sqrt(theta)) + c0,
     "1" = c1 * 0.5 / sqrt(theta * (1 - theta)),
     "2" = (c1 / 4) * (2 * theta - 1) / (theta *
           (1 - theta))^(3/2),
     "3" = (c1 / 4) * 1.5 * ((2 * theta - 1)^2) / (
           theta * (1 - theta))^(5/2) +
           (c1 / 4) * 2 / (theta * (1 - theta))^(3/2),
     "4" = c1 * (15 / 16) * ((2 * theta - 1)^3) / (
           theta * (1 - theta))^(7/2) +
           (c1 / 4) * 9 * (2 * theta - 1) / (theta *
           (1 - theta))^(5/2),
     stop("argument 'deriv' unmatched"))
  }
}  # asinlink







 lcalogitlink <-
  function(theta,
           bvalue = NULL,   # .Machine$double.eps,
           pmix.logit = 0.01,
           tol = 1e-13, nmax = 99,
           inverse = FALSE, deriv = 0,
           short = TRUE, tag = FALSE,
           c10 = c(4, -pi)) {  # === asinlink()
  T <- TRUE; F <- FALSE
  if (!is.Numeric(c10, length.arg = 2))
    stop("bad input for 'c10'")
  if ((c1 <- c10[1]) <= 0)
    stop("c10[1] must be positive")
  c0 <- c10[2]
  plain <- c1 == 1 && c0 == 0
  fc1 <- format(c1, digits = 4)
  fc0 <- format(c0, digits = 4)
  if (!is.Numeric(pmix.logit, length.arg = 1) ||
     pmix.logit < 0 || pmix.logit > 1)
   stop("bad input for argument 'pmix.logit'")
 
  if (is.character(theta)) {
    string <- if (short)
      paste0("lcalogitlink(", theta, ")") else {
      Theta <- as.char.expression(theta)
      paste0(1 - pmix.logit, "*asinlink(",
             Theta, ")+",
             pmix.logit, "*logitlink(",
             Theta, ")")
    }
    if (tag)
      string <- paste("LC-asin-logit:", string)
    return(string)
  }

  if (!inverse && length(bvalue)) {
    theta[theta <= 0.0] <- bvalue
    theta[theta >= 1.0] <- 1.0 - bvalue
  }

lcalogitmix <-  # pmix.logit local/global
  function(p, deriv = 0, Value = 0, c10 = c(4, -pi))
(1 - pmix.logit) *
asinlink(p, deriv = deriv, c10 = c10) +
pmix.logit *
logitlink(p, deriv = deriv) - Value
 
  if (inverse) {
    P <- pmix.logit
    switch(as.character(deriv),
    "0" = {  eps <- 0  #  1e-13
       Lo <- numeric(length(theta)) + eps
       Up <- numeric(length(theta)) + 1 - eps
       bisection.basic(lcalogitmix,
                       tol = tol, nmax = nmax,
                       Lo, Up, Value = theta,
                       c10 = c10)
    },
    "1" =  1 / Recall(theta, pm = P, deriv = 1, c10 = c10),
    "2" = {
  iD1 <- Recall(theta, de = 1, pm=P, inv = T, c10 = c10)
  DD2 <- Recall(theta, de = 2, pm=P, inv = F, c10 = c10)
  -(iD1^3) * DD2
          },
    "3" = {
  iD1 <- Recall(theta, de = 1, pm=P, inv = T, c10 = c10)
  DD2 <- Recall(theta, de = 2, pm=P, inv = F, c10 = c10)
  DD3 <- Recall(theta, de = 3, pm=P, inv = F, c10 = c10)
  (iD1^4) * (3 * iD1 * DD2^2 - DD3)
        },
        "4" = {
  iD1 <- Recall(theta, de = 1, pm=P, inv = T, c10 = c10)
  iD2 <- Recall(theta, de = 2, pm=P, inv = T, c10 = c10)
  iD3 <- Recall(theta, de = 3, pm=P, inv = T, c10 = c10)
  DD1 <- Recall(theta, de = 1, pm=P, inv = F, c10 = c10)
  DD2 <- Recall(theta, de = 2, pm=P, inv = F, c10 = c10)
  DD3 <- Recall(theta, de = 3, pm=P, inv = F, c10 = c10)
  DD4 <- Recall(theta, de = 4, pm=P, inv = F, c10 = c10)
  (iD1^3) * (15 * iD1 * iD2 * (DD2^2) +
              6 * (iD1^3) * DD2 * DD3 -
              4 * iD2 * DD3 - (iD1^2) * DD4)
        },
           stop("argument 'deriv' unmatched"))
  } else {
    switch(as.character(deriv),
    "0" = lcalogitmix(theta, deriv = 0, c10 = c10),
    "1" = lcalogitmix(theta, deriv = 1, c10 = c10),
    "2" = lcalogitmix(theta, deriv = 2, c10 = c10),
    "3" = lcalogitmix(theta, deriv = 3, c10 = c10),
    "4" = lcalogitmix(theta, deriv = 4, c10 = c10),
      stop("argument 'deriv' unmatched"))
  }
}  # lcalogitlink





 sqrtlink <-
  function(theta,
           inverse = FALSE, deriv = 0,
           short = TRUE, tag = FALSE,
           c10 = c(2, -2)) {  # poissonff-sloglink
  T <- TRUE; F <- FALSE
  if (!is.Numeric(c10, length.arg = 2))
    stop("bad input for 'c10'")
  if ((c1 <- c10[1]) <= 0)
    stop("c10[1] must be positive")
  c0 <- c10[2]
  plain <- c1 == 1 && c0 == 0
  fc1 <- format(c1, digits = 4)
  fc0 <- format(c0, digits = 4)
  if (is.character(theta)) {
    string <- if (short) {
      if (plain)
        paste0("sqrtlink(", theta, ")") else
        paste0(fc1, "*sqrtlink(", theta, ")",
           ifelse(c0 < 0, fc0,
           ifelse(c0 > 0, paste0("+", fc0), "")))
      } else {
        theta <- as.char.expression(theta)
        if (plain)
        paste0("sqrt(", theta, ")") else
        paste0(fc1, "*sqrt(", theta, ")",
           ifelse(c0 < 0, fc0,
           ifelse(c0 > 0, paste0("+", fc0), "")))
    }
    if (tag)
      string <- paste("Square root:", string)
    return(string)
  }

  if (inverse) {
    switch(as.character(deriv),
      "0" = ((theta - c0) / c1)^2,
      "1" = 2 * sqrt(theta) / c1,
      "2" = 2 / c1^2,
      "3" = 0, "4" = 0, "5" = 0, "6" = 0,
      stop("argument 'deriv' unmatched"))
  } else {
    switch(as.character(deriv),
       "0" =  c1 * sqrt(theta) + c0,
       "1" =  c1 * 0.5 / sqrt(theta),
       "2" = -c1 * 0.25 / theta^1.5,
       "3" =  c1 * (3 / 8) / theta^2.5,
       "4" = -c1 * (15 / 16) / theta^3.5,
       "5" = stop("5th deriv not implemented"),
    stop("argument 'deriv' unmatched"))
  }
}  # sqrtlink






 lcsloglink <-
  function(theta,
           bvalue = NULL,   # .Machine$double.eps,
           pmix.log = 0.01,
           tol = 1e-13, nmax = 99,
           inverse = FALSE, deriv = 0,
           short = TRUE, tag = FALSE,
           c10 = c(2, -2)) {  # poissonff-sloglink
  T <- TRUE; F <- FALSE
  if (!is.Numeric(c10, length.arg = 2))
    stop("bad input for 'c10'")
  if ((c1 <- c10[1]) <= 0)
    stop("c10[1] must be positive")
  c0 <- c10[2]
  plain <- c1 == 1 && c0 == 0
  fc1 <- format(c1, digits = 4)
  fc0 <- format(c0, digits = 4)
  if (!is.Numeric(pmix.log, length.arg = 1) ||
      pmix.log < 0 || pmix.log > 1)
    stop("bad input for argument 'pmix.log'")
 
  if (is.character(theta)) {
    string <- if (short)
      paste0("lcsloglink(", theta, ")") else {
      Theta <- as.char.expression(theta)
      paste0(1- pmix.log, "*sqrtlink(",
             Theta, ")+",
             pmix.log, "*loglink(",
             Theta, ")")
    }
    if (tag)
      string <- paste("LC-sqrt-log:", string)
    return(string)
  }

  if (!inverse && length(bvalue)) {
    theta[theta <= 0.0] <- bvalue
    theta[theta >= 1.0] <- 1.0 - bvalue
  }

lcslogmix <-  # pmix.log local/global
  function(mu, deriv = 0, Value = 0, c10 = c(2, -2))
(1 - pmix.log) *
sqrtlink(mu, deriv = deriv, c10 = c10) +
pmix.log *
loglink(mu, deriv = deriv) - Value

  if (inverse) {
    P <- pmix.log
    switch(as.character(deriv),
    "0" = {  smallno <- 1e-10  # 1e-12
       Lo <- pmin(((theta - c0) / c1)^2,
                  exp(theta)) - 2.5
       Up <- pmax(((theta - c0) / c1)^2,
                  exp(theta)) + 2.5
       Up <- pmax(Up, max(theta, na.rm = TRUE))
       bisection.basic(lcslogmix,
           tol = tol, nmax = nmax,
           Lo, Up, Value = theta, c10 = c10)
    },
    "1" =  1 / Recall(theta, pm = P, deriv = 1, c10 = c10),
    "2" = {
  iD1 <- Recall(theta, der = 1, inv = T, c10 = c10)
  DD2 <- Recall(theta, der = 2, inv = F, c10 = c10)
  -(iD1^3) * DD2
          },
    "3" = {
  iD1 <- Recall(theta, de = 1, pm=P, inv = T, c10 = c10)
  DD2 <- Recall(theta, de = 2, pm=P, inv = F, c10 = c10)
  DD3 <- Recall(theta, de = 3, pm=P, inv = F, c10 = c10)
  (iD1^4) * (3 * iD1 * DD2^2 - DD3)
        },
        "4" = {
  iD1 <- Recall(theta, de = 1, pm=P, inv = T, c10 = c10)
  iD2 <- Recall(theta, de = 2, pm=P, inv = T, c10 = c10)
  iD3 <- Recall(theta, de = 3, pm=P, inv = T, c10 = c10)
  DD1 <- Recall(theta, de = 1, pm=P, inv = F, c10 = c10)
  DD2 <- Recall(theta, de = 2, pm=P, inv = F, c10 = c10)
  DD3 <- Recall(theta, de = 3, pm=P, inv = F, c10 = c10)
  DD4 <- Recall(theta, de = 4, pm=P, inv = F, c10 = c10)
  (iD1^3) * (15 * iD1 * iD2 * (DD2^2) +
              6 * (iD1^3) * DD2 * DD3 -
              4 * iD2 * DD3 - (iD1^2) * DD4)
        },
           stop("argument 'deriv' unmatched"))
  } else {
    switch(as.character(deriv),
    "0" = lcslogmix(theta, deriv = 0, c10 = c10),
    "1" = lcslogmix(theta, deriv = 1, c10 = c10),
    "2" = lcslogmix(theta, deriv = 2, c10 = c10),
    "3" = lcslogmix(theta, deriv = 3, c10 = c10),
    "4" = lcslogmix(theta, deriv = 4, c10 = c10),
      stop("argument 'deriv' unmatched"))
  }
}  # lcsloglink







   sloglink <-
  function(theta,
           bvalue = NULL,   # .Machine$double.eps,
           taumix.log = 1,  # [0, Inf)
           tol = 1e-13, nmax = 99,
           inverse = FALSE, deriv = 0,
           short = TRUE, tag = FALSE,
           c10 = c(2, -2)) {  # poissonff-sloglink
  T <- TRUE; F <- FALSE
  if (!is.Numeric(c10, length.arg = 2))
    stop("bad input for 'c10'")
  if ((c1 <- c10[1]) <= 0)
    stop("c10[1] must be positive")
  c0 <- c10[2]
  plain <- c1 == 1 && c0 == 0
  fc1 <- format(c1, digits = 4)
  fc0 <- format(c0, digits = 4)
  if (!is.Numeric(taumix.log, length.arg = 1) ||
      taumix.log < 0)
    stop("bad input for argument 'taumix.log'")
 
  if (is.character(theta)) {
    string <- if (short)
      paste0("sloglink(", theta, ")") else {
      Theta <- as.char.expression(theta)
      paste0("-expm1(-p*", taumix.log,
             ")*sqrtlink(",
             Theta, ")+exp(-p*",
             taumix.log, ")*loglink(",
             Theta, ")")
    }
    if (tag)
      string <- paste("EW-sqrt-log:", string)
    return(string)
  }

  if (!inverse && length(bvalue)) {
    theta[theta <= 0.0] <- bvalue
    theta[theta >= 1.0] <- 1.0 - bvalue
  }

ewslogmix <-
  function(mu, deriv = 0, Value = 0,
           c10 = c(2, -2)) {
  eta <- (-taumix.log * mu)
  if (deriv == 0) return(
    -expm1(eta) *
    sqrtlink(mu, deriv = deriv, c10 = c10) +
    exp(eta) *
    loglink(mu, deriv = deriv) - Value)
  if (deriv == 1) {
    dd1 <- (-taumix.log)
    return(
    sqrtlink(mu, deriv = deriv, c10 = c10) -
    exp(eta) * (
    sqrtlink(mu, deriv = deriv, c10 = c10) -
     loglink(mu, deriv = deriv) + dd1 * (
    sqrtlink(mu, c10 = c10) - loglink(mu))))
    }
  if (deriv == 2) {
    dd1 <- (-taumix.log)
    dd2 <- 0
    return(
    sqrtlink(mu, deriv = deriv, c10 = c10) -
    exp(eta) * (
    sqrtlink(mu, deriv = deriv, c10 = c10) -
     loglink(mu, deriv = deriv) +
    (dd2 + dd1^2) * ( 
    sqrtlink(mu, c10 = c10) - loglink(mu)) +
    2 * dd1 * ( 
    sqrtlink(mu, deriv = 1, c10 = c10) -
    loglink(mu, deriv = 1))))
  }
  if (deriv == 3) {
    dd1 <- (-taumix.log)
    dd2 <- dd3 <- 0
    return(
    sqrtlink(mu, deriv = deriv, c10 = c10) -
    exp(eta) * (
    sqrtlink(mu, deriv = deriv, c10 = c10) -
     loglink(mu, deriv = deriv) +
    (dd3 + dd1 * (3 * dd2 + dd1^2)) * ( 
        sqrtlink(mu, c10 = c10) -
        loglink(mu)) + 3 * dd1 * ( 
    sqrtlink(mu, deriv = 2, c10 = c10) -
     loglink(mu, deriv = 2)) +
    3 * (dd2 + dd1^2) * ( 
    sqrtlink(mu, deriv = 1, c10 = c10) -
     loglink(mu, deriv = 1))))
  }
  if (deriv == 4) {
    stop("yettodo")
  }
  }  # ewslogmix

     
  if (inverse) {
    P <- taumix.log      
    switch(as.character(deriv),
    "0" = { eps <- 0   # 1e-14
       Lo <- pmin(((theta - c0) / c1)^2,
                  exp(theta)) - 2.5
       Lo <- pmax(Lo, eps)  # 20240110
       Up <- pmax(((theta - c0) / c1)^2,
                  exp(theta)) + 2.5
       Up <- pmax(Up, max(theta, na.rm = TRUE))
       bisection.basic(ewslogmix, c10 = c10,
                       tol = tol, nmax = nmax,
                       Lo, Up, Value = theta)
    },
    "1" =  1 / Recall(theta, tau=P, der = 1, c10 = c10),
    "2" = {
  iD1 <- Recall(theta, de = 1, tau=P, inv = T, c10 = c10)
  DD2 <- Recall(theta, de = 2, tau=P, inv = F, c10 = c10)
  -(iD1^3) * DD2
          },
    "3" = {
  iD1 <- Recall(theta, de = 1, tau=P, inv = T, c10 = c10)
  DD2 <- Recall(theta, de = 2, tau=P, inv = F, c10 = c10)
  DD3 <- Recall(theta, de = 3, tau=P, inv = F, c10 = c10)
  (iD1^4) * (3 * iD1 * DD2^2 - DD3)
        },
        "4" = {
  iD1 <- Recall(theta, de = 1, tau=P, inv = T, c10 = c10)
  iD2 <- Recall(theta, de = 2, tau=P, inv = T, c10 = c10)
  iD3 <- Recall(theta, de = 3, tau=P, inv = T, c10 = c10)
  DD1 <- Recall(theta, de = 1, tau=P, inv = F, c10 = c10)
  DD2 <- Recall(theta, de = 2, tau=P, inv = F, c10 = c10)
  DD3 <- Recall(theta, de = 3, tau=P, inv = F, c10 = c10)
  DD4 <- Recall(theta, de = 4, tau=P, inv = F, c10 = c10)
  (iD1^3) * (15 * iD1 * iD2 * (DD2^2) +
              6 * (iD1^3) * DD2 * DD3 -
              4 * iD2 * DD3 - (iD1^2) * DD4)
        },
           stop("argument 'deriv' unmatched"))
  } else {
    switch(as.character(deriv),
    "0" = ewslogmix(theta, deriv = 0, c10 = c10),
    "1" = ewslogmix(theta, deriv = 1, c10 = c10),
    "2" = ewslogmix(theta, deriv = 2, c10 = c10),
    "3" = ewslogmix(theta, deriv = 3, c10 = c10),
    "4" = ewslogmix(theta, deriv = 4, c10 = c10),
      stop("argument 'deriv' unmatched"))
  }
}  # sloglink





   alogitlink <-
  function(theta,
           bvalue = NULL,   # .Machine$double.eps,
           taumix.logit = 1,  # [0, Inf)
           tol = 1e-13, nmax = 99,
           inverse = FALSE, deriv = 0,
           short = TRUE, tag = FALSE,
           c10 = c(4, -pi)) {  # === asinlink()
  T <- TRUE; F <- FALSE
  if (!is.Numeric(c10, length.arg = 2))
    stop("bad input for 'c10'")
  if ((c1 <- c10[1]) <= 0)
    stop("c10[1] must be positive")
  c0 <- c10[2]
  plain <- c1 == 1 && c0 == 0
  fc1 <- format(c1, digits = 4)
  fc0 <- format(c0, digits = 4)
  ftm <- format(taumix.logit, digits = 4)
  if (!is.Numeric(taumix.logit, length.arg = 1) ||
      taumix.logit < 0)
    stop("bad input for argument 'taumix.logit'")
 
  if (is.character(theta)) {
    string <- if (short)
      paste0("alogitlink(", theta, ")") else {
      Theta <- as.char.expression(theta)
      paste0("-expm1(-p*(1-p)*", ftm,
             ")*asinlink(",
             Theta, ")+exp(-p*(1-p)*",
             ftm, ")*logitlink(",
             Theta, ")")
    }
    if (tag)
      string <- paste("EW-asin-logit:", string)
    return(string)
  }

  if (!inverse && length(bvalue)) {
    theta[theta <= 0.0] <- bvalue
    theta[theta >= 1.0] <- 1.0 - bvalue
  }

ewalogitmix <-
  function(mu, deriv = 0, Val = 0,
           c10 = c(4, -pi)) {
  eta <- (-taumix.logit * mu * (1 - mu))
  if (deriv == 0) return(
    -expm1(eta) *
    asinlink(mu, deriv = deriv, c10 = c10) +
    exp(eta) *
    logitlink(mu, deriv = deriv) - Val)
  if (deriv == 1) {
    dd1 <- (-taumix.logit) * (1 - 2 * mu)
    return(
     asinlink(mu, deriv = deriv, c10 = c10) -
     exp(eta) * (
     asinlink(mu, deriv = deriv, c10 = c10) -
     logitlink(mu, deriv = deriv) +
     dd1 * (
         asinlink(mu, c10 = c10) -
         logitlink(mu))))
    }
  if (deriv == 2) {
    dd1 <- (-taumix.logit) * (1 - 2 * mu)
    dd2 <-   taumix.logit * 2
    return(
     asinlink(mu, deriv = deriv, c10 = c10) -
     exp(eta) * (
     asinlink(mu, deriv = deriv, c10 = c10) -
    logitlink(mu, deriv = deriv) +
    (dd2 + dd1^2) * ( 
     asinlink(mu, c10 = c10) -
    logitlink(mu)) + 2 * dd1 * ( 
     asinlink(mu, deriv = 1, c10 = c10) -
    logitlink(mu, deriv = 1))))
  }
  if (deriv == 3) {
    dd1 <- (-taumix.logit) * (1 - 2 * mu)
    dd2 <-   taumix.logit * 2
    dd3 <- 0
    return(
     asinlink(mu, deriv = deriv, c10 = c10) -
     exp(eta) * (
     asinlink(mu, deriv = deriv, c10 = c10) -
    logitlink(mu, deriv = deriv) +
    (dd3 + dd1 * (3 * dd2 + dd1^2)) * ( 
     asinlink(mu, c10 = c10) -
     logitlink(mu)) + 3 * dd1 * ( 
     asinlink(mu, deriv = 2, c10 = c10) -
    logitlink(mu, deriv = 2)) +
    3 * (dd2 + dd1^2) * ( 
     asinlink(mu, deriv = 1, c10 = c10) -
    logitlink(mu, deriv = 1))))
  }
  if (deriv == 4) {
    stop("yettodo")
  }
  }  # ewalogitmix
      
  if (inverse) {
    P <- taumix.logit
    switch(as.character(deriv),
    "0" = {  eps <- 0  # 1e-13
       Lo <- numeric(length(theta)) + eps
       Up <- numeric(length(theta)) + 1 - eps
       bisection.basic(ewalogitmix,
                       tol = tol, nmax = nmax,
                       Lo, Up, Val = theta,
                       c10 = c10)
    },
    "1" =  1 / Recall(theta, tau=P, der = 1, c10 = c10),
    "2" = {
  iD1 <- Recall(theta, de = 1, tau=P, inv = T, c10 = c10)
  DD2 <- Recall(theta, de = 2, tau=P, inv = F, c10 = c10)
  -(iD1^3) * DD2
          },
    "3" = {
  iD1 <- Recall(theta, de = 1, tau=P, inv = T, c10 = c10)
  DD2 <- Recall(theta, de = 2, tau=P, inv = F, c10 = c10)
  DD3 <- Recall(theta, de = 3, tau=P, inv = F, c10 = c10)
  (iD1^4) * (3 * iD1 * DD2^2 - DD3)
        },
        "4" = {
  iD1 <- Recall(theta, de = 1, tau=P, inv = T, c10 = c10)
  iD2 <- Recall(theta, de = 2, tau=P, inv = T, c10 = c10)
  iD3 <- Recall(theta, de = 3, tau=P, inv = T, c10 = c10)
  DD1 <- Recall(theta, de = 1, tau=P, inv = F, c10 = c10)
  DD2 <- Recall(theta, de = 2, tau=P, inv = F, c10 = c10)
  DD3 <- Recall(theta, de = 3, tau=P, inv = F, c10 = c10)
  DD4 <- Recall(theta, de = 4, tau=P, inv = F, c10 = c10)
  (iD1^3) * (15 * iD1 * iD2 * (DD2^2) +
              6 * (iD1^3) * DD2 * DD3 -
              4 * iD2 * DD3 - (iD1^2) * DD4)
        },
           stop("argument 'deriv' unmatched"))
  } else {
    switch(as.character(deriv),
    "0" = ewalogitmix(theta, deriv = 0, c10 = c10),
    "1" = ewalogitmix(theta, deriv = 1, c10 = c10),
    "2" = ewalogitmix(theta, deriv = 2, c10 = c10),
    "3" = ewalogitmix(theta, deriv = 3, c10 = c10),
    "4" = ewalogitmix(theta, deriv = 4, c10 = c10),
      stop("argument 'deriv' unmatched"))
  }
}  # alogitlink




























# These functions are
# Copyright (C) 1998-2018 T.W. Yee, University of Auckland.
# All rights reserved.











if (FALSE)
log1pexp <- function(x) {

  ans <- log1p(exp(x))
  big <- (x > 10)
  ans[big] <- x[big] + log1p(exp(-x[big]))
  ans
}







erf <- function(x, inverse = FALSE) {
  if (inverse) {
    ans <- qnorm((x+1)/2) / sqrt(2)
    ans[x <  -1] <- NA
    ans[x >  +1] <- NA
    ans[x == -1] <- -Inf
    ans[x == +1] <-  Inf
    ans
  } else {
    2 * pnorm(x * sqrt(2)) - 1
  }
}



erfc <- function(x, inverse = FALSE) {
  if (inverse) {
    ans <- qnorm(x/2, lower.tail = FALSE) / sqrt(2)
    ans[x <  0] <- NA
    ans[x >  2] <- NA
    ans[x == 0] <-  Inf
    ans[x == 2] <- -Inf
    ans
  } else {
    2 * pnorm(x * sqrt(2), lower.tail = FALSE)
  }
}







lambertW <- function(x, tolerance = 1.0e-10, maxit = 50) {
  if (any(Im(x) != 0.0))
    stop("argument 'x' must be real, not complex!")

  ans <- x
  ans[!is.na(x) & x <  -exp(-1)] <- NA
  ans[!is.na(x) & x >= -exp(-1)] <- log1p(x[!is.na(x) & x >= -exp(-1)])
  ans[!is.na(x) & x >= 0       ] <-  sqrt(x[!is.na(x) & x >= 0   ]) / 2

  cutpt <- 3.0
  if (any(myTF <- !is.na(x) & x > cutpt)) {
    L1 <- log(x[!is.na(x) & x > cutpt])  # log(as.complex(x))
    L2 <- log(L1)  # log(as.complex(L1))
    wzinit <- L1 - L2 +
          (L2 +
          (L2*( -2 + L2)/(2) +
          (L2*(  6 + L2*(-9 + L2*   2)) / (6) +
           L2*(-12 + L2*(36 + L2*(-22 + L2*3))) / (12*L1)) / L1) / L1) / L1

    ans[myTF] <- wzinit
  }

  for (ii in 1:maxit) {
    exp1 <- exp(ans)
    exp2 <- ans * exp1
    delta <- (exp2 - x) / (exp2 + exp1 -
                ((ans + 2) * (exp2 - x) / (2 * (ans + 1.0))))
    ans <- ans - delta
    if (all(is.na(delta) ||
        max(abs(delta), na.rm = TRUE) < tolerance)) break
    if (ii == maxit)
      warning("did not converge")
  }
  ans[x == Inf] <- Inf
  ans
}






 pgamma.deriv <- function(q, shape, tmax = 100) {

  nnn <- max(length(q), length(shape))
  if (length(q)     != nnn) q     <- rep_len(q,     nnn)
  if (length(shape) != nnn) shape <- rep_len(shape, nnn)

  if (!is.Numeric(q, positive = TRUE))
    stop("bad input for argument 'q'")
  if (!is.Numeric(shape, positive = TRUE))
    stop("bad input for argument 'shape'")

  if (!is.Numeric(tmax, length.arg = 1, positive = TRUE))
    stop("bad input for argument 'tmax'")
  if (tmax < 10)
    warning("probably argument 'tmax' is too small")


  gplog  <- lgamma(shape)
  gp1log <- gplog + log(shape)
  psip   <- digamma(shape)
  psip1  <- psip + 1 / shape
  psidp  <- trigamma(shape)
  psidp1 <- psidp - 1 / shape^2

  fred <-
    .C("VGAM_C_vdigami",
         d = as.double(matrix(0, 6, nnn)),
         x = as.double(q), p = as.double(shape),
         as.double(gplog), as.double(gp1log), as.double(psip),
         as.double(psip1), as.double(psidp), as.double(psidp1),
         ifault = integer(nnn),
         tmax = as.double(tmax),
         as.integer(nnn))
  answer <- matrix(fred$d, nnn, 6, byrow = TRUE)
  dimnames(answer) <- list(names(q),
                           c("q", "q^2", "shape", "shape^2",
                             "q.shape", "pgamma(q, shape)"))

  if (any(fred$ifault != 0)) {
    indices <- which(fred$ifault != 0)
    warning("convergence problems with elements ",
             indices)
  }

  answer
}








expint <- function (x, deriv = 0) {
  if (deriv == 0) {
    LLL <- length(x)
    answer <- .C("sf_C_expint", x = as.double(x), size = as.integer(LLL),
                 ans = double(LLL))$ans
    answer[x < 0] <- NA
    answer[x == 0] <- NA
    answer
  } else {
    if (!is.Numeric(deriv, integer.valued = TRUE, positive = TRUE) ||
        deriv > 3)
      stop("Bad input for argument 'deriv'")
    answer <- rep_len(0, length(x))
    if (deriv == 1) {
      answer <- exp(x) / x
    }
    if (deriv == 2) {
      answer <- exp(x) / x - exp(x) / x^2
    }
    if (deriv == 3) {
      answer <- exp(x) / x - 2 * exp(x) / x^2 +
        2 * exp(x) / x^3
    }
    answer
  }
}


expexpint <- function (x, deriv = 0) {
  LLL <- length(x)
  answer <- .C("sf_C_expexpint", x = as.double(x),
                size = as.integer(LLL), ans = double(LLL))$ans
  answer[x <  0] <- NA
  answer[x == 0] <- NA
  if (deriv > 0) {
    if (!is.Numeric(deriv, integer.valued = TRUE, positive = TRUE) ||
        deriv > 3)
      stop("Bad input for argument 'deriv'")
    if (deriv >= 1) {
      answer <- -answer + 1 / x
    }
    if (deriv >= 2) {
      answer <- -answer - 1 / x^2
    }
    if (deriv == 3) {
      answer <- -answer + 2 / x^3
    }
  }
  answer
}


expint.E1 <- function (x, deriv = 0) {
  if (deriv == 0) {
    LLL <- length(x)
    answer <- .C("sf_C_expint_e1", x = as.double(x),
                  size = as.integer(LLL), ans = double(LLL))$ans
    answer[x < 0] <- NA
    answer[x == 0] <- NA
  } else {
    if (!is.Numeric(deriv, integer.valued = TRUE, positive = TRUE) ||
        deriv > 3)
      stop("Bad input for argument 'deriv'")
    answer <- rep_len(0, length(x))
    if (deriv == 1) {
      answer <- exp(-x) / x
    }
    if (deriv == 2) {
      answer <- exp(-x) / x + exp(-x) / x^2
    }
    if (deriv == 3) {
      answer <- exp(-x) / x + 2 * exp(-x) / x^2 +
        2 * exp(-x) / x^3
    }
    answer <- (-1)^deriv * answer
  }
  answer
}











if (FALSE)
expint <- function(x) {


  LLL <- length(x)
  answer <- .C("sf_C_expint",
                 x = as.double(x),
                 size = as.integer(LLL),
                 ans = double(LLL))$ans

  answer[x  < 0] <- NA
  answer[x == 0] <- NA

  answer
}



if (FALSE)
expexpint <- function(x) {




  LLL <- length(x)
  answer <- .C("sf_C_expexpint",
                 x = as.double(x),
                 size = as.integer(LLL),
                 ans = double(LLL))$ans

  answer[x  < 0] <- NA
  answer[x == 0] <- NA

  answer
}


if (FALSE)
pochhammer <- function (x, n) {
  exp(lgamma(x+n) - lgamma(x))
}







if (FALSE)
expint.E1 <- function(x) {




  LLL <- length(x)
  answer <- .C("sf_C_expint_e1",
                 x = as.double(x),
                 size = as.integer(LLL),
                 ans = double(LLL))$ans

  answer[x  < 0] <- NA
  answer[x == 0] <- NA

  answer
}




 Zeta.aux <- function(shape, qq, shift = 1) {



  LLL <- max(length(shape), length(qq))
  if (length(shape) != LLL) shape <- rep_len(shape, LLL)
  if (length(qq   ) != LLL) qq    <- rep_len(qq,    LLL)

  if (any(qq < 12-1))
    warning("all values of argument 'q' should be 12 or more")
  aa <- qq


  B2 <- c(1/6, -1/30, 1/42, -1/30, 5/66, -691/2730, 7/6, -3617/510)
  kk <- length(B2)  # 8
  ans <- 1 / ((shape-1) * (shift + aa)^(shape-1)) +
         0.5 / (shift + aa)^shape

  term <- (shape/2) / (shift + aa)^(shape+1)
  ans <- ans + term * B2[1]

  for (mm in 2:kk) {
    term <- term * (shape+2*mm-3) *
            (shape+2*mm-2) / ((2*mm-1) * 2 * mm * (shift + aa)^2)
    ans <- ans + term * B2[mm]
  }
  ifelse(aa - 1 <= qq, ans, rep(0, length(ans)))  # Handled above
}







 zeta <- function(x, deriv = 0, shift = 1) {




  deriv.arg <- deriv
  rm(deriv)
  if (!is.Numeric(deriv.arg, length.arg = 1,
                  integer.valued = TRUE))
    stop("'deriv' must be a single non-negative integer")
  if (deriv.arg < 0 || deriv.arg > 2)
    stop("'deriv' must be 0, 1, or 2")


  if (deriv.arg > 0)
    return(Zeta.derivative(x, deriv.arg = deriv.arg, shift = shift))



  if (any(special <- Re(x) <= 1)) {
    ans <- x
    ans[special] <- Inf  # For Re(x) == 1

    special3 <- Re(x) < 1
    ans[special3] <- NA  # For 0 < Re(x) < 1

    special4 <- (0 < Re(x)) & (Re(x) < 1) & (Im(x) == 0)
    ans[special4] <-
      Zeta.derivative(x[special4], deriv.arg = deriv.arg, shift = shift)


    special2 <- Re(x) < 0
    if (any(special2)) {
      x2 <- x[special2]
      cx <- 1 - x2
      ans[special2] <- 2^(x2) * pi^(x2-1) * sin(pi*x2/2) *
                       gamma(cx) * Recall(cx)
    }  # special2

    if (any(!special)) {
      ans[!special] <- Recall(x[!special])
    }
    return(ans)
  }  # special

  aa <- 12
  ans <- 0
  for (ii in 0:(aa-1))
    ans <- ans + 1 / (shift + ii)^x

  ans <- ans + Zeta.aux(shape = x, aa, shift = shift)
  ans[shift <= 0] <- NaN
  ans
}  # zeta



 Zeta.derivative <- function(x, deriv.arg = 0, shift = 1) {



  if (!all(shift == 1))
    stop("currently 'shift' must all be 1")


  if (!is.Numeric(deriv.arg, length.arg = 1,
                  integer.valued = TRUE))
    stop("'deriv.arg' must be a single non-negative integer")
  if (deriv.arg < 0 || deriv.arg > 2)
    stop("'deriv.arg' must be 0, 1, or 2")

  if (any(Im(x) != 0))
    stop("Sorry, currently can only handle x real, not complex")
  if (any(x < 0))
    stop("Sorry, currently cannot handle x < 0")

  ok <- is.finite(x) & x > 0 & x != 1  # Handles NAs
  ans <- rep_len(NA_real_, length(x))
  nn <- sum(ok)  # Effective length (excludes x < 0 and x = 1 values)
  if (nn)
    ans[ok] <- .C("vzetawr", as.double(x[ok]), ans = double(nn),
                  as.integer(deriv.arg), as.integer(nn))$ans



  if (deriv.arg == 0)
    ans[is.finite(x) & abs(x) < 1.0e-12] <- -0.5

  ans
}





mills.ratio <- function(x) {
  ans <- exp(dnorm(x, log = TRUE) - pnorm(x, log.p = TRUE))
  if (any(vecTF <- (x < -1e2))) {
    xvneg <- x[vecTF]
    ans[vecTF] <- -xvneg / (1 - 1/xvneg^2 + 3 / xvneg^4)
  }
  ans
}  # mills.ratio()



mills.ratio2 <- function(x) {
  ans <- exp(2 * dnorm(x, log = TRUE) - pnorm(x, log.p = TRUE))
  ans[x < -40] <- 0
  ans
}  # mills.ratio2()











ghn100 <-
    c(-13.4064873381449, -12.8237997494878, -12.3429642228597,
      -11.9150619431142, 
-11.521415400787, -11.1524043855851, -10.8022607536847, -10.4671854213428, 
-10.1445099412928, -9.83226980777795, -9.5289658233901, -9.23342089021916, 
-8.94468921732547, -8.66199616813451, -8.38469694041627, -8.11224731116279, 
-7.84418238446082, -7.58010080785749, -7.31965282230454, -7.06253106024886, 
-6.80846335285879, -6.55720703192154, -6.30854436111214, -6.0622788326143, 
-5.81823213520352, -5.57624164932992, -5.33615836013836, -5.09784510508914, 
-4.86117509179121, -4.62603063578716, -4.39230207868269, -4.15988685513103, 
-3.92868868342767, -3.69861685931849, -3.46958563641859, -3.24151367963101, 
-3.01432358033115, -2.78794142398199, -2.56229640237261, -2.33732046390688, 
-2.11294799637119, -1.88911553742701, -1.66576150874151, -1.44282597021593, 
-1.22025039121895, -0.997977436098106, -0.775950761540146,
-0.554114823591618, 
-0.332414692342232, -0.11079587242244, 0.110795872422439,
0.332414692342232, 
0.554114823591617, 0.775950761540145, 0.997977436098105, 1.22025039121895, 
1.44282597021593, 1.66576150874151, 1.88911553742701, 2.11294799637119, 
2.33732046390688, 2.5622964023726, 2.78794142398199, 3.01432358033115, 
3.24151367963101, 3.46958563641859, 3.69861685931849, 3.92868868342767, 
4.15988685513103, 4.39230207868269, 4.62603063578716, 4.86117509179121, 
5.09784510508914, 5.33615836013836, 5.57624164932993, 5.81823213520351, 
6.0622788326143, 6.30854436111214, 6.55720703192153, 6.80846335285879, 
7.06253106024886, 7.31965282230453, 7.58010080785749, 7.84418238446082, 
8.11224731116279, 8.38469694041626, 8.66199616813451, 8.94468921732548, 
9.23342089021915, 9.52896582339012, 9.83226980777795, 10.1445099412928, 
10.4671854213428, 10.8022607536847, 11.1524043855851, 11.521415400787, 
11.9150619431142, 12.3429642228597, 12.8237997494878, 13.4064873381449
)



ghw100 <-
c(5.90806786503149e-79, 1.97286057487953e-72, 3.08302899000321e-67, 
9.01922230369242e-63, 8.51888308176111e-59, 3.45947793647603e-55, 
7.19152946346349e-52, 8.59756395482676e-49, 6.42072520534849e-46, 
3.18521787783596e-43, 1.10047068271428e-40, 2.74878488435709e-38, 
5.11623260438594e-36, 7.27457259688812e-34, 8.06743427870884e-32, 
7.10181222638517e-30, 5.03779116621273e-28, 2.91735007262926e-26, 
1.39484152606877e-24, 5.56102696165936e-23, 1.86499767513029e-21, 
5.30231618313167e-20, 1.28683292112113e-18, 2.68249216476057e-17, 
4.82983532170314e-16, 7.5488968779154e-15, 1.02887493735098e-13, 
1.22787851441009e-12, 1.28790382573158e-11, 1.19130063492903e-10, 
9.74792125387112e-10, 7.07585728388942e-09, 4.568127508485e-08, 
2.62909748375372e-07, 1.35179715911036e-06, 6.22152481777778e-06, 
2.56761593845487e-05, 9.51716277855096e-05, 0.000317291971043304, 
0.000952692188548621, 0.00257927326005907, 0.00630300028560806, 
0.0139156652202317, 0.0277791273859335, 0.0501758126774289,
0.0820518273912242, 
0.121537986844105, 0.163130030502782, 0.198462850254188,
0.218892629587438, 
0.21889262958744, 0.198462850254186, 0.163130030502783,
0.121537986844104, 
0.082051827391225, 0.0501758126774289, 0.0277791273859336,
0.0139156652202318, 
0.00630300028560809, 0.00257927326005912, 0.000952692188548612, 
0.000317291971043303, 9.51716277855086e-05, 2.5676159384549e-05, 
6.22152481777782e-06, 1.35179715911039e-06, 2.62909748375376e-07, 
4.56812750848495e-08, 7.07585728388942e-09, 9.74792125387167e-10, 
1.19130063492907e-10, 1.28790382573154e-11, 1.22787851441012e-12, 
1.02887493735101e-13, 7.5488968779154e-15, 4.82983532170362e-16, 
2.68249216476036e-17, 1.28683292112121e-18, 5.30231618313197e-20, 
1.86499767513026e-21, 5.56102696165912e-23, 1.39484152606877e-24, 
2.91735007262916e-26, 5.03779116621305e-28, 7.10181222638506e-30, 
8.06743427870919e-32, 7.2745725968875e-34, 5.1162326043855e-36, 
2.74878488435732e-38, 1.10047068271418e-40, 3.18521787783605e-43, 
6.42072520534922e-46, 8.59756395482676e-49, 7.1915294634638e-52, 
3.45947793647628e-55, 8.51888308176039e-59, 9.01922230369063e-63, 
3.08302899000303e-67, 1.97286057487992e-72, 5.90806786503182e-79
)




bell <- function(n) {





bell218 <- c(
  1.000000000000000e+00,  1.000000000000000e+00,  2.000000000000000e+00,
  5.000000000000000e+00,  1.500000000000000e+01,  5.200000000000000e+01,
  2.030000000000000e+02,  8.770000000000000e+02,  4.140000000000000e+03,
  2.114700000000000e+04,  1.159750000000000e+05,  6.785700000000000e+05,
  4.213597000000000e+06,  2.764443700000000e+07,  1.908993220000000e+08,
  1.382958545000000e+09,  1.048014214700000e+10,  8.286486980400000e+10,
  6.820768061590000e+11,  5.832742205057000e+12,  5.172415823537200e+13,
  4.748698161567510e+14,  4.506715738447323e+15,  4.415200585508434e+16,
  4.459588692948053e+17,  4.638590332230000e+18,  4.963124652361875e+19,
  5.457170479360600e+20,  6.160539404599935e+21,  7.133980193886027e+22,
  8.467490145118094e+23,  1.029335894622638e+25,  1.280646700499087e+26,
  1.629595892846008e+27,  2.119503938864036e+28,  2.816002030195603e+29,
  3.819714729894818e+30,  5.286836620855045e+31,  7.462898920956253e+32,
  1.073882333077469e+34,  1.574505883912049e+35,  2.351152507740618e+36,
  3.574254919887262e+37,  5.529501187971655e+38,  8.701963427387055e+39,
  1.392585052662637e+41,  2.265418219334494e+42,  3.745005950246151e+43,
  6.289197963031185e+44,  1.072613715457336e+46,  1.857242687710783e+47,
  3.263983870004112e+48,  5.820533802419587e+49,  1.052928518014714e+51,
  1.931728758914562e+52,  3.593340859686228e+53,  6.775685320645825e+54,
  1.294826619475070e+56,  2.507136358984296e+57,  4.917674333630963e+58,
  9.769393074670076e+59,  1.965236447154794e+61,  4.002373048214548e+62,
  8.250771700405626e+63,  1.721341433573589e+65,  3.633778785457900e+66,
  7.760590723884368e+67,  1.676501284301524e+69,  3.662822420669614e+70,
  8.092127683879479e+71,  1.807500389834051e+73,  4.081300934104643e+74,
  9.314528182092654e+75,  2.148346235684789e+77,  5.006908024247926e+78,
  1.178960269208583e+80,  2.804379077740745e+81,  6.737944959525486e+82,
  1.635000770532737e+84,  4.006416684408436e+85,  9.912679888084249e+86,
  2.476128871846587e+88,  6.243874544294799e+89,  1.589229281329695e+91,
  4.082481412918058e+92,  1.058332187322824e+94,  2.768444430541609e+95,
  7.306720755827531e+96,  1.945538974039657e+98,  5.225728505358478e+99,
 1.415803181233929e+101, 3.868731362280703e+102, 1.066117978927398e+104,
 2.962614388531219e+105, 8.301204355096730e+106, 2.345129936856330e+108,
 6.679085342279742e+109, 1.917593350464113e+111, 5.549467792774635e+112,
 1.618706027446069e+114, 4.758539127676484e+115, 1.409730628836818e+117,
 4.208466654083319e+118, 1.265919065795175e+120, 3.836647504164687e+121,
 1.171472088078324e+123, 3.603435930172301e+124, 1.116548875515524e+126,
 3.484869565157084e+127, 1.095507758559136e+129, 3.468464920150717e+130,
 1.105924120760088e+132, 3.551021092739916e+133, 1.148141058308286e+135,
 3.737885009905923e+136, 1.225237748576812e+138, 4.043468020910462e+139,
 1.343391970369164e+141, 4.493066567595946e+142, 1.512693687409340e+144,
 5.126302378999328e+145, 1.748557611821212e+147, 6.002828727345698e+148,
 2.074007725804645e+150, 7.211434692028923e+151, 2.523298790135570e+153,
 8.884449525134851e+154, 3.147651706923194e+156, 1.122062101168634e+158,
 4.024401240170187e+159, 1.452179133414943e+161, 5.271744094072983e+162,
 1.925238851994869e+164, 7.072815377340572e+165, 2.613719303956171e+167,
 9.715524863383836e+168, 3.632422266133662e+170, 1.365939384725138e+172,
 5.165996528087503e+173, 1.964930390893073e+175, 7.516118869070622e+176,
 2.891191035313989e+178, 1.118356108235045e+180, 4.349980802230563e+181,
 1.701305914243053e+183, 6.690360034164158e+184, 2.645287905991050e+186,
 1.051568003439830e+188, 4.202691619499678e+189, 1.688606328843198e+191,
 6.820641270431348e+192, 2.769512201318386e+194, 1.130441765592415e+196,
 4.638159591084169e+197, 1.912853398932712e+199, 7.929435924143752e+200,
 3.303800276765733e+202, 1.383510971580284e+204, 5.822846957871418e+205,
 2.462968564095297e+207, 1.046983958386285e+209, 4.472660677663019e+210,
 1.920100208581896e+212, 8.283259472755864e+213, 3.590753745804391e+215,
 1.564100551321700e+217, 6.845832306150363e+218, 3.010636875670279e+220,
 1.330298684661114e+222, 5.905910761818993e+223, 2.634267885372520e+225,
 1.180475203722732e+227, 5.314548104307013e+228, 2.403682782369313e+230,
 1.092139766309350e+232, 4.984924571688584e+233, 2.285637911063909e+235,
 1.052722792076366e+237, 4.870435313406448e+238, 2.263383846911113e+240,
 1.056513271054947e+242, 4.953449716133311e+243, 2.332633314212639e+245,
 1.103268112653746e+247, 5.240848508276247e+248, 2.500335400094731e+250,
 1.198012461909184e+252, 5.764759380024167e+253, 2.785789501038456e+255,
 1.351927034914996e+257, 6.588502688121199e+258, 3.224330064591281e+260,
 1.584537108486101e+262, 7.819274340696755e+263, 3.874562322516795e+265,
 1.927800791834669e+267, 9.631106570740573e+268, 4.831211697148529e+270,
 2.433286203557306e+272, 1.230492725107423e+274, 6.247484776193702e+275,
 3.184661755035448e+277, 1.629840477536268e+279, 8.374181674001929e+280,
 4.319634877400987e+282, 2.236923067904868e+284, 1.162910920454423e+286,
 6.069114899917592e+287, 3.179654797737542e+289, 1.672255197855375e+291,
 8.828469928638078e+292, 4.678655241500936e+294, 2.488867863122317e+296,
 1.328985911975485e+298, 7.123103780254685e+299, 3.832138592196165e+301,
 2.069326088380941e+303, 1.121567107165201e+305, 6.101309833875322e+306)





  maxn <- max(n[is.finite(n)], na.rm = TRUE)
  lbell218 <- length(bell218)
  ans <- numeric(length(n))
  nok <- !is.na(n) & n == round(n) & n >= 0 & is.finite(n)
  ans[ nok] <- bell218[n[nok] + 1]
  ans[!nok] <- NaN
  ans[is.na(n)] <- NA
  ans[!is.na(n) & n >= lbell218 & n == round(n)] <- Inf
  ans
}










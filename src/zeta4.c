/* 20091207   Thomas W Yee                                                  */
/* This is my attempt to convert zeta.pas into zeta.r                       */
/* It seems that zeta.pas implements 2 algorithms; Ive chosen               */
/* the latest algorithm.                                                    */
/*                                                                          */
/* 20020802; Ive added a "8" to the function names: just in case            */ 
/*                                                                          */
/* Adapted from zeta.c                                                      */
/* Adapted from Garry Tees Pascal program Gram: # 1999-8-11                 */
/* It all works, since Ive compared the results to Garrys for               */
/* s=1.1,...,4.9                                                            */
/*                                                                          */
/* Evaluate zeta(s) by Grams method, 1903                                   */
/*                                                                          */
/*       gn = 20;    # Stieltjes coefficients gg(0..gn)  of length(gn+1)    */
/*                   # Stieltjes for zeta(s) O.K. for s<4                   */
/*                                                                          */
/* Note:                                                                    */
/* Subroutine vstcoef(gg)                                                   */
/* Subroutine vzsdzs(s, zs, dzs, gg)                                        */
/*  commented out (not needed)                                              */

/* ------------------------------------------------------------------------ */

#include <math.h>
#include <Rmath.h>
#include <R.h>

/* Function prototype                                                       */
void aaaa_vzetawr(double x[], double *ans, long *deriv, long *nn);
double aaaa_zeta8(double s, double B2[]);
double aaaa_dzeta8(double s, double B2[]);
double Upsilon8(double s, double B2[]);
double aaaa_ddzeta8(double s, double B2[]);
double duds(double s, double B2[]);
void vbecoef(double B2[]);

/* ------------------------------------------------------------------------ */

void aaaa_vzetawr(double x[], double *ans, long *deriv, long *nn) {

/* Wrapper function for the zeta function and first 2 derivs                */

/* double precision x(nn), ans(nn),                                         */
/* double           x[nn-1], ans[nn-1]                                      */
  
  int    ilocal;

/* Bernoulli numbers B(2k)                                                  */
  double B2[12];
 
  vbecoef(B2);
  
/* Generate the Bernoulli coefficients into B2                              */
  
  if (*deriv == 0) {
      for (ilocal = 0; ilocal < *nn; ilocal++) {
          ans[ilocal] = aaaa_zeta8(x[ilocal], B2);
      }
  }
  else if (*deriv == 1) {
      for (ilocal = 0; ilocal < *nn; ilocal++) {
          ans[ilocal] = aaaa_dzeta8(x[ilocal], B2);
      }
  }
  else if (*deriv == 2) {
      for (ilocal = 0; ilocal < *nn; ilocal++) {
          ans[ilocal] = aaaa_ddzeta8(x[ilocal], B2);
      }
  }
  else {
      Rprintf("Error: *deriv must equal 0, 1 or 2\n");  
  }
  
}

/* ------------------------------------------------------------------------ */

double aaaa_zeta8(double s, double B2[]) {

/* Grams method for zeta(s), s>0,s<>1                                       */
/* cf Bengt Markman, BIT 5 (1965), 138-141                                  */

/* Bernoulli numbers B(2k) assumed set                                      */
  
  int    alocal, klocal;

/* Local var                                                                */
  int    mlocal, nlocal, m2local;
  double sum, plocal, a2local, fred;

  alocal = 12;
/* For 62 significant binary digits in zeta(s)                              */
  klocal = 8;       

  a2local = pow(alocal, 2.0);
  plocal = s / 2.000 / a2local;
  sum = 1.000 / (s - 1.000) + 0.500 / alocal + B2[0] * plocal;
  
  for (mlocal = 2; mlocal <= klocal; mlocal++){
    m2local = mlocal + mlocal;
    plocal *= (s + m2local - 3.000) * 
             (s + m2local - 2.000) / (m2local - 1.000) / m2local / a2local;
    sum += plocal * B2[mlocal-1];
  }
  fred = exp((s - 1.000) * log( (double) alocal)); 
  sum = 1.000 + sum / fred; 
  
  for (nlocal = 2; nlocal < alocal; nlocal++){
/* nlocal ** (-s) # power(nlocal, -s)                                       */
    sum += exp(-s * log( (double) nlocal ));
  }

  return sum;
}

/* ------------------------------------------------------------------------ */

double aaaa_dzeta8(double s, double B2[]) {

/* Grams method for d zeta(s)/ds, s>0,s<>1                                  */                       
/* cf Bengt Markman, BIT 5 (1965), 138-141                                  */ 

/* Local var                                                                */
  int    alocal, klocal;
  int    mlocal, nlocal, m2local;
  double sum, plocal, qlocal, a2local, loga, logn;
  double fred;

  alocal = 12;
/* For 62 significant binary digits in zeta8(s)                             */ 
  klocal = 8;     

  loga = log( (double) alocal );
  a2local = alocal * alocal;
  plocal = s / 2.000 / a2local;
  qlocal = 1.000 / s - loga;
  sum = B2[0] * plocal * qlocal;
  
  for (mlocal = 2; mlocal <= klocal; mlocal++){
    m2local = mlocal + mlocal;
    plocal *= (s + m2local - 3.000) * 
             (s + m2local - 2.000) / (m2local - 1.000) / m2local / a2local;
    qlocal += 1.000 / (s + m2local - 3.000) + 1.000 / (s + m2local - 2.000);
    sum += B2[mlocal-1] * plocal * qlocal;
  }
  fred = exp((1.000 - s) * loga);
/* nb. sqr in pascal is square, not square root                             */
  sum = (sum - 1.000/ pow((s - 1.000), 2.0) - loga * (1.000/(s - 1.000) +
         0.5000/alocal)) * fred;
  
  for (nlocal = 2; nlocal < alocal; nlocal++){	 
    logn = log( (double) nlocal );
    sum -= logn / exp(logn * s);
  }

  return sum;
}

/* ------------------------------------------------------------------------ */

double Upsilon8(double s, double B2[]) {
/* double precision dzeta8, zeta8;                                          */
    double Upsilon8; 
    
    Upsilon8 = -aaaa_dzeta8(s, B2) / aaaa_zeta8(s, B2);
    
    return Upsilon8;
}

/* ------------------------------------------------------------------------ */

double aaaa_ddzeta8(double s, double B2[]) {

/* Grams method for zeta"(s), s>0,s<>1                                      */
/* cf Bengt Markman, BIT 5 (1965), 138-141                                  */


/* Local var                                                                */
    int      alocal, klocal;
    int      mlocal, nlocal, m2local;
    double   sum, plocal, qlocal, rlocal, a2local, loga, logn;
    double   fred, fred2;

    alocal = 12;
/* For 62 significant binary digits in zeta8(s)                             */
    klocal = 8;     

    loga = log( (double) alocal );
    a2local = alocal * alocal;
    plocal = s / 2.000 / a2local;
    qlocal = 1.000 / s - loga;
    rlocal = 1.000 / s / s;
    sum = B2[0] * plocal * (pow(qlocal, 2.0) - rlocal);

    for (mlocal = 2; mlocal < klocal; mlocal++){ 
        m2local = mlocal + mlocal;
        plocal *= (s + m2local - 3.000) *
                (s + m2local - 2.000) / (m2local - 1.000) / m2local / a2local;
/* plocal  = s(s+1)...(s+2mlocal-2)/(2mlocal)!/a^(2mlocal)                  */
        qlocal += 1.000 / (s + m2local - 3.000) +
                  1.000 / (s + m2local - 2.000);
/* qlocal = 1/s + 1/(s+1) + ... + 1/(s+2mlocal-2)  - log alocal             */
        rlocal += 1.000 / pow((s + m2local - 3.000), 2.0) +
                  1.000 / pow((s + m2local - 2.000), 2.0);
/* rlocal =  -dq/ds = 1/s^2 + 1/(s+1)^2 + ... + 1/(s+2mlocal-2)^2           */
        sum += B2[mlocal-1] * plocal * (qlocal * qlocal - rlocal); 
    }
    fred = exp((1.000 - s) * loga);
    fred2 = pow(loga, 2.0) * (1.000/(s - 1.000) + 0.500/alocal); 
    sum = (sum + 2.000 / pow((s - 1.000), 3.0) +
          2.000 * loga / pow((s - 1.000), 2.0) + fred2) * fred;
	
    for (nlocal = 2; nlocal < alocal; nlocal++) {
        logn = log( (double) nlocal );
        sum += pow(logn, 2.0) / exp(logn * s);
    }
     
    return sum;
}

/* ------------------------------------------------------------------------ */

double duds(double s, double B2[]) {
/* d Upsilon / ds                                                           */
/* double zs, zeta8, dzeta8, ddzeta8                                        */

   double zs, duds;

   zs = aaaa_zeta8(s, B2);
   duds = pow(aaaa_dzeta8(s, B2) / zs, 2.0) - aaaa_ddzeta8(s, B2) / zs;
   
   return duds;
}

/* ------------------------------------------------------------------------ */

void vbecoef(double B2[]) {
/* Bernoulli numbers B(2k)                                                  */

/* Generate (set) the Bernoulli numbers                                     */
/* Reference: p.810 Abramowitz and Stegun                                   */

    B2[0] = 1.000 / 6.000;
    B2[1] = -1.000 / 30.000;
    B2[2] = 1.000 / 42.000;
    B2[3] = -1.000 / 30.000;
    B2[4] = 5.000 / 66.000;
    B2[5] = -691.000 / 2730.000;
    B2[6] = 7.000 / 6.000;
    B2[7] = -3617.000 / 510.000;
    B2[8] = 4386.700 / 79.800;
    B2[9] = -1746.1100 / 3.3000;
    B2[10] = 8545.1300 / 1.3800;
    B2[11] = -2363.6409100 / 0.0273000;
}

/* ------------------------------------------------------------------------ */

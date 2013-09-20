

#include <math.h>
#include <Rmath.h>
#include <R.h>

void vzetawr(double sjwyig9t[], double *bqelz3cy, int *kpzavbj3, int *f8yswcat);
double fvlmz9iyzeta8(double , double kxae8glp[]);
double fvlmz9iydzeta8(double , double kxae8glp[]);
double fvlmz9iyddzeta8(double , double kxae8glp[]);
void vbecoef(double kxae8glp[]);


void vzetawr(double sjwyig9t[], double *bqelz3cy, int *kpzavbj3, int *f8yswcat) {



  int    ayfnwr1v;
  double *qnwamo0e1, *qnwamo0e2;

  double kxae8glp[12];

  vbecoef(kxae8glp);

  qnwamo0e1 = bqelz3cy;
  qnwamo0e2 = sjwyig9t;
  if (*kpzavbj3 == 0) {
    for (ayfnwr1v = 0; ayfnwr1v < *f8yswcat; ayfnwr1v++) {
      *qnwamo0e1++ = fvlmz9iyzeta8(*qnwamo0e2++, kxae8glp);
    }
  } else
  if (*kpzavbj3 == 1) {
    for (ayfnwr1v = 0; ayfnwr1v < *f8yswcat; ayfnwr1v++) {
      *qnwamo0e1++ = fvlmz9iydzeta8(*qnwamo0e2++, kxae8glp);
    }
  } else
  if (*kpzavbj3 == 2) {
    for (ayfnwr1v = 0; ayfnwr1v < *f8yswcat; ayfnwr1v++) {
      *qnwamo0e1++ = fvlmz9iyddzeta8(*qnwamo0e2++, kxae8glp);
    }
  } else {
    Rprintf("Error: *kpzavbj3 must equal 0, 1 or 2 in C function vzetawr\n");
  }
}



double fvlmz9iyzeta8(double ghz9vuba, double kxae8glp[]) {



  int    ayfnwr1v, gp1jxzuh, uw3favmo, nsvdbx3tk, m2svdbx3tk;
  double q6zdcwxk, xvr7bonh, a2svdbx3tk, fred;

  ayfnwr1v = 12;
  gp1jxzuh = 8;

  a2svdbx3tk = pow((double) ayfnwr1v, (double) 2.0);
  xvr7bonh = ghz9vuba / 2.000 / a2svdbx3tk;
  q6zdcwxk = 1.000 / (ghz9vuba - 1.000) + 0.500 / ayfnwr1v + kxae8glp[0] * xvr7bonh;

  for (uw3favmo = 2; uw3favmo <= gp1jxzuh; uw3favmo++) {
    m2svdbx3tk = uw3favmo + uw3favmo;
    xvr7bonh *= (ghz9vuba + m2svdbx3tk - 3.000) * 
              (ghz9vuba + m2svdbx3tk - 2.000) / (m2svdbx3tk -
              1.000) / m2svdbx3tk / a2svdbx3tk;
    q6zdcwxk += xvr7bonh * kxae8glp[uw3favmo-1];
  }
  fred = pow((double) ayfnwr1v, (double) 1.0 - ghz9vuba);
  q6zdcwxk = 1.000 + q6zdcwxk * fred;

  for (nsvdbx3tk = 2; nsvdbx3tk < ayfnwr1v; nsvdbx3tk++) {
    q6zdcwxk += pow((double) nsvdbx3tk, (double) -ghz9vuba);
  }

  return q6zdcwxk;
}



double fvlmz9iydzeta8(double ghz9vuba, double kxae8glp[]) {


  int    ayfnwr1v, gp1jxzuh, uw3favmo, nsvdbx3tk, m2svdbx3tk;
  double q6zdcwxk, xvr7bonh, dh9mgvze, a2svdbx3tk, ugqvjoe5a, ugqvjoe5n, fred;

  ayfnwr1v = 12;
  gp1jxzuh = 8;

  ugqvjoe5a = log( (double) ayfnwr1v );
  a2svdbx3tk = ayfnwr1v * ayfnwr1v;
  xvr7bonh = ghz9vuba / 2.000 / a2svdbx3tk;
  dh9mgvze = 1.000 / ghz9vuba - ugqvjoe5a;
  q6zdcwxk = kxae8glp[0] * xvr7bonh * dh9mgvze;

  for (uw3favmo = 2; uw3favmo <= gp1jxzuh; uw3favmo++) {
    m2svdbx3tk = uw3favmo + uw3favmo;
    xvr7bonh *= (ghz9vuba + m2svdbx3tk - 3.0) * 
              (ghz9vuba + m2svdbx3tk - 2.0) / (m2svdbx3tk - 1.0) / m2svdbx3tk / a2svdbx3tk;
    dh9mgvze += 1.0 / (ghz9vuba + m2svdbx3tk - 3.0) + 1.0 / (ghz9vuba + m2svdbx3tk - 2.0);
    q6zdcwxk   += kxae8glp[uw3favmo-1] * xvr7bonh * dh9mgvze;
  }
  fred = pow((double) ayfnwr1v, (double) 1.0 - ghz9vuba);
  q6zdcwxk = (q6zdcwxk - 1.000 / pow(ghz9vuba - 1.000, (double) 2.0) -
         ugqvjoe5a * (1.000 /    (ghz9vuba - 1.000) + 0.5000 / ayfnwr1v)) * fred;

  for (nsvdbx3tk = 2; nsvdbx3tk < ayfnwr1v; nsvdbx3tk++) {
    ugqvjoe5n = log( (double) nsvdbx3tk );
    q6zdcwxk -= ugqvjoe5n / exp(ugqvjoe5n * ghz9vuba);
  }

  return q6zdcwxk;
}








double fvlmz9iyddzeta8(double ghz9vuba, double kxae8glp[]) {


  int      ayfnwr1v, gp1jxzuh, uw3favmo, nsvdbx3tk, m2svdbx3tk;
  double   q6zdcwxk, xvr7bonh, dh9mgvze, hpmwnav2, a2svdbx3tk, ugqvjoe5a, ugqvjoe5n,
           fred1, fred2;

  ayfnwr1v = 12;
  gp1jxzuh = 8;

  ugqvjoe5a = log( (double) ayfnwr1v );
  a2svdbx3tk = ayfnwr1v * ayfnwr1v;
  xvr7bonh = ghz9vuba / 2.000 / a2svdbx3tk;
  dh9mgvze = 1.000 / ghz9vuba - ugqvjoe5a;
  hpmwnav2 = 1.000 / ghz9vuba / ghz9vuba;
  q6zdcwxk = kxae8glp[0] * xvr7bonh * (pow(dh9mgvze, (double) 2.0) - hpmwnav2);

  for (uw3favmo = 2; uw3favmo < gp1jxzuh; uw3favmo++) {
    m2svdbx3tk = uw3favmo + uw3favmo;
    xvr7bonh *= (ghz9vuba + m2svdbx3tk - 3.000) *
              (ghz9vuba + m2svdbx3tk - 2.000) / (m2svdbx3tk -
              1.0) / m2svdbx3tk / a2svdbx3tk;
    dh9mgvze += 1.000 / (ghz9vuba + m2svdbx3tk - 3.000) +
              1.000 / (ghz9vuba + m2svdbx3tk - 2.000);
    hpmwnav2 += 1.000 / pow(ghz9vuba + m2svdbx3tk - 3.000, (double) 2.0) +
              1.000 / pow(ghz9vuba + m2svdbx3tk - 2.000, (double) 2.0);
    q6zdcwxk += kxae8glp[uw3favmo-1] * xvr7bonh * (dh9mgvze * dh9mgvze - hpmwnav2);
  }
  fred1 = pow((double) ayfnwr1v, (double) 1.0 - ghz9vuba);
  fred2 = pow(ugqvjoe5a, (double) 2.0) * (1.0 / (ghz9vuba - 1.0) + 0.50 / ayfnwr1v);
  q6zdcwxk = (q6zdcwxk + 2.0 / pow(ghz9vuba - 1.0, (double) 3.0) +
          2.0 * ugqvjoe5a / pow(ghz9vuba - 1.0, (double) 2.0) + fred2) * fred1;

  for (nsvdbx3tk = 2; nsvdbx3tk < ayfnwr1v; nsvdbx3tk++) {
    ugqvjoe5n = log( (double) nsvdbx3tk );
    q6zdcwxk += pow(ugqvjoe5n, (double) 2.0) / exp(ugqvjoe5n * ghz9vuba);
  }

  return q6zdcwxk;
}









void vbecoef(double kxae8glp[]) {



  kxae8glp[0] =  1.000 / 6.000;
  kxae8glp[1] = -1.000 / 30.000;
  kxae8glp[2] =  1.000 / 42.000;
  kxae8glp[3] = -1.000 / 30.000;
  kxae8glp[4] =  5.000 / 66.000;
  kxae8glp[5] = -691.000 / 2730.000;
  kxae8glp[6] =  7.000 / 6.000;
  kxae8glp[7] = -3617.000 / 510.000;
  kxae8glp[8] = 4386.700 / 79.800;
  kxae8glp[9] = -1746.1100 / 3.3000;
  kxae8glp[10] = 8545.1300 / 1.3800;
  kxae8glp[11] = -2363.6409100 / 0.0273000;
}




void conmax_Z(double *lamvec, double *nuvec, double *bqelz3cy,
              int    *nlength, int   *kpzavbj3,
              double *qaltf0nz) {


  double *pq6zdcwxk, denom = 0.0, yq6lorbx, prevterm;
  int    ayfnwr1v;

  *qaltf0nz = 1.0e-6;

  if (*kpzavbj3 == 0) {
    pq6zdcwxk = bqelz3cy;
    for (ayfnwr1v = 0; ayfnwr1v < *nlength; ayfnwr1v++) {
      prevterm = 1.0 + *lamvec;
      denom = 1.0;
      *pq6zdcwxk = prevterm;
      yq6lorbx = 2.0;

      if (*nuvec == 0.0 && *lamvec >= 1.0) {
        Rprintf("Error: series will not converge. Returning 0.0\n");
        *pq6zdcwxk = 0.0;
      } else {
        while (prevterm > *qaltf0nz) {
          denom = denom * pow(yq6lorbx, *lamvec);
          prevterm = prevterm * *lamvec / denom;
          *pq6zdcwxk += prevterm;
          yq6lorbx += 1.0;
        }
      }
      lamvec++;
      nuvec++;
      pq6zdcwxk++;
    }
  } else if (*kpzavbj3 == 1) {

    } else if (*kpzavbj3 == 2) {

    }
}








#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<R.h>
#include<Rmath.h>



void sf_C_expint(double *x, int *size, double *bzmd6ftv);
void sf_C_expexpint(double *x, int *size, double *bzmd6ftv);
void sf_C_expint_e1(double *x, int *size, double *bzmd6ftv);
void VGAM_C_kend_tau(double *x, double *y, int *f8yswcat, double *bqelz3cy);


void F77_NAME(einlib)(double*, double*);
void F77_NAME(expeinl)(double*, double*);
void F77_NAME(eonenl)(double*, double*);







void sf_C_expint(double *x,
                 int *size,
                 double *bzmd6ftv) {
  int ayfnwr1v;

  for (ayfnwr1v = 0; ayfnwr1v < *size; ayfnwr1v++)
    F77_NAME(einlib)(x + ayfnwr1v, bzmd6ftv + ayfnwr1v);
}



void sf_C_expexpint(double *x,
                    int *size,
                    double *bzmd6ftv) {
  int ayfnwr1v;

  for (ayfnwr1v = 0; ayfnwr1v < *size; ayfnwr1v++)
    F77_NAME(expeinl)(x + ayfnwr1v, bzmd6ftv + ayfnwr1v);
}



void sf_C_expint_e1(double *x,
                    int *size,
                    double *bzmd6ftv) {
  int ayfnwr1v;

  for (ayfnwr1v = 0; ayfnwr1v < *size; ayfnwr1v++)
    F77_NAME(eonenl)(x + ayfnwr1v, bzmd6ftv + ayfnwr1v);
}




void VGAM_C_kend_tau(double *x, double *y, int *f8yswcat, double *bqelz3cy) {


  int ayfnwr1v, yq6lorbx, gp1jxzuh = *f8yswcat    ;
  double q6zdcwxk1, q6zdcwxk2;

  for (ayfnwr1v = 0; ayfnwr1v < 3; ayfnwr1v++)
    bqelz3cy[ayfnwr1v] = 0.0;

  for (ayfnwr1v = 0; ayfnwr1v < gp1jxzuh; ayfnwr1v++) {
    for (yq6lorbx = ayfnwr1v + 1; yq6lorbx < *f8yswcat; yq6lorbx++) {
      q6zdcwxk1 = x[ayfnwr1v] - x[yq6lorbx];
      q6zdcwxk2 = y[ayfnwr1v] - y[yq6lorbx];

      if (q6zdcwxk1 == 0.0 || q6zdcwxk2 == 0.0) {
        bqelz3cy[1] += 1.0;
      } else if ((q6zdcwxk1 < 0.0 && q6zdcwxk2 < 0.0) ||
                 (q6zdcwxk1 > 0.0 && q6zdcwxk2 > 0.0)) {
        bqelz3cy[0] += 1.0;
      } else {
        bqelz3cy[2] += 1.0;
      }
    }
  }
}






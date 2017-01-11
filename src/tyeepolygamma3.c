


#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<R.h>
#include<Rmath.h>

void tyee_C_vdgam1(double *xval, double *lfu2qhid, int *dvhw1ulq);
void tyee_C_vtgam1(double *xval, double *lfu2qhid, int *dvhw1ulq);
void tyee_C_dgam1w(double sjwyig9t[], double lfu2qhid[], int *f8yswcat, int *dvhw1ulq);
void tyee_C_tgam1w(double sjwyig9t[], double lfu2qhid[], int *f8yswcat, int *dvhw1ulq);
void tyee_C_cum8sum(double ci1oyxas[], double lfu2qhid[], int *nlfu2qhid,
                    double valong[], int *ntot, int *notdvhw1ulq);
void eimpnbinomspecialp(int *interceptonly, double *nrows,
			double *ncols, double *sizevec,
			double *pnbinommat,
                        double *rowsums);


void tyee_C_vdgam1(double *xval, double *lfu2qhid, int *dvhw1ulq) {


  double wval, series, obr6tcex = 0.0, tmp1;

  *dvhw1ulq = 1;
  if (*xval <= 0.0e0) {
      *dvhw1ulq = 0;
      return;
  }

  if (*xval < 6.0e0) {
      tmp1 = *xval + 6.0e0;
      tyee_C_vdgam1(&tmp1, &obr6tcex, dvhw1ulq);
      *lfu2qhid = obr6tcex - 1.0e0 /  *xval          - 1.0e0 / (*xval + 1.0e0) -
                       1.0e0 / (*xval + 2.0e0) - 1.0e0 / (*xval + 3.0e0) -
                       1.0e0 / (*xval + 4.0e0) - 1.0e0 / (*xval + 5.0e0);
      return;
  }
  wval = 1.0e0 / (*xval * *xval);
  series = ((wval * ( -1.0e0 /   12.0e0 +
           ((wval * (  1.0e0 /  120.0e0 +
           ((wval * ( -1.0e0 /  252.0e0 +
           ((wval * (  1.0e0 /  240.0e0 +
           ((wval * ( -1.0e0 /  132.0e0 +
           ((wval * (691.0e0 /32760.0e0 +
           ((wval * ( -1.0e0 /   12.0e0 +
            (wval * 3617.0e0)/ 8160.0e0)))))))))))))))))))));
  *lfu2qhid = log(*xval) - 0.5e0 / *xval + series;
}



void tyee_C_vtgam1(double *xval, double *lfu2qhid, int *dvhw1ulq) {


  double wval, series, obr6tcex = 0.0, tmp1;
  *dvhw1ulq = 1;

  if (*xval <= 0.0e0) {
      *dvhw1ulq = 0;
      return;
  }

  if (*xval < 6.0e0) {
      tmp1 = *xval + 6.0e0;
      tyee_C_vtgam1(&tmp1, &obr6tcex, dvhw1ulq);
      *lfu2qhid = obr6tcex +
                1.0e0 / pow( (double)  *xval,          (double) 2.0) +
                1.0e0 / pow( (double) (*xval + 1.0e0), (double) 2.0) +
                1.0e0 / pow( (double) (*xval + 2.0e0), (double) 2.0) +
                1.0e0 / pow( (double) (*xval + 3.0e0), (double) 2.0) +
                1.0e0 / pow( (double) (*xval + 4.0e0), (double) 2.0) +
                1.0e0 / pow( (double) (*xval + 5.0e0), (double) 2.0);
      return;
  }
  wval = 1.0e0 / (*xval * *xval);
  series = 1.0e0 +
           (wval * (   1.0e0 /   6.0e0 +
           (wval * (  -1.0e0 /  30.0e0 +
           (wval * (   1.0e0 /  42.0e0 +
           (wval * (  -1.0e0 /  30.0e0 +
           (wval * (   5.0e0 /  66.0e0 +
           (wval * (-691.0e0 /2370.0e0 +
           (wval * (   7.0e0 /   6.0e0 -
           (wval *  3617.0e0)/ 510.0e0))))))))))))));
  *lfu2qhid = 0.5e0 * wval + series / *xval;
}



void tyee_C_dgam1w(double sjwyig9t[], double lfu2qhid[], int *f8yswcat, int *dvhw1ulq) {

  int    ayfnwr1v, okobr6tcex;
  double *qnwamo0e1, *qnwamo0e2;

  *dvhw1ulq = 1;

  qnwamo0e1 = sjwyig9t; qnwamo0e2 = lfu2qhid;
  for (ayfnwr1v = 1; ayfnwr1v <= *f8yswcat; ayfnwr1v++) {
      tyee_C_vdgam1(qnwamo0e1++, qnwamo0e2++, &okobr6tcex);
      if (okobr6tcex != 1) *dvhw1ulq = okobr6tcex;
  }
}


void tyee_C_tgam1w(double sjwyig9t[], double lfu2qhid[], int *f8yswcat, int *dvhw1ulq) {

  int    ayfnwr1v, okobr6tcex;
  double *qnwamo0e1, *qnwamo0e2;

  *dvhw1ulq = 1;

  qnwamo0e1 = sjwyig9t; qnwamo0e2 = lfu2qhid;
  for (ayfnwr1v = 1; ayfnwr1v <= *f8yswcat; ayfnwr1v++) {
      tyee_C_vtgam1(qnwamo0e1++, qnwamo0e2++, &okobr6tcex);
      if (okobr6tcex != 1) *dvhw1ulq = okobr6tcex;
  }
}



void tyee_C_cum8sum(double ci1oyxas[], double lfu2qhid[], int *nlfu2qhid,
                    double valong[], int *ntot, int *notdvhw1ulq) {


  int    ayfnwr1v, iii = 1;

  lfu2qhid[iii-1] = ci1oyxas[iii-1];
  for (ayfnwr1v = 2; ayfnwr1v <= *ntot; ayfnwr1v++) {
      if (valong[ayfnwr1v-1] > valong[ayfnwr1v-2]) {
          lfu2qhid[iii-1] += ci1oyxas[ayfnwr1v-1];
      } else {
          iii++;
          lfu2qhid[iii-1]  = ci1oyxas[ayfnwr1v-1];
      }
  }

  *notdvhw1ulq = (iii == *nlfu2qhid) ? 0 : 1;
}






void eimpnbinomspecialp(int *interceptonly,
			double *nrows,
			double *ncols,
			double *sizevec,     /* length is nrows */
			double *pnbinommat,
                        double *rowsums) {


  double ayfnwr1v, yq6lorbx, tmp1 = 0.0, tmp2;
  double *fpdlcqk9rowsums, *fpdlcqk9sizevec;


  if (*interceptonly == 1) {
    for (yq6lorbx = 0; yq6lorbx < *ncols; yq6lorbx++) {
      tmp2 = (*sizevec + yq6lorbx);
      tmp1 += *pnbinommat++ / (tmp2 * tmp2);
    }
    *rowsums = tmp1;
    return;
  }



  fpdlcqk9rowsums = rowsums;
  for (ayfnwr1v = 0; ayfnwr1v < *nrows; ayfnwr1v++)
    *fpdlcqk9rowsums++ = 0.0;

  for (yq6lorbx = 0; yq6lorbx < *ncols; yq6lorbx++) {
    fpdlcqk9rowsums = rowsums;
    fpdlcqk9sizevec = sizevec;
    for (ayfnwr1v = 0; ayfnwr1v < *nrows; ayfnwr1v++) {
      tmp2 = (yq6lorbx + *fpdlcqk9sizevec++);
      tmp1 = *pnbinommat++ / (tmp2 * tmp2);
      *fpdlcqk9rowsums++ += tmp1;
    }
  }
}






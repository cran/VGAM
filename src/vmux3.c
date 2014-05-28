

#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<R.h>
#include<Rmath.h>

void fvlmz9iyC_qpsedg8x(int tgiyxdw1[], int dufozmt7[], int *wy1vqfzu);
int  fvlmz9iyC_VIAM(int *cz8qdfyj, int *rvy1fpli, int *wy1vqfzu);
void fvlmz9iyC_vm2a(double mtlgduey8[], double bzmd6ftvmat[], int *dim1m, int *f8yswcat,
         int *wy1vqfzu, int *irb1onzwu, int tgiyxdw1[], int dufozmt7[], int *oey3ckps);
void fvlmz9iyC_mux22(double wpuarq2m[], double tlgduey8[], double bzmd6ftvmat[],
                  int *npjlv3mr, int *f8yswcat, int *wy1vqfzu);
void fvlmz9iyC_vbks(double wpuarq2m[], double unvxka0m[],
                 int *wy1vqfzu, int *f8yswcat, int *dimu);
void fvlmz9iyjdbomp0g(double rbne6ouj[], double unvxka0m[],
                  int *wy1vqfzu, int *dvhw1ulq, int *i_solve);
void fvlmz9iyC_mux17(double wpuarq2m[], double he7mqnvy[],
                  int *wy1vqfzu, int *xjc4ywlh, int *f8yswcat, int *dimu, int *rutyk8mg);
void fvlmz9iyC_lkhnw9yq(double wpuarq2m[], double ks3wejcv[],
                   int *npjlv3mr, int *wy1vqfzu, int *dvhw1ulq);
double fvlmz9iyC_tldz5ion(double xx);
void fvlmz9iyC_enbin9(double bzmd6ftv[], double hdqsx7bk[], double nm0eljqk[],
                   double *n2kersmx, int *f8yswcat, int *dvhw1ulq, int *zy1mchbf,
                   double *ux3nadiw, double *rsynp1go, int *sguwj9ty);
void fvlmz9iyC_enbin8(double bzmd6ftv[], double hdqsx7bk[], double hsj9bzaq[],
                   double *n2kersmx, int *f8yswcat, int *dvhw1ulq, int *zy1mchbf,
                   double *ux3nadiw, double *rsynp1go);
void fvlmz9iyC_mbessI0(double unvxka0m[], int *f8yswcat, int *kpzavbj3,
                    double dvector0[], double dvector1[], double dvector2[],
                    int *zjkrtol8, double *qaltf0nz);
void VGAM_C_mux34(double he7mqnvy[], double Dmat[], int *vnc1izfy, int *e0nmabdk,
                  int *ui4ntmvd, double bqelz3cy[]);


void fvlmz9iyC_qpsedg8x(int tgiyxdw1[], int dufozmt7[], int *wy1vqfzu) {

  int  urohxe6t, bpvaqm5z, *ptri;

  ptri = tgiyxdw1;
  for (urohxe6t = *wy1vqfzu; urohxe6t >= 1; urohxe6t--) {
      for (bpvaqm5z = 1; bpvaqm5z <= urohxe6t; bpvaqm5z++) {
          *ptri++ = bpvaqm5z;
      }
  }

  ptri = dufozmt7;
  for (urohxe6t = 1; urohxe6t <= *wy1vqfzu; urohxe6t++) {
      for (bpvaqm5z = urohxe6t; bpvaqm5z <= *wy1vqfzu; bpvaqm5z++) {
          *ptri++ = bpvaqm5z;
      }
  }

}


int fvlmz9iyC_VIAM(int *cz8qdfyj, int *rvy1fpli, int *wy1vqfzu) {


  int urohxe6t;

  int    *wkumc9idtgiyxdw1, *wkumc9iddufozmt7;
  int    imk5wjxg = *wy1vqfzu * (*wy1vqfzu + 1) / 2;
  wkumc9idtgiyxdw1  = Calloc(imk5wjxg, int);
  wkumc9iddufozmt7  = Calloc(imk5wjxg, int);
  fvlmz9iyC_qpsedg8x(wkumc9idtgiyxdw1, wkumc9iddufozmt7, wy1vqfzu);

  for (urohxe6t = 1; urohxe6t <= imk5wjxg; urohxe6t++) {
      if ((wkumc9idtgiyxdw1[urohxe6t-1]== *cz8qdfyj && wkumc9iddufozmt7[urohxe6t-1] == *rvy1fpli) ||
          (wkumc9idtgiyxdw1[urohxe6t-1]== *rvy1fpli && wkumc9iddufozmt7[urohxe6t-1] == *cz8qdfyj)) {
          Free(wkumc9idtgiyxdw1);    Free(wkumc9iddufozmt7);
          return urohxe6t;
      }
  }

  Free(wkumc9idtgiyxdw1);    Free(wkumc9iddufozmt7);
  return 0;
}


void fvlmz9iyC_vm2a(double mtlgduey8[], double bzmd6ftvmat[], int *dim1m, int *f8yswcat,
         int *wy1vqfzu, int *irb1onzwu, int tgiyxdw1[], int dufozmt7[], int *oey3ckps) {



  int      ayfnwr1v, yq6lorbx, gp1jxzuh, urohxe6t;
  int      bpvaqm5z, usvdbx3tk, i_size_bzmd6ftvmat, imk5wjxg = *wy1vqfzu * (*wy1vqfzu + 1) / 2,
           zyojx5hw   =  *wy1vqfzu * *wy1vqfzu;
  double   *qnwamo0e;

  if (*oey3ckps == 1) {
      if (*irb1onzwu == 1 || *dim1m != imk5wjxg) {
          i_size_bzmd6ftvmat = zyojx5hw * *f8yswcat;
          qnwamo0e = bzmd6ftvmat;
          for (ayfnwr1v = 0; ayfnwr1v < i_size_bzmd6ftvmat; ayfnwr1v++) {
              *qnwamo0e++ = 0.0e0;
          }
      }
  }

  if (irb1onzwu == 0) {
      for (ayfnwr1v = 1; ayfnwr1v <= *f8yswcat; ayfnwr1v++) {
          urohxe6t = (ayfnwr1v-1) *  zyojx5hw;
          for (yq6lorbx = 1; yq6lorbx <= *dim1m; yq6lorbx++) {
               bpvaqm5z =   tgiyxdw1[yq6lorbx-1] - 1  +
                         (dufozmt7[yq6lorbx-1] - 1) * *wy1vqfzu + urohxe6t;
               usvdbx3tk =   dufozmt7[yq6lorbx-1] - 1  +
                         (tgiyxdw1[yq6lorbx-1] - 1) * *wy1vqfzu + urohxe6t;
               gp1jxzuh = (yq6lorbx-1) + (ayfnwr1v-1) * *dim1m;
               bzmd6ftvmat[usvdbx3tk] =
               bzmd6ftvmat[bpvaqm5z] = mtlgduey8[gp1jxzuh];
          }
      }
  } else {
      for (ayfnwr1v = 1; ayfnwr1v <= *f8yswcat; ayfnwr1v++) {
          urohxe6t = (ayfnwr1v-1) *  zyojx5hw;
          for (yq6lorbx = 1; yq6lorbx <= *dim1m; yq6lorbx++) {
               bpvaqm5z =  tgiyxdw1[yq6lorbx-1] - 1  +
                         (dufozmt7[yq6lorbx-1] - 1) * *wy1vqfzu + urohxe6t;
               gp1jxzuh = (ayfnwr1v-1) * *dim1m + (yq6lorbx-1);
               bzmd6ftvmat[bpvaqm5z] = mtlgduey8[gp1jxzuh];
          }
      }
  }
}


void fvlmz9iyC_mux22(double wpuarq2m[], double tlgduey8[], double bzmd6ftvmat[],
                  int *npjlv3mr, int *f8yswcat, int *wy1vqfzu) {



  int    ayfnwr1v, yq6lorbx, bpvaqm5z, pqneb2ra = 1, djaq7ckz = 1, oey3ckps = 0;
  int    zyojx5hw  = *wy1vqfzu *  *wy1vqfzu,
         imk5wjxg = *wy1vqfzu * (*wy1vqfzu + 1) / 2;
  int    *wkumc9idtgiyxdw1, *wkumc9iddufozmt7;
  double q6zdcwxk;

  double *wkumc9idwk12;
  wkumc9idwk12 = Calloc(zyojx5hw, double);

  wkumc9idtgiyxdw1  = Calloc(imk5wjxg, int);
  wkumc9iddufozmt7  = Calloc(imk5wjxg, int);
  fvlmz9iyC_qpsedg8x(wkumc9idtgiyxdw1, wkumc9iddufozmt7, wy1vqfzu);

  for (ayfnwr1v = 1; ayfnwr1v <= *f8yswcat; ayfnwr1v++) {
      fvlmz9iyC_vm2a(wpuarq2m + (ayfnwr1v - 1) * *npjlv3mr, wkumc9idwk12, npjlv3mr, &pqneb2ra,
                  wy1vqfzu, &djaq7ckz, wkumc9idtgiyxdw1, wkumc9iddufozmt7, &oey3ckps);

      for (yq6lorbx = 1; yq6lorbx <= *wy1vqfzu; yq6lorbx++) {
          q6zdcwxk = 0.0e0;
          for (bpvaqm5z = yq6lorbx; bpvaqm5z <= *wy1vqfzu; bpvaqm5z++) {
            q6zdcwxk += wkumc9idwk12[yq6lorbx-1 + (bpvaqm5z-1) * *wy1vqfzu] *
                         tlgduey8[ayfnwr1v-1 + (bpvaqm5z-1) * *f8yswcat];
          }
          bzmd6ftvmat[yq6lorbx-1 + (ayfnwr1v-1) * *wy1vqfzu] = q6zdcwxk;
      }
  }

  Free(wkumc9idwk12);
  Free(wkumc9idtgiyxdw1);    Free(wkumc9iddufozmt7);
}


void fvlmz9iyC_vbks(double wpuarq2m[], double unvxka0m[],
                 int *wy1vqfzu, int *f8yswcat, int *npjlv3mr) {





  int    ayfnwr1v, yq6lorbx, gp1jxzuh, pqneb2ra = 1, djaq7ckz = 1, oey3ckps = 0,
         zyojx5hw = *wy1vqfzu * *wy1vqfzu,
         imk5wjxg = *wy1vqfzu * (*wy1vqfzu + 1) / 2;
  int    *wkumc9idtgiyxdw1, *wkumc9iddufozmt7;
  double q6zdcwxk;

  double *wkumc9idwk12;
  wkumc9idwk12 = Calloc(zyojx5hw , double);

  wkumc9idtgiyxdw1  = Calloc(imk5wjxg, int);
  wkumc9iddufozmt7  = Calloc(imk5wjxg, int);
  fvlmz9iyC_qpsedg8x(wkumc9idtgiyxdw1, wkumc9iddufozmt7, wy1vqfzu);

  for (ayfnwr1v = 1; ayfnwr1v <= *f8yswcat; ayfnwr1v++) {
      fvlmz9iyC_vm2a(wpuarq2m + (ayfnwr1v - 1) * *npjlv3mr, wkumc9idwk12, npjlv3mr, &pqneb2ra,
                  wy1vqfzu, &djaq7ckz, wkumc9idtgiyxdw1, wkumc9iddufozmt7, &oey3ckps);

      for (yq6lorbx = *wy1vqfzu; yq6lorbx >= 1; yq6lorbx--) {
          q6zdcwxk = unvxka0m[yq6lorbx-1 + (ayfnwr1v-1) * *wy1vqfzu];
          for (gp1jxzuh = yq6lorbx+1; gp1jxzuh <= *wy1vqfzu; gp1jxzuh++) {
              q6zdcwxk -= wkumc9idwk12[yq6lorbx-1 + (gp1jxzuh-1) * *wy1vqfzu] *
                        unvxka0m[gp1jxzuh-1 + (ayfnwr1v-1) * *wy1vqfzu];
          }
                   unvxka0m[yq6lorbx-1 + (ayfnwr1v-1) * *wy1vqfzu] =
          q6zdcwxk / wkumc9idwk12[yq6lorbx-1 + (yq6lorbx-1) * *wy1vqfzu];
      }
  }

  Free(wkumc9idwk12);
  Free(wkumc9idtgiyxdw1);   Free(wkumc9iddufozmt7);
}


void fvlmz9iyjdbomp0g(double rbne6ouj[], double unvxka0m[],
                  int *wy1vqfzu, int *dvhw1ulq, int *i_solve) {




  double q6zdcwxk;
  int    ayfnwr1v, yq6lorbx, gp1jxzuh;

  *dvhw1ulq = 1;

  for (ayfnwr1v = 1; ayfnwr1v <= *wy1vqfzu; ayfnwr1v++) {
      q6zdcwxk = 0.0e0;
      for (gp1jxzuh = 1; gp1jxzuh <= ayfnwr1v-1; gp1jxzuh++) {
          q6zdcwxk += pow(rbne6ouj[gp1jxzuh-1 + (ayfnwr1v-1) * *wy1vqfzu], (double) 2.0);
      }
      rbne6ouj[ayfnwr1v-1 + (ayfnwr1v-1) * *wy1vqfzu] -= q6zdcwxk;

      if (rbne6ouj[ayfnwr1v-1 + (ayfnwr1v-1) * *wy1vqfzu] <= 0.0e0) {
          Rprintf("Error in fvlmz9iyjdbomp0g: not pos-def.\n");
          *dvhw1ulq = 0;
          return;
      }
           rbne6ouj[ayfnwr1v-1 + (ayfnwr1v-1) * *wy1vqfzu] =
      sqrt(rbne6ouj[ayfnwr1v-1 + (ayfnwr1v-1) * *wy1vqfzu]);

      for (yq6lorbx = ayfnwr1v+1; yq6lorbx <= *wy1vqfzu; yq6lorbx++) {
          q6zdcwxk = 0.0e0;
          for (gp1jxzuh = 1; gp1jxzuh <= ayfnwr1v-1; gp1jxzuh++) {
              q6zdcwxk += rbne6ouj[gp1jxzuh-1 + (ayfnwr1v-1) * *wy1vqfzu] *
                      rbne6ouj[gp1jxzuh-1 + (yq6lorbx-1) * *wy1vqfzu];
          }
                  rbne6ouj[ayfnwr1v-1 + (yq6lorbx-1) * *wy1vqfzu] =
                 (rbne6ouj[ayfnwr1v-1 + (yq6lorbx-1) * *wy1vqfzu] -
          q6zdcwxk) / rbne6ouj[ayfnwr1v-1 + (ayfnwr1v-1) * *wy1vqfzu];
      }
  }

  if (*i_solve == 0) {
      for (ayfnwr1v = 2; ayfnwr1v <= *wy1vqfzu; ayfnwr1v++) {
          for (yq6lorbx = 1; yq6lorbx <= ayfnwr1v-1; yq6lorbx++) {
              rbne6ouj[ayfnwr1v-1 + (yq6lorbx-1) * *wy1vqfzu] = 0.0e0;
          }
          return;
      }
  }

  for (yq6lorbx = 1; yq6lorbx <= *wy1vqfzu; yq6lorbx++) {
      q6zdcwxk = unvxka0m[yq6lorbx-1];
      for (gp1jxzuh = 1; gp1jxzuh <= yq6lorbx-1; gp1jxzuh++) {
          q6zdcwxk -= rbne6ouj[gp1jxzuh-1 + (yq6lorbx-1) * *wy1vqfzu] * unvxka0m[gp1jxzuh-1];
      }
      unvxka0m[yq6lorbx-1] = q6zdcwxk / rbne6ouj[yq6lorbx-1 + (yq6lorbx-1) * *wy1vqfzu];
  }

  for(yq6lorbx = *wy1vqfzu; yq6lorbx >= 1; yq6lorbx--) {
      q6zdcwxk = unvxka0m[yq6lorbx-1];
      for(gp1jxzuh = yq6lorbx+1; gp1jxzuh <= *wy1vqfzu; gp1jxzuh++) {
          q6zdcwxk -= rbne6ouj[yq6lorbx-1 + (gp1jxzuh-1) * *wy1vqfzu] * unvxka0m[gp1jxzuh-1];
      }
      unvxka0m[yq6lorbx-1] = q6zdcwxk / rbne6ouj[yq6lorbx-1 + (yq6lorbx-1) * *wy1vqfzu];
  }
}


void fvlmz9iyC_mux17(double wpuarq2m[], double he7mqnvy[],
                  int *wy1vqfzu, int *xjc4ywlh, int *f8yswcat,
                  int *npjlv3mr, int *rutyk8mg) {




  double q6zdcwxk;
  int    ayfnwr1v, yq6lorbx, gp1jxzuh, bpvaqm5z;

  double *wkumc9idwk12, *wkumc9idwk34;
  int    *wkumc9idtgiyxdw1, *wkumc9iddufozmt7,
         imk5wjxg = *wy1vqfzu * (*wy1vqfzu + 1) / 2,
         zyojx5hw  = *wy1vqfzu *  *wy1vqfzu,
         dz1lbtph  = *wy1vqfzu *  *xjc4ywlh;
  wkumc9idtgiyxdw1  = Calloc(imk5wjxg, int);
  wkumc9iddufozmt7  = Calloc(imk5wjxg, int);
  fvlmz9iyC_qpsedg8x(wkumc9idtgiyxdw1, wkumc9iddufozmt7, wy1vqfzu);

  wkumc9idwk12  = Calloc(zyojx5hw, double);
  wkumc9idwk34  = Calloc(dz1lbtph, double);

  for (ayfnwr1v = 1; ayfnwr1v <= *f8yswcat; ayfnwr1v++) {
      for (bpvaqm5z = 1; bpvaqm5z <= *npjlv3mr; bpvaqm5z++) {
          yq6lorbx =  wkumc9idtgiyxdw1[bpvaqm5z-1] - 1  +
                   (wkumc9iddufozmt7[bpvaqm5z-1] - 1) * *wy1vqfzu;
          wkumc9idwk12[yq6lorbx] = wpuarq2m[bpvaqm5z-1 + (ayfnwr1v-1) * *npjlv3mr];
      }

      for (gp1jxzuh = 1; gp1jxzuh <= *xjc4ywlh; gp1jxzuh++) {
          for (yq6lorbx = 1; yq6lorbx <= *wy1vqfzu; yq6lorbx++) {
                             wkumc9idwk34[yq6lorbx-1 + (gp1jxzuh-1) * *wy1vqfzu] =
              he7mqnvy[(ayfnwr1v-1) * *wy1vqfzu + yq6lorbx-1 + (gp1jxzuh-1) * *rutyk8mg];
          }
      }

    for (gp1jxzuh = 1; gp1jxzuh <= *xjc4ywlh; gp1jxzuh++) {
       for (yq6lorbx = 1; yq6lorbx <= *wy1vqfzu; yq6lorbx++) {
            q6zdcwxk = 0.0e0;
            for (bpvaqm5z = yq6lorbx; bpvaqm5z <= *wy1vqfzu; bpvaqm5z++) {
                   q6zdcwxk += wkumc9idwk12[yq6lorbx-1 + (bpvaqm5z-1) * *wy1vqfzu] *
                           wkumc9idwk34[bpvaqm5z-1 + (gp1jxzuh-1) * *wy1vqfzu];
            }
            he7mqnvy[(ayfnwr1v-1) * *wy1vqfzu + yq6lorbx-1 + (gp1jxzuh-1) * *rutyk8mg] = q6zdcwxk;
       }
    }
  }

  Free(wkumc9idwk12);
  Free(wkumc9idwk34);
  Free(wkumc9idtgiyxdw1);    Free(wkumc9iddufozmt7);
}


void fvlmz9iyC_lkhnw9yq(double wpuarq2m[], double ks3wejcv[],
                   int *npjlv3mr, int *wy1vqfzu, int *dvhw1ulq) {

 



  int     ayfnwr1v, yq6lorbx, gp1jxzuh, uaoynef0,
          zyojx5hw = *wy1vqfzu * *wy1vqfzu;
  double  q6zdcwxk, vn3iasxugno = 1.0e-14;
  double  *wkumc9idwrk;
  wkumc9idwrk = Calloc(zyojx5hw, double);

  *dvhw1ulq = 1;

  for (ayfnwr1v = 1; ayfnwr1v <= *wy1vqfzu; ayfnwr1v++) {
      for (yq6lorbx = ayfnwr1v; yq6lorbx >= 1; yq6lorbx--) {
          q6zdcwxk = (yq6lorbx == ayfnwr1v) ? 1.0e0 : 0.0e0;
          for (gp1jxzuh = yq6lorbx+1; gp1jxzuh <= ayfnwr1v; gp1jxzuh++) {
              q6zdcwxk -=     wpuarq2m[yq6lorbx-1 + (gp1jxzuh-1) * *npjlv3mr] *
                      wkumc9idwrk[gp1jxzuh-1 + (ayfnwr1v-1) * *wy1vqfzu];
          }
          if (fabs(wpuarq2m[yq6lorbx-1 + (yq6lorbx-1) * *npjlv3mr]) < vn3iasxugno) {
              Rprintf("Error in fvlmz9iyC_lkhnw9yq: U(cz8qdfyj,cz8qdfyj) is zero.\n");
              *dvhw1ulq = 0;
          } else {
                 wkumc9idwrk[yq6lorbx-1 + (ayfnwr1v-1) * *wy1vqfzu] =
              q6zdcwxk / wpuarq2m[yq6lorbx-1 + (yq6lorbx-1) * *npjlv3mr];
          }
      }
  }

  for (yq6lorbx = 1; yq6lorbx <= *wy1vqfzu; yq6lorbx++) {
      for (ayfnwr1v = yq6lorbx; ayfnwr1v <= *wy1vqfzu; ayfnwr1v++) {
          uaoynef0 = (yq6lorbx < ayfnwr1v) ? ayfnwr1v : yq6lorbx;
          q6zdcwxk = 0.0e0;
          for(gp1jxzuh = uaoynef0; gp1jxzuh <= *wy1vqfzu; gp1jxzuh++) {
              q6zdcwxk += wkumc9idwrk[yq6lorbx-1 + (gp1jxzuh-1) * *wy1vqfzu] *
                      wkumc9idwrk[ayfnwr1v-1 + (gp1jxzuh-1) * *wy1vqfzu];
          }
          ks3wejcv[yq6lorbx-1 + (ayfnwr1v-1) * *wy1vqfzu] =
          ks3wejcv[ayfnwr1v-1 + (yq6lorbx-1) * *wy1vqfzu] = q6zdcwxk;
      }
  }
  Free(wkumc9idwrk);
}


double fvlmz9iyC_tldz5ion(double xval) {

  double hofjnx2e, xd4mybgj[6], q6zdcwxk = 1.000000000190015, tmp_y = xval;
  int    yq6lorbx;

  xd4mybgj[0]=  76.18009172947146e0;
  xd4mybgj[1]= -86.50532032941677e0;
  xd4mybgj[2]=  24.01409824083091e0;
  xd4mybgj[3]=  -1.231739572450155e0;
  xd4mybgj[4]=   0.1208650973866179e-2;
  xd4mybgj[5]=  -0.5395239384953e-5;
  hofjnx2e  =  xval + 5.50;
  hofjnx2e -= (xval + 0.50) * log(hofjnx2e);
  for (yq6lorbx = 0; yq6lorbx < 6; yq6lorbx++) {
      tmp_y += 1.0e0;
      q6zdcwxk  += xd4mybgj[yq6lorbx] / tmp_y;
  }
  return -hofjnx2e + log(2.5066282746310005e0 * q6zdcwxk / xval);
}


void fvlmz9iyC_enbin9(double bzmd6ftvmat[], double hdqsx7bk[], double nm0eljqk[],
                   double *n2kersmx, int *f8yswcat, int *dvhw1ulq, int *zy1mchbf,
                   double *ux3nadiw, double *rsynp1go, int *sguwj9ty) {






  int    ayfnwr1v, kij0gwer, esql7umk;
  double vjz5sxty, pvcjl2na, mwuvskg1, btiehdm2 = 100.0e0 * *rsynp1go,
         ydb, ft3ijqmy, q6zdcwxk, plo6hkdr, csi9ydge, oxjgzv0e = 0.001e0;
  double bk3ymcih = -1.0;




  csi9ydge   = bk3ymcih;
  bk3ymcih += bk3ymcih;
  bk3ymcih += csi9ydge;






  if (*n2kersmx <= 0.80e0 || *n2kersmx >= 1.0e0) {
      Rprintf("Error in fvlmz9iyC_enbin9: bad n2kersmx value.\n");
      *dvhw1ulq = 0;
      return;
  }

  *dvhw1ulq = 1;
  for (kij0gwer = 1; kij0gwer <= *zy1mchbf; kij0gwer++) {
      for (ayfnwr1v = 1; ayfnwr1v <= *f8yswcat; ayfnwr1v++) {
          vjz5sxty =   nm0eljqk[ayfnwr1v-1 + (kij0gwer-1) * *f8yswcat]
                 / hdqsx7bk[ayfnwr1v-1 + (kij0gwer-1) * *f8yswcat];

          if ((vjz5sxty < oxjgzv0e) ||
              ( nm0eljqk[ayfnwr1v-1 + (kij0gwer-1) * *f8yswcat] > 1.0e5)) {
             bzmd6ftvmat[ayfnwr1v-1 + (kij0gwer-1) * *f8yswcat] =
              -nm0eljqk[ayfnwr1v-1 + (kij0gwer-1) * *f8yswcat] * (1.0e0 +
             hdqsx7bk[ayfnwr1v-1 + (kij0gwer-1) * *f8yswcat]
          / (hdqsx7bk[ayfnwr1v-1 + (kij0gwer-1) * *f8yswcat] +
               nm0eljqk[ayfnwr1v-1 + (kij0gwer-1) * *f8yswcat]))
       / pow(hdqsx7bk[ayfnwr1v-1 + (kij0gwer-1) * *f8yswcat], (double) 2.0);
            if (bzmd6ftvmat[ayfnwr1v-1 + (kij0gwer-1) * *f8yswcat] > -btiehdm2)
                bzmd6ftvmat[ayfnwr1v-1 + (kij0gwer-1) * *f8yswcat] = -btiehdm2;
              goto ceqzd1hi20;
          }

          q6zdcwxk = 0.0e0; 
          pvcjl2na =  hdqsx7bk[ayfnwr1v-1 + (kij0gwer-1) * *f8yswcat]
               / (hdqsx7bk[ayfnwr1v-1 + (kij0gwer-1) * *f8yswcat] +
                    nm0eljqk[ayfnwr1v-1 + (kij0gwer-1) * *f8yswcat]);
          mwuvskg1 = 1.0e0 - pvcjl2na;
          csi9ydge = hdqsx7bk[ayfnwr1v-1 + (kij0gwer-1) * *f8yswcat];
          if (pvcjl2na < btiehdm2)
            pvcjl2na = btiehdm2;
          if (mwuvskg1 < btiehdm2)
            mwuvskg1 = btiehdm2;
          esql7umk = 100 + 15 * floor(nm0eljqk[ayfnwr1v-1 + (kij0gwer-1) * *f8yswcat]);
          if (esql7umk < *sguwj9ty) {
            esql7umk = *sguwj9ty;
          }

          ft3ijqmy = pow(pvcjl2na, csi9ydge);
          *ux3nadiw = ft3ijqmy;
          plo6hkdr = (1.0e0 - *ux3nadiw) / pow(hdqsx7bk[ayfnwr1v-1 +
                   (kij0gwer-1) * *f8yswcat], (double) 2.0);
          q6zdcwxk += plo6hkdr;

          ydb = 1.0e0;
          ft3ijqmy = hdqsx7bk[ayfnwr1v-1 + (kij0gwer-1) * *f8yswcat] * mwuvskg1 * ft3ijqmy;
          *ux3nadiw += ft3ijqmy;
          plo6hkdr = (1.0e0 - *ux3nadiw) / pow((hdqsx7bk[ayfnwr1v-1 +
                    (kij0gwer-1) * *f8yswcat] + ydb), (double) 2.0);
          q6zdcwxk += plo6hkdr;

          ydb = 2.0e0;
          while (((*ux3nadiw <= *n2kersmx) || (plo6hkdr > 1.0e-4))
                && (ydb < esql7umk)) {
              ft3ijqmy = (hdqsx7bk[ayfnwr1v-1 + (kij0gwer-1) * *f8yswcat] - 1.0 + ydb) *
                       mwuvskg1 * ft3ijqmy / ydb;
              *ux3nadiw += ft3ijqmy;
              plo6hkdr =  (1.0e0 - *ux3nadiw) / pow((hdqsx7bk[ayfnwr1v-1 +
                        (kij0gwer-1) * *f8yswcat] + ydb), (double) 2.0);
              q6zdcwxk += plo6hkdr;
              ydb  += 1.0e0;
          }
          bzmd6ftvmat[ayfnwr1v-1 + (kij0gwer-1) * *f8yswcat] = -q6zdcwxk;

          ceqzd1hi20: bk3ymcih = 0.0e0;
      }
  }
}



void fvlmz9iyC_enbin8(double bzmd6ftvmat[], double hdqsx7bk[], double hsj9bzaq[],
                   double *n2kersmx, int *f8yswcat, int *dvhw1ulq, int *zy1mchbf,
                   double *ux3nadiw, double *rsynp1go) {




  int    ayfnwr1v, kij0gwer;
  double ft3ijqmy, tad5vhsu, o3jyipdf, pq0hfucn, q6zdcwxk,
         plo6hkdr, qtce8hzo1 = 0.0e0, qtce8hzo2 = 0.0e0;
  int    fw2rodat, rx8qfndg, mqudbv4y;
  double onemse, nm0eljqk, ydb, btiehdm2 = -100.0 * *rsynp1go,
         kbig = 1.0e4, oxjgzv0e = 0.0010;

  if (*n2kersmx <= 0.80e0 || *n2kersmx >= 1.0e0) {
      Rprintf("returning since n2kersmx <= 0.8 or >= 1\n");
      *dvhw1ulq = 0;
      return;
  }

  onemse = 1.0e0 / (1.0e0 + oxjgzv0e);
  *dvhw1ulq = 1;

  for (kij0gwer = 1; kij0gwer <= *zy1mchbf; kij0gwer++) {
      for (ayfnwr1v = 1; ayfnwr1v <= *f8yswcat; ayfnwr1v++) {

        if ( hdqsx7bk[ayfnwr1v-1 + (kij0gwer-1) * *f8yswcat] > kbig)
             hdqsx7bk[ayfnwr1v-1 + (kij0gwer-1) * *f8yswcat] = kbig;
        if (hsj9bzaq[ayfnwr1v-1 + (kij0gwer-1) * *f8yswcat] < oxjgzv0e)
            hsj9bzaq[ayfnwr1v-1 + (kij0gwer-1) * *f8yswcat] = oxjgzv0e;

        if (hsj9bzaq[ayfnwr1v-1 + (kij0gwer-1) * *f8yswcat] > onemse) {
            nm0eljqk =       hdqsx7bk[ayfnwr1v-1 + (kij0gwer-1) * *f8yswcat] *
               (1.0e0 / hsj9bzaq[ayfnwr1v-1 + (kij0gwer-1) * *f8yswcat] - 1.0e0);
            bzmd6ftvmat[ayfnwr1v-1 + (kij0gwer-1) * *f8yswcat] = -nm0eljqk * (1.0e0 +
            hdqsx7bk[ayfnwr1v-1 + (kij0gwer-1) * *f8yswcat] /
           (hdqsx7bk[ayfnwr1v-1 + (kij0gwer-1) * *f8yswcat] + nm0eljqk))
      / pow(hdqsx7bk[ayfnwr1v-1 + (kij0gwer-1) * *f8yswcat], (double) 2.0);
            if (bzmd6ftvmat[ayfnwr1v-1 + (kij0gwer-1) * *f8yswcat] > btiehdm2)
                bzmd6ftvmat[ayfnwr1v-1 + (kij0gwer-1) * *f8yswcat] = btiehdm2;
            goto ceqzd1hi20;
        }

        q6zdcwxk = 0.0e0; 
        fw2rodat = 1;
        rx8qfndg = hsj9bzaq[ayfnwr1v-1 + (kij0gwer-1)**f8yswcat] < (1.0 - *rsynp1go)
                   ? 1 : 0;
        mqudbv4y = fw2rodat && rx8qfndg ? 1 : 0;

        if (mqudbv4y) {
            qtce8hzo2 = hdqsx7bk[ayfnwr1v-1 + (kij0gwer-1) * *f8yswcat] *
                 log(hsj9bzaq[ayfnwr1v-1 + (kij0gwer-1) * *f8yswcat]);
            *ux3nadiw = exp(qtce8hzo2);   
        } else {
            *ux3nadiw = 0.0e0;
        }
        plo6hkdr = (1.0e0 - *ux3nadiw) / pow(hdqsx7bk[ayfnwr1v-1 +
                 (kij0gwer-1) * *f8yswcat], (double) 2.0);
        q6zdcwxk += plo6hkdr;
        o3jyipdf = fvlmz9iyC_tldz5ion(hdqsx7bk[ayfnwr1v-1 + (kij0gwer-1) * *f8yswcat]);

        ydb = 1.0e0;
        tad5vhsu = fvlmz9iyC_tldz5ion(ydb + hdqsx7bk[ayfnwr1v-1 + (kij0gwer-1) * *f8yswcat]);
        pq0hfucn = 0.0e0;
        if (mqudbv4y) {
            qtce8hzo1 = log(1.0e0 - hsj9bzaq[ayfnwr1v-1 + (kij0gwer-1) * *f8yswcat]);
            ft3ijqmy = exp(ydb * qtce8hzo1 + qtce8hzo2 + tad5vhsu - o3jyipdf - pq0hfucn);
        } else {
            ft3ijqmy = 0.0e0;
        }
        *ux3nadiw += ft3ijqmy;
        plo6hkdr = (1.0e0 - *ux3nadiw) / pow(hdqsx7bk[ayfnwr1v-1 +
                 (kij0gwer-1) * *f8yswcat] + ydb, (double) 2.0);
        q6zdcwxk  += plo6hkdr;

        ydb = 2.0e0;
        while((*ux3nadiw <= *n2kersmx) || (plo6hkdr > 1.0e-4)) {

            tad5vhsu += log(ydb + hdqsx7bk[ayfnwr1v-1+(kij0gwer-1) * *f8yswcat] - 1.0);
            pq0hfucn += log(ydb);
            if (mqudbv4y) {
              ft3ijqmy = exp(ydb * qtce8hzo1 + qtce8hzo2 + tad5vhsu - o3jyipdf - pq0hfucn);
            } else {
                ft3ijqmy = 0.0e0;
            }
            *ux3nadiw += ft3ijqmy;
            plo6hkdr = (1.0e0 - *ux3nadiw) / pow(hdqsx7bk[ayfnwr1v-1 +
                     (kij0gwer-1) * *f8yswcat] + ydb, (double) 2.0);
            q6zdcwxk += plo6hkdr;
            ydb += 1.0e0;

            if (ydb > 1.0e3) goto ceqzd1hi21;
        }

        ceqzd1hi21: bzmd6ftvmat[ayfnwr1v-1 + (kij0gwer-1) * *f8yswcat] = -q6zdcwxk;

        ceqzd1hi20: tad5vhsu = 0.0e0;
      }
  }

}





void fvlmz9iyC_mbessI0(double unvxka0m[], int *f8yswcat, int *kpzavbj3,
                    double dvector0[], double dvector1[], double dvector2[],
                    int *zjkrtol8, double *qaltf0nz) {





  int    ayfnwr1v, gp1jxzuh, c5aesxkus;
  double f0, t0, m0, f1, t1, m1, f2, t2, m2, Toobig = 20.0e0;

  *zjkrtol8 = 0;
  if (!(*kpzavbj3 == 0 || *kpzavbj3 == 1 || *kpzavbj3 == 2)) {
      Rprintf("Error in fvlmz9iyC_mbessI0: kpzavbj3 not in 0:2. Returning.\n");
      *zjkrtol8 = 1;
      return;
  }

  for (gp1jxzuh = 1; gp1jxzuh <= *f8yswcat; gp1jxzuh++) {
      if (fabs(unvxka0m[gp1jxzuh-1]) > Toobig) {
          Rprintf("Error in fvlmz9iyC_mbessI0: unvxka0m[] value > too big.\n");
          *zjkrtol8 = 1;
          return;
      }
      t1 = unvxka0m[gp1jxzuh-1] / 2.0e0;
      f1 = t1;
      t0 = t1 * t1;
      f0 = 1.0e0 + t0;
      t2 = 0.50e0;
      f2 = t2; 
      c5aesxkus = 15;
      if (fabs(unvxka0m[gp1jxzuh-1]) > 10.0) c5aesxkus = 25;
      if (fabs(unvxka0m[gp1jxzuh-1]) > 15.0) c5aesxkus = 35;
      if (fabs(unvxka0m[gp1jxzuh-1]) > 20.0) c5aesxkus = 40;
      if (fabs(unvxka0m[gp1jxzuh-1]) > 30.0) c5aesxkus = 55;

      for (ayfnwr1v = 1; ayfnwr1v <= c5aesxkus; ayfnwr1v++) {  
          m0 = pow(unvxka0m[gp1jxzuh-1] / (2.0 * (ayfnwr1v + 1.0)), (double) 2);
          m1 = m0 * (1.0e0 + 1.0e0 / ayfnwr1v);
          m2 = m1 * (2.0e0 * ayfnwr1v + 1.0e0) / (2.0e0 * ayfnwr1v - 1.0e0);
          t0 = t0 * m0;
          t1 = t1 * m1;
          t2 = t2 * m2;
          f0 = f0 + t0;
          f1 = f1 + t1;
          f2 = f2 + t2;
          if ((fabs(t0) < *qaltf0nz) &&
              (fabs(t1) < *qaltf0nz) &&
              (fabs(t2) < *qaltf0nz))
              break;
      }
      if (0 <= *kpzavbj3) dvector0[gp1jxzuh-1] = f0;
      if (1 <= *kpzavbj3) dvector1[gp1jxzuh-1] = f1;
      if (2 <= *kpzavbj3) dvector2[gp1jxzuh-1] = f2;
  }

}


void VGAM_C_mux34(double he7mqnvy[], double Dmat[], int *vnc1izfy, int *e0nmabdk,
                  int *ui4ntmvd, double bqelz3cy[]) {




  int    ayfnwr1v, yq6lorbx, gp1jxzuh;
  double *qnwamo0e1, *qnwamo0e2;

  if (*e0nmabdk == 1) {
      qnwamo0e1 = bqelz3cy;  qnwamo0e2 = he7mqnvy;
      for (ayfnwr1v = 0; ayfnwr1v < *vnc1izfy; ayfnwr1v++) {
          *qnwamo0e1++ = *Dmat * pow(*qnwamo0e2++, (double) 2.0);
      }
      return;
  }

  if (*ui4ntmvd == 1) {
      for (ayfnwr1v = 1; ayfnwr1v <= *vnc1izfy; ayfnwr1v++) {
          bqelz3cy[ayfnwr1v-1] = 0.0e0;
          for (yq6lorbx = 1; yq6lorbx <= *e0nmabdk; yq6lorbx++) {
              bqelz3cy[ayfnwr1v-1] += Dmat[yq6lorbx-1 + (yq6lorbx-1) * *e0nmabdk] *
                              pow(he7mqnvy[ayfnwr1v-1 + (yq6lorbx-1) * *vnc1izfy],
                                  (double) 2.0);
          }
          if (*e0nmabdk > 1) {
            for (yq6lorbx = 1; yq6lorbx <= *e0nmabdk; yq6lorbx++) {
              for (gp1jxzuh = yq6lorbx+1; gp1jxzuh <= *e0nmabdk; gp1jxzuh++) {
                bqelz3cy[ayfnwr1v-1] += Dmat[yq6lorbx-1 + (gp1jxzuh-1) * *e0nmabdk] *
                                    he7mqnvy[ayfnwr1v-1 + (yq6lorbx-1) * *vnc1izfy] *
                                    he7mqnvy[ayfnwr1v-1 + (gp1jxzuh-1) * *vnc1izfy] *
                                    2.0;
              }
            }
          }
      }
  } else {
      for (ayfnwr1v = 1; ayfnwr1v <= *vnc1izfy; ayfnwr1v++) {
          bqelz3cy[ayfnwr1v-1] = 0.0e0;
          for (yq6lorbx = 1; yq6lorbx <= *e0nmabdk; yq6lorbx++) {
              for (gp1jxzuh = 1; gp1jxzuh <= *e0nmabdk; gp1jxzuh++) {
                  bqelz3cy[ayfnwr1v-1] += Dmat[yq6lorbx-1 + (gp1jxzuh-1) * *e0nmabdk] *
                                      he7mqnvy[ayfnwr1v-1 + (yq6lorbx-1) * *vnc1izfy] *
                                      he7mqnvy[ayfnwr1v-1 + (gp1jxzuh-1) * *vnc1izfy];
              }
          }
      }
  }
}



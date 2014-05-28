




#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<R.h>
#include<Rmath.h>


void n5aioudkdnaoqj0l(double *qgnl3toc,
                    double sjwyig9t[], double bhcji9gl[], double po8rwsmy[],
                    int *kuzxj1lo, int *acpios9q,
                    double gkdx5jal[], double rpyis2kc[], double imdvf4hx[],
                    double ifys6woa[],
                    double *wbkq9zyi, double jstx4uwe[4],
                    double xecbg0pf[], double z4grbpiq[], double d7glzhbj[], double v2eydbxs[],
                    double *tt2,
                    int *cvnjhg2u, int l3zpbstu[3],
                    int *xtov9rbf, int *wep0oibc, int *fbd5yktj);
void n5aioudkhbzuprs6(double *qgnl3toc,
            double sjwyig9t[], double bhcji9gl[], double po8rwsmy[],
            int *kuzxj1lo, int *acpios9q, double gkdx5jal[],
            double *rpyis2kc, double *imdvf4hx, double *ifys6woa,
            double *i9mwnvqt, int *pn9eowxc, int *ic5aesxku,
            double *mynl7uaq, double *zustx4fw, double *nbe4gvpq, double *qaltf0nz,
            int *cvnjhg2u,
            double xwy[],
            double zvau2lct[], double f6lsuzax[], double fvh2rwtc[], double dcfir2no[],
            double xecbg0pf[], double z4grbpiq[], double d7glzhbj[], double v2eydbxs[],
            double *tt2,
            double buhyalv4[],
            double fulcp8wa[], double plj0trqx[],
            int *xtov9rbf, int *wep0oibc, int *fbd5yktj);
void n5aioudkzosq7hub(double xecbg0pf[], double z4grbpiq[], double d7glzhbj[], double v2eydbxs[],
                   double gkdx5jal[], int *acpios9q);
void n5aioudkvmnweiy2(double buhyalv4[], double fulcp8wa[], double plj0trqx[], int *xtov9rbf,
                    int *acpios9q, int *wep0oibc, int *iflag);
void n5aioudkwmhctl9x( double *qgnl3toc, double sjwyig9t[],
                    double po8rwsmy[], int    *kuzxj1lo, int *acpios9q,
                    int *pn9eowxc, // int *icrit,
                    double gkdx5jal[], double rpyis2kc[], double imdvf4hx[],
                    double ifys6woa[], double *i9mwnvqt, double xwy[],
                    double *qcpiaj7f,
                    double zvau2lct[], double f6lsuzax[], double fvh2rwtc[], double dcfir2no[],
                    double xecbg0pf[], double z4grbpiq[], double d7glzhbj[], double v2eydbxs[],
                    double buhyalv4[], double fulcp8wa[], double plj0trqx[],
                    int    *xtov9rbf, int *wep0oibc, int *algpft4y);
void n5aioudkgt9iulbf(double sjwyig9t[], double ghz9vuba[], double po8rwsmy[],
                   double gkdx5jal[], int *rvy1fpli, int *kuzxj1lo, double zyupcmk6[],
                   double zvau2lct[], double f6lsuzax[], double fvh2rwtc[], double dcfir2no[]);

void F77_NAME(vinterv)(double*, int*, double*, int*, int*);
void F77_NAME(vbsplvd)(double*, int*, double*, int*, double*, double*,
                       int*);
void F77_NAME(dpbfa8)(double*, int*, int*, int*, int*);
void F77_NAME(dpbsl8)(double*, int*, int*, int*, double*);
void F77_NAME(wbvalue)(double*, double*, int*, int*, double*, int*,
                       double*);


void n5aioudkdnaoqj0l(double *qgnl3toc,
                    double sjwyig9t[], double bhcji9gl[], double po8rwsmy[],
                    int *kuzxj1lo, int *acpios9q,
                    double gkdx5jal[], double rpyis2kc[], double imdvf4hx[],
                    double ifys6woa[],
                    double *wbkq9zyi, double jstx4uwe[4],
                    double xecbg0pf[], double z4grbpiq[], double d7glzhbj[], double v2eydbxs[],
                    double *tt2,
                    int *cvnjhg2u, int l3zpbstu[3],
                    int *xtov9rbf, int *wep0oibc, int *fbd5yktj) {





  double   *wkumc9idxwy,  *wkumc9idbuhyalv4,
           *wkumc9idzvau2lct,  *wkumc9idf6lsuzax,  *wkumc9idfvh2rwtc, *wkumc9iddcfir2no,
           *wkumc9idfulcp8wa, *wkumc9idplj0trqx;

  wkumc9idxwy      = Calloc(*acpios9q,            double);
  wkumc9idzvau2lct      = Calloc(*acpios9q,            double);
  wkumc9idf6lsuzax      = Calloc(*acpios9q,            double);
  wkumc9idfvh2rwtc      = Calloc(*acpios9q,            double);
  wkumc9iddcfir2no      = Calloc(*acpios9q,            double);
  wkumc9idbuhyalv4      = Calloc(*xtov9rbf  * *acpios9q,  double);
  wkumc9idfulcp8wa     = Calloc(*xtov9rbf  * *acpios9q,  double);
  wkumc9idplj0trqx     = Calloc( (int)  1      ,    double);

  n5aioudkhbzuprs6(qgnl3toc, sjwyig9t, bhcji9gl,
                po8rwsmy, kuzxj1lo, acpios9q, gkdx5jal,
                rpyis2kc, imdvf4hx, ifys6woa,
                wbkq9zyi, l3zpbstu + 1, l3zpbstu + 2,
                jstx4uwe, jstx4uwe + 1, jstx4uwe + 2, jstx4uwe + 3,
                cvnjhg2u,
                wkumc9idxwy,
                wkumc9idzvau2lct, wkumc9idf6lsuzax, wkumc9idfvh2rwtc, wkumc9iddcfir2no,
                xecbg0pf, z4grbpiq, d7glzhbj, v2eydbxs,
                tt2,
                wkumc9idbuhyalv4,
                wkumc9idfulcp8wa, wkumc9idplj0trqx,
                xtov9rbf, wep0oibc, fbd5yktj);

  Free(wkumc9idxwy);  Free(wkumc9idbuhyalv4);
  Free(wkumc9idzvau2lct);  Free(wkumc9idf6lsuzax);  Free(wkumc9idfvh2rwtc);  Free(wkumc9iddcfir2no);
  Free(wkumc9idfulcp8wa); Free(wkumc9idplj0trqx);
}




void n5aioudkhbzuprs6(double *qgnl3toc,
            double sjwyig9t[], double bhcji9gl[], double po8rwsmy[],
            int *kuzxj1lo, int *acpios9q, double gkdx5jal[],
            double *rpyis2kc, double *imdvf4hx, double *ifys6woa,
            double *wbkq9zyi, int *pn9eowxc, int *ic5aesxku,
            double *mynl7uaq, double *zustx4fw, double *nbe4gvpq, double *qaltf0nz,
            int *cvnjhg2u,
            double xwy[],
            double zvau2lct[], double f6lsuzax[], double fvh2rwtc[], double dcfir2no[],
            double xecbg0pf[], double z4grbpiq[], double d7glzhbj[], double v2eydbxs[],
            double *tt2,
            double buhyalv4[],
            double fulcp8wa[], double plj0trqx[],
            int *xtov9rbf, int *wep0oibc, int *fbd5yktj) {



  static const double c_Gold = 0.381966011250105151795413165634;
  double tt1 = 0.0, g2dnwteb,
         wkumc9ida, wkumc9idb,       wkumc9idd, wkumc9ide,
         wkumc9idxm, wkumc9idp, wkumc9idq, wkumc9idr, // qaltf0nz,
         Tol1, Tol2, wkumc9idu, wkumc9idv, wkumc9idw,
         wkumc9idfu, wkumc9idfv, wkumc9idfw, wkumc9idfx, wkumc9idx, wkumc9idax, wkumc9idbx;
  int    ayfnwr1v, viter = 0;
  double yjpnro8d = 8.0e88, bk3ymcih = 0.0e0,
         *qcpiaj7f,  qcpiaj7f0 = 0.0;


  qcpiaj7f = &qcpiaj7f0;



  g2dnwteb  = bk3ymcih;
  bk3ymcih += bk3ymcih;
  bk3ymcih *= bk3ymcih;
  bk3ymcih += g2dnwteb;































  wkumc9idd = 0.0; wkumc9idfu = 0.0e0; wkumc9idu = 0.0e0;


  if (*cvnjhg2u == 0) {
      n5aioudkzosq7hub(xecbg0pf, z4grbpiq, d7glzhbj, v2eydbxs, gkdx5jal, acpios9q);

      *tt2 = 0.0;
      for (ayfnwr1v = 3; ayfnwr1v <= (*acpios9q - 3); ayfnwr1v++) {
          *tt2 += xecbg0pf[ayfnwr1v-1];
      }

      *cvnjhg2u = 1;
  } else {
  }

  n5aioudkgt9iulbf(sjwyig9t, bhcji9gl, po8rwsmy,
                gkdx5jal, kuzxj1lo, acpios9q, xwy,
                zvau2lct, f6lsuzax, fvh2rwtc, dcfir2no);

  for (ayfnwr1v = 3; ayfnwr1v <= (*acpios9q - 3); ayfnwr1v++) {
    tt1 += zvau2lct[ayfnwr1v-1];
  }
  g2dnwteb = tt1 / *tt2;

  if (*pn9eowxc == 1) {

    *mynl7uaq = g2dnwteb * pow(16.0, *wbkq9zyi * 6.0 - 2.0);
    n5aioudkwmhctl9x(qgnl3toc, sjwyig9t,
                   po8rwsmy, kuzxj1lo, acpios9q,
                   pn9eowxc, // icrit, (icrit used to be used solely)
                   gkdx5jal, rpyis2kc, imdvf4hx,
                   ifys6woa, mynl7uaq, xwy,
                   qcpiaj7f,  // Not used here
                   zvau2lct, f6lsuzax, fvh2rwtc, dcfir2no,
                   xecbg0pf, z4grbpiq, d7glzhbj, v2eydbxs,
                   buhyalv4, fulcp8wa, plj0trqx,
                   xtov9rbf, wep0oibc, fbd5yktj);
    return;
  }



      wkumc9idax = *mynl7uaq;
      wkumc9idbx = *zustx4fw;












 /* Initialization.                                                       */
      wkumc9ida = wkumc9idax;
      wkumc9idb = wkumc9idbx;
      wkumc9idv = wkumc9ida + c_Gold * (wkumc9idb - wkumc9ida);
      wkumc9idw =
      wkumc9idx = wkumc9idv;
      wkumc9ide = 0.0e0;

      *wbkq9zyi = wkumc9idx;
      *mynl7uaq = g2dnwteb * pow((double) 16.0, (double) *wbkq9zyi * 6.0 - 2.0);
      n5aioudkwmhctl9x(qgnl3toc, sjwyig9t,
                     po8rwsmy, kuzxj1lo, acpios9q,
                     pn9eowxc, // icrit,
                     gkdx5jal, rpyis2kc, imdvf4hx,
                     ifys6woa, mynl7uaq, xwy,
                     qcpiaj7f,
                     zvau2lct, f6lsuzax, fvh2rwtc, dcfir2no,
                     xecbg0pf, z4grbpiq, d7glzhbj, v2eydbxs,
                     buhyalv4, fulcp8wa, plj0trqx,
                     xtov9rbf, wep0oibc, fbd5yktj);


      wkumc9idfx = *qcpiaj7f;
      wkumc9idfv =
      wkumc9idfw = wkumc9idfx;


      while (*fbd5yktj == 0) {
        viter++;
        wkumc9idxm = 0.5e0 * (wkumc9ida + wkumc9idb);
        Tol1 = *qaltf0nz * fabs(wkumc9idx) + *nbe4gvpq / 3.0e0;
        Tol2 = 2.0e0 * Tol1;

          if ((fabs(wkumc9idx - wkumc9idxm) <= (Tol2 - 0.5 * (wkumc9idb - wkumc9ida)))
           || (viter > *ic5aesxku))
             goto L_End;

          if ((fabs(wkumc9ide) <= Tol1)  ||
              (wkumc9idfx      >= yjpnro8d) ||
              (wkumc9idfv      >= yjpnro8d) ||
              (wkumc9idfw      >= yjpnro8d)) goto a3bdsirf;

        wkumc9idr = (wkumc9idx - wkumc9idw) * (wkumc9idfx - wkumc9idfv);
        wkumc9idq = (wkumc9idx - wkumc9idv) * (wkumc9idfx - wkumc9idfw);
        wkumc9idp = (wkumc9idx - wkumc9idv) * wkumc9idq - (wkumc9idx - wkumc9idw) * wkumc9idr;
        wkumc9idq = 2.0e0 * (wkumc9idq - wkumc9idr);
        if (wkumc9idq > 0.0e0) wkumc9idp = -wkumc9idp;
        wkumc9idq = fabs(wkumc9idq);
        wkumc9idr = wkumc9ide;
        wkumc9ide = wkumc9idd;

        if (fabs(wkumc9idp) >= fabs(0.5 * wkumc9idq * wkumc9idr) ||
            wkumc9idq == 0.0e0) {
          goto a3bdsirf;
        }
        if (wkumc9idp <= wkumc9idq * (wkumc9ida - wkumc9idx) ||
            wkumc9idp >= wkumc9idq * (wkumc9idb - wkumc9idx))
          goto a3bdsirf;

        wkumc9idd = wkumc9idp / wkumc9idq;




        wkumc9idu = wkumc9idx + wkumc9idd;

        if (wkumc9idu - wkumc9ida < Tol2 ||
            wkumc9idb - wkumc9idu < Tol2)
          wkumc9idd = fsign(Tol1, wkumc9idxm - wkumc9idx);

        goto ceqzd1hi50;

a3bdsirf:
         wkumc9ide = (wkumc9idx >= wkumc9idxm) ? wkumc9ida - wkumc9idx : wkumc9idb - wkumc9idx;
         wkumc9idd = c_Gold * wkumc9ide;

ceqzd1hi50: wkumc9idu = wkumc9idx +
                  ((fabs(wkumc9idd) >= Tol1) ? wkumc9idd : fsign(Tol1, wkumc9idd));

        *wbkq9zyi = wkumc9idu;
        *mynl7uaq = g2dnwteb * pow((double) 16.0, (double) *wbkq9zyi * 6.0 - 2.0);
        n5aioudkwmhctl9x(qgnl3toc, sjwyig9t,
                       po8rwsmy, kuzxj1lo, acpios9q,
                       pn9eowxc, // icrit,
                       gkdx5jal, rpyis2kc, imdvf4hx,
                       ifys6woa, mynl7uaq, xwy,
                       qcpiaj7f,
                       zvau2lct, f6lsuzax, fvh2rwtc, dcfir2no,
                       xecbg0pf, z4grbpiq, d7glzhbj, v2eydbxs,
                       buhyalv4, fulcp8wa, plj0trqx,
                       xtov9rbf, wep0oibc, fbd5yktj);

        wkumc9idfu = *qcpiaj7f;

        if (wkumc9idfu > yjpnro8d)
            wkumc9idfu = 2.0e0 * yjpnro8d;

        if (wkumc9idfu <= wkumc9idfx) {
          if (wkumc9idu >= wkumc9idx) wkumc9ida = wkumc9idx; else wkumc9idb = wkumc9idx;
          wkumc9idv = wkumc9idw; wkumc9idfv = wkumc9idfw;
          wkumc9idw = wkumc9idx; wkumc9idfw = wkumc9idfx;
          wkumc9idx = wkumc9idu; wkumc9idfx = wkumc9idfu;
        } else {
        if (wkumc9idu < wkumc9idx) wkumc9ida = wkumc9idu; else wkumc9idb = wkumc9idu;
          if (wkumc9idfu <= wkumc9idfw || wkumc9idw == wkumc9idx) {
            wkumc9idv = wkumc9idw; wkumc9idfv = wkumc9idfw;
            wkumc9idw = wkumc9idu; wkumc9idfw = wkumc9idfu;
          } else
            if (wkumc9idfu <= wkumc9idfv ||
                wkumc9idv  == wkumc9idx  ||
                wkumc9idv  == wkumc9idw) {
              wkumc9idv  = wkumc9idu;
              wkumc9idfv = wkumc9idfu;
            }
        }
    }
    L_End: bk3ymcih = 0.0e0;

    *wbkq9zyi = wkumc9idx;
    *qcpiaj7f = wkumc9idfx;
    return;
}







void n5aioudkzosq7hub(double xecbg0pf[], double z4grbpiq[], double d7glzhbj[], double v2eydbxs[],
                   double gkdx5jal[], int *acpios9q) {


  int    dqlr5bse, pqzfxw4i, bvsquk3z = 3, h2dpsbkr = 4, nkplus1 = *acpios9q + 1;
  int    ayfnwr1v, gp1jxzuh, yq6lorbx;
  int    urohxe6t;
  double g9fvdrbw[12], ms0qypiw[16], yw1[4], yw2[4], wrk1, othird = 1.0 / 3.0,
         *qnwamo0e0, *qnwamo0e1,  *qnwamo0e2, *qnwamo0e3;





  qnwamo0e0 = xecbg0pf; qnwamo0e1 = z4grbpiq; qnwamo0e2 = d7glzhbj; qnwamo0e3 = v2eydbxs;
  for (ayfnwr1v = 0; ayfnwr1v < *acpios9q; ayfnwr1v++) {
      *qnwamo0e0++ = *qnwamo0e1++ = *qnwamo0e2++ = *qnwamo0e3++ = 0.0e0;
  }


  for (ayfnwr1v = 1; ayfnwr1v <= *acpios9q; ayfnwr1v++) {


    F77_CALL(vinterv)(gkdx5jal, &nkplus1, gkdx5jal + ayfnwr1v-1, &dqlr5bse, &pqzfxw4i);

    F77_CALL(vbsplvd)(gkdx5jal, &h2dpsbkr, gkdx5jal + ayfnwr1v - 1, &dqlr5bse, ms0qypiw,
                      g9fvdrbw, &bvsquk3z);

    for (gp1jxzuh = 1; gp1jxzuh <= 4; gp1jxzuh++) {
      yw1[gp1jxzuh-1] = g9fvdrbw[gp1jxzuh-1 + 2*4];
    }

    F77_CALL(vbsplvd)(gkdx5jal, &h2dpsbkr, gkdx5jal + ayfnwr1v, &dqlr5bse, ms0qypiw,
                      g9fvdrbw, &bvsquk3z);

    for (gp1jxzuh = 1; gp1jxzuh <= 4; gp1jxzuh++) { 
      yw2[gp1jxzuh-1] = g9fvdrbw[gp1jxzuh-1 + 2*4] - yw1[gp1jxzuh-1];
    }
    wrk1 = gkdx5jal[ayfnwr1v] - gkdx5jal[ayfnwr1v-1];

    if (dqlr5bse >= 4) {
      for (gp1jxzuh = 1; gp1jxzuh <= 4; gp1jxzuh++) {
        yq6lorbx = gp1jxzuh;
        urohxe6t = dqlr5bse - 4 + gp1jxzuh;
        xecbg0pf[urohxe6t-1] +=
           wrk1 * (yw1[gp1jxzuh-1]*yw1[yq6lorbx-1] +
                  (yw2[gp1jxzuh-1]*yw1[yq6lorbx-1] +
                   yw2[yq6lorbx-1]*yw1[gp1jxzuh-1]) * 0.50 +
                   yw2[gp1jxzuh-1]*yw2[yq6lorbx-1]  * othird);
        yq6lorbx = gp1jxzuh + 1;
        if (yq6lorbx <= 4) {
          z4grbpiq[urohxe6t-1] +=
            wrk1 * (yw1[gp1jxzuh-1]*yw1[yq6lorbx-1] +
                   (yw2[gp1jxzuh-1]*yw1[yq6lorbx-1] +
                    yw2[yq6lorbx-1]*yw1[gp1jxzuh-1]) * 0.50  +
                    yw2[gp1jxzuh-1]*yw2[yq6lorbx-1]  * othird);
        }
        yq6lorbx = gp1jxzuh + 2;
        if (yq6lorbx <= 4) {
          d7glzhbj[urohxe6t-1] +=
          wrk1 * (yw1[gp1jxzuh-1]*yw1[yq6lorbx-1] +
                 (yw2[gp1jxzuh-1]*yw1[yq6lorbx-1] +
                  yw2[yq6lorbx-1]*yw1[gp1jxzuh-1]) * 0.50  +
                  yw2[gp1jxzuh-1]*yw2[yq6lorbx-1]  * othird);
        }
        yq6lorbx = gp1jxzuh + 3;
        if (yq6lorbx <= 4) {
          v2eydbxs[urohxe6t-1] +=
          wrk1 * (yw1[gp1jxzuh-1]*yw1[yq6lorbx-1] +
                 (yw2[gp1jxzuh-1]*yw1[yq6lorbx-1] +
                  yw2[yq6lorbx-1]*yw1[gp1jxzuh-1]) * 0.50  +
                  yw2[gp1jxzuh-1]*yw2[yq6lorbx-1]  * othird);
        }
      }
    } else if (dqlr5bse == 3) {
      for (gp1jxzuh = 1; gp1jxzuh <= 3; gp1jxzuh++) {
        yq6lorbx = gp1jxzuh;
        urohxe6t = dqlr5bse - 3 + gp1jxzuh;
        xecbg0pf[urohxe6t-1] +=
           wrk1 * (yw1[gp1jxzuh-1]*yw1[yq6lorbx-1] +
                  (yw2[gp1jxzuh-1]*yw1[yq6lorbx-1] +
                   yw2[yq6lorbx-1]*yw1[gp1jxzuh-1]) * 0.50  +
                   yw2[gp1jxzuh-1]*yw2[yq6lorbx-1]  * othird);
        yq6lorbx = gp1jxzuh + 1;
        if (yq6lorbx <= 3) {
          z4grbpiq[urohxe6t-1] +=
              wrk1 * (yw1[gp1jxzuh-1]*yw1[yq6lorbx-1] +
                     (yw2[gp1jxzuh-1]*yw1[yq6lorbx-1] +
                      yw2[yq6lorbx-1]*yw1[gp1jxzuh-1]) * 0.50  +
                      yw2[gp1jxzuh-1]*yw2[yq6lorbx-1]  * othird);
        }
        yq6lorbx = gp1jxzuh + 2;
        if (yq6lorbx <= 3) {
          d7glzhbj[urohxe6t-1] +=
             wrk1 * (yw1[gp1jxzuh-1]*yw1[yq6lorbx-1] +
                    (yw2[gp1jxzuh-1]*yw1[yq6lorbx-1] +
                     yw2[yq6lorbx-1]*yw1[gp1jxzuh-1]) * 0.50  +
                     yw2[gp1jxzuh-1]*yw2[yq6lorbx-1]  * othird);
        }
      }
    } else if (dqlr5bse == 2) {
      for (gp1jxzuh = 1; gp1jxzuh <= 2; gp1jxzuh++) {
        yq6lorbx = gp1jxzuh;
        urohxe6t = dqlr5bse - 2 + gp1jxzuh;
        xecbg0pf[urohxe6t-1] +=
          wrk1 * (yw1[gp1jxzuh-1]*yw1[yq6lorbx-1] +
                 (yw2[gp1jxzuh-1]*yw1[yq6lorbx-1] +
                  yw2[yq6lorbx-1]*yw1[gp1jxzuh-1]) * 0.50  +
                  yw2[gp1jxzuh-1]*yw2[yq6lorbx-1]  * othird);
        yq6lorbx = gp1jxzuh + 1;
        if (yq6lorbx <= 2) {
          z4grbpiq[urohxe6t-1] +=
            wrk1 * (yw1[gp1jxzuh-1]*yw1[yq6lorbx-1] +
                   (yw2[gp1jxzuh-1]*yw1[yq6lorbx-1] +
                    yw2[yq6lorbx-1]*yw1[gp1jxzuh-1]) * 0.50  +
                    yw2[gp1jxzuh-1]*yw2[yq6lorbx-1]  * othird);
        }
      }
    } else if (dqlr5bse == 1) {
      for (gp1jxzuh = 1; gp1jxzuh <= 1; gp1jxzuh++) {
        yq6lorbx = gp1jxzuh;
        urohxe6t = dqlr5bse - 1 + gp1jxzuh;
        xecbg0pf[urohxe6t-1] +=
          wrk1 * (yw1[gp1jxzuh-1]*yw1[yq6lorbx-1] +
                 (yw2[gp1jxzuh-1]*yw1[yq6lorbx-1] +
                  yw2[yq6lorbx-1]*yw1[gp1jxzuh-1]) * 0.50  +
                  yw2[gp1jxzuh-1]*yw2[yq6lorbx-1]  * othird);
      }
    }
  }
}




void n5aioudkvmnweiy2(double buhyalv4[], double fulcp8wa[], double plj0trqx[], int *xtov9rbf,
                    int *acpios9q, int *wep0oibc, int *iflag) {


  int    ayfnwr1v, yq6lorbx, gp1jxzuh;
  double wjm3[3], wjm2[2], wjm1[1], c0, c1, c2, c3;
  double pcsuow9k, qdbgu6oi, upwkh5xz, rul5fnyd, ueydbrg6,
         plce2srm, k3yvomnh, bfdjhu7l, ctfvwdu0;

  c1 = c2 = c3 = 0.0e0;







  wjm3[0] = wjm3[1] = wjm3[2] =
  wjm2[0] = wjm2[1] =
  wjm1[0] = 0.0e0;

  for (ayfnwr1v = 1; ayfnwr1v <= *acpios9q; ayfnwr1v++) {
    yq6lorbx = *acpios9q - ayfnwr1v + 1;
    c0 = 1.0e0 / buhyalv4[3 + (yq6lorbx-1) * *xtov9rbf];
    if (yq6lorbx <= (*acpios9q-3)) {
      c1 = buhyalv4[0 + (yq6lorbx+2) * *xtov9rbf] * c0;
      c2 = buhyalv4[1 + (yq6lorbx+1) * *xtov9rbf] * c0;
      c3 = buhyalv4[2 + (yq6lorbx+0) * *xtov9rbf] * c0;
    } else if (yq6lorbx == (*acpios9q - 2)) {
      c1 = 0.0e0;
      c2 = buhyalv4[1 + (yq6lorbx+1) * *xtov9rbf] * c0;
      c3 = buhyalv4[2 +  yq6lorbx    * *xtov9rbf] * c0;
    } else if (yq6lorbx == (*acpios9q - 1)) {
      c1 =
      c2 = 0.0e0;
      c3 = buhyalv4[2 +  yq6lorbx    * *xtov9rbf] * c0;
    } else if (yq6lorbx ==  *acpios9q) {
      c1 =
      c2 =
      c3 = 0.0e0;
    }

    pcsuow9k = c1 * wjm3[0];
    qdbgu6oi = c2 * wjm3[1];
    upwkh5xz = c3 * wjm3[2];
    rul5fnyd = c1 * wjm3[1];
    ueydbrg6 = c2 * wjm2[0];
    plce2srm = c3 * wjm2[1];
    k3yvomnh = c1 * wjm3[2];
    bfdjhu7l = c2 * wjm2[1];
    ctfvwdu0 = c3 * wjm1[0];
    fulcp8wa[0 + (yq6lorbx-1) * *xtov9rbf] = 0.0 - (pcsuow9k+qdbgu6oi+upwkh5xz);
    fulcp8wa[1 + (yq6lorbx-1) * *xtov9rbf] = 0.0 - (rul5fnyd+ueydbrg6+plce2srm);
    fulcp8wa[2 + (yq6lorbx-1) * *xtov9rbf] = 0.0 - (k3yvomnh+bfdjhu7l+ctfvwdu0);

    fulcp8wa[3 + (yq6lorbx-1) * *xtov9rbf] = pow(c0, (double) 2.0) +
                c1 * (pcsuow9k + 2.0e0 * (qdbgu6oi + upwkh5xz)) +
                c2 * (ueydbrg6 + 2.0e0 *  plce2srm) +
                c3 *  ctfvwdu0;


    wjm3[0] = wjm2[0];
    wjm3[1] = wjm2[1];
    wjm3[2] = fulcp8wa[1 + (yq6lorbx-1) * *xtov9rbf];
    wjm2[0] = wjm1[0];
    wjm2[1] = fulcp8wa[2 + (yq6lorbx-1) * *xtov9rbf];
    wjm1[0] = fulcp8wa[3 + (yq6lorbx-1) * *xtov9rbf];
  }


  if (*iflag == 0) {
    return;
  }
  Rprintf("plj0trqx must not be a double of length one!\n");

    for (ayfnwr1v = 1; ayfnwr1v <= *acpios9q; ayfnwr1v++) {
      yq6lorbx = *acpios9q - ayfnwr1v + 1;
      for (gp1jxzuh = 1; gp1jxzuh <= 4 &&
                       yq6lorbx + gp1jxzuh-1 <= *acpios9q; gp1jxzuh++) {
           plj0trqx[yq6lorbx-1 + (yq6lorbx+gp1jxzuh-2) * *wep0oibc] =
           fulcp8wa[4-gp1jxzuh + (yq6lorbx-1)        * *xtov9rbf];
      }
    }

    for (ayfnwr1v = 1; ayfnwr1v <= *acpios9q; ayfnwr1v++) {
      yq6lorbx = *acpios9q - ayfnwr1v + 1;
      for (gp1jxzuh = yq6lorbx-4; gp1jxzuh >= 1; gp1jxzuh--) {
        c0 = 1.0 / buhyalv4[3 + (gp1jxzuh-1) * *xtov9rbf];
        c1 = buhyalv4[0 + (gp1jxzuh+2) * *xtov9rbf] * c0;
        c2 = buhyalv4[1 + (gp1jxzuh+1) * *xtov9rbf] * c0;
        c3 = buhyalv4[2 +  gp1jxzuh    * *xtov9rbf] * c0;
                 plj0trqx[gp1jxzuh-1 + (yq6lorbx-1) * *wep0oibc] = 0.0e0 -
          ( c1 * plj0trqx[gp1jxzuh+2 + (yq6lorbx-1) * *wep0oibc] +
            c2 * plj0trqx[gp1jxzuh+1 + (yq6lorbx-1) * *wep0oibc] +
            c3 * plj0trqx[gp1jxzuh   + (yq6lorbx-1) * *wep0oibc] );
      }
    }
}


void n5aioudkwmhctl9x(double *qgnl3toc, double sjwyig9t[],
                    double po8rwsmy[], int    *kuzxj1lo, int *acpios9q,
                    int *pn9eowxc, // int *icrit,
                    double gkdx5jal[], double rpyis2kc[], double imdvf4hx[],
                    double ifys6woa[], double *i9mwnvqt, double xwy[],
                    double *qcpiaj7f,
                    double zvau2lct[], double f6lsuzax[], double fvh2rwtc[], double dcfir2no[],
                    double xecbg0pf[], double z4grbpiq[], double d7glzhbj[], double v2eydbxs[],
                    double buhyalv4[], double fulcp8wa[], double plj0trqx[],
                    int    *xtov9rbf, int *wep0oibc, int *algpft4y) {


  double ms0qypiw[16], b0, b1, b2, b3, qaltf0nz = 0.1e-10,
         g9fvdrbw[4], qtce8hzo, *chw8lzty, egwbdua212 = 0.0e0;
  int    yu6izdrc = 0, pqneb2ra = 1, bvsquk3z = 3, h2dpsbkr = 4,
         pqzfxw4i, ayfnwr1v, yq6lorbx, dqlr5bse, nkp1 = *acpios9q + 1;
  double *qnwamo0e1, *qnwamo0e2;





  qnwamo0e1 = rpyis2kc; qnwamo0e2 = xwy;
  for (ayfnwr1v = 0; ayfnwr1v <   *acpios9q; ayfnwr1v++) {
    *qnwamo0e1++ = *qnwamo0e2++;
  }

  qnwamo0e1 = zvau2lct; qnwamo0e2 = xecbg0pf;
  for (ayfnwr1v = 0; ayfnwr1v <   *acpios9q;    ayfnwr1v++) {
    buhyalv4[3 +  ayfnwr1v    * *xtov9rbf] = *qnwamo0e1++ + *i9mwnvqt * *qnwamo0e2++;
  }

  qnwamo0e1 = f6lsuzax; qnwamo0e2 = z4grbpiq;
  for (ayfnwr1v = 1; ayfnwr1v <= (*acpios9q-1); ayfnwr1v++) {
    buhyalv4[2 +  ayfnwr1v    * *xtov9rbf] = *qnwamo0e1++ + *i9mwnvqt * *qnwamo0e2++;
  }

  qnwamo0e1 = fvh2rwtc; qnwamo0e2 = d7glzhbj;
  for (ayfnwr1v = 1; ayfnwr1v <= (*acpios9q-2); ayfnwr1v++) {
    buhyalv4[1 + (ayfnwr1v+1) * *xtov9rbf] = *qnwamo0e1++ + *i9mwnvqt * *qnwamo0e2++;
  }

  qnwamo0e1 = dcfir2no; qnwamo0e2 = v2eydbxs;
  for (ayfnwr1v = 1; ayfnwr1v <= (*acpios9q-3); ayfnwr1v++) {
    buhyalv4[    (ayfnwr1v+2) * *xtov9rbf] = *qnwamo0e1++ + *i9mwnvqt * *qnwamo0e2++;
  }

  F77_CALL(dpbfa8)(buhyalv4, xtov9rbf, acpios9q, &bvsquk3z, algpft4y);
  if (*algpft4y != 0) {
    Rprintf("In C function wmhctl9x; Error:\n");
    Rprintf("Leading minor of order %d is not pos-def\n", *algpft4y);
    return;
  }
  F77_CALL(dpbsl8)(buhyalv4, xtov9rbf, acpios9q, &bvsquk3z, rpyis2kc);

  chw8lzty = sjwyig9t;  qnwamo0e1 = imdvf4hx;
  for (ayfnwr1v = 1; ayfnwr1v <= *kuzxj1lo; ayfnwr1v++) {
    F77_CALL(wbvalue)(gkdx5jal, rpyis2kc, acpios9q, &h2dpsbkr,
                      chw8lzty++, &yu6izdrc, qnwamo0e1++);
  }


    n5aioudkvmnweiy2(buhyalv4, fulcp8wa, plj0trqx, xtov9rbf, acpios9q, wep0oibc, &yu6izdrc);

  //Rprintf("first one n5aioudkwmhctl9x pow(po8rwsmy[0], (double) 1.0) = ");
  //Rprintf("%9.5e\n", pow(po8rwsmy[0], (double) 1.0));

    chw8lzty = sjwyig9t;
    for (ayfnwr1v = 1; ayfnwr1v <= *kuzxj1lo; ayfnwr1v++) {

      F77_CALL(vinterv)(gkdx5jal, &nkp1, chw8lzty, &dqlr5bse, &pqzfxw4i);

      if (pqzfxw4i == -1) {
        dqlr5bse = 4;
        *chw8lzty = gkdx5jal[3]       + qaltf0nz;
      } else
      if (pqzfxw4i ==  1) {
        dqlr5bse = *acpios9q;
        *chw8lzty = gkdx5jal[*acpios9q] - qaltf0nz;
      }
      yq6lorbx = dqlr5bse-3;

      F77_CALL(vbsplvd)(gkdx5jal, &h2dpsbkr, chw8lzty++, &dqlr5bse,
                        ms0qypiw, g9fvdrbw, &pqneb2ra);

      b0 = g9fvdrbw[0]; b1 = g9fvdrbw[1]; b2 = g9fvdrbw[2]; b3 = g9fvdrbw[3];

      qtce8hzo = (b0   * (fulcp8wa[3 + (yq6lorbx-1) * *xtov9rbf] * b0 +
                2.0e0* (fulcp8wa[2 + (yq6lorbx-1) * *xtov9rbf] * b1 +
                        fulcp8wa[1 + (yq6lorbx-1) * *xtov9rbf] * b2 +
                        fulcp8wa[0 + (yq6lorbx-1) * *xtov9rbf] * b3)) +
                b1   * (fulcp8wa[3 +  yq6lorbx    * *xtov9rbf] * b1 +
                2.0e0* (fulcp8wa[2 +  yq6lorbx    * *xtov9rbf] * b2 +
                        fulcp8wa[1 +  yq6lorbx    * *xtov9rbf] * b3)) +
                b2   * (fulcp8wa[3 + (yq6lorbx+1) * *xtov9rbf] * b2 +
                2.0e0*  fulcp8wa[2 + (yq6lorbx+1) * *xtov9rbf] * b3) +
                        fulcp8wa[3 + (yq6lorbx+2) * *xtov9rbf] *
                        pow(b3, (double) 2.0)) *
                    po8rwsmy[ayfnwr1v-1];
      ifys6woa[ayfnwr1v-1] = qtce8hzo;
    }

    if (*pn9eowxc == 1) {
      return;
    }


    for (ayfnwr1v = 1; ayfnwr1v <= *kuzxj1lo; ayfnwr1v++) {
      egwbdua212 += ifys6woa[ayfnwr1v-1];
    }
    *qcpiaj7f = pow(*qgnl3toc - egwbdua212, (double) 2.0);
}



void n5aioudkgt9iulbf(double sjwyig9t[], double ghz9vuba[], double po8rwsmy[],
                   double gkdx5jal[], int *rvy1fpli, int *kuzxj1lo, double zyupcmk6[],
                   double zvau2lct[], double f6lsuzax[], double fvh2rwtc[], double dcfir2no[]) {




  double g9fvdrbw[12];  /* 20140522 Effectively g9fvdrbw(4,3), just in case     */

  double ms0qypiw[16], wsvdbx3tk, wv2svdbx3tk, qaltf0nz = 0.1e-9;
  int    ayfnwr1v, yq6lorbx, dqlr5bse, pqzfxw4i, nhnpt1zym1 = *kuzxj1lo + 1,
         pqneb2ra = 1, h2dpsbkr = 4;
  double *qnwamo0e0, *qnwamo0e1,  *qnwamo0e2, *qnwamo0e3, *qnwamo0e4;


  qnwamo0e0 = zvau2lct; qnwamo0e1 = f6lsuzax;  qnwamo0e2 = fvh2rwtc; qnwamo0e3 = dcfir2no; qnwamo0e4 = zyupcmk6;
  for (ayfnwr1v = 0; ayfnwr1v < *kuzxj1lo; ayfnwr1v++) {
    *qnwamo0e0++ = *qnwamo0e1++ = *qnwamo0e2++ = *qnwamo0e3++ = *qnwamo0e4++ = 0.0e0;
  }

  //Rprintf("first one n5aioudkgt9iulbf pow(po8rwsmy[0], (double) 1.0) = ");
  //Rprintf("%9.5e\n", pow(po8rwsmy[0], (double) 1.0));


  for (ayfnwr1v = 1; ayfnwr1v <= *rvy1fpli; ayfnwr1v++) {

    F77_CALL(vinterv)(gkdx5jal, &nhnpt1zym1, sjwyig9t + ayfnwr1v - 1, &dqlr5bse, &pqzfxw4i);

    if (pqzfxw4i == 1) {
      if (sjwyig9t[ayfnwr1v-1] <= (gkdx5jal[dqlr5bse-1] + qaltf0nz)) {
        dqlr5bse--;
      } else {
        return;
      }
    }

    F77_CALL(vbsplvd)(gkdx5jal, &h2dpsbkr, sjwyig9t + ayfnwr1v - 1, &dqlr5bse,
                      ms0qypiw, g9fvdrbw, &pqneb2ra);


    yq6lorbx = dqlr5bse - 4 + 1;
    wsvdbx3tk =     po8rwsmy[ayfnwr1v-1];
    wv2svdbx3tk = wsvdbx3tk * g9fvdrbw[0];




    zyupcmk6[yq6lorbx-1] += wv2svdbx3tk * ghz9vuba[ayfnwr1v-1];


    zvau2lct[yq6lorbx-1]     += wv2svdbx3tk * g9fvdrbw[0];
    f6lsuzax[yq6lorbx-1]     += wv2svdbx3tk * g9fvdrbw[1];
    fvh2rwtc[yq6lorbx-1]     += wv2svdbx3tk * g9fvdrbw[2];
    dcfir2no[yq6lorbx-1]     += wv2svdbx3tk * g9fvdrbw[3];

    yq6lorbx = dqlr5bse - 4 + 2;
    wv2svdbx3tk = wsvdbx3tk * g9fvdrbw[1];
    zyupcmk6[yq6lorbx-1] += wv2svdbx3tk * ghz9vuba[ayfnwr1v-1];
    zvau2lct[yq6lorbx-1]     += wv2svdbx3tk * g9fvdrbw[1];
    f6lsuzax[yq6lorbx-1]     += wv2svdbx3tk * g9fvdrbw[2];
    fvh2rwtc[yq6lorbx-1]     += wv2svdbx3tk * g9fvdrbw[3];

    yq6lorbx = dqlr5bse - 4 + 3;
    wv2svdbx3tk = wsvdbx3tk * g9fvdrbw[2];
    zyupcmk6[yq6lorbx-1] += wv2svdbx3tk * ghz9vuba[ayfnwr1v-1];
    zvau2lct[yq6lorbx-1]     += wv2svdbx3tk * g9fvdrbw[2];
    f6lsuzax[yq6lorbx-1]     += wv2svdbx3tk * g9fvdrbw[3];
    yq6lorbx = dqlr5bse;
    wv2svdbx3tk = wsvdbx3tk * g9fvdrbw[3];
    zyupcmk6[yq6lorbx-1] += wv2svdbx3tk * ghz9vuba[ayfnwr1v-1];
    zvau2lct[yq6lorbx-1]     += wv2svdbx3tk * g9fvdrbw[3];
  }
}


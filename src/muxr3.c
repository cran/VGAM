

#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<R.h>
#include<Rmath.h>


void vdec(int *hqipzx3n, int *exz2jrhq, int *dimm);
void m2a(double *m, double *a, int *dimm, int *hqipzx3n,
       int *exz2jrhq, int *n, int *M, int *rb1onzwu);
void a2m(double *a, double *m, int *dimm, int *hqipzx3n,
       int *exz2jrhq, int *n, int *M);
void mux2(double *cc, double *tlgduey8,
          double *bzmd6ftv, int *p, int *n, int *M);
void mux22(double *cc, double *tlgduey8, double *bzmd6ftv, int *dimm,
       int *hqipzx3n, int *exz2jrhq, 
       int *n, int *M, double *wk, int *rb1onzwu);
void mux5(double *cc, double *x,
          double *bzmd6ftv, int *M, int *n, int *r,
          int *dimm,
          int *dimr,
          int *matrix,
          double *wk, double *wk2,
          int *hqipzx3n_M, int *exz2jrhq_M,
          int *hqipzx3n_r, int *exz2jrhq_r);
void mux55(double *evects, double *evals, double *bzmd6ftv, double *wk, 
           double *wk2, int *hqipzx3n, int *exz2jrhq,
           int *M, int *n);
void mux7(double *cc, double *x,
          double *bzmd6ftv, int *M, int *q, int *n, int *r);
void mux111(double *cc, double *the7mqnvy, int *M, int *R, int *n,
        double *wk, double *wk2, int *hqipzx3n, int *exz2jrhq,
	    int *dimm, int *rb1onzwu);
void mux15(double *cc, double *x,
           double *bzmd6ftv, int *M, int *n);
void vchol(double *cc, int *M, int *n, int *ok, double *wk,
           int *hqipzx3n, int *exz2jrhq, int *dimm);
void vforsub(double *cc, double *b, int *M, int *n,
             double *wk, int *hqipzx3n,
             int *exz2jrhq, int *dimm);
void vbacksub(double *cc, double *b, int *M, int *n,
              double *wk, int *hqipzx3n,
              int *exz2jrhq, int *dimm);
void tapply_mat1(double *mat, int *nr, int *nc, int *type);
















void vdec(int *hqipzx3n, int *exz2jrhq, int *dimm) {
  int ayfnwr1v;

  for(ayfnwr1v = 0; ayfnwr1v < *dimm; ayfnwr1v++) {
    hqipzx3n[ayfnwr1v] -= 1;
    exz2jrhq[ayfnwr1v] -= 1;
  }
}


void m2a(double *m, double *a, int *dimm, int *hqipzx3n,
         int *exz2jrhq, int *n, int *M, int *rb1onzwu) {
  int ayfnwr1v, gp1jxzuh, MM = *M * *M, MMn = *M * *M * *n;

  if(*rb1onzwu == 1 || *dimm != *M * (*M + 1) / 2)
    for(gp1jxzuh = 0; gp1jxzuh < MMn; gp1jxzuh++)
        a[gp1jxzuh] = 0.0;

  for(gp1jxzuh = 0; gp1jxzuh < *n; gp1jxzuh++) {
    for(ayfnwr1v = 0; ayfnwr1v < *dimm; ayfnwr1v++) {
      a[hqipzx3n[ayfnwr1v] + exz2jrhq[ayfnwr1v] * *M] = m[ayfnwr1v];
      if(*rb1onzwu == 0)
        a[exz2jrhq[ayfnwr1v] + hqipzx3n[ayfnwr1v] * *M] = m[ayfnwr1v];
    }
    a += MM;
    m += *dimm;
  }
}


void a2m(double *a, double *m, int *dimm, int *hqipzx3n,
         int *exz2jrhq, int *n, int *M) {
  int ayfnwr1v, gp1jxzuh, MM= *M * *M;

  for(gp1jxzuh = 0; gp1jxzuh < *n; gp1jxzuh++) {
    for(ayfnwr1v = 0; ayfnwr1v < *dimm; ayfnwr1v++)
      m[ayfnwr1v] = a[hqipzx3n[ayfnwr1v] + exz2jrhq[ayfnwr1v] * *M];
    a += MM;
    m += *dimm;
  }
}




void mux2(double *cc, double *tlgduey8,
          double *bzmd6ftv, int *p, int *n, int *M) {
  double urohxe6t;
  int ayfnwr1v, yq6lorbx, bpvaqm5z, Mp = *M * *p;

  for(ayfnwr1v = 0; ayfnwr1v < *n; ayfnwr1v++) {
    for(yq6lorbx = 0; yq6lorbx < *M; yq6lorbx++) {
      urohxe6t = 0.0;
      for(bpvaqm5z = 0; bpvaqm5z < *p; bpvaqm5z++)
        urohxe6t += cc[yq6lorbx + bpvaqm5z * *M] * tlgduey8[bpvaqm5z];
      *bzmd6ftv++ = urohxe6t;
    }
    tlgduey8 += *p;
    cc += Mp;
  }
}



void mux22(double *cc, double *tlgduey8, double *bzmd6ftv, int *dimm,
       int *hqipzx3n, int *exz2jrhq, 
       int *n, int *M, double *wk, int *rb1onzwu) {
  double urohxe6t;
  int yq6lorbx, bpvaqm5z, gp1jxzuh, one = 1, nzqklc9x;

  vdec(hqipzx3n, exz2jrhq, dimm);
  for(gp1jxzuh = 0; gp1jxzuh < *n; gp1jxzuh++) {
    m2a(cc, wk, dimm, hqipzx3n, exz2jrhq, &one, M, rb1onzwu);

    for(yq6lorbx = 0; yq6lorbx < *M; yq6lorbx++) {
      urohxe6t = 0.0;
      nzqklc9x = *rb1onzwu == 0 ? 0 : yq6lorbx; 
      for(bpvaqm5z = nzqklc9x; bpvaqm5z < *M; bpvaqm5z++)
          urohxe6t += wk[yq6lorbx + bpvaqm5z * *M] * tlgduey8[bpvaqm5z];
      *bzmd6ftv++ = urohxe6t;
    }
    tlgduey8 += *M;
    cc += *dimm;
  }
}


void mux5(double *cc, double *x,
          double *bzmd6ftv, int *M, int *n, int *r,
          int *dimm,
          int *dimr,
          int *matrix,
          double *wk, double *wk2,
          int *hqipzx3n_M, int *exz2jrhq_M,
          int *hqipzx3n_r, int *exz2jrhq_r) {
  double urohxe6t, *pd, *pd2;
  int ayfnwr1v, yq6lorbx, gp1jxzuh, bpvaqm5z, Mr = *M * *r,
      rr = *r * *r, MM = *M * *M, usvdbx3tk,
      jM, jr, kM, kr, one=1, rb1onzwu=0;

  if(*matrix == 1) {
    vdec(hqipzx3n_M, exz2jrhq_M, dimm);
    vdec(hqipzx3n_r, exz2jrhq_r, dimr);
    pd = wk;
    pd2 = wk2;
  } else {

    pd = pd2 = wk;

  }

  for(ayfnwr1v = 0; ayfnwr1v < *n; ayfnwr1v++) {
    if(*matrix == 1)
      m2a(cc, pd, dimm, hqipzx3n_M, exz2jrhq_M, &one, M, &rb1onzwu);
    else {
      pd = cc;
      pd2 = bzmd6ftv;
    }

    for(yq6lorbx = 0; yq6lorbx < *r; yq6lorbx++) {
      jM = yq6lorbx * *M;
      jr = yq6lorbx * *r;
      for(gp1jxzuh = yq6lorbx; gp1jxzuh < *r; gp1jxzuh++) {
        kM = gp1jxzuh * *M;
        kr = gp1jxzuh * *r;
        urohxe6t = 0.0;
        for(bpvaqm5z = 0; bpvaqm5z < *M; bpvaqm5z++)
          for(usvdbx3tk = 0; usvdbx3tk < *M; usvdbx3tk++)
            urohxe6t +=  x[bpvaqm5z + jM] * pd[bpvaqm5z + usvdbx3tk * *M] *
                       x[usvdbx3tk + kM];
        pd2[yq6lorbx + kr] =
        pd2[gp1jxzuh + jr] = urohxe6t;
      }
    }

    if(*matrix == 1)
      a2m(pd2, bzmd6ftv, dimr, hqipzx3n_r, exz2jrhq_r, &one, r);

    cc += (*matrix == 1 ? *dimm : MM);
    x += Mr;
    bzmd6ftv += (*matrix == 1 ? *dimr : rr);
  }
}



void mux55(double *evects, double *evals, double *bzmd6ftv, double *wk, 
           double *wk2, int *hqipzx3n, int *exz2jrhq,
           int *M, int *n) {
  double *pd, *pd2, bpvaqm5z;
  int ayfnwr1v, yq6lorbx, gp1jxzuh, urohxe6t, MM = *M * *M, one = 1,
      imk5wjxg = *M * (*M + 1)/2;

  vdec(hqipzx3n, exz2jrhq, &imk5wjxg);

  for(ayfnwr1v = 0; ayfnwr1v < *n; ayfnwr1v++) {
    pd = evects;
    pd2 = wk2;
    for(yq6lorbx = 0; yq6lorbx < *M; yq6lorbx++)
      for(gp1jxzuh = 0; gp1jxzuh < *M; gp1jxzuh++)
        *pd2++ = *pd++ * evals[yq6lorbx];

    for(yq6lorbx = 0; yq6lorbx < *M; yq6lorbx++)
      for(gp1jxzuh = yq6lorbx; gp1jxzuh < *M; gp1jxzuh++) {
        bpvaqm5z = 0.0; 
        for(urohxe6t = 0; urohxe6t < *M; urohxe6t++)
            bpvaqm5z +=    wk2[yq6lorbx + urohxe6t * *M] *
                      evects[gp1jxzuh + urohxe6t * *M];
        wk[yq6lorbx + gp1jxzuh * *M] =
        wk[gp1jxzuh + yq6lorbx * *M] = bpvaqm5z;
      }

    a2m(wk, bzmd6ftv, &imk5wjxg, hqipzx3n, exz2jrhq, &one, M);

    bzmd6ftv += imk5wjxg;
    evals += *M;
    evects += MM;
  }
}





void mux7(double *cc, double *x,
          double *bzmd6ftv, int *M, int *q, int *n, int *r) {
  double urohxe6t;
  int ayfnwr1v, yq6lorbx, gp1jxzuh, bpvaqm5z,
      Mq = *M * *q, qr = *q * *r, Mr = *M * *r,
      kq, kM;

  for(ayfnwr1v = 0; ayfnwr1v < *n; ayfnwr1v++) {
    for(yq6lorbx = 0; yq6lorbx < *M; yq6lorbx++) {
      for(gp1jxzuh = 0; gp1jxzuh < *r; gp1jxzuh++) {
        kq = gp1jxzuh * *q;
        kM = gp1jxzuh * *M;
        urohxe6t = 0.0;
        for(bpvaqm5z = 0; bpvaqm5z < *q; bpvaqm5z++)
          urohxe6t += cc[yq6lorbx + bpvaqm5z * *M] * x[bpvaqm5z + kq];
        bzmd6ftv[yq6lorbx + kM] = urohxe6t;
      }
    }
    cc += Mq;
    bzmd6ftv += Mr;
    x += qr;
  }
}



void mux111(double *cc, double *the7mqnvy, int *M, int *R, int *n,
            double *wk, double *wk2, int *hqipzx3n, int *exz2jrhq,
            int *dimm, int *rb1onzwu) {
  double urohxe6t, *pd2;
  int ayfnwr1v, yq6lorbx, gp1jxzuh, bpvaqm5z,
      MM = *M * *M, MR = *M * *R, nzqklc9x;

  vdec(hqipzx3n, exz2jrhq, dimm);

  for(ayfnwr1v = 0; ayfnwr1v < MM; ayfnwr1v++)
    wk[ayfnwr1v] = 0.0;

  for(bpvaqm5z = 0; bpvaqm5z < *n; bpvaqm5z++) {
    for(ayfnwr1v = 0; ayfnwr1v < *dimm; ayfnwr1v++) {
      if(*rb1onzwu == 0)
        wk[hqipzx3n[ayfnwr1v] + exz2jrhq[ayfnwr1v] * *M] =
        wk[exz2jrhq[ayfnwr1v] + hqipzx3n[ayfnwr1v] * *M] = *cc++;
      else
        wk[hqipzx3n[ayfnwr1v] + exz2jrhq[ayfnwr1v] * *M] = *cc++;
    }

    pd2 = the7mqnvy;
    for(ayfnwr1v = 0; ayfnwr1v < *M; ayfnwr1v++)
      for(yq6lorbx = 0; yq6lorbx < *R; yq6lorbx++)
        wk2[ayfnwr1v + yq6lorbx * *M] = *pd2++;

    for(ayfnwr1v = 0; ayfnwr1v < *M; ayfnwr1v++)
      for(yq6lorbx = 0; yq6lorbx < *R; yq6lorbx++) {
        urohxe6t = 0.0;
        nzqklc9x = *rb1onzwu == 0 ? 0 : ayfnwr1v;
        for(gp1jxzuh = nzqklc9x; gp1jxzuh < *M; gp1jxzuh++)
          urohxe6t += wk[ayfnwr1v + gp1jxzuh * *M] * wk2[gp1jxzuh + yq6lorbx * *M];
        the7mqnvy[yq6lorbx + ayfnwr1v * *R] = urohxe6t;
      }
    the7mqnvy += MR;
  }
}




void mux15(double *cc, double *x,
           double *bzmd6ftv, int *M, int *n) {
  double *pd, *pd2;
  int ayfnwr1v, yq6lorbx, gp1jxzuh, MM = *M * *M;

  for(ayfnwr1v = 0; ayfnwr1v < *n; ayfnwr1v++) {
    pd = cc;
    pd2 = bzmd6ftv;
    for(yq6lorbx = 0; yq6lorbx < *M; yq6lorbx++)
      for(gp1jxzuh = 0; gp1jxzuh < *M; gp1jxzuh++)
        *pd2++ = *pd++ * x[yq6lorbx];

    pd2 = bzmd6ftv;
    for(yq6lorbx = 0; yq6lorbx < *M; yq6lorbx++)
      for(gp1jxzuh = 0; gp1jxzuh < *M; gp1jxzuh++) {
        *pd2 *= x[gp1jxzuh];
        pd2++;
      }

    bzmd6ftv += MM;
    x += *M;
  }
}




void vchol(double *cc, int *M, int *n, int *ok, double *wk,
           int *hqipzx3n, int *exz2jrhq, int *dimm) {
  double urohxe6t, *pd;
  int bpvaqm5z, ayfnwr1v, yq6lorbx, gp1jxzuh, iM, iiM, rb1onzwu = 0, one = 1;

  vdec(hqipzx3n, exz2jrhq, dimm);
  pd = wk;

  for(bpvaqm5z = 0; bpvaqm5z < *n; bpvaqm5z++) {
    *ok = 1; 

    m2a(cc, wk, dimm, hqipzx3n, exz2jrhq, &one, M, &rb1onzwu);

    for(ayfnwr1v = 0; ayfnwr1v < *M; ayfnwr1v++) {
      urohxe6t = 0.0;
      iM = ayfnwr1v * *M;
      iiM = ayfnwr1v + iM;
      for(gp1jxzuh = 0; gp1jxzuh < ayfnwr1v; gp1jxzuh++)
        urohxe6t += pd[gp1jxzuh + iM] * pd[gp1jxzuh + iM];

      pd[iiM] -= urohxe6t;
      if(pd[iiM] < 0.0) {
        *ok = 0;
        break;
      }
      pd[iiM] = sqrt(pd[iiM]);

      for(yq6lorbx = ayfnwr1v+1; yq6lorbx < *M; yq6lorbx++) {
        urohxe6t = 0.0;
        for(gp1jxzuh = 0; gp1jxzuh < ayfnwr1v; gp1jxzuh++)
          urohxe6t += pd[gp1jxzuh + iM] * pd[gp1jxzuh + yq6lorbx * *M];
        pd[ayfnwr1v + yq6lorbx * *M] = (pd[ayfnwr1v + yq6lorbx * *M] -
                                    urohxe6t) / pd[iiM];
      }
    }

    a2m(wk, cc, dimm, hqipzx3n, exz2jrhq, &one, M);

    cc += *dimm;
    ok++;
  }
}



void vforsub(double *cc, double *b, int *M, int *n,
             double *wk, int *hqipzx3n,
             int *exz2jrhq, int *dimm) {
  double urohxe6t, *pd;
  int yq6lorbx, gp1jxzuh, bpvaqm5z, rb1onzwu = 1, one = 1;

  pd = wk;
  vdec(hqipzx3n, exz2jrhq, dimm);

  for(bpvaqm5z = 0; bpvaqm5z < *n; bpvaqm5z++) {
    m2a(cc, wk, dimm, hqipzx3n, exz2jrhq, &one, M, &rb1onzwu);

    for(yq6lorbx = 0; yq6lorbx < *M; yq6lorbx++) {
      urohxe6t = b[yq6lorbx];
      for(gp1jxzuh = 0; gp1jxzuh < yq6lorbx; gp1jxzuh++)
        urohxe6t -= pd[gp1jxzuh + yq6lorbx * *M] * b[gp1jxzuh];
      b[yq6lorbx] = urohxe6t / pd[yq6lorbx + yq6lorbx * *M];
    }
    cc += *dimm;
    b += *M;
  }
}




void vbacksub(double *cc, double *b, int *M, int *n,
              double *wk, int *hqipzx3n,
              int *exz2jrhq, int *dimm) {
  double urohxe6t, *pd;
  int yq6lorbx, gp1jxzuh, bpvaqm5z, rb1onzwu = 1, one = 1;

  pd = wk;
  vdec(hqipzx3n, exz2jrhq, dimm);

  for(bpvaqm5z = 0; bpvaqm5z < *n; bpvaqm5z++) {
    m2a(cc, wk, dimm, hqipzx3n, exz2jrhq, &one, M, &rb1onzwu);

    for(yq6lorbx = *M - 1; yq6lorbx >= 0; yq6lorbx--) {
      urohxe6t = b[yq6lorbx];
      for(gp1jxzuh = yq6lorbx + 1; gp1jxzuh < *M; gp1jxzuh++)
        urohxe6t -= pd[yq6lorbx + gp1jxzuh * *M] * b[gp1jxzuh];
      b[yq6lorbx] = urohxe6t / pd[yq6lorbx + yq6lorbx * *M];
    }
    cc += *dimm;
    b += *M;
  }
}



void tapply_mat1(double *mat, int *nr, int *nc, int *type) {
  double *pd = mat, *pd2 = mat + *nr;
  int ayfnwr1v, yq6lorbx;

  if(*type == 1)
    for(yq6lorbx = 2; yq6lorbx <= *nc; yq6lorbx++)
      for(ayfnwr1v = 0; ayfnwr1v < *nr; ayfnwr1v++, pd2++)
        *pd2 += *pd++;

  if(*type == 2) {
    pd2 = mat + *nr * *nc - 1;
    pd = pd2 - *nr;
    for(yq6lorbx = *nc; yq6lorbx >= 2; yq6lorbx--)
      for(ayfnwr1v = 0; ayfnwr1v < *nr; ayfnwr1v++, pd2--)
        *pd2 -= *pd--;
  }

  if(*type == 3)
    for(yq6lorbx = 2; yq6lorbx <= *nc; yq6lorbx++)
      for(ayfnwr1v = 0; ayfnwr1v < *nr; ayfnwr1v++, pd2++)
        *pd2 *= *pd++;

  if(*type < 1 || *type > 3)
    Rprintf("Error: *type not ezlgm2uped\n");
}





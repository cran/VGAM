/*
This code is
Copyright (C) 1998-2005 T. W. Yee, University of Auckland. All rights reserved.
*/


#include <stdio.h>
#include <math.h>


void vdec(int *row_index, int *col_index, int *dimm)
{
    int i;

    for(i = 0; i < *dimm; i++) {
        row_index[i] -= 1;
        col_index[i] -= 1;
    }
}


void m2a(double *m, double *a, int *dimm, int *row_index,
       int *col_index, int *n, int *M, int *upper)
{
    int i, k, MM = *M * *M, MMn = *M * *M * *n;

    if(*upper == 1 || *dimm != *M * (*M + 1) / 2)
        for(k = 0; k < MMn; k++)
            a[k] = 0.0;

    for(k = 0; k < *n; k++)
    {
        for(i = 0; i < *dimm; i++)
        {
            a[row_index[i] + col_index[i] * *M] = m[i];
            if(*upper == 0)
                a[col_index[i] + row_index[i] * *M] = m[i];
        }
        a += MM;
        m += *dimm;
    }
}


void a2m(double *a, double *m, int *dimm, int *row_index,
       int *col_index, int *n, int *M)
{
    int i, k, MM= *M * *M;

    for(k = 0; k < *n; k++)
    {
        for(i = 0; i < *dimm; i++)
            m[i] = a[row_index[i] + col_index[i] * *M];
        a += MM;
        m += *dimm;
    }
}




void mux2(double *cc, double *ymat,
          double *ans, int *p, int *n, int *M)
{
    double s;
    int i, j, t, Mp = *M * *p;

    for(i = 0; i < *n; i++)
    {
        for(j = 0; j < *M; j++)
        {
            s = 0.0;
            for(t = 0; t < *p; t++)
                s += cc[j + t * *M] * ymat[t];
            *ans++ = s;
        }
        ymat += *p;
        cc += Mp;
    }
}



void mux22(double *cc, double *ymat, double *ans, int *dimm,
       int *row_index, int *col_index, 
       int *n, int *M, double *wk, int *upper)
{
    double s;
    int j, t, k, one = 1, lower;

    vdec(row_index, col_index, dimm);
    for(k = 0; k < *n; k++)
    {
        m2a(cc, wk, dimm, row_index, col_index, &one, M, upper);

        for(j = 0; j < *M; j++)
        {
            s = 0.0;
            lower = *upper == 0 ? 0 : j; 
            for(t = lower; t < *M; t++)
                s += wk[j + t * *M] * ymat[t];
            *ans++ = s;
        }
        ymat += *M;
        cc += *dimm;
    }
}


void mux5(double *cc, double *x,
          double *ans, int *M, int *n, int *r,
          int *dimm,
          int *dimr,
          int *matrix,
          double *wk, double *wk2,
          int *row_index_M, int *col_index_M,
          int *row_index_r, int *col_index_r)
{
    double s, *pd, *pd2;
    int i, j, k, t, Mr = *M * *r, rr = *r * *r, MM = *M * *M, u,
         jM, jr, kM, kr, one=1, upper=0;

    if(*matrix == 1) 
    {
        vdec(row_index_M, col_index_M, dimm);
        vdec(row_index_r, col_index_r, dimr);
        pd = wk;
        pd2 = wk2;
    } else {
/* Commented out on 2/5/06. Need to fix this up more cleanly.
        printf("Error: can only handle matrix.arg == 1\n");
        exit(-1); 
*/

/*
26/9/07:
The following line was added only to avoid a warning message from the compiler
*/
        pd = pd2 = wk;

    }

    for(i = 0; i < *n; i++)
    {
        if(*matrix == 1)
            m2a(cc, pd, dimm, row_index_M, col_index_M, &one, M, &upper);
        else
        {
            pd = cc;
            pd2 = ans;
        }

        for(j = 0; j < *r; j++)
        {
            jM = j * *M;
            jr = j * *r;
            for(k = j; k < *r; k++)
            {
                kM = k * *M;
                kr = k * *r;
                s = 0.0;
                for(t = 0; t < *M; t++)
                    for(u = 0; u < *M; u++)
                        s +=  x[t + jM] * pd[t + u * *M] * x[u + kM];
                pd2[j + kr] =
                pd2[k + jr] = s;
            }

        }

        if(*matrix == 1)
            a2m(pd2, ans, dimr, row_index_r, col_index_r, &one, r);

        cc += (*matrix == 1 ? *dimm : MM);
        x += Mr;
        ans += (*matrix == 1 ? *dimr : rr);
    }
}



void mux55(double *evects, double *evals, double *ans, double *wk, 
           double *wk2, int *row_index, int *col_index,
           int *M, int *n)
{
    double *pd, *pd2, t;
    int i, j, k, s, MM = *M * *M, one=1,
         MM12 = *M * (*M + 1)/2;

    vdec(row_index, col_index, &MM12);

    for(i = 0; i < *n; i++)
    {
        pd = evects;
        pd2 = wk2;
        for(j = 0; j < *M; j++)
            for(k = 0; k < *M; k++)
                *pd2++ = *pd++ * evals[j];

        for(j = 0; j < *M; j++)
            for(k = j; k < *M; k++)
            {
                t = 0.0; 
                for(s = 0; s < *M; s++)
                    t += wk2[j + s * *M] * evects[k + s * *M];
                wk[j + k * *M] =
                wk[k + j * *M] = t; 
            }

        a2m(wk, ans, &MM12, row_index, col_index, &one, M);

        ans += MM12;
        evals += *M;
        evects += MM;
    }
}





void mux7(double *cc, double *x,
          double *ans, int *M, int *q, int *n, int *r)
{
    double s;
    int i, j, k, t, Mq = *M * *q, qr = *q * *r, Mr = *M * *r,
         kq, kM;

    for(i = 0; i < *n; i++)
    {
        for(j = 0; j < *M; j++)
        {
            for(k = 0; k < *r; k++)
            {
                kq = k * *q;
                kM = k * *M;
                s = 0.0;
                for(t = 0; t < *q; t++)
                    s += cc[j + t * *M] * x[t + kq];
                ans[j + kM] = s;
            }
        }
        cc += Mq;
        ans += Mr;
        x += qr;
    }
}



void mux111(double *cc, double *txmat, int *M, int *R, int *n,
        double *wk, double *wk2, int *row_index, int *col_index,
        int *dimm, int *upper)
{
    double s, *pd2;
    int i, j, k, t, MM = *M * *M, MR = *M * *R, lower;

    vdec(row_index, col_index, dimm);

    for(i = 0; i < MM; i++)
        wk[i] = 0.0;

    for(t = 0; t < *n; t++)
    {
        for(i = 0; i < *dimm; i++) 
        {
            if(*upper == 0)
                wk[row_index[i] + col_index[i] * *M] =
                wk[col_index[i] + row_index[i] * *M] = *cc++;
            else
                wk[row_index[i] + col_index[i] * *M] = *cc++;
        }

        pd2 = txmat;
        for(i = 0; i < *M; i++)
            for(j = 0; j < *R; j++)
                wk2[i + j * *M] = *pd2++;

        for(i = 0; i < *M; i++)
            for(j = 0; j < *R; j++)
            {
                s = 0.0;
                lower = *upper == 0 ? 0 : i;
                for(k = lower; k < *M; k++)
                    s += wk[i + k * *M] * wk2[k + j * *M];
                txmat[j + i * *R] = s;
            }
        txmat += MR;
    }
}




void mux15(double *cc, double *x,
           double *ans, int *M, int *n)
{
    double *pd, *pd2;
    int i, j, k, MM = *M * *M;

    for(i = 0; i < *n; i++)
    {
        pd = cc;
        pd2 = ans;
        for(j = 0; j < *M; j++)
            for(k = 0; k < *M; k++)
                *pd2++ = *pd++ * x[j];

        pd2 = ans;
        for(j = 0; j < *M; j++)
            for(k = 0; k < *M; k++)
            {
                *pd2 *= x[k];
                pd2++;
            }

        ans += MM;
        x += *M;
    }
}




void vchol(double *cc, int *M, int *n, int *ok, double *wk,
           int *row_index, int *col_index, int *dimm)
 
{
    double s, *pd;
    int t, i, j, k, iM, iiM, upper = 0, one = 1;

    vdec(row_index, col_index, dimm);
    pd = wk;

    for(t = 0; t < *n; t++)
    {
        *ok = 1; 

        m2a(cc, wk, dimm, row_index, col_index, &one, M, &upper);

        for(i = 0; i < *M; i++)
        {
            s = 0.0;
            iM = i * *M;
            iiM = i + iM;
            for(k = 0; k < i; k++)
                s += pd[k + iM] * pd[k + iM];

            pd[iiM] -= s;
            if(pd[iiM] < 0.0) 
            {
	        *ok = 0;
                break;
            }
            pd[iiM] = sqrt(pd[iiM]);

            for(j = i+1; j < *M; j++)
            {
                s = 0.0;
                for(k = 0; k < i; k++)
                    s += pd[k + iM] * pd[k + j * *M];
                pd[i + j * *M] = (pd[i + j * *M] - s) / pd[iiM];
            }

        }

        a2m(wk, cc, dimm, row_index, col_index, &one, M);

        cc += *dimm;
        ok++;
    }
}



void vforsub(double *cc, double *b, int *M, int *n,
             double *wk, int *row_index,
             int *col_index, int *dimm)
{
    double s, *pd;
    int j, k, t, upper = 1, one = 1;

    pd = wk;
    vdec(row_index, col_index, dimm);

    for(t = 0; t < *n; t++)
    {
        m2a(cc, wk, dimm, row_index, col_index, &one, M, &upper);

        for(j = 0; j < *M; j++)
        {
            s = b[j];
            for(k = 0; k < j; k++)
                s -= pd[k + j * *M] * b[k];
            b[j] = s / pd[j + j * *M];
        }
        cc += *dimm;
        b += *M;
    }
}




void vbacksub(double *cc, double *b, int *M, int *n,
              double *wk, int *row_index,
              int *col_index, int *dimm)
{
    double s, *pd;
    int j, k, t, upper = 1, one = 1;

    pd = wk;
    vdec(row_index, col_index, dimm);

    for(t = 0; t < *n; t++)
    {
        m2a(cc, wk, dimm, row_index, col_index, &one, M, &upper);

        for(j = *M - 1; j >= 0; j--)
        {
            s = b[j];
            for(k = j + 1; k < *M; k++)
                s -= pd[j + k * *M] * b[k];
            b[j] = s / pd[j + j * *M];
        }
        cc += *dimm;
        b += *M;
    }
}



void tapplymat1(double *mat, int *nr, int *nc, int *type)
{
    double *pd = mat, *pd2 = mat + *nr;
    int i, j;

    if(*type==1) 
        for(j = 2; j <= *nc; j++)
            for(i = 0; i < *nr; i++, pd2++)
                *pd2 += *pd++;

    if(*type==2) 
    {
        pd2 = mat + *nr * *nc - 1;
        pd = pd2 - *nr;
        for(j = *nc; j >= 2; j--)
            for(i = 0; i < *nr; i++, pd2--)
                *pd2 -= *pd--;
    }

    if(*type==3) 
        for(j = 2; j <= *nc; j++)
            for(i = 0; i < *nr; i++, pd2++)
                *pd2 *= *pd++;

    if(*type < 1 || *type > 3)
        printf("Error: *type not matched\n");
}


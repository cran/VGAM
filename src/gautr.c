#include "math.h"

/* Frequently used numerical constants: */
#define OneUponSqrt2Pi .39894228040143267794
#define twopi 6.283195307179587
#define LnSqrt2Pi -0.9189385332046727417803296 
#define SQRT2 1.414213562373095049
#define SQRTPI 1.772453850905516027

/* ---------------------------------------------------------------------------

   UNIVARIATE NORMAL PROBABILITY

   ---------------------------------------------------------------------------*/

#define UPPERLIMIT 20.0 
/* I won't return either of univariate normal density or 
	probability when x < -UPPERLIMIT  or x > UPPERLIMIT. */

#define P10 242.66795523053175
#define P11 21.979261618294152
#define P12 6.9963834886191355
#define P13 -.035609843701815385
#define Q10 215.05887586986120
#define Q11 91.164905404514901
#define Q12 15.082797630407787
#define Q13 1.0

#define P20 300.4592610201616005
#define P21 451.9189537118729422
#define P22 339.3208167343436870
#define P23 152.9892850469404039
#define P24 43.16222722205673530
#define P25 7.211758250883093659
#define P26 .5641955174789739711
#define P27 -.0000001368648573827167067
#define Q20 300.4592609569832933
#define Q21 790.9509253278980272
#define Q22 931.3540948506096211
#define Q23 638.9802644656311665
#define Q24 277.5854447439876434
#define Q25 77.00015293522947295
#define Q26 12.78272731962942351
#define Q27 1.0

#define P30 -.00299610707703542174
#define P31 -.0494730910623250734
#define P32 -.226956593539686930
#define P33 -.278661308609647788
#define P34 -.0223192459734184686
#define Q30 .0106209230528467918
#define Q31 .191308926107829841
#define Q32 1.05167510706793207
#define Q33 1.98733201817135256
#define Q34 1.0

double pnorm1(double x)
{
	int    sn;
	double R1, R2, y, y2, y3, y4, y5, y6, y7;
	double erf, erfc, z, z2, z3, z4;
	double phi;

	if (x < -UPPERLIMIT) return 0.0;
	if (x >  UPPERLIMIT) return 1.0;

	y = x / SQRT2;
	if (y < 0) 
        {
		y = -y;
		sn = -1;
	}
	else
		sn = 1;

	y2 = y * y;
	y4 = y2 * y2;
	y6 = y4 * y2;

	if(y < 0.46875) 
        {
		R1 = P10 + P11 * y2 + P12 * y4 + P13 * y6;
		R2 = Q10 + Q11 * y2 + Q12 * y4 + Q13 * y6;
		erf = y * R1 / R2;
		if (sn == 1)
			phi = 0.5 + 0.5*erf;
		else 
			phi = 0.5 - 0.5*erf;
	}
	else 
		if (y < 4.0) 
                {
			y3 = y2 * y;
			y5 = y4 * y;
			y7 = y6 * y;
			R1 = P20 + P21 * y + P22 * y2 + P23 * y3 + 
			    P24 * y4 + P25 * y5 + P26 * y6 + P27 * y7;
			R2 = Q20 + Q21 * y + Q22 * y2 + Q23 * y3 + 
			    Q24 * y4 + Q25 * y5 + Q26 * y6 + Q27 * y7;
			erfc = exp(-y2) * R1 / R2;
			if (sn == 1)
				phi = 1.0 - 0.5*erfc;
			else
				phi = 0.5*erfc;
		}
		else 
                {
			z = y4;
			z2 = z * z;
			z3 = z2 * z;
			z4 = z2 * z2;
			R1 = P30 + P31 * z + P32 * z2 + P33 * z3 + P34 * z4;
			R2 = Q30 + Q31 * z + Q32 * z2 + Q33 * z3 + Q34 * z4;
			erfc = (exp(-y2)/y) * (1.0 / SQRTPI + R1 / (R2 * y2));
			if (sn == 1)
				phi = 1.0 - 0.5*erfc;
			else 
				phi = 0.5*erfc;
		}

	return phi;
}



/* ---------------------------------------------------------------------------

   UNIVARIATE NORMAL DENSITY

   ---------------------------------------------------------------------------*/

double dnorm1(double x)
{
	if (x < -UPPERLIMIT) return 0.0;
	if (x >  UPPERLIMIT) return 0.0;
	return OneUponSqrt2Pi * exp(-0.5 * x * x);
}



/* ---------------------------------------------------------------------------

   LN OF UNIVARIATE NORMAL DENSITY

   ---------------------------------------------------------------------------*/

double lndnorm1(double x)
{
	return LnSqrt2Pi - (0.5*x*x);
}



/*---------------------------------------------------------------------------

  BIVARIATE NORMAL PROBABILITY

  ---------------------------------------------------------------------------*/

#define con (twopi / 2.0) * 10.0e-10

double bivnor(double ah, double ak, double r)
{
/*
    based on alg 462 comm. acm oct 73
    gives the probability that a bivariate normal exceeds (ah,ak).
    gh and gk are .5 times the right tail areas of ah, ak under a n(0,1)

    Tranlated from FORTRAN to ratfor using struct; from ratfor to C by hand.
*/
	double a2, ap, b, cn, conex, ex, g2, gh, gk, gw, h2, h4, rr, s1, s2, 
	sgn, sn, sp, sqr, t, temp, w2, wh, wk;
	int is;

	temp = -ah;
	gh = pnorm1(temp);
	gh = gh / 2.0;
	temp = -ak;
	gk = pnorm1(temp);
	gk = gk / 2.0;

	b = 0;

	if (r==0)
		b = 4*gh*gk;
	else {
		rr = 1-r*r;
		if (rr<0)
			return 0;  /* zz; 29/6/02; was originally return; not sure */
		if (rr!=0) {
			sqr = sqrt(rr);
			if (ah!=0) {
				b = gh;
				if (ah*ak<0)
					b = b-.5;
				else if (ah*ak==0)
					goto label10;
			}
			else if (ak==0) {
				b = atan(r/sqr)/twopi+.25;
				goto label50;
			}
			b = b+gk;
			if (ah==0)
				goto label20;
label10:
			wh = -ah;
			wk = (ak/ah-r)/sqr;
			gw = 2*gh;
			is = -1;
			goto label30;
label20:
			do {
				wh = -ak;
				wk = (ah/ak-r)/sqr;
				gw = 2*gk;
				is = 1;
label30:
				sgn = -1;
				t = 0;
				if (wk!=0) {
					if (fabs(wk)>=1) {
                                         /* this brace added 28/6/02 by tyee */
						if (fabs(wk)==1) {
							t = wk*gw*(1-gw)/2;
							goto label40;
						}
						else {
							sgn = -sgn;
							wh = wh*wk;
							g2 = pnorm1(wh);
							wk = 1/wk;
							if (wk<0)
								b = b+.5;
							b = b-(gw+g2)/2+gw*g2;
						}
					}
					h2 = wh*wh;
					a2 = wk*wk;
					h4 = h2*.5;
					ex = 0;
					if (h4<150.0)
						ex = exp(-h4);
					w2 = h4*ex;
					ap = 1;
					s2 = ap-ex;
					sp = ap;
					s1 = 0;
					sn = s1;
					conex = fabs(con/wk);
					do {
						cn = ap*s2/(sn+sp);
						s1 = s1+cn;
						if (fabs(cn)<=conex)
							break;
						sn = sp;
						sp = sp+1;
						s2 = s2-w2;
						w2 = w2*h4/sp;
						ap = -ap*a2;
					} while (1);
					t = (atan(wk)-wk*s1)/twopi;
label40:
					b = b+sgn*t;
				}
				if (is>=0)
					break;
			} while(ak!=0);
		}
		else if (r>=0)
			if (ah>=ak)
				b = 2*gh;
			else
				b = 2*gk;
		else if (ah+ak<0)
			b = 2*(gh+gk)-1;
	}
label50:
	if (b<0)
		b = 0;
	if (b>1)
		b = 1;

        return(b);  
}

/* in the following function
size measures the dimension of x
singler == 1 if r is a scalar; otherwise r is same size as x & y
*/

/* This is called by S */

void pnorm2(double *x, double *y, double *r, 
            int *size, int *singler, double *ans)
{
    int i;

    if(*singler == 1)
    {
        for(i = 0; i < *size; i++)
            ans[i] = bivnor(x[i], y[i], *r);
    }
    else
    {
        for(i = 0; i < *size; i++)
            ans[i] = bivnor(x[i], y[i], r[i]);
    }
}



/*
main()
{
    int i;
    double x,y,r;

    x = 0.0;
    y = 0.0;

    for(i = -9; i<=9; i++)
    {
        r = i / 10.0;
        Rprintf("%10.2f   %10.6f  \n",r,bivnor(x,y,r));

    }

}
*/





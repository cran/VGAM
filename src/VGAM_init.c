#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void a2mccc(void *, void *, void *, void *, void *, void *, void *);
extern void cqo_1(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void cqo_2(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void dcqo1(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void eimpnbinomspecialp(void *, void *, void *, void *, void *, void *);
extern void lerchphi123(void *, void *, void *, void *, void *, void *, void *, void *);
extern void m2accc(void *, void *, void *, void *, void *, void *, void *, void *);
extern void mux111ccc(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void mux111ddd(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void mux15ccc(void *, void *, void *, void *, void *);
extern void mux2ccc(void *, void *, void *, void *, void *, void *);
extern void mux22ccc(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void mux5ccc(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void mux55ccc(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void mux7ccc(void *, void *, void *, void *, void *, void *, void *);
extern void pnorm2ccc(void *, void *, void *, void *, void *, void *);
extern void sf_C_expexpint(void *, void *, void *);
extern void sf_C_expint(void *, void *, void *);
extern void sf_C_expint_e1(void *, void *, void *);
extern void tapply_mat1(void *, void *, void *, void *);
extern void tyee_C_cum8sum(void *, void *, void *, void *, void *, void *);
extern void vbacksubccc(void *, void *, void *, void *, void *, void *, void *, void *);
extern void vcao6(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void vcholccc(void *, void *, void *, void *, void *, void *, void *, void *);
extern void vdcao6(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void vforsubccc(void *, void *, void *, void *, void *, void *, void *, void *);
extern void VGAM_C_kend_tau(void *, void *, void *, void *);
extern void VGAM_C_mux34(void *, void *, void *, void *, void *, void *);
extern void VGAM_C_vdigami(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void vknootl2(void *, void *, void *, void *, void *);
extern void vsuff9(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void vzetawr(void *, void *, void *, void *);
extern void Yee_pknootl2(void *, void *, void *, void *);
extern void Yee_spline(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void Yee_vbfa(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void Yee_vbvs(void *, void *, void *, void *, void *, void *, void *, void *);

/* .Fortran calls */
extern void F77_NAME(veigenf)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(yjngintf)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"a2mccc",                (DL_FUNC) &a2mccc,                 7},
    {"cqo_1",              (DL_FUNC) &cqo_1,              24},
    {"cqo_2",              (DL_FUNC) &cqo_2,              24},
    {"dcqo1",              (DL_FUNC) &dcqo1,              29},
    {"eimpnbinomspecialp", (DL_FUNC) &eimpnbinomspecialp,  6},
    {"lerchphi123",        (DL_FUNC) &lerchphi123,         8},
    {"m2accc",                (DL_FUNC) &m2accc,                 8},
    {"mux111ccc",             (DL_FUNC) &mux111ccc,             11},
    {"mux111ddd",             (DL_FUNC) &mux111ddd,             12},
    {"mux15ccc",              (DL_FUNC) &mux15ccc,               5},
    {"mux2ccc",               (DL_FUNC) &mux2ccc,                6},
    {"mux22ccc",              (DL_FUNC) &mux22ccc,              10},
    {"mux5ccc",               (DL_FUNC) &mux5ccc,               16},
    {"mux55ccc",              (DL_FUNC) &mux55ccc,               9},
    {"mux7ccc",               (DL_FUNC) &mux7ccc,                7},
    {"pnorm2ccc",             (DL_FUNC) &pnorm2ccc,              6},
    {"sf_C_expexpint",     (DL_FUNC) &sf_C_expexpint,      3},
    {"sf_C_expint",        (DL_FUNC) &sf_C_expint,         3},
    {"sf_C_expint_e1",     (DL_FUNC) &sf_C_expint_e1,      3},
    {"tapply_mat1",        (DL_FUNC) &tapply_mat1,         4},
    {"tyee_C_cum8sum",     (DL_FUNC) &tyee_C_cum8sum,      6},
    {"vbacksubccc",           (DL_FUNC) &vbacksubccc,            8},
    {"vcao6",              (DL_FUNC) &vcao6,              42},
    {"vcholccc",              (DL_FUNC) &vcholccc,               8},
    {"vdcao6",             (DL_FUNC) &vdcao6,             47},
    {"vforsubccc",            (DL_FUNC) &vforsubccc,             8},
    {"VGAM_C_kend_tau",    (DL_FUNC) &VGAM_C_kend_tau,     4},
    {"VGAM_C_mux34",       (DL_FUNC) &VGAM_C_mux34,        6},
    {"VGAM_C_vdigami",     (DL_FUNC) &VGAM_C_vdigami,     12},
    {"vknootl2",           (DL_FUNC) &vknootl2,            5},
    {"vsuff9",             (DL_FUNC) &vsuff9,             21},
    {"vzetawr",            (DL_FUNC) &vzetawr,             4},
    {"Yee_pknootl2",       (DL_FUNC) &Yee_pknootl2,        4},
    {"Yee_spline",         (DL_FUNC) &Yee_spline,         28},
    {"Yee_vbfa",           (DL_FUNC) &Yee_vbfa,           30},
    {"Yee_vbvs",           (DL_FUNC) &Yee_vbvs,            8},
    {NULL, NULL, 0}
};

static const R_FortranMethodDef FortranEntries[] = {
    {"veigenf",         (DL_FUNC) &F77_NAME(veigenf),         13},
    {"yjngintf",       (DL_FUNC) &F77_NAME(yjngintf),       11},
    {NULL, NULL, 0}
};

void R_init_VGAM(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}





#ifndef TRICUBIC_H
#define TRICUBIC_H
extern const char *tricubic_version(void);
extern void tricubic_get_coeff(double a[64], double f[8], double dfdx[8], double dfdy[8], double dfdz[8], double d2fdxdy[8], double d2fdxdz[8], double d2fdydz[8], double d3fdxdydz[8]);
extern double tricubic_eval(double a[64], double x, double y, double z);
extern double tricubic_eval(double a[64], double x, double y, double z, int derx, int dery, int derz);

extern void tricubic_pointID2xyz(int id, int *x, int *y, int *z);
extern void tricubic_pointID2xyz(int id, double *x, double *y, double *z);
#endif

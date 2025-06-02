#ifndef INTERPOLATION_H
#define INTERPOLATION_H

#include "zheevd.h"

class Interpolate{
public:
   Interpolate(int, int, int, int, doublecomplex **);
   ~Interpolate();
 
   void set_method();
   void execute(double *, doublecomplex *);
   void reset_gamma();
 
   int UseGamma;

   class UserInput *input;

private:
   void tricubic_init();
   void tricubic(double *, doublecomplex *);
   void trilinear(double *, doublecomplex *);
   class Memory *memory;
 
   int which;
   int Nx, Ny, Nz, Npt, ndim;
   int flag_reset_gamma, flag_allocated_dfs;
 
   doublecomplex **data;
   doublecomplex **Dfdx, **Dfdy, **Dfdz, **D2fdxdy, **D2fdxdz, **D2fdydz, **D3fdxdydz;
   double a[64], f[8], dfdx[8], dfdy[8], dfdz[8], d2fdxdy[8], d2fdxdz[8], d2fdydz[8], d3fdxdydz[8];
   int vidx[8];
};

#endif

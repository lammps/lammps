#ifndef LMP_MANIFOLD_SPHEROCYLINDER_H
#define LMP_MANIFOLD_SPHEROCYLINDER_H

#include "manifold.h"

namespace LAMMPS_NS {

namespace user_manifold {


  // A sphere:
  class manifold_spherocylinder : public manifold {
   public:
    enum { NPARAMS = 2 };
    manifold_spherocylinder( LAMMPS *lmp, int, char ** ) : manifold(lmp){}

    virtual ~manifold_spherocylinder(){}
    virtual double g( const double *x )
    {
      double xy2 = x[0]*x[0] + x[1]*x[1];

      if( fabs(x[2]) < Lh ){
        // Cylinder part:
        return xy2 - R2;
      }else{
        // Sphere part:

        return xy2 + x[2]*x[2] - R2;
      }
    }

    virtual double g_and_n( const double *x, double *nn )
    {
      double xy2 = x[0]*x[0] + x[1]*x[1];

      nn[0] = -2*x[0];
      nn[1] = -2*x[1];
      nn[2] = 0.0;
      if( fabs(x[2]) > Lh ){
        nn[2] = -2*x[2];
        return xy2 + x[2]*x[2] - R2;
      }else{
        return xy2 - R2;
      }
    }

    virtual void   n( const double *x, double *nn )
    {
      double xy2 = x[0]*x[0] + x[1]*x[1];

      nn[0] = -2*x[0];
      nn[1] = -2*x[1];
      nn[2] = 0.0;
      if( fabs(x[2]) > Lh ){
        nn[2] = -2*x[2];
      }
    }

    virtual void post_param_init()
    {
      double R = params[0];
      double L = params[1];
      R2 = R*R;
      Lh = 0.5*L;
    }

    static const char* type(){ return "spherocylinder"; }
    virtual const char *id(){ return type(); }
    static int expected_argc(){ return NPARAMS; }
    virtual int nparams(){ return NPARAMS; }

   private:
    double R2, Lh;
  };
}

}


#endif // LMP_MANIFOLD_SPHEROCYLINDER_H

#ifndef LMP_MANIFOLD_ELLIPSOID_H
#define LMP_MANIFOLD_ELLIPSOID_H

#include "manifold.h"


namespace LAMMPS_NS {

namespace user_manifold {
  // An ellipsoid:
  class manifold_ellipsoid : public manifold {
   public:
    enum { NPARAMS = 3 };
    manifold_ellipsoid( LAMMPS *lmp, int, char ** );
    virtual ~manifold_ellipsoid(){}
    virtual double g( const double *x );
    virtual void   n( const double *x, double *n );

    static const char* type(){ return "ellipsoid"; }
    virtual const char *id(){ return type(); }
    static int expected_argc(){ return NPARAMS; }
    virtual int nparams(){ return NPARAMS; }
  };
}

}

#endif // LMP_MANIFOLD_ELLIPSOID_H

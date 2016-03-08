#ifndef LMP_MANIFOLD_PLANE_H
#define LMP_MANIFOLD_PLANE_H

#include "manifold.h"


namespace LAMMPS_NS {

namespace user_manifold {


  // A 2D plane
  class manifold_plane : public manifold {
   public:
    enum { NPARAMS = 6 }; // Number of parameters.
    manifold_plane( LAMMPS *lmp, int, char ** );
    virtual ~manifold_plane(){}
    virtual double g( const double *x );
    virtual void   n( const double *x, double *n );
    static const char *type(){ return "plane"; }
    virtual const char *id(){ return type(); }
    static int expected_argc(){ return NPARAMS; }
    virtual int nparams(){ return NPARAMS; }
  };
}


}

#endif // LMP_MANIFOLD_PLANE_H

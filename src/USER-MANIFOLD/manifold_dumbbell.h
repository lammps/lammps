#ifndef LMP_MANIFOLD_DUMBBELL_H
#define LMP_MANIFOLD_DUMBBELL_H

#include "manifold.h"



namespace LAMMPS_NS {

namespace user_manifold {

  // A dendritic dumbbell approximation:
  class manifold_dumbbell : public manifold {
   public:
    enum { NPARAMS = 4 }; // Number of parameters.
    manifold_dumbbell( LAMMPS *lmp, int, char ** );
    virtual ~manifold_dumbbell(){}
    virtual double g      ( const double *x );
    virtual void   n      ( const double *x, double *nn );

    static const char* type(){ return "dumbbell"; }
    virtual const char *id(){ return type(); }

    static int expected_argc(){ return NPARAMS; }
    virtual int nparams(){ return NPARAMS; }
  };
}

}

#endif // LMP_MANIFOLD_DUMBBELL_H

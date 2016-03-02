#ifndef LMP_MANIFOLD_SPINE_H
#define LMP_MANIFOLD_SPINE_H

#include "manifold.h"



namespace LAMMPS_NS {
  // A dendritic spine approximation:
  class manifold_spine : public manifold, protected Pointers {
   public:
    enum { NPARAMS = 5 }; // Number of parameters.
    manifold_spine( LAMMPS *lmp, int, char ** );
    virtual ~manifold_spine(){}
    virtual double g      ( const double *x );
    virtual void   n      ( const double *x, double *nn );

    static const char* type(){ return "spine"; }
    virtual const char *id(){ return type(); }

    static int expected_argc(){ return NPARAMS; }
    virtual int nparams(){ return NPARAMS; }
  };
}



#endif // LMP_MANIFOLD_SPINE_H

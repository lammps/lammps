#ifndef LMP_MANIFOLD_SPINE_H
#define LMP_MANIFOLD_SPINE_H

#include "manifold.h"


namespace LAMMPS_NS {

namespace user_manifold {

  // A dendritic spine approximation:
  class manifold_spine : public manifold {
   public:
    enum { NPARAMS = 5 }; // Number of parameters.
    manifold_spine( LAMMPS *lmp, int, char ** );
    virtual ~manifold_spine(){}
    virtual double g      ( const double *x );
    virtual void   n      ( const double *x, double *nn );
    virtual double g_and_n( const double *x, double *nn );

    static const char* type(){ return "spine"; }
    virtual const char *id(){ return type(); }

    static int expected_argc(){ return NPARAMS; }
    virtual int nparams(){ return NPARAMS; }
   protected:
    int power;
  };

  class manifold_spine_two : public manifold_spine {
   public:
    manifold_spine_two( LAMMPS *lmp, int, char **);

    static const char* type(){ return "spine/two"; }
    virtual const char *id(){ return type(); }

  };
}

}


#endif // LMP_MANIFOLD_SPINE_H

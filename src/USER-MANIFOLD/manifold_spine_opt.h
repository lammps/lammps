#ifndef LMP_MANIFOLD_SPINE_OPT_H
#define LMP_MANIFOLD_SPINE_OPT_H

#include "manifold.h"
#include "sincos_lookup.h"


namespace LAMMPS_NS {
  // A dendritic spine approximation:
  class manifold_spine_opt : public manifold, protected Pointers {
   public:
    enum { NPARAMS = 5 }; // Number of parameters.
    manifold_spine_opt( LAMMPS *lmp, int, char ** );
    virtual ~manifold_spine_opt();
    virtual double g      ( const double *x );
    virtual void   n      ( const double *x, double *nn  );
    virtual double g_and_n( const double *x, double *nn );


    static const char* type(){ return "spine/opt"; }
    virtual const char *id(){ return type(); }

    static int expected_argc(){ return NPARAMS; }
    virtual int nparams(){ return NPARAMS; }
	    
    virtual void post_param_init();
   private:
    double c2, inv_c2, a2, A4;
    
    sincos_lookup lu_table;
    
  };
}

#endif // LMP_MANIFOLD_SPINE_OPT_H

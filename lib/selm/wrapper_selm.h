/* ----------------------------------------------------------------------

  Stochastic Eulerian Lagrangian Methods (SELMs) Package

  Paul J. Atzberger
  http://atzberger.org/
  
  Please cite the follow paper when referencing this package
  
  "Fluctuating Hydrodynamics Methods for Dynamic Coarse-Grained Implicit-Solvent Simulations in LAMMPS," 
  Wang, Y. and Sigurdsson, J. K. and Atzberger, P. J., SIAM Journal on Scientific Computing, 38(5), 2016.
  
  @article{atz_selm_lammps_fluct_hydro,
    title = {Fluctuating Hydrodynamics Methods for Dynamic 
    Coarse-Grained Implicit-Solvent Simulations in LAMMPS},
    author = {Wang, Y. and Sigurdsson, J. K. and Atzberger, P. J.},
    journal = {SIAM Journal on Scientific Computing},
    volume = {38},
    number = {5},
    pages = {S62-S77},
    year = {2016},
    doi = {10.1137/15M1026390},
    URL = {https://doi.org/10.1137/15M1026390},
  }  
    
  For latest releases, examples, and additional information see 
  http://mango-selm.org/
 
------------------------------------------------------------------------- */

#ifndef LMP_WRAPPER_SELM_H
#define LMP_WRAPPER_SELM_H

//#ifndef FFT_FFTW3 /* determine if FFTW package specified */
//  #warning "No FFTW package specified for SELM codes.  The SELM functionality will be disabled."
//#else
//  #define USE_PACKAGE_FFTW3
//#endif

/*
#include "Atz_XML_Package.h"
#include "SELM_CouplingOperator.h"
#include "SELM_Eulerian.h"
#include "SELM_Integrator.h"
#include "SELM_Interaction.h"
#include "SELM_Lagrangian.h"
*/

// forward declaration
/*
class SELM_Eulerian;
class SELM_Lagrangian;
class SELM_CouplingOperator;
class SELM_Interaction;
class SELM_Integrator;
*/


//#ifdef USE_PACKAGE_FFTW3
//  #include <fftw3.h>
//#endif

//#include "fix_selm.h"

//#include "mpi.h"
//#include "fix.h"
#include "lammps.h"
//#include "random_mars.h"

#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cstddef>

using namespace std;

namespace LAMMPS_NS {

// forward declaration 
class DriverSELM;
class FixSELM;
class LAMMPS;

struct WrapperSELM {

 public:

  WrapperSELM(class FixSELM *, class LAMMPS *, int, char **);
  WrapperSELM();
  ~WrapperSELM();

  int          setmask();
  virtual void init();
  void         setup(int vflag);
  virtual void initial_integrate(int);
  virtual void final_integrate();
  void         reset_dt();
  void         post_force(int vflag);
  
  // additional methods/hooks
  void pre_exchange();
  void end_of_step();

  void init_from_fix();

  DriverSELM *driver_selm;
  
};

}

#endif


/* ----------------------------------------------------------------------

  Stochastic Eulerian Lagrangian Methods (SELMs) Package (library version)

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
    
  For latest version of the codes, examples, and additional information see 
  http://mango-selm.org/
 
------------------------------------------------------------------------- 
*/

//#include <cstdlib>
//#include <cstddef>
//#include <cstdio>
//#include <cstring>
//#include <iostream>
//#include <sstream>

/* LAMMPS includes */
//#include "atom.h"
//#include "force.h"
//#include "update.h"
//#include "respa.h"
#include "error.h"
//#include "comm.h"
//#include "universe.h"
//#include "version.h" 
//#include "random_mars.h"
#include "citeme.h"
#include "lammps.h"

/* SELM_includes */
#include "fix_selm.h"
#include "wrapper/wrapper_selm.h"

//using namespace SELM;
using namespace LAMMPS_NS;
using namespace FixConst;
using namespace std;

static const char cite_selm_str[] =
  "USER-SELM Package: Fluctuating Hydrodynamics \n\n"
  "@article{atz_selm_lammps_fluct_hydro,\n"
  "title = {Fluctuating Hydrodynamics Methods for Dynamic\n"
  "Coarse-Grained Implicit-Solvent Simulations in LAMMPS},\n"
  "author = {Wang, Y. and Sigurdsson, J. K. and Atzberger, P. J.},\n"
  "journal = {SIAM Journal on Scientific Computing},\n"
  "volume = {38},\n"
  "number = {5},\n"
  "pages = {S62-S77},\n"
  "year = {2016},\n"
  "doi = {10.1137/15M1026390},\n"
  "URL = {https://doi.org/10.1137/15M1026390},\n"
  "}\n\n";

/* =========================== Class definitions =========================== */
FixSELM::FixSELM() : Fix(NULL, 0, NULL) {
  wrapper_selm = new WrapperSELM(); 
}

/* =========================== Class definitions =========================== */
FixSELM::FixSELM(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg) 
{  

  /* add to citation collection */
  if (lmp->citeme) lmp->citeme->add(cite_selm_str); 
 
  /* set to 1 for fix performing integration, 0 if fix does not */
  time_integrate = 1;
  
  /* setup the wrapper for SELM library */
  wrapper_selm = new WrapperSELM(this,lmp,narg,arg);  

}

/* destructor */
FixSELM::~FixSELM() {
  delete wrapper_selm;
}


/* ---------------------------------------------------------------------- */
void FixSELM::setup(int vflag) { wrapper_selm->setup(vflag); }

/* ---------------------------------------------------------------------- */
int FixSELM::setmask()
{
  /*
    mask |= INITIAL_INTEGRATE;
    mask |= FINAL_INTEGRATE;
   */

  // pass value back from the wrapper
  SELM_integrator_mask = wrapper_selm->setmask();

  return SELM_integrator_mask; /* syncronize the SELM mask with that returned to LAMMPS */
}

/* ---------------------------------------------------------------------- */
void FixSELM::pre_exchange() { wrapper_selm->pre_exchange(); }

/* ---------------------------------------------------------------------- */
void FixSELM::end_of_step() { wrapper_selm->end_of_step(); }

/* ---------------------------------------------------------------------- */
void FixSELM::init() { wrapper_selm->init_from_fix(); }

/* ---------------------------------------------------------------------- */
void FixSELM::initial_integrate(int vflag) { wrapper_selm->initial_integrate(vflag); }

/* ---------------------------------------------------------------------- */
void FixSELM::final_integrate() { wrapper_selm->final_integrate(); }

/* ---------------------------------------------------------------------- */
void FixSELM::reset_dt() { wrapper_selm->reset_dt(); }

/* ---------------------------------------------------------------------- */
void FixSELM::post_force(int vflag) { wrapper_selm->post_force(vflag); }

/*****************************************************************************************/
/* Supporting SELM codes */
/*****************************************************************************************/
void FixSELM::packageError(int code, void *extras) {
  std::exit(code); // exit, may want also to notify lammps 
}


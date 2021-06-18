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
    
  For latest releases, examples, and additional information see 
  http://mango-selm.org/
 
------------------------------------------------------------------------- 
*/

/* SELM_includes */
#include "wrapper_selm.h"
#include "driver_selm.h"

#include "fix_selm.h"
#include "lammps.h"

//#include "driver_selm.cpp"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */
WrapperSELM::WrapperSELM() { driver_selm = new DriverSELM(); }

/* ---------------------------------------------------------------------- */
WrapperSELM::WrapperSELM(FixSELM *fixSELM, LAMMPS *lmp, int narg, char **arg) { driver_selm = new DriverSELM(fixSELM,lmp,narg,arg); }
     
/* ---------------------------------------------------------------------- */
WrapperSELM::~WrapperSELM() { delete driver_selm; }

/* ---------------------------------------------------------------------- */
void WrapperSELM::setup(int vflag) { driver_selm->setup(vflag); }

/* ---------------------------------------------------------------------- */
int WrapperSELM::setmask() { return driver_selm->setmask(); }

/* ---------------------------------------------------------------------- */
void WrapperSELM::pre_exchange() { driver_selm->pre_exchange(); }

/* ---------------------------------------------------------------------- */
void WrapperSELM::end_of_step() { driver_selm->end_of_step(); }

/* ---------------------------------------------------------------------- */
void WrapperSELM::init() { driver_selm->init(); }

/* ---------------------------------------------------------------------- */
void WrapperSELM::initial_integrate(int vflag) { driver_selm->initial_integrate(vflag); }

/* ---------------------------------------------------------------------- */
void WrapperSELM::final_integrate() { driver_selm->final_integrate(); }

/* ---------------------------------------------------------------------- */
void WrapperSELM::reset_dt() { driver_selm->reset_dt(); }

/* ---------------------------------------------------------------------- */
void WrapperSELM::post_force(int vflag) { driver_selm->post_force(vflag); }

/* ---------------------------------------------------------------------- */
void WrapperSELM::init_from_fix() { driver_selm->init_from_fix(); }

/* ---------------------------------------------------------------------- */


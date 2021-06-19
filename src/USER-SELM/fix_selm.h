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

  For latest version of the codes, examples, and additional information see
  http://mango-selm.org/

------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(selm,FixSELM)

#else

#ifndef LMP_FIX_SELM_H
#define LMP_FIX_SELM_H

//#include "wrapper_selm.h"

#include "fix.h"
#include "lammps.h"

#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cstddef>

using namespace std;

namespace LAMMPS_NS {

class WrapperSELM;

class FixSELM : public Fix {

 public:

  int SELM_integrator_mask;

  /* ================= Function prototypes ================= */
  FixSELM(class LAMMPS *, int, char **);
  FixSELM();
  ~FixSELM();

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

  /* =========================== Function Calls =========================== */
  void packageError(int code, void *extras);

  /* =========================== Variables =========================== */
  LAMMPS                                                    *lammps;
  WrapperSELM                                               *wrapper_selm;

};

}

#endif
#endif



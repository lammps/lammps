/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef KSPACE_CLASS

KSpaceStyle(pppm/tip4p/omp,PPPMTIP4POMP)

#else

#ifndef LMP_PPPM_TIP4P_OMP_H
#define LMP_PPPM_TIP4P_OMP_H

#include "pppm_omp.h"

namespace LAMMPS_NS {

class PPPMTIP4POMP : public PPPMOMP {
 public:
  PPPMTIP4POMP(class LAMMPS *, int, char **);
  virtual ~PPPMTIP4POMP () {};
  virtual void init();
    
 protected:
  virtual void particle_map();
  virtual void fieldforce();
  virtual void make_rho();

 private:
  void find_M(int, int &, int &, double *); 

//  void slabcorr(int);

};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Kspace style pppm/tip4p/omp requires newton on

Self-explanatory.

*/

/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef INTEGRATE_CLASS

IntegrateStyle(verlet,Verlet)

#else

#ifndef LMP_VERLET_H
#define LMP_VERLET_H

#include "integrate.h"

namespace LAMMPS_NS {

class Verlet : public Integrate {
 public:
  Verlet(class LAMMPS *, int, char **);
  virtual ~Verlet() {}
  virtual void init();
  virtual void setup();
  virtual void setup_minimal(int);
  virtual void run(int);
  void cleanup();

 protected:
  int triclinic;                    // 0 if domain is orthog, 1 if triclinic
  int torqueflag,extraflag;

  virtual void force_clear();
};

}

#endif
#endif

/* ERROR/WARNING messages:

W: No fixes defined, atoms won't move

If you are not using a fix like nve, nvt, npt then atom velocities and
coordinates will not be updated during timestepping.

*/

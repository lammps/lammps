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

#ifdef FIX_CLASS

FixStyle(press/berendsen,FixPressBerendsen)

#else

#ifndef LMP_FIX_PRESS_BERENDSEN_H
#define LMP_FIX_PRESS_BERENDSEN_H

#include "fix.h"

namespace LAMMPS_NS {

class FixPressBerendsen : public Fix {
 public:
  FixPressBerendsen(class LAMMPS *, int, char **);
  ~FixPressBerendsen();
  int setmask();
  void init();
  void setup(int);
  void end_of_step();
  int modify_param(int, char **);

 protected:
  int dimension,which;
  double bulkmodulus;

  int pstyle,pcouple,allremap;
  int p_flag[3];                   // 1 if control P on this dim, 0 if not
  double p_start[3],p_stop[3];
  double p_period[3],p_target[3];
  double p_current[3],dilation[3];
  double factor[3];
  int kspace_flag;                 // 1 if KSpace invoked, 0 if not
  int nrigid;                      // number of rigid fixes
  int *rfix;                       // indices of rigid fixes

  char *id_temp,*id_press;
  class Compute *temperature,*pressure;
  int tflag,pflag;

  void couple();
  void remap();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Invalid fix press/berendsen for a 2d simulation

The z component of pressure cannot be controlled for a 2d model.

E: Invalid fix press/berendsen pressure settings

Settings for coupled dimensions must be the same.

E: Cannot use fix press/berendsen on a non-periodic dimension

Self-explanatory.

E: Fix press/berendsen damping parameters must be > 0.0

Self-explanatory.

E: Cannot use fix press/berendsen with triclinic box

Self-explanatory.

E: Cannot use fix press/berendsen and fix deform on same component of stress tensor

These commands both change the box size/shape, so you cannot use both
together.

E: Temperature ID for fix press/berendsen does not exist

Self-explanatory.

E: Pressure ID for fix press/berendsen does not exist

The compute ID needed to compute pressure for the fix does not
exist.

E: Could not find fix_modify temperature ID

The compute ID for computing temperature does not exist.

E: Fix_modify temperature ID does not compute temperature

The compute ID assigned to the fix must compute temperature.

W: Temperature for NPT is not for group all

User-assigned temperature to NPT fix does not compute temperature for
all atoms.  Since NPT computes a global pressure, the kinetic energy
contribution from the temperature is assumed to also be for all atoms.
Thus the pressure used by NPT could be inaccurate.

E: Could not find fix_modify pressure ID

The compute ID for computing pressure does not exist.

E: Fix_modify pressure ID does not compute pressure

The compute ID assigned to the fix must compute pressure.

*/

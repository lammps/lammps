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

FixStyle(pour,FixPour)

#else

#ifndef LMP_FIX_POUR_H
#define LMP_FIX_POUR_H

#include "fix.h"

namespace LAMMPS_NS {

class FixPour : public Fix {
  friend class PairGranHertzHistory;
  friend class PairGranHertzHistoryOMP;
  friend class PairGranHooke;
  friend class PairGranHookeOMP;
  friend class PairGranHookeHistory;
  friend class PairGranHookeHistoryOMP;
  friend class PairGranHookeCuda;

 public:
  FixPour(class LAMMPS *, int, char **);
  ~FixPour();
  int setmask();
  void init();
  void pre_exchange();
  void reset_dt();

 private:
  int ninsert,ntype,seed;
  double radius_lo,radius_hi;
  double density_lo,density_hi;
  double volfrac;
  int maxattempt;
  int region_style;
  double rate;
  double vxlo,vxhi,vylo,vyhi,vy,vz;
  double xlo,xhi,ylo,yhi,zlo,zhi;
  double xc,yc,rc;
  double grav;

  int me,nprocs;
  int *recvcounts,*displs;
  int nfreq,nfirst,ninserted,nper;
  double lo_current,hi_current;
  class FixShearHistory *fix_history;
  class RanPark *random;

  int overlap(int);
  void xyz_random(double, double *);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Fix pour requires atom attributes radius, rmass

The atom style defined does not have these attributes.

E: Fix pour region ID does not exist

Self-explanatory.

E: Must specify a region in fix pour

The region keyword must be specified with this fix.

E: Fix pour region does not support a bounding box

Not all regions represent bounded volumes.  You cannot use
such a region with the fix pour command.

E: Fix pour region cannot be dynamic

Only static regions can be used with fix pour.

E: Insertion region extends outside simulation box

Region specified with fix pour command extends outside the global
simulation box.

E: Must use a z-axis cylinder with fix pour

The axis of the cylinder region used with the fix pour command must
be oriented along the z dimension.

E: Must use a block or cylinder region with fix pour

Self-explanatory.

E: Must use a block region with fix pour for 2d simulations

Self-explanatory.

E: No fix gravity defined for fix pour

Cannot add poured particles without gravity to move them.

E: Cannot use fix pour with triclinic box

This feature is not yet supported.

E: Gravity must point in -z to use with fix pour in 3d

Gravity must be pointing "down" in a 3d box, i.e. theta = 180.0.

E: Gravity must point in -y to use with fix pour in 2d

Gravity must be pointing "down" in a 2d box.

E: Gravity changed since fix pour was created

Gravity must be static and not dynamic for use with fix pour.

W: Less insertions than requested

Less atom insertions occurred on this timestep due to the fix pour
command than were scheduled.  This is probably because there were too
many overlaps detected.

E: Cannot change timestep with fix pour

This fix pre-computes some values based on the timestep, so it cannot
be changed during a simulation run.

*/

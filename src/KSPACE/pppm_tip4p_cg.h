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

KSpaceStyle(pppm/tip4p/cg,PPPMTIP4PCG)

#else

#ifndef LMP_PPPM_TIP4P_CG_H
#define LMP_PPPM_TIP4P_CG_H

#include "pppm_tip4p.h"

namespace LAMMPS_NS {

class PPPMTIP4PCG : public PPPMTIP4P {
 public:
  PPPMTIP4PCG(class LAMMPS *, int, char **);
  virtual ~PPPMTIP4PCG ();
  virtual void compute(int,int);

 protected:
  int num_charged;
  int *is_charged;
  double smallq;

 protected:
  virtual void particle_map();
  virtual void make_rho();
  virtual void fieldforce_ik();
  virtual void fieldforce_ad();
  virtual void fieldforce_peratom();
  void find_M_cg(int, int &, int &, double *);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Kspace style pppm/tip4p requires newton on

Self-explanatory.

E: Out of range atoms - cannot compute PPPM

One or more atoms are attempting to map their charge to a PPPM grid
point that is not owned by a processor.  This is likely for one of two
reasons, both of them bad.  First, it may mean that an atom near the
boundary of a processor's sub-domain has moved more than 1/2 the
"neighbor skin distance"_neighbor.html without neighbor lists being
rebuilt and atoms being migrated to new processors.  This also means
you may be missing pairwise interactions that need to be computed.
The solution is to change the re-neighboring criteria via the
"neigh_modify"_neigh_modify command.  The safest settings are "delay 0
every 1 check yes".  Second, it may mean that an atom has moved far
outside a processor's sub-domain or even the entire simulation box.
This indicates bad physics, e.g. due to highly overlapping atoms, too
large a timestep, etc.

E: TIP4P hydrogen is missing

The TIP4P pairwise computation failed to find the correct H atom
within a water molecule.

E: TIP4P hydrogen has incorrect atom type

The TIP4P pairwise computation found an H atom whose type does not
agree with the specified H type.

*/

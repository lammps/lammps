/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef KSPACE_CLASS
// clang-format off
KSpaceStyle(pppm/disp/tip4p,PPPMDispTIP4P);
// clang-format on
#else

#ifndef LMP_PPPM_DISP_TIP4P_H
#define LMP_PPPM_DISP_TIP4P_H

#include "pppm_disp.h"

namespace LAMMPS_NS {

class PPPMDispTIP4P : public PPPMDisp {
 public:
  PPPMDispTIP4P(class LAMMPS *);
  virtual ~PPPMDispTIP4P(){};
  void init();

 protected:
  virtual void particle_map_c(double, double, double, double, int **, int, int, int, int, int, int,
                              int, int);
  virtual void make_rho_c();
  virtual void fieldforce_c_ik();
  virtual void fieldforce_c_ad();
  virtual void fieldforce_c_peratom();

 private:
  void find_M(int, int &, int &, double *);
};

}    // namespace LAMMPS_NS

#endif
#endif

/* ERROR/WARNING messages:

E: Kspace style pppm/disp/tip4p requires newton on

Self-explanatory.

E: Non-numeric box dimensions - simulation unstable

The box size has apparently blown up.

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

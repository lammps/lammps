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

KSpaceStyle(pppm/dipole/spin,PPPMDipoleSpin)

#else

#ifndef LMP_PPPM_DIPOLE_SPIN_H
#define LMP_PPPM_DIPOLE_SPIN_H

#include "pppm_dipole.h"

namespace LAMMPS_NS {

class PPPMDipoleSpin : public PPPMDipole {
 public:
  PPPMDipoleSpin(class LAMMPS *);
  virtual ~PPPMDipoleSpin();
  void init();
  void compute(int, int);

 protected:
  double hbar;                  // reduced Planck's constant      
  double mub;                   // Bohr's magneton                
  double mu_0;                  // vacuum permeability
  double mub2mu0;               // prefactor for mech force
  double mub2mu0hbinv;          // prefactor for mag force

  void slabcorr();

  // spin

  void make_rho_spin();
  void fieldforce_ik_spin();
  void fieldforce_peratom_spin();
  void spsum_spsq();

};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Cannot (yet) use charges with Kspace style PPPMDipoleSpin

Charge-spin interactions are not yet implemented in PPPMDipoleSpin so this
feature is not yet supported.

E: Must redefine kspace_style after changing to triclinic box

Self-explanatory.

E: Kspace style requires atom attribute mu

The atom style defined does not have this attribute.

E: Cannot (yet) use kspace_modify diff ad with spins

This feature is not yet supported.

E: Cannot (yet) use 'electron' units with spins

This feature is not yet supported.

E: Cannot yet use triclinic cells with PPPMDipoleSpin

This feature is not yet supported.

E: Cannot yet use TIP4P with PPPMDipoleSpin

This feature is not yet supported.

E: Cannot use nonperiodic boundaries with PPPM

For kspace style pppm, all 3 dimensions must have periodic boundaries
unless you use the kspace_modify command to define a 2d slab with a
non-periodic z dimension.

E: Incorrect boundaries with slab PPPM

Must have periodic x,y dimensions and non-periodic z dimension to use
2d slab option with PPPM.

E: PPPM order cannot be < 2 or > than %d

This is a limitation of the PPPM implementation in LAMMPS.

E: KSpace style is incompatible with Pair style

Setting a kspace style requires that a pair style with matching
long-range spin components be used.

W: Reducing PPPM order b/c stencil extends beyond nearest neighbor processor

This may lead to a larger grid than desired. See the kspace_modify overlap
command to prevent changing of the PPPM order.

E: PPPM order < minimum allowed order

The default minimum order is 2. This can be reset by the
kspace_modify minorder command.

E: PPPM grid stencil extends beyond nearest neighbor processor

This is not allowed if the kspace_modify overlap setting is no.

E: Cannot (yet) compute per-atom virial with kspace style pppm/dipole/spin

This feature is not yet supported.

E: KSpace accuracy must be > 0

The kspace accuracy designated in the input must be greater than zero.

E: Could not compute grid size

The code is unable to compute a grid size consistent with the desired
accuracy. This error should not occur for typical problems. Please
send an email to the developers.

E: PPPM grid is too large

The global PPPM grid is larger than OFFSET in one or more dimensions.
OFFSET is currently set to 4096. You likely need to decrease the
requested accuracy.

E: Could not compute g_ewald

The Newton-Raphson solver failed to converge to a good value for
g_ewald. This error should not occur for typical problems. Please
send an email to the developers.

E: Non-numeric box dimensions - simulation unstable

The box size has apparently blown up.

E: Out of range atoms - cannot compute PPPM

One or more atoms are attempting to map their charge to a PPPM grid
point that is not owned by a processor. This is likely for one of two
reasons, both of them bad. First, it may mean that an atom near the
boundary of a processor's sub-domain has moved more than 1/2 the
"neighbor skin distance"_neighbor.html without neighbor lists being
rebuilt and atoms being migrated to new processors. This also means
you may be missing pairwise interactions that need to be computed.
The solution is to change the re-neighboring criteria via the
"neigh_modify"_neigh_modify command. The safest settings are "delay 0
every 1 check yes". Second, it may mean that an atom has moved far
outside a processor's sub-domain or even the entire simulation box.
This indicates bad physics, e.g. due to highly overlapping atoms, too
large a timestep, etc.

E: Using kspace solver PPPMDipoleSpin on system with no spins

Must have non-zero spins with PPPMDipoleSpin.

E: Must use kspace_modify gewald for system with no spins

Self-explanatory.

*/

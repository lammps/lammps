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

KSpaceStyle(msm/cg/omp,MSMCGOMP)

#else

#ifndef LMP_MSM_CG_OMP_H
#define LMP_MSM_CG_OMP_H

#include "msm_omp.h"

namespace LAMMPS_NS {

class MSMCGOMP : public MSMOMP {
 public:
  MSMCGOMP(class LAMMPS *);
  virtual ~MSMCGOMP();
  virtual void settings(int, char **);
  virtual void compute(int, int);
  virtual double memory_usage();

 protected:
  int num_charged;
  int *is_charged;
  double smallq;

 protected:
  virtual void particle_map();
  virtual void make_rho();
  virtual void fieldforce();
  virtual void fieldforce_peratom();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Must use 'kspace_modify pressure/scalar no' with kspace_style msm/cg/omp

The kspace scalar pressure option is not compatible with kspace_style msm/cg/omp.

E: Cannot (yet) use MSM with triclinic box

This feature is not yet supported.

E: Cannot (yet) use MSM with 2d simulation

This feature is not yet supported.

E: Kspace style requires atom attribute q

The atom style defined does not have these attributes.

E: Cannot use slab correction with MSM

Slab correction can only be used with Ewald and PPPM, not MSM.

E: MSM order must be 4, 6, 8, or 10

This is a limitation of the MSM implementation in LAMMPS:
the MSM order can only be 4, 6, 8, or 10.

E: Cannot (yet) use single precision with MSM (remove -DFFT_SINGLE from Makefile and recompile)

Single precision cannot be used with MSM.

E: KSpace style is incompatible with Pair style

Setting a kspace style requires that a pair style with a long-range
Coulombic component be selected that is compatible with MSM.  Note
that TIP4P is not (yet) supported by MSM.

E: Cannot use kspace solver on system with no charge

No atoms in system have a non-zero charge.

E: System is not charge neutral, net charge = %g

The total charge on all atoms on the system is not 0.0, which
is not valid for MSM.

E: MSM grid is too large

The global MSM grid is larger than OFFSET in one or more dimensions.
OFFSET is currently set to 16384.  You likely need to decrease the
requested accuracy.

W: MSM mesh too small, increasing to 2 points in each direction

The global MSM grid is too small, so the number of grid points has been
increased

E: KSpace accuracy must be > 0

The kspace accuracy designated in the input must be greater than zero.

W: Number of MSM mesh points increased to be a multiple of 2

MSM requires that the number of grid points in each direction be a multiple
of two and the number of grid points in one or more directions have been
adjusted to meet this requirement.

W: Adjusting Coulombic cutoff for MSM, new cutoff = %g

The adjust/cutoff command is turned on and the Coulombic cutoff has been
adjusted to match the user-specified accuracy.

E: Out of range atoms - cannot compute MSM

One or more atoms are attempting to map their charge to a MSM grid point
that is not owned by a processor.  This is likely for one of two
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

*/

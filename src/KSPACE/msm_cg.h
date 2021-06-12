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
KSpaceStyle(msm/cg,MSMCG);
// clang-format on
#else

#ifndef LMP_MSM_CG_H
#define LMP_MSM_CG_H

#include "msm.h"

namespace LAMMPS_NS {

class MSMCG : public MSM {
 public:
  MSMCG(class LAMMPS *);
  virtual ~MSMCG();
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

}    // namespace LAMMPS_NS

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Must use 'kspace_modify pressure/scalar no' with kspace_style msm/cg

The kspace scalar pressure option is not compatible with kspace_style msm/cg.

E: Non-numeric box dimensions - simulation unstable

The box size has apparently blown up.

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

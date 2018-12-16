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

/*------------------------------------------------------------------------
  Contributing Authors : Romain Vermorel (LFCR), Laurent Joly (ULyon)
  --------------------------------------------------------------------------*/

#ifdef COMPUTE_CLASS

ComputeStyle(stress/mop/profile,ComputeStressMopProfile)

#else

#ifndef LMP_COMPUTE_STRESS_MOP_PROFILE_H
#define LMP_COMPUTE_STRESS_MOP_PROFILE_H

#include "compute.h"

namespace LAMMPS_NS {

  class ComputeStressMopProfile : public Compute {
  public:
    ComputeStressMopProfile(class LAMMPS *, int, char **);
    virtual ~ComputeStressMopProfile();
    void init();
    void init_list(int, class NeighList *);
    void compute_array();

  private:

    void compute_pairs();
    void setup_bins();

    int me,nvalues,dir;
    int *which;

    int originflag;
    double origin,delta,offset,invdelta;
    int nbins;
    double **coord,**coordp;
    double **values_local,**values_global;

    double dt,nktv2p,ftm2v;
    double area;
    class NeighList *list;

  };

}

#endif
#endif

/* ERROR/WARNING messages:

   E: Illegal ... command

   Self-explanatory.  Check the input script syntax and compare to the
   documentation for the command.  You can use -echo screen as a
   command-line option when running LAMMPS to see the offending line.

   E: Compute stress/mop/profile incompatible with simulation dimension

   Compute stress/mop/profile only works with 3D simulations.

   E: Compute stress/mop/profile incompatible with triclinic simulation box

   Self-explanatory.

   E: Compute stress/mop/profile requires a fixed simulation box

   Compute stress/mop/profile is not compatible with any change of volume or shape
   or boundary conditions of the simulation box.

   E: No pair style is defined for compute stress/mop/profile

   Self-explanatory. Compute stress/mop/profile requires the definition of a pair style.

   E: Pair style does not support compute stress/mop/profile

   The pair style does not have a single() function, so it can
   not be invoked by compute stress/mop/profile.

   E: Origin of bins for compute stress/mop/profile is out of bounds

   Self-explanatory.

   W: compute stress/mop/profile does not account for bond potentials

   W: compute stress/mop/profile does not account for angle potentials

   W: compute stress/mop/profile does not account for dihedral potentials

   W: compute stress/mop/profile does not account for improper potentials

   W: compute stress/mop/profile does not account for kspace contributions

   Compute stress/mop/profile only accounts for pairwise additive interactions for
   the computation of local stress tensor components.

*/


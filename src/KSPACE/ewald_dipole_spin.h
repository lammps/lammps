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

KSpaceStyle(ewald/dipole/spin,EwaldDipoleSpin)

#else

#ifndef LMP_EWALD_DIPOLE_SPIN_H
#define LMP_EWALD_DIPOLE_SPIN_H

#include "ewald_dipole.h"

namespace LAMMPS_NS {

class EwaldDipoleSpin : public EwaldDipole {
 public:
  //EwaldDipoleSpin(class LAMMPS *, int, char **);
  EwaldDipoleSpin(class LAMMPS *);
  virtual ~EwaldDipoleSpin();
  void init();
  void setup();
  void compute(int, int);

 protected:
  double hbar;                  // reduced Planck's constant      
  double mub;                   // Bohr's magneton                
  double mu_0;                  // vacuum permeability
  double mub2mu0;               // prefactor for mech force
  double mub2mu0hbinv;          // prefactor for mag force

  void spsum_musq(); 
  virtual void eik_dot_r();
  void slabcorr();

};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Cannot use EwaldDipoleSpin with 2d simulation

The kspace style ewald cannot be used in 2d simulations.  You can use
2d EwaldDipoleSpin in a 3d simulation; see the kspace_modify command.

E: Kspace style requires atom attribute q

The atom style defined does not have these attributes.

E: Cannot use nonperiodic boundaries with EwaldDipoleSpin

For kspace style ewald, all 3 dimensions must have periodic boundaries
unless you use the kspace_modify command to define a 2d slab with a
non-periodic z dimension.

E: Incorrect boundaries with slab EwaldDipoleSpin

Must have periodic x,y dimensions and non-periodic z dimension to use
2d slab option with EwaldDipoleSpin.

E: Cannot (yet) use EwaldDipoleSpin with triclinic box and slab correction

This feature is not yet supported.

E: KSpace style is incompatible with Pair style

Setting a kspace style requires that a pair style with matching
long-range Coulombic or dispersion components be used.

E: KSpace accuracy must be > 0

The kspace accuracy designated in the input must be greater than zero.

E: Must use 'kspace_modify gewald' for uncharged system

UNDOCUMENTED

E: Cannot (yet) use K-space slab correction with compute group/group for triclinic systems

This option is not yet supported.

*/

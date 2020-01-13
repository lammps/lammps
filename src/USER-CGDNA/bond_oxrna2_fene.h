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

#ifdef BOND_CLASS

BondStyle(oxrna2/fene,BondOxrna2Fene)

#else

#ifndef LMP_BOND_OXRNA2_FENE_H
#define LMP_BOND_OXRNA2_FENE_H

#include "bond_oxdna_fene.h"

namespace LAMMPS_NS {

class BondOxrna2Fene : public BondOxdnaFene {
 public:
  BondOxrna2Fene(class LAMMPS *);
  virtual ~BondOxrna2Fene();
  virtual void compute_interaction_sites(double *, double *, double *,
    double *);
};

}

#endif
#endif

/* ERROR/WARNING messages:

W: FENE bond too long: %ld %d %d %g

A FENE bond has stretched dangerously far.  It's interaction strength
will be truncated to attempt to prevent the bond from blowing up.

E: Bad FENE bond

Two atoms in a FENE bond have become so far apart that the bond cannot
be computed.

E: Incorrect args for bond coefficients

Self-explanatory.  Check the input script or data file.

W: Use special bonds = 0,1,1 with bond style oxrna

Most FENE models need this setting for the special_bonds command.

W: FENE bond too long: %ld %g

A FENE bond has stretched dangerously far.  It's interaction strength
will be truncated to attempt to prevent the bond from blowing up.

*/

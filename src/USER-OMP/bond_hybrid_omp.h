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

/* ----------------------------------------------------------------------
   Contributing author: Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#ifdef BOND_CLASS

BondStyle(hybrid/omp,BondHybridOMP)

#else

#ifndef LMP_BOND_HYBRID_OMP_H
#define LMP_BOND_HYBRID_OMP_H

#include "bond_hybrid.h"

namespace LAMMPS_NS {

class BondHybridOMP : public BondHybrid {

 public:
  BondHybridOMP(class LAMMPS *);
};

}

#endif
#endif

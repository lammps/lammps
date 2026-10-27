/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef LMP_ACE_IMPL_H
#define LMP_ACE_IMPL_H

#include "ace-evaluator/ace_c_basis.h"
#include "ace-evaluator/ace_recursive.h"

namespace LAMMPS_NS {
struct ACEImpl {
  ACEImpl() : basis_set(nullptr), ace(nullptr) {}
  ~ACEImpl()
  {
    delete basis_set;
    delete ace;
  }
  ACECTildeBasisSet *basis_set;
  ACERecursiveEvaluator *ace;
};
}    // namespace LAMMPS_NS
#endif

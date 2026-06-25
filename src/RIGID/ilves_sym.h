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

/* ----------------------------------------------------------------------
   ILVES-FAST: symmetric (quasi-Newton) variant, banded Cholesky/LDLT factored
   once per step.  Adapted from GROMACS 2021 ILVES (LGPL-2.1),
   src/gromacs/mdlib/ilves_sym.{cpp,h}.  See ilves_graph.h for attribution.
------------------------------------------------------------------------- */

#ifndef LMP_ILVES_SYM_H
#define LMP_ILVES_SYM_H

#include "ilves.h"

namespace LAMMPS_NS {
namespace ILVES {

class IlvesSym : public Ilves {
 public:
  IlvesSym(LAMMPS *lmp, int nbonds, const int *catom1, const int *catom2, const real *cdist,
           const real *invmass, int nthreads);

  real prepare(double **x, double **xprime) override;
  void step(double **dx) override;
};

}    // namespace ILVES
}    // namespace LAMMPS_NS

#endif

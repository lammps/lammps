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
   ILVES (full / exact-Newton variant): structurally-symmetric Jacobian,
   re-assembled and LU-factored each Newton iteration.  Adapted from GROMACS
   2021 ILVES (LGPL-2.1), src/gromacs/mdlib/ilves_asym.{cpp,h}.  See
   ilves_graph.h for attribution.
------------------------------------------------------------------------- */

#ifndef LMP_ILVES_ASYM_H
#define LMP_ILVES_ASYM_H

#include "ilves.h"

namespace LAMMPS_NS {
namespace ILVES {

class IlvesAsym : public Ilves {
 public:
  IlvesAsym(LAMMPS *lmp, int nbonds, const int *catom1, const int *catom2, const double *cdist,
            const double *invmass, int nthreads);

  double prepare(double **x, double **xprime) override;
  void step(double **dx) override;
};

}    // namespace ILVES
}    // namespace LAMMPS_NS

#endif

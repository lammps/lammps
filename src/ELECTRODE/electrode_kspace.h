/* ----------------------------------------------------------------------
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
   Contributing authors: Ludwig Ahrens-Iwers (TUHH), Shern Tee (UQ), Robert Meißner (TUHH)
------------------------------------------------------------------------- */

#ifndef LMP_ELECTRODE_KSPACE_H
#define LMP_ELECTRODE_KSPACE_H

#include "lmptype.h"

namespace LAMMPS_NS {
class ElectrodeKSpace {
 public:
  virtual void compute_vector(double *, int, int, bool) = 0;
  virtual void compute_vector_corr(double *, int, int, bool) = 0;
  virtual void compute_matrix(bigint *, double **, bool) = 0;
  virtual void compute_matrix_corr(bigint *, double **) = 0;
};
}    // namespace LAMMPS_NS

#endif

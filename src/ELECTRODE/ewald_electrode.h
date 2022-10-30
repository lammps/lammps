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

#ifdef KSPACE_CLASS
// clang-format off
KSpaceStyle(ewald/electrode, EwaldElectrode);
// clang-format on
#else

#ifndef LMP_EWALD_ELECTRODE_H
#define LMP_EWALD_ELECTRODE_H

#include "electrode_kspace.h"
#include "ewald.h"

namespace LAMMPS_NS {

class EwaldElectrode : public Ewald, public ElectrodeKSpace {
 public:
  EwaldElectrode(class LAMMPS *);
  ~EwaldElectrode() override;
  void init() override;
  void setup() override;
  void compute(int, int) override;
  void compute_group_group(int, int, int) override;

  // k-space part of coulomb matrix computation
  void compute_vector(double *, int, int, bool) override;
  void compute_vector_corr(double *, int, int, bool) override;
  void compute_matrix(bigint *, double **, bool) override;
  void compute_matrix_corr(bigint *, double **) override;

 protected:
  class BoundaryCorrection *boundcorr;
  double area;
  void coeffs() override;
  void eik_dot_r() override;

 private:
  int eikr_step;
  void update_eikr(bool);
};

}    // namespace LAMMPS_NS

#endif
#endif

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
   Contributing authors: Ludwig Ahrens-Iwers (TUHH), Shern Tee (UQ), Robert Meissner (TUHH)
------------------------------------------------------------------------- */

#ifdef KSPACE_CLASS
// clang-format off
KSpaceStyle(pppm/electrode, PPPMElectrode);
// clang-format on
#else

#ifndef LMP_PPPM_ELECTRODE_H
#define LMP_PPPM_ELECTRODE_H

#include "electrode_kspace.h"
#include "pppm.h"

namespace LAMMPS_NS {

class PPPMElectrode : public PPPM, public ElectrodeKSpace {
 public:
  PPPMElectrode(class LAMMPS *);
  ~PPPMElectrode() override;
  void init() override;
  void setup() override;
  void reset_grid() override;
  void compute(int, int) override;

  void compute_vector(double *, int, int, bool) override;
  void compute_vector_corr(double *, int, int, bool) override;
  void compute_matrix(bigint *, double **, bool) override;
  void compute_matrix_corr(bigint *, double **) override;

  void compute_group_group(int, int, int) override;

 protected:
  FFT_SCALAR ***electrolyte_density_brick;
  FFT_SCALAR *electrolyte_density_fft;
  class BoundaryCorrection *boundcorr;

  void set_grid_global() override;
  void set_grid_local() override;

  void allocate() override;
  void deallocate() override;
  void allocate_peratom() override;
  double compute_df_kspace() override;
  double compute_qopt() override;
  void compute_gf_ik() override;
  void compute_gf_ad() override;

 private:
  int compute_step;
  int last_source_grpbit;
  bool last_invert_source;
  void start_compute();
  void make_rho_in_brick(int, FFT_SCALAR ***, bool);
  void project_psi(double *, int);
  void one_step_multiplication(bigint *, double *, double **, double **, const int, bool);
  void two_step_multiplication(bigint *, double *, double **, double **, const int, bool);
  void build_amesh(int, int, int, double *, double *);
  bool compute_vector_called;
};

}    // namespace LAMMPS_NS

#endif
#endif

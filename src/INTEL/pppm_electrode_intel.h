/* *- c++ -*- -----------------------------------------------------------
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
   Contributing authors: William McDoniel (RWTH Aachen University)
                         Rodrigo Canales (RWTH Aachen University)
                         Markus Hoehnerbach (RWTH Aachen University)
                         W. Michael Brown (Intel)
------------------------------------------------------------------------- */

#ifdef KSPACE_CLASS

// clang-format off
KSpaceStyle(pppm/electrode/intel,PPPMElectrodeIntel)
// clang-format on

#else

#ifndef LMP_PPPM_ELECTRODE_INTEL_H
#define LMP_PPPM_ELECTRODE_INTEL_H

#include "pppm.h"
#include "fix_intel.h"
#include <algorithm>
#include "electrode_kspace.h"
#include "boundary_correction.h"
#include "pppm_intel.h"

namespace LAMMPS_NS {

class PPPMElectrodeIntel : public PPPMIntel, public ElectrodeKSpace {
 public:
  PPPMElectrodeIntel(class LAMMPS *);
  virtual ~PPPMElectrodeIntel();
  virtual void init();
  virtual void setup();
  virtual void compute(int, int);

  void compute_vector(bigint *, double *);
  void compute_vector_corr(bigint *, double *);
  void compute_matrix(bigint *, double **);
  void compute_matrix_corr(bigint *, double **);

  virtual void compute_group_group(int, int, int);

 protected:
  FFT_SCALAR ***electrolyte_density_brick;
  FFT_SCALAR *electrolyte_density_fft;
  class BoundaryCorrection *boundcorr;

  // virtual void set_grid_global();
  // void set_grid_local();

  virtual void allocate();
  // virtual void deallocate();
  // double compute_df_kspace();

 private:
  int compute_step;
  void start_compute();
  template <class flt_t, class acc_t, int use_table>
  void make_rho_in_brick(IntelBuffers<flt_t, acc_t> *buffers, bigint *,
                         FFT_SCALAR ***, bool);
  template <class flt_t, class acc_t>
  void make_rho_in_brick(IntelBuffers<flt_t, acc_t> *buffers, bigint *imat,
                         FFT_SCALAR ***scratch_brick, bool which_particles) {
    if (_use_table == 1)
      make_rho_in_brick<flt_t, acc_t, 1>(buffers, imat, scratch_brick,
                                         which_particles);
    else
      make_rho_in_brick<flt_t, acc_t, 0>(buffers, imat, scratch_brick,
                                         which_particles);
  }
  template <class flt_t, class acc_t, int use_table>
  void project_psi(IntelBuffers<flt_t, acc_t> *buffers, bigint *,
                         double *, double *);
  template <class flt_t, class acc_t>
  void project_psi(IntelBuffers<flt_t, acc_t> *buffers, bigint *imat,
                         double *vec, double *brick_psi) {
    if (_use_table == 1)
      project_psi<flt_t, acc_t, 1>(buffers, imat, vec,
                                         brick_psi);
    else
      project_psi<flt_t, acc_t, 0>(buffers, imat, vec,
                                         brick_psi);
  }


  void one_step_multiplication(bigint *, std::vector<double>,
                               double **, double **, int const);
  void two_step_multiplication(bigint *, std::vector<double>,
                               double **, double **, int const);
  bool compute_vector_called;
  bigint *imat_cached;
};

}  // namespace LAMMPS_NS

#endif
#endif


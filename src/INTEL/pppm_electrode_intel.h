/* *- c++ -*- -----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
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

#include "boundary_correction.h"
#include "electrode_kspace.h"
#include "fix_intel.h"
#include "pppm.h"
#include "pppm_intel.h"
#include <algorithm>

namespace LAMMPS_NS {

class PPPMElectrodeIntel : public PPPMIntel, public ElectrodeKSpace {
 public:
  PPPMElectrodeIntel(class LAMMPS *);
  ~PPPMElectrodeIntel();
  void init() override;
  void setup() override;
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

  void allocate() override;
  void deallocate() override;
  void allocate_peratom() override;

 private:
  int compute_step;
  int last_source_grpbit;
  bool last_invert_source;
  void start_compute();
  template <class flt_t, class acc_t, int use_table>
  void make_rho_in_brick(IntelBuffers<flt_t, acc_t> *buffers, int, FFT_SCALAR ***, bool);
  template <class flt_t, class acc_t>
  void make_rho_in_brick(IntelBuffers<flt_t, acc_t> *buffers, int source_grpbit,
                         FFT_SCALAR ***scratch_brick, bool invert_source)
  {
    if (_use_table == 1)
      make_rho_in_brick<flt_t, acc_t, 1>(buffers, source_grpbit, scratch_brick, invert_source);
    else
      make_rho_in_brick<flt_t, acc_t, 0>(buffers, source_grpbit, scratch_brick, invert_source);
  }
  template <class flt_t, class acc_t, int use_table>
  void project_psi(IntelBuffers<flt_t, acc_t> *buffers, double *, int);
  template <class flt_t, class acc_t>
  void project_psi(IntelBuffers<flt_t, acc_t> *buffers, double *vec, int sensor_grpbit)
  {
    if (_use_table == 1)
      project_psi<flt_t, acc_t, 1>(buffers, vec, sensor_grpbit);
    else
      project_psi<flt_t, acc_t, 0>(buffers, vec, sensor_grpbit);
  }

  void one_step_multiplication(bigint *, std::vector<double>, double **, double **, int const,
                               bool);
  void two_step_multiplication(bigint *, std::vector<double>, double **, double **, int const,
                               bool);
  bool compute_vector_called;
};

}    // namespace LAMMPS_NS

#endif
#endif

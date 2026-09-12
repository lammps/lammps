
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
   Contributing authors: Ludwig Ahrens-Iwers (TUHH), Robert Meissner (Hereon, TUHH),
   Shern Tee (GU) with LLM (GLM-5.1)
------------------------------------------------------------------------- */

#ifndef LMP_SLAB_INTEL_DIPOLE_H
#define LMP_SLAB_INTEL_DIPOLE_H

#include "boundary_correction.h"
#include "intel_buffers.h"
#include "fix_intel.h"

namespace LAMMPS_NS {

class SlabDipoleIntel : public BoundaryCorrection {
 public:
  SlabDipoleIntel(LAMMPS *, FixIntel *);
  void vector_corr(double *, int, int, bool) override;
  void matrix_corr(bigint *, double **) override;
  void compute_corr(double, int, int, double &, double *) override;

 private:
  template <class flt_t, class acc_t>
  void vector_corr(IntelBuffers<flt_t, acc_t> *buffers,
                  double *, int, int, bool);
  template <class flt_t, class acc_t>
  void matrix_corr(IntelBuffers<flt_t, acc_t> *buffers,
                  bigint *, double **);
  template <class flt_t, class acc_t>
  void compute_corr(IntelBuffers<flt_t, acc_t> *buffers,
                  double, int, int, double &, double *);

  FixIntel* fix;
  int _use_lrt;
};

}    // namespace LAMMPS_NS
#endif

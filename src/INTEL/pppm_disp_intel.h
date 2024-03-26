// clang-format off
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
   Contributing authors: William McDoniel (RWTH Aachen University)
------------------------------------------------------------------------- */

#ifdef KSPACE_CLASS
// clang-format off
KSpaceStyle(pppm/disp/intel,PPPMDispIntel);
// clang-format on
#else

#ifndef LMP_PPPMINTEL_DISP_H
#define LMP_PPPMINTEL_DISP_H

#include "fix_intel.h"
#include "pppm_disp.h"

namespace LAMMPS_NS {

class PPPMDispIntel : public PPPMDisp {
 public:
  PPPMDispIntel(class LAMMPS *);
  ~PPPMDispIntel() override;
  void init() override;
  void compute(int, int) override;

#ifdef _LMP_INTEL_OFFLOAD
  int use_base();
#endif

 protected:
  FixIntel *fix;

  int _use_lrt;
  FFT_SCALAR **perthread_density;
  FFT_SCALAR *particle_ekx;
  FFT_SCALAR *particle_eky;
  FFT_SCALAR *particle_ekz;
  FFT_SCALAR *particle_ekx0;
  FFT_SCALAR *particle_eky0;
  FFT_SCALAR *particle_ekz0;
  FFT_SCALAR *particle_ekx1;
  FFT_SCALAR *particle_eky1;
  FFT_SCALAR *particle_ekz1;
  FFT_SCALAR *particle_ekx2;
  FFT_SCALAR *particle_eky2;
  FFT_SCALAR *particle_ekz2;
  FFT_SCALAR *particle_ekx3;
  FFT_SCALAR *particle_eky3;
  FFT_SCALAR *particle_ekz3;
  FFT_SCALAR *particle_ekx4;
  FFT_SCALAR *particle_eky4;
  FFT_SCALAR *particle_ekz4;
  FFT_SCALAR *particle_ekx5;
  FFT_SCALAR *particle_eky5;
  FFT_SCALAR *particle_ekz5;
  FFT_SCALAR *particle_ekx6;
  FFT_SCALAR *particle_eky6;
  FFT_SCALAR *particle_ekz6;

  int _use_table;
  int rho_points;
  FFT_SCALAR **rho_lookup;
  FFT_SCALAR **rho6_lookup;
  FFT_SCALAR **drho_lookup;
  FFT_SCALAR **drho6_lookup;
  FFT_SCALAR half_rho_scale, half_rho_scale_plus;

  int _use_packing;

#ifdef _LMP_INTEL_OFFLOAD
  int _use_base;
#endif

  template <class flt_t, class acc_t>
  void particle_map_intel(double, double, double, double, int **, int, int, int, int, int, int, int, int,
                          IntelBuffers<flt_t, acc_t> *buffers);

  template <class flt_t, class acc_t, int use_table>
  void make_rho_c_intel(IntelBuffers<flt_t, acc_t> *buffers);
  template <class flt_t, class acc_t> void make_rho_c_intel(IntelBuffers<flt_t, acc_t> *buffers)
  {
    if (_use_table == 1) {
      make_rho_c_intel<flt_t, acc_t, 1>(buffers);
    } else {
      make_rho_c_intel<flt_t, acc_t, 0>(buffers);
    }
  }

  template <class flt_t, class acc_t, int use_table>
  void make_rho_g_intel(IntelBuffers<flt_t, acc_t> *buffers);
  template <class flt_t, class acc_t> void make_rho_g_intel(IntelBuffers<flt_t, acc_t> *buffers)
  {
    if (_use_table == 1) {
      make_rho_g_intel<flt_t, acc_t, 1>(buffers);
    } else {
      make_rho_g_intel<flt_t, acc_t, 0>(buffers);
    }
  }

  template <class flt_t, class acc_t, int use_table>
  void make_rho_a_intel(IntelBuffers<flt_t, acc_t> *buffers);
  template <class flt_t, class acc_t> void make_rho_a_intel(IntelBuffers<flt_t, acc_t> *buffers)
  {
    if (_use_table == 1) {
      make_rho_a_intel<flt_t, acc_t, 1>(buffers);
    } else {
      make_rho_a_intel<flt_t, acc_t, 0>(buffers);
    }
  }

  template <class flt_t, class acc_t, int use_table>
  void make_rho_none_intel(IntelBuffers<flt_t, acc_t> *buffers);
  template <class flt_t, class acc_t> void make_rho_none_intel(IntelBuffers<flt_t, acc_t> *buffers)
  {
    if (_use_table == 1) {
      make_rho_none_intel<flt_t, acc_t, 1>(buffers);
    } else {
      make_rho_none_intel<flt_t, acc_t, 0>(buffers);
    }
  }

  template <class flt_t, class acc_t, int use_table>
  void fieldforce_c_ik_intel(IntelBuffers<flt_t, acc_t> *buffers);
  template <class flt_t, class acc_t> void fieldforce_c_ik_intel(IntelBuffers<flt_t, acc_t> *buffers)
  {
    if (_use_table == 1) {
      fieldforce_c_ik_intel<flt_t, acc_t, 1>(buffers);
    } else {
      fieldforce_c_ik_intel<flt_t, acc_t, 0>(buffers);
    }
  }

  template <class flt_t, class acc_t, int use_table>
  void fieldforce_c_ad_intel(IntelBuffers<flt_t, acc_t> *buffers);
  template <class flt_t, class acc_t> void fieldforce_c_ad_intel(IntelBuffers<flt_t, acc_t> *buffers)
  {
    if (_use_table == 1) {
      fieldforce_c_ad_intel<flt_t, acc_t, 1>(buffers);
    } else {
      fieldforce_c_ad_intel<flt_t, acc_t, 0>(buffers);
    }
  }

  template <class flt_t, class acc_t, int use_table>
  void fieldforce_g_ik_intel(IntelBuffers<flt_t, acc_t> *buffers);
  template <class flt_t, class acc_t> void fieldforce_g_ik_intel(IntelBuffers<flt_t, acc_t> *buffers)
  {
    if (_use_table == 1) {
      fieldforce_g_ik_intel<flt_t, acc_t, 1>(buffers);
    } else {
      fieldforce_g_ik_intel<flt_t, acc_t, 0>(buffers);
    }
  }

  template <class flt_t, class acc_t, int use_table>
  void fieldforce_g_ad_intel(IntelBuffers<flt_t, acc_t> *buffers);
  template <class flt_t, class acc_t> void fieldforce_g_ad_intel(IntelBuffers<flt_t, acc_t> *buffers)
  {
    if (_use_table == 1) {
      fieldforce_g_ad_intel<flt_t, acc_t, 1>(buffers);
    } else {
      fieldforce_g_ad_intel<flt_t, acc_t, 0>(buffers);
    }
  }

  template <class flt_t, class acc_t, int use_table>
  void fieldforce_a_ik_intel(IntelBuffers<flt_t, acc_t> *buffers);
  template <class flt_t, class acc_t> void fieldforce_a_ik_intel(IntelBuffers<flt_t, acc_t> *buffers)
  {
    if (_use_table == 1) {
      fieldforce_a_ik_intel<flt_t, acc_t, 1>(buffers);
    } else {
      fieldforce_a_ik_intel<flt_t, acc_t, 0>(buffers);
    }
  }

  template <class flt_t, class acc_t, int use_table>
  void fieldforce_a_ad_intel(IntelBuffers<flt_t, acc_t> *buffers);
  template <class flt_t, class acc_t> void fieldforce_a_ad_intel(IntelBuffers<flt_t, acc_t> *buffers)
  {
    if (_use_table == 1) {
      fieldforce_a_ad_intel<flt_t, acc_t, 1>(buffers);
    } else {
      fieldforce_a_ad_intel<flt_t, acc_t, 0>(buffers);
    }
  }
  template <class flt_t, class acc_t, int use_table>
  void fieldforce_none_ik_intel(IntelBuffers<flt_t, acc_t> *buffers);
  template <class flt_t, class acc_t> void fieldforce_none_ik_intel(IntelBuffers<flt_t, acc_t> *buffers)
  {
    if (_use_table == 1) {
      fieldforce_none_ik_intel<flt_t, acc_t, 1>(buffers);
    } else {
      fieldforce_none_ik_intel<flt_t, acc_t, 0>(buffers);
    }
  }

  template <class flt_t, class acc_t, int use_table>
  void fieldforce_none_ad_intel(IntelBuffers<flt_t, acc_t> *buffers);
  template <class flt_t, class acc_t> void fieldforce_none_ad_intel(IntelBuffers<flt_t, acc_t> *buffers)
  {
    if (_use_table == 1) {
      fieldforce_none_ad_intel<flt_t, acc_t, 1>(buffers);
    } else {
      fieldforce_none_ad_intel<flt_t, acc_t, 0>(buffers);
    }
  }

  void precompute_rho();
};

}    // namespace LAMMPS_NS
#endif
#endif

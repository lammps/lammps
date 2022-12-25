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
   Contributing author: Markus Hohnerbach (RWTH)
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(airebo/intel,PairAIREBOIntel);
// clang-format on
#else

#ifndef LMP_PAIR_AIREBO_INTEL_H
#define LMP_PAIR_AIREBO_INTEL_H

#include "fix_intel.h"
#include "pair.h"
#include "pair_airebo.h"

namespace LAMMPS_NS {

template <class flt_t, class acc_t> struct PairAIREBOIntelParam;

class PairAIREBOIntel : public PairAIREBO {
 public:
  PairAIREBOIntel(class LAMMPS *);
  ~PairAIREBOIntel() override;
  void compute(int, int) override;
  void init_style() override;

 protected:
  template <class flt_t, class acc_t>
  void compute(int eflag, int vflag, IntelBuffers<flt_t, acc_t> *buffers);

  template <int EVFLAG, int EFLAG, class flt_t, class acc_t>
  void eval(const int offload, const int vflag, IntelBuffers<flt_t, acc_t> *buffers,
            const int astart, const int aend);

  template <class flt_t, class acc_t> void pack_force_const(IntelBuffers<flt_t, acc_t> *buffers);

  template <class flt_t, class acc_t> PairAIREBOIntelParam<flt_t, acc_t> get_param();

  FixIntel *fix;
  int _cop;

  int *REBO_cnumneigh;
  int *REBO_num_skin;
  int *REBO_list_data;
};

}    // namespace LAMMPS_NS

#endif
#endif

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
   Contributing author: W. Michael Brown (Intel)
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(lj/cut/intel,PairLJCutIntel);
// clang-format on
#else

#ifndef LMP_PAIR_LJ_CUT_INTEL_H
#define LMP_PAIR_LJ_CUT_INTEL_H

#include "fix_intel.h"
#include "pair_lj_cut.h"

namespace LAMMPS_NS {

class PairLJCutIntel : public PairLJCut {

 public:
  PairLJCutIntel(class LAMMPS *);

  void compute(int, int) override;
  void init_style() override;

 private:
  FixIntel *fix;
  int _cop, _onetype;

  template <class flt_t> class ForceConst;
  template <class flt_t, class acc_t>
  void compute(int eflag, int vflag, IntelBuffers<flt_t, acc_t> *buffers,
               const ForceConst<flt_t> &fc);
  template <int ONETYPE, int EFLAG, int NEWTON_PAIR, class flt_t, class acc_t>
  void eval(const int offload, const int vflag, IntelBuffers<flt_t, acc_t> *buffers,
            const ForceConst<flt_t> &fc, const int astart, const int aend);

  template <class flt_t, class acc_t>
  void pack_force_const(ForceConst<flt_t> &fc, IntelBuffers<flt_t, acc_t> *buffers);

  // ----------------------------------------------------------------------

  template <class flt_t> class ForceConst {
   public:
    typedef struct {
      flt_t cutsq, lj1, lj2, offset;
    } fc_packed1;
    typedef struct {
      flt_t lj3, lj4;
    } fc_packed2;

    _alignvar(flt_t special_lj[4], 64);
    fc_packed1 **ljc12o;
    fc_packed2 **lj34;

    ForceConst() : _ntypes(0) {}
    ~ForceConst() noexcept(false) { set_ntypes(0, nullptr, _cop); }

    void set_ntypes(const int ntypes, Memory *memory, const int cop);

   private:
    int _ntypes, _cop;
    Memory *_memory;
  };
  ForceConst<float> force_const_single;
  ForceConst<double> force_const_double;
};

}    // namespace LAMMPS_NS

#endif
#endif

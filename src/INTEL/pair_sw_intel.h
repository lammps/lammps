/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

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
PairStyle(sw/intel,PairSWIntel);
// clang-format on
#else

#ifndef LMP_PAIR_SW_INTEL_H
#define LMP_PAIR_SW_INTEL_H

#include "fix_intel.h"
#include "pair_sw.h"

namespace LAMMPS_NS {

class PairSWIntel : public PairSW {
 public:
  PairSWIntel(class LAMMPS *);

  void compute(int, int) override;
  void init_style() override;

 protected:
  FixIntel *fix;
  int _cop;
  template <class flt_t> class ForceConst;

  void allocate() override;

  template <class flt_t, class acc_t>
  void compute(int eflag, int vflag, IntelBuffers<flt_t, acc_t> *buffers,
               const ForceConst<flt_t> &fc);
  template <int SPQ, int ONETYPE, int EFLAG, class flt_t, class acc_t>
  void eval(const int offload, const int vflag, IntelBuffers<flt_t, acc_t> *buffers,
            const ForceConst<flt_t> &fc, const int astart, const int aend);

  template <class flt_t, class acc_t>
  void pack_force_const(ForceConst<flt_t> &fc, IntelBuffers<flt_t, acc_t> *buffers);

  int _ccache_stride, _spq, _onetype, _onetype3;
#ifdef LMP_USE_AVXCD
  int _ccache_stride3;
#endif

  // ----------------------------------------------------------------------

  template <class flt_t> class ForceConst {
   public:
    typedef struct {
      flt_t cutsq, cut, sigma_gamma, pad;
    } fc_packed0;
    typedef struct {
      flt_t powerp, powerq, cut, sigma;
    } fc_packed1;
    typedef struct {
      flt_t c1, c2, c3, c4;
    } fc_packed1p2;
    typedef struct {
      flt_t c5, c6, d1, d2;
    } fc_packed2;
    typedef struct {
      flt_t costheta, lambda_epsilon, lambda_epsilon2, pad;
    } fc_packed3;

    fc_packed0 **p2;
    fc_packed1 **p2f;
    fc_packed1p2 **p2f2;
    fc_packed2 **p2e;
    fc_packed3 ***p3;

    ForceConst() : p2(0), p2f(0), p2f2(0), p2e(0), p3(0), _ntypes(0) {}
    ~ForceConst() { set_ntypes(0, nullptr, _cop); }

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

// clang-format off
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
                        Shun Xu (Computer Network Information Center, CAS)
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(dpd/intel,PairDPDIntel);
// clang-format on
#else

#ifndef LMP_PAIR_DPD_INTEL_H
#define LMP_PAIR_DPD_INTEL_H

#include "fix_intel.h"
#include "pair_dpd.h"

#ifdef LMP_USE_MKL_RNG
#include "mkl_vsl.h"
#else
#include "random_mars.h"
#endif

namespace LAMMPS_NS {

class PairDPDIntel : public PairDPD {

 public:
  PairDPDIntel(class LAMMPS *);
  ~PairDPDIntel() override;

  void compute(int, int) override;
  void settings(int, char **) override;
  void init_style() override;
  void read_restart_settings(FILE *) override;

 private:
  FixIntel *fix;
  int _cop, _onetype, _nrandom_thread;

#ifdef LMP_USE_MKL_RNG
  VSLStreamStatePtr *random_thread;
#else
  RanMars **random_thread;
#endif

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
      flt_t icut, a0, gamma, sigma;
    } fc_packed1;

    _alignvar(flt_t special_lj[4], 64);
    fc_packed1 **param;
    flt_t **rand_buffer_thread;
    int *rngi;

    ForceConst() : _ntypes(0) {}
    ~ForceConst() { set_ntypes(0, 0, 0, nullptr, _cop); }

    void set_ntypes(const int ntypes, const int nthreads, const int max_nbors, Memory *memory,
                    const int cop);

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

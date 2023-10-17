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

#ifdef PAIR_CLASS
// clang-format off
// Currently the Intel compilers are required for this pair style.
#ifdef __INTEL_COMPILER
PairStyle(tersoff/intel,PairTersoffIntel);
#endif
// clang-format on
#else

#ifndef LMP_PAIR_TERSOFF_INTEL_H
#define LMP_PAIR_TERSOFF_INTEL_H

#include "fix_intel.h"
#include "pair.h"
#include "pair_tersoff.h"

#ifdef __INTEL_COMPILER

namespace LAMMPS_NS {

class PairTersoffIntel : public PairTersoff {
 public:
  PairTersoffIntel(class LAMMPS *);
  virtual void compute(int, int);
  void init_style();

 protected:
  typedef struct {
    float x, y, z;
    int w;
  } sng4_t;

 private:
  FixIntel *fix;
  int _cop;

 public:    // wo needs secrets?
  // ----------------------------------------------------------------------
  //
  template <class flt_t> class ForceConst {
   public:
    typedef struct {
      flt_t cutsq;
    } c_cutoff_t;
    typedef struct {
      flt_t bigr, bigd, lam1, biga;
    } c_first_loop_t;
    typedef struct {
      flt_t lam2, beta, bigb, powern, c1, c2, c3, c4;
    } c_second_loop_t;
    typedef struct {
      flt_t lam3, bigr, bigd, c2, d2, h, gamma, powermint;
    } c_inner_loop_t;
    typedef struct {
      flt_t cutsq, pad[3];
      flt_t bigr, bigd, lam1, biga;
      flt_t lam2, beta, bigb, powern;
      flt_t c1, c2, c3, c4;
    } c_outer_t;
    typedef struct {
      flt_t cutsq, pad[7];
      flt_t lam3, powermint, bigr, bigd;
      flt_t c2, d2, h, gamma;
    } c_inner_t;
    c_cutoff_t **c_cutoff_outer;
    c_cutoff_t ***c_cutoff_inner;
    c_first_loop_t **c_first_loop;
    c_second_loop_t **c_second_loop;
    c_inner_loop_t ***c_inner_loop;
    c_outer_t **c_outer;
    c_inner_t ***c_inner;
    ForceConst() : _ntypes(0) {}
    ~ForceConst() noexcept(false) { set_ntypes(0, nullptr, _cop); }

    void set_ntypes(const int ntypes, Memory *memory, const int cop);

   private:
    int _ntypes, _cop;
    Memory *_memory;
  };
  ForceConst<float> force_const_single;
  ForceConst<double> force_const_double;

  template <class flt_t, class acc_t>
  void compute(int eflag, int vflag, IntelBuffers<flt_t, acc_t> *buffers,
               const ForceConst<flt_t> &fc);
  template <int EFLAG, class flt_t, class acc_t>
  void eval(const int offload, const int vflag, IntelBuffers<flt_t, acc_t> *buffers,
            const ForceConst<flt_t> &fc, const int astart, const int aend);

  template <class flt_t, class acc_t>
  void pack_force_const(ForceConst<flt_t> &fc, IntelBuffers<flt_t, acc_t> *buffers);
};

}    // namespace LAMMPS_NS
#endif    // __INTEL_COMPILER

#endif
#endif

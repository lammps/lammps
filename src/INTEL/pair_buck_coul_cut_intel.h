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
   Contributing author: Rodrigo Canales (RWTH Aachen University)
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(buck/coul/cut/intel,PairBuckCoulCutIntel);
// clang-format on
#else

#ifndef LMP_PAIR_BUCK_COUL_CUT_INTEL_H
#define LMP_PAIR_BUCK_COUL_CUT_INTEL_H

#include "fix_intel.h"
#include "pair_buck_coul_cut.h"

namespace LAMMPS_NS {

class PairBuckCoulCutIntel : public PairBuckCoulCut {

 public:
  PairBuckCoulCutIntel(class LAMMPS *);

  void compute(int, int) override;
  void init_style() override;
  typedef struct {
    float x, y, z;
    int w;
  } sng4_t;

 private:
  FixIntel *fix;
  int _cop;

  template <class flt_t> class ForceConst;

  template <class flt_t, class acc_t>
  void compute(int eflag, int vflag, IntelBuffers<flt_t, acc_t> *buffers,
               const ForceConst<flt_t> &fc);

  template <int EFLAG, int NEWTON_PAIR, class flt_t, class acc_t>
  void eval(const int offload, const int vflag, IntelBuffers<flt_t, acc_t> *buffers,
            const ForceConst<flt_t> &fc, const int astart, const int aend);

  template <class flt_t, class acc_t>
  void pack_force_const(ForceConst<flt_t> &fc, IntelBuffers<flt_t, acc_t> *buffers);

  template <class flt_t> class ForceConst {
   public:
    typedef struct {
      flt_t buck1, buck2, rhoinv, pad;
    } c_force_t;
    typedef struct {
      flt_t cutsq, cut_ljsq, cut_coulsq, pad;
    } c_cut_t;
    typedef struct {
      flt_t a, c, offset, pad;
    } c_energy_t;
    _alignvar(flt_t special_coul[4], 64);
    _alignvar(flt_t special_lj[4], 64);

    c_force_t **c_force;
    c_energy_t **c_energy;
    c_cut_t **c_cut;

    ForceConst() : _ntypes(0), _ntable(0) {}
    ~ForceConst() noexcept(false) { set_ntypes(0, 0, nullptr, _cop); }

    void set_ntypes(const int ntypes, const int ntable, Memory *memory, const int cop);

   private:
    int _ntypes, _ntable, _cop;
    Memory *_memory;
  };
  ForceConst<float> force_const_single;
  ForceConst<double> force_const_double;
};

}    // namespace LAMMPS_NS

#endif    // LMP_PAIR_BUCK_COUL_CUT_INTEL_H
#endif

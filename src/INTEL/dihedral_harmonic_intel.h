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

#ifdef DIHEDRAL_CLASS
// clang-format off
DihedralStyle(harmonic/intel,DihedralHarmonicIntel);
// clang-format on
#else

#ifndef LMP_DIHEDRAL_HARMONIC_INTEL_H
#define LMP_DIHEDRAL_HARMONIC_INTEL_H

#include "dihedral_harmonic.h"
#include "fix_intel.h"

namespace LAMMPS_NS {

class DihedralHarmonicIntel : public DihedralHarmonic {

 public:
  DihedralHarmonicIntel(class LAMMPS *lmp);
  void compute(int, int) override;
  void init_style() override;

 private:
  FixIntel *fix;

  template <class flt_t> class ForceConst;
  template <class flt_t, class acc_t>
  void compute(int eflag, int vflag, IntelBuffers<flt_t, acc_t> *buffers,
               const ForceConst<flt_t> &fc);
  template <int EVFLAG, int EFLAG, int NEWTON_BOND, class flt_t, class acc_t>
  void eval(const int vflag, IntelBuffers<flt_t, acc_t> *buffers, const ForceConst<flt_t> &fc);
  template <class flt_t, class acc_t>
  void pack_force_const(ForceConst<flt_t> &fc, IntelBuffers<flt_t, acc_t> *buffers);

#ifdef _LMP_INTEL_OFFLOAD
  int _use_base;
#endif

  template <class flt_t> class ForceConst {
   public:
    typedef struct {
      flt_t cos_shift, sin_shift, k;
      int multiplicity;
    } fc_packed1;

    fc_packed1 *fc;

    ForceConst() : fc(nullptr), _ndihderaltypes(0) {}
    ~ForceConst() noexcept(false) { set_ntypes(0, nullptr); }

    void set_ntypes(const int ndihderaltypes, Memory *memory);

   private:
    int _ndihderaltypes;
    Memory *_memory;
  };
  ForceConst<float> force_const_single;
  ForceConst<double> force_const_double;
};

}    // namespace LAMMPS_NS

#endif
#endif

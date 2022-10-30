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

#ifdef IMPROPER_CLASS
// clang-format off
ImproperStyle(harmonic/intel,ImproperHarmonicIntel);
// clang-format on
#else

#ifndef LMP_IMPROPER_HARMONIC_INTEL_H
#define LMP_IMPROPER_HARMONIC_INTEL_H

#include "fix_intel.h"
#include "improper_harmonic.h"

namespace LAMMPS_NS {

class ImproperHarmonicIntel : public ImproperHarmonic {
 public:
  ImproperHarmonicIntel(class LAMMPS *);
  void compute(int, int) override;
  void init_style() override;

 protected:
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
      flt_t k, chi;
    } fc_packed1;

    fc_packed1 *fc;

    ForceConst() : fc(nullptr), _nimpropertypes(0) {}
    ~ForceConst() { set_ntypes(0, nullptr); }

    void set_ntypes(const int nimpropertypes, Memory *memory);

   private:
    int _nimpropertypes;
    Memory *_memory;
  };
  ForceConst<float> force_const_single;
  ForceConst<double> force_const_double;
};

}    // namespace LAMMPS_NS

#endif
#endif

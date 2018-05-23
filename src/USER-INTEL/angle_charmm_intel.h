/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
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

#ifdef ANGLE_CLASS

AngleStyle(charmm/intel,AngleCharmmIntel)

#else

#ifndef LMP_ANGLE_CHARMM_INTEL_H
#define LMP_ANGLE_CHARMM_INTEL_H

#include <cstdio>
#include "angle_charmm.h"
#include "fix_intel.h"

namespace LAMMPS_NS {

class AngleCharmmIntel : public AngleCharmm {
 public:
  AngleCharmmIntel(class LAMMPS *);
  virtual ~AngleCharmmIntel();
  virtual void compute(int, int);
  virtual void init_style();

 protected:
  FixIntel *fix;

  template <class flt_t> class ForceConst;
  template <class flt_t, class acc_t>
  void compute(int eflag, int vflag, IntelBuffers<flt_t,acc_t> *buffers,
               const ForceConst<flt_t> &fc);
  template <int EVFLAG, int EFLAG, int NEWTON_BOND, class flt_t, class acc_t>
  void eval(const int vflag, IntelBuffers<flt_t,acc_t> * buffers,
            const ForceConst<flt_t> &fc);
  template <class flt_t, class acc_t>
  void pack_force_const(ForceConst<flt_t> &fc,
                        IntelBuffers<flt_t, acc_t> *buffers);

  #ifdef _LMP_INTEL_OFFLOAD
  int _use_base;
  #endif

  template <class flt_t>
  class ForceConst {
   public:
    typedef struct { flt_t k, theta0, k_ub, r_ub; } fc_packed1;

    fc_packed1 *fc;
    ForceConst() : _nangletypes(0)  {}
    ~ForceConst() { set_ntypes(0, NULL); }

    void set_ntypes(const int nangletypes, Memory *memory);

   private:
    int _nangletypes;
    Memory *_memory;
  };
  ForceConst<float> force_const_single;
  ForceConst<double> force_const_double;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Incorrect args for angle coefficients

Self-explanatory.  Check the input script or data file.

*/

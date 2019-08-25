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

#ifdef DIHEDRAL_CLASS

DihedralStyle(charmm/intel,DihedralCharmmIntel)

#else

#ifndef LMP_DIHEDRAL_CHARMM_INTEL_H
#define LMP_DIHEDRAL_CHARMM_INTEL_H

#include "dihedral_charmm.h"
#include "fix_intel.h"

namespace LAMMPS_NS {

class DihedralCharmmIntel : public DihedralCharmm {

 public:
  DihedralCharmmIntel(class LAMMPS *lmp);
  virtual void compute(int, int);
  void init_style();

 private:
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
    typedef struct { flt_t lj1, lj2, lj3, lj4; } fc_packed1;
    typedef struct { flt_t cos_shift, sin_shift, k;
                     int multiplicity; } fc_packed3;

    fc_packed1 **ljp;
    fc_packed3 *bp;
    flt_t *weight;

    ForceConst() : _npairtypes(0), _nbondtypes(0)  {}
    ~ForceConst() { set_ntypes(0, 0, NULL); }

    void set_ntypes(const int npairtypes, const int nbondtypes, Memory *memory);

   private:
    int _npairtypes, _nbondtypes;
    Memory *_memory;
  };
  ForceConst<float> force_const_single;
  ForceConst<double> force_const_double;
};

}

#endif
#endif

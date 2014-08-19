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

#ifdef PAIR_CLASS

PairStyle(gayberne/intel,PairGayBerneIntel)

#else

#ifndef LMP_PAIR_GAYBERNE_INTEL_H
#define LMP_PAIR_GAYBERNE_INTEL_H

#include "pair_gayberne.h"
#include "fix_intel.h"

namespace LAMMPS_NS {

class PairGayBerneIntel : public PairGayBerne {

 public:
  PairGayBerneIntel(class LAMMPS *);

  virtual void compute(int, int);
  void init_style();

 private:
  template <class flt_t> class ForceConst;

  template <class flt_t, class acc_t>
  void compute(int eflag, int vflag, IntelBuffers<flt_t,acc_t> *buffers,
               const ForceConst<flt_t> &fc);
  template <int EVFLAG, int EFLAG, int NEWTON_PAIR, class flt_t, class acc_t>
  void eval(const int offload, const int vflag,
            IntelBuffers<flt_t,acc_t> * buffers,
            const ForceConst<flt_t> &fc, const int astart, const int aend);

  template <class flt_t, class acc_t>
  void pack_force_const(ForceConst<flt_t> &fc,
                        IntelBuffers<flt_t, acc_t> *buffers);

  template <class flt_t>
  class ForceConst {
   public:
    typedef struct { 
      flt_t cutsq, lj1, lj2, offset, sigma, epsilon, lshape; 
      int form; 
    } fc_packed1;
    typedef struct { flt_t lj3, lj4; } fc_packed2;
    typedef struct { flt_t shape2[4], well[4]; } fc_packed3;

    __declspec(align(64)) flt_t special_lj[4], gamma, upsilon, mu;
    fc_packed1 **ijc;
    fc_packed2 **lj34;
    fc_packed3 *ic;

    flt_t **rsq_form, **delx_form, **dely_form, **delz_form;
    int **jtype_form, **jlist_form;

    ForceConst() : _ntypes(0)  {}
    ~ForceConst() { set_ntypes(0, 0, 0, NULL, _cop); }

    void set_ntypes(const int ntypes, const int one_length,
                    const int nthreads, Memory *memory, const int cop);

   private:
    int _ntypes, _cop;
    Memory *_memory;
  };

  ForceConst<float> force_const_single;
  ForceConst<double> force_const_double;
  int _max_nbors;

  double gayberne_lj(const int i, const int j, double a1[3][3],
                     double b1[3][3], double g1[3][3], double *r12,
                     const double rsq, double *fforce, double *ttor);

  FixIntel *fix;
  int _cop;
};

}

#endif
#endif

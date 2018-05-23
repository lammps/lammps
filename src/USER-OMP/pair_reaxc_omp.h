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

#ifdef PAIR_CLASS

PairStyle(reax/c/omp,PairReaxCOMP)

#else

#ifndef LMP_PAIR_REAXC_OMP_H
#define LMP_PAIR_REAXC_OMP_H

#include "pair_reaxc.h"
#include "thr_omp.h"
#include "suffix.h"

namespace LAMMPS_NS {

class PairReaxCOMP : public PairReaxC, public ThrOMP {
 public:
  PairReaxCOMP(class LAMMPS *);
  ~PairReaxCOMP();
  virtual void compute(int, int);
  virtual void init_style();

  inline FixOMP *getFixOMP() {
        return fix;
  };

  inline void ev_setup_thr_proxy(int eflagparm, int vflagparm, int nallparm,
                                 double *eatomparm, double **vatomparm, ThrData *thrparm) {
    ev_setup_thr(eflagparm, vflagparm, nallparm, eatomparm, vatomparm, thrparm);
  };

  // reduce per thread data as needed
  inline void reduce_thr_proxy(void * const styleparm, const int eflagparm,
                               const int vflagparm, ThrData * const thrparm) {
    reduce_thr(styleparm, eflagparm, vflagparm, thrparm);
  }

  inline void ev_tally_thr_proxy(Pair * const pairparm, const int iparm, const int jparm,
                                 const int nlocalparm, const int newton_pairparm,
                                 const double evdwlparm, const double ecoulparm,
                                 const double fpairparm, const double delxparm,
                                 const double delyparm, const double delzparm,
                                 ThrData * const thrparm) {
    ev_tally_thr(pairparm, iparm, jparm, nlocalparm, newton_pairparm,
                 evdwlparm, ecoulparm, fpairparm, delxparm, delyparm, delzparm, thrparm);
  }

  inline void ev_tally_xyz_thr_proxy(Pair * const pairparm, const int iparm, const int jparm,
                                     const int nlocalparm, const int newton_pairparm,
                                     const double evdwlparm, const double ecoulparm,
                                     const double fxparm, const double fyparm, const double fzparm,
                                     const double delxparm, const double delyparm,
                                     const double delzparm, ThrData * const thrparm) {
    ev_tally_xyz_thr(pairparm, iparm, jparm, nlocalparm, newton_pairparm,
                     evdwlparm, ecoulparm, fxparm, fyparm, fzparm,
                     delxparm, delyparm, delzparm, thrparm);
  }

  inline void ev_tally3_thr_proxy(Pair * const pairparm,int i, int j, int k,
                                  double evdwl, double ecoul, double *fj, double *fk,
                                  double *drji, double *drki, ThrData * const thrparm) {
    ev_tally3_thr(pairparm, i, j, k, evdwl, ecoul, fj, fk, drji, drki, thrparm);
  }

 protected:
  virtual void setup();
  virtual void write_reax_atoms();
  virtual int estimate_reax_lists();
  virtual int write_reax_lists();
  virtual void read_reax_forces(int);
  virtual void FindBond();

  // work array used in write_reax_lists()
  int * num_nbrs_offset;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Too many ghost atoms

Number of ghost atoms has increased too much during simulation and has exceeded
the size of reax/c arrays.  Increase safe_zone and min_cap in pair_style reax/c
command

*/

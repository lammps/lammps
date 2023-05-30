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
PairStyle(reaxff/omp,PairReaxFFOMP);
PairStyle(reax/c/omp,PairReaxFFOMP);
// clang-format on
#else

#ifndef LMP_PAIR_REAXFF_OMP_H
#define LMP_PAIR_REAXFF_OMP_H

#include "pair_reaxff.h"
#include "thr_omp.h"    // IWYU pragma: export

namespace LAMMPS_NS {

class PairReaxFFOMP : public PairReaxFF, public ThrOMP {
 public:
  PairReaxFFOMP(class LAMMPS *);
  ~PairReaxFFOMP() override;
  void compute(int, int) override;
  void init_style() override;

  inline FixOMP *getFixOMP() { return fix; };

  inline void ev_setup_thr_proxy(int eflagparm, int vflagparm, int nallparm, double *eatomparm,
                                 double **vatomparm, double **cvatomparm, ThrData *thrparm)
  {
    ev_setup_thr(eflagparm, vflagparm, nallparm, eatomparm, vatomparm, cvatomparm, thrparm);
  };

  // reduce per thread data as needed
  inline void reduce_thr_proxy(void *const styleparm, const int eflagparm, const int vflagparm,
                               ThrData *const thrparm)
  {
    reduce_thr(styleparm, eflagparm, vflagparm, thrparm);
  }

  inline void ev_tally_thr_proxy(const int iparm, const int jparm, const int nlocalparm,
                                 const int newton_pairparm, const double evdwlparm,
                                 const double ecoulparm, const double fpairparm,
                                 const double delxparm, const double delyparm,
                                 const double delzparm, ThrData *const thrparm)
  {
    ev_tally_thr(this, iparm, jparm, nlocalparm, newton_pairparm, evdwlparm, ecoulparm, fpairparm,
                 delxparm, delyparm, delzparm, thrparm);
  }

  inline void ev_tally_xyz_thr_proxy(const int iparm, const int jparm, const int nlocalparm,
                                     const int newton_pairparm, const double evdwlparm,
                                     const double ecoulparm, const double fxparm,
                                     const double fyparm, const double fzparm,
                                     const double delxparm, const double delyparm,
                                     const double delzparm, ThrData *const thrparm)
  {
    ev_tally_xyz_thr(this, iparm, jparm, nlocalparm, newton_pairparm, evdwlparm, ecoulparm, fxparm,
                     fyparm, fzparm, delxparm, delyparm, delzparm, thrparm);
  }

  inline void ev_tally3_thr_proxy(int i, int j, int k, double evdwl, double ecoul, double *fj,
                                  double *fk, double *drji, double *drki, ThrData *const thrparm)
  {
    ev_tally3_thr(this, i, j, k, evdwl, ecoul, fj, fk, drji, drki, thrparm);
  }

  inline void v_tally2_newton_thr_proxy(const int i, const double *const fi,
                                        const double *const deli, ThrData *const thrparm)
  {
    v_tally2_newton_thr(this, i, fi, deli, thrparm);
  }

  inline void v_tally3_thr_proxy(const int i, const int j, const int k, const double *const fi,
                                 const double *const fk, const double *const drij,
                                 const double *const drkj, ThrData *const thrparm)
  {
    v_tally3_thr(this, i, j, k, fi, fk, drij, drkj, thrparm);
  }

  inline void v_tally4_thr_proxy(const int i, const int j, const int k, const int l,
                                 const double *const fi, const double *const fj,
                                 const double *const fk, const double *const dril,
                                 const double *const drjl, const double *const drkl,
                                 ThrData *const thrparm)
  {
    v_tally4_thr(this, i, j, k, l, fi, fj, fk, dril, drjl, drkl, thrparm);
  }

 protected:
  void setup() override;
  virtual void write_reax_atoms();
  virtual int estimate_reax_lists();
  virtual int write_reax_lists();
  virtual void read_reax_forces(int);
  virtual void FindBond();

  // work array used in write_reax_lists()
  int *num_nbrs_offset;
};

}    // namespace LAMMPS_NS

#endif
#endif

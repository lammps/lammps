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
   Contributing author: David Immel (d.immel@fz-juelich.de, FZJ, Germany)
------------------------------------------------------------------------- */

#ifdef FIX_CLASS
// clang-format off
FixStyle(lambda/la/csp/apip,FixLambdaLACSPAPIP);
// clang-format on
#else

#ifndef LMP_FIX_LAMBDA_LA_CSP_H
#define LMP_FIX_LAMBDA_LA_CSP_H

#include "fix.h"

namespace LAMMPS_NS {

class CSPpairAPIP {

 public:
  CSPpairAPIP() { tag_smaller = tag_larger = 0; };
  CSPpairAPIP(tagint, tagint);

  tagint tag_smaller;
  tagint tag_larger;
};

inline bool operator==(const CSPpairAPIP &pl, const CSPpairAPIP &pr)
{
  return (pl.tag_smaller == pr.tag_smaller && pl.tag_larger == pr.tag_larger ? true : false);
}

inline bool operator!=(const CSPpairAPIP &pl, const CSPpairAPIP &pr)
{
  return !(pl == pr);
}

inline bool operator<(const CSPpairAPIP &pl, const CSPpairAPIP &pr)
{
  return ((pl.tag_smaller < pr.tag_smaller) ||
                  ((pl.tag_smaller == pr.tag_smaller) && (pl.tag_larger < pr.tag_larger))
              ? true
              : false);
}

inline bool operator>(const CSPpairAPIP &pl, const CSPpairAPIP &pr)
{
  return pr < pl;
}

class FixLambdaLACSPAPIP : public Fix {
 public:
  FixLambdaLACSPAPIP(class LAMMPS *, int, char **);
  ~FixLambdaLACSPAPIP() override;
  void init() override;
  void init_list(int, NeighList *) override;
  void post_constructor() override;
  void setup_post_neighbor() override;
  void post_neighbor() override;
  int setmask() override;
  void setup_pre_force(int) override;
  void pre_force(int) override;
  void *extract(const char *, int &) override;
  int pack_reverse_comm(int, int, double *) override;
  void unpack_reverse_comm(int, int *, double *) override;
  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;
  double compute_scalar() override;
  void write_restart(FILE *) override;
  void restart(char *) override;
  void setup_pre_reverse(int, int) override;
  void pre_reverse(int, int) override;

 private:
  double cut_lo;    ///< distance at which the cutoff function of the transition zone decays from 1
  double cut_hi;    ///< distance at which the cutoff function of the transition zone is 0
  double cut_width;           ///< cut_hi - cut_lo
  double cut_hi_sq;           ///< cut_hi_sq * cut_hi_sq
  double threshold_lo;        ///< threshold above which the fast potential starts to be turned off
  double threshold_hi;        ///< threshold above which the fast potential is turned off completely
  double threshold_width;     ///< threshold_hi - threshold_lo
  int nnn;                    ///< number of nearest neighbors used for the CSP calculation
  double csp_cutsq;           ///< cutof in which neighbors are considered for the CSP calculation
  double lambda_non_group;    ///< lambda for atoms that are not in the group of this fix
  double cutsq_combined;      ///< max of cut_hi_sq and csp_cutsq

  int **ngh_pairs;       ///< used neighbor pairs of local atoms
  int ngh_pairs_size;    ///< number of allocated neighbor pairs

  double **f_lambda;    ///< force due to nabla lambda on atoms
  int size_f_lambda;    ///< size of f_lambda

  int maxneigh;      ///< size of distsq, nearest
  double *distsq;    ///< squared distance to nearest neighbors
  int *nearest;      ///< index of nearest neighbors

  class NeighList *list;

  class FixStoreAtom *fixstore_pairs;      ///< ptr to stored neighbour pairs of last timestep
  class FixStoreAtom *fixstore_la_avg;     ///< ptr to stored local averaging data
  class FixStoreAtom *fixstore_la_inp;     ///< ptr to stored local averaging data
  class FixStoreAtom *fixstore_la_norm;    ///< ptr to stored local averaging data
  bool tags_stored;                        ///< CSP-pairs are stored/not stored for true/false
  int counter_changed_csp_nghs;    ///< scalar return value that indicates lost conservativity
  bool
      const_ngh_flag;    ///< use/do not use constant CSP-pairs to get a conservative APIP for true/false
  bool calculate_forces_flag;    ///< calculate/do not calculate forces for true/false
  bool store_stats;              ///< store per-atom stats/no stats for true/false

  enum { FORWARD_INP_LAMBDA, FORWARD_PREFACTOR };
  int comm_forward_flag;    // flag that determines which variables are communicated in comm forward

  double *prefactor1;     // per atom array for the force calculation
  int prefactor1_size;    // own + ghosts
  double *prefactor2;     // per atom array for the force calculation
  int prefactor2_size;    // own

  void calculate_forces(int);
  void store_f_lambda_before();
  void store_f_lambda_after();
  void store_la();
  void pre_force_const_pairs();
  void pre_force_dyn_pairs();
  double switching_function_poly(double);
  double der_switching_function_poly(double);
  double weighting_function_poly(double);
  double der_weighting_function_poly(double);
  void select2(int, int, double *, int *);
  void ev_tally2(int, int, double, double, double, double);
  void ev_tally3(int, int, int, double *, double, double, double, double, double, double);
};

}    // namespace LAMMPS_NS

#endif
#endif

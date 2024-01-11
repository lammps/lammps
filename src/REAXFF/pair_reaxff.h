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
   Contributing author: Hasan Metin Aktulga, Purdue University
   (now at Lawrence Berkeley National Laboratory, hmaktulga@lbl.gov)

   Please cite the related publication:
   H. M. Aktulga, J. C. Fogarty, S. A. Pandit, A. Y. Grama,
   "Parallel Reactive Molecular Dynamics: Numerical Methods and
   Algorithmic Techniques", Parallel Computing, in press.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(reaxff,PairReaxFF);
// clang-format on
#else

#ifndef LMP_PAIR_REAXFF_H
#define LMP_PAIR_REAXFF_H

#include "pair.h"

namespace ReaxFF {
struct API;
struct far_neighbor_data;
}    // namespace ReaxFF

namespace LAMMPS_NS {

class PairReaxFF : public Pair {
 public:
  PairReaxFF(class LAMMPS *);
  ~PairReaxFF() override;
  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  void init_style() override;
  double init_one(int, int) override;
  void *extract(const char *, int &) override;
  int fixbond_flag, fixspecies_flag;
  int **tmpid;
  double **tmpbo, **tmpr;
  std::vector<std::string> eletype;

  ReaxFF::API *api;
  typedef double rvec[3];

 protected:
  char *fix_id;
  double cutmax;
  class FixReaxFF *fix_reaxff;

  double *chi, *eta, *gamma, *bcut_acks2;
  int qeqflag;
  int setup_flag;
  int firstwarn;
  int list_blocking_flag;

  void allocate();
  void setup() override;
  void create_compute();
  void create_fix();
  void write_reax_atoms();
  void get_distance(rvec, rvec, double *, rvec *);
  void set_far_nbr(ReaxFF::far_neighbor_data *, int, double, rvec);
  int estimate_reax_lists();
  int write_reax_lists();
  void read_reax_forces(int);

  int nmax;
  void FindBond();
  double memory_usage() override;
};

}    // namespace LAMMPS_NS

#endif
#endif

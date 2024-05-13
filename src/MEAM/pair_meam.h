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
PairStyle(meam,PairMEAM);
// clang-format on
#else

#ifndef LMP_PAIR_MEAM_H
#define LMP_PAIR_MEAM_H

#include "pair.h"

namespace LAMMPS_NS {

class PairMEAM : public Pair {
 public:
  PairMEAM(class LAMMPS *);
  ~PairMEAM() override;
  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  void init_style() override;
  void init_list(int, class NeighList *) override;
  double init_one(int, int) override;
  void *extract(const char *, int &) override;

  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;
  int pack_reverse_comm(int, int, double *) override;
  void unpack_reverse_comm(int, int *, double *) override;
  double memory_usage() override;

 protected:
  class MEAM *meam_inst;
  double cutmax;                           // max cutoff for all elements
  int nlibelements;                        // # of library elements
  int msmeamflag;                          // 0 (default) for normal MEAM, 1 for MS-MEAM
  std::string myname;                      // name of the pair style
  std::vector<std::string> libelements;    // names of library elements
  std::vector<double> mass;                // mass of library element

  double **scale;    // scaling factor for adapt

  void allocate();
  void read_files(const std::string &, const std::string &);
  void read_global_meam_file(const std::string &);
  void read_user_meam_file(const std::string &);
  void neigh_strip(int, int *, int *, int **);
};

}    // namespace LAMMPS_NS

#endif
#endif

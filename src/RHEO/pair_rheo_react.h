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
PairStyle(rheo/react,PairRHEOReact)
// clang-format on
#else

#ifndef LMP_PAIR_RHEO_REACT_H
#define LMP_PAIR_RHEO_REACT_H

#include "pair.h"

namespace LAMMPS_NS {

class PairRHEOReact : public Pair {
 public:
  PairRHEOReact(class LAMMPS *);
  ~PairRHEOReact() override;
  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  void init_style() override;
  void setup() override;
  double init_one(int, int) override;
  void write_restart(FILE *) override;
  void read_restart(FILE *) override;
  int pack_reverse_comm(int, int, double *) override;
  void unpack_reverse_comm(int, int *, double *) override;

 protected:
  double **cut, **cutbond, **cutbsq, **k, **eps, **gamma, **t_form, **rlimit;

  void allocate();
  void transfer_history(double*, double*);

  int size_history;
  int *dbond, *nbond;

  int index_nb, nmax_store;
  char *id_fix;

  class FixDummy *fix_dummy;
  class FixNeighHistory *fix_history;
  class FixRHEO *fix_rheo;
  class ComputeRHEOSurface *compute_surface;
};

}    // namespace LAMMPS_NS

#endif
#endif

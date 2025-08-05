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

#ifdef FIX_CLASS
// clang-format off
FixStyle(continuum/chunk,FixContinuumChunk);
// clang-format on
#else

#ifndef LMP_FIX_CONTINUUM_CHUNK_H
#define LMP_FIX_CONTINUUM_CHUNK_H

#include "fix.h"

namespace LAMMPS_NS {

class FixContinuumChunk : public Fix {
 public:
  FixContinuumChunk(class LAMMPS *, int, char **);
  ~FixContinuumChunk() override;
  int setmask() override;
  void init() override;
  void init_list(int, class NeighList *) override;
  void setup(int) override;
  void end_of_step() override;
  double compute_array(int, int) override;
  double memory_usage() override;

 private:
  std::vector<std::pair<int, int>> values;

  int dim, pstyle, single_needed;
  double w_cut, w_cut_sq, w_sd, w_sd_sq, w_scale, w_offset;

  int nvalues, nrepeat, nfreq, irepeat;
  int borderflag, overwrite, colextra;
  bigint nvalid, nvalid_last;
  char *format, *format_user;
  FILE *fp;

  class NeighList *list;

  int volflag;        // SCALAR/VECTOR for density normalization by volume
  double chunk_volume_scalar;
  double *chunk_volume_vec;

  int ave, nwindow;
  int normcount, iwindow, window_limit;

  int nchunk, maxchunk;
  char *idchunk;
  class ComputeChunkAtom *cchunk;
  int lockforever;

  bigint filepos;

  int maxvar;
  double *varatom;

  // one,many,sum vecs/arrays are used with a single Nfreq epoch
  // total,list vecs/arrays are used across epochs

  double *count_one, *count_many, *count_sum;
  double **values_one, **values_many, **values_sum;
  double *count_total, **count_list;
  double **values_total, ***values_list;

  void allocate();
  bigint nextvalid();
  inline double calc_w(const double) const ;
  inline double calc_w_int(const double, const double, const double) const ;
  void get_chunk_center(int, int*, double*, double*, double*);
};

}    // namespace LAMMPS_NS

#endif
#endif

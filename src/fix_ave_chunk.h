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
FixStyle(ave/chunk,FixAveChunk);
// clang-format on
#else

#ifndef LMP_FIX_AVE_CHUNK_H
#define LMP_FIX_AVE_CHUNK_H

#include "fix.h"

namespace LAMMPS_NS {

class FixAveChunk : public Fix {
 public:
  FixAveChunk(class LAMMPS *, int, char **);
  ~FixAveChunk() override;
  int setmask() override;
  void init() override;
  void setup(int) override;
  void end_of_step() override;
  double compute_array(int, int) override;
  double memory_usage() override;

 private:
  struct value_t {
    int which;         // type of data: COMPUTE, FIX, VARIABLE
    int argindex;      // 1-based index if data is vector, else 0
    std::string id;    // compute/fix/variable ID
    union {
      class Compute *c;
      class Fix *f;
      int v;
    } val;
  };
  std::vector<value_t> values;

  int nvalues, nrepeat, nfreq, irepeat;
  int normflag, scaleflag, overwrite, biasflag, colextra;
  bigint nvalid, nvalid_last;
  double adof, cdof;
  char *format, *format_user;
  char *tstring, *sstring, *id_bias;
  class Compute *tbias;    // ptr to additional bias compute
  FILE *fp;

  int densityflag;    // 1 if density/number or density/mass requested
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
};

}    // namespace LAMMPS_NS

#endif
#endif

/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(meam,PairMEAM);
PairStyle(meam/c,PairMEAM);
// clang-format on
#else

#ifndef LMP_PAIR_MEAM_H
#define LMP_PAIR_MEAM_H

#include "pair.h"

namespace LAMMPS_NS {

class PairMEAM : public Pair {
 public:
  PairMEAM(class LAMMPS *);
  ~PairMEAM();
  void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  void init_style();
  void init_list(int, class NeighList *);
  double init_one(int, int);
  virtual void *extract(const char *, int &);

  int pack_forward_comm(int, int *, double *, int, int *);
  void unpack_forward_comm(int, int, double *);
  int pack_reverse_comm(int, int, double *);
  void unpack_reverse_comm(int, int *, double *);
  double memory_usage();

 private:
  class MEAM *meam_inst;
  double cutmax;                           // max cutoff for all elements
  int nlibelements;                        // # of library elements
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

/* ERROR/WARNING messages:

E: MEAM library error %d

A call to the MEAM Fortran library returned an error.

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

E: Pair style MEAM requires newton pair on

See the newton command.  This is a restriction to use the MEAM
potential.

E: Cannot open MEAM potential file %s

The specified MEAM potential file cannot be opened.  Check that the
path and name are correct.

E: Incorrect format in MEAM library file

Incorrect number of words per line in the potential file.

E: Too many elements extracted from MEAM library.

Increase 'maxelt' in meam.h and recompile.

E: Unrecognized lattice type in MEAM library/parameter file

The lattice type in an entry of the MEAM library/parameter file is not
valid.

E: Unsupported parameter in MEAM library file: ...

Self-explanatory.

E: Mismatched parameter in MEAM library file: z!=lat

The coordination number and lattice do not match, check that consistent values are given.

E: Did not find all elements in MEAM library file

Some requested elements were not found in the MEAM file. Check spelling etc.

E: Keyword %s in MEAM parameter file not recognized

Self-explanatory.

E: Error in MEAM parameter file: keyword %s (further information)

Self-explanatory. Check the parameter file.

*/

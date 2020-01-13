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

PairStyle(meam/c,PairMEAMC)
PairStyle(meam,PairMEAMC)

#else

#ifndef LMP_PAIR_MEAMC_H
#define LMP_PAIR_MEAMC_H

#include "pair.h"

namespace LAMMPS_NS {

class PairMEAMC : public Pair {
 public:
  PairMEAMC(class LAMMPS *);
  ~PairMEAMC();
  void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  void init_style();
  void init_list(int, class NeighList *);
  double init_one(int, int);

  int pack_forward_comm(int, int *, double *, int, int *);
  void unpack_forward_comm(int, int, double *);
  int pack_reverse_comm(int, int, double *);
  void unpack_reverse_comm(int, int *, double *);
  double memory_usage();

 private:
  class MEAM *meam_inst;
  double cutmax;                // max cutoff for all elements
  int nelements;                // # of unique elements
  char **elements;              // names of unique elements
  double *mass;                 // mass of each element

  int *map;                     // mapping from atom types (1-indexed) to elements (1-indexed)

  void allocate();
  void read_files(char *, char *);
  void neigh_strip(int, int *, int *, int **);
};

}

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

E: Incorrect format in MEAM potential file

Incorrect number of words per line in the potential file.

E: Unrecognized lattice type in MEAM file 1

The lattice type in an entry of the MEAM library file is not
valid.

E: Did not find all elements in MEAM library file

The requested elements were not found in the MEAM file.

E: Keyword %s in MEAM parameter file not recognized

Self-explanatory.

E: Unrecognized lattice type in MEAM file 2

The lattice type in an entry of the MEAM parameter file is not
valid.

*/

/* ----------------------------------------------------------------------
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

PairStyle(sw/table,PairSWTable)

#else

#ifndef LMP_PAIR_SW_TABLE_H
#define LMP_PAIR_SW_TABLE_H

#include "pair_sw.h"

namespace LAMMPS_NS {

class PairSWTable : public PairSW {
 public:
  PairSWTable(class LAMMPS *);
  ~PairSWTable();
  void compute(int, int);
  void settings(int, char **);
  double memory_usage();

 protected:
  int ntable;
  double deltaR2;
  double oneOverDeltaR2;
  double ***forceTable;         // table of forces per element pair
  double ***potentialTable;     // table of potential energies

  int neigh3BodyMax;            // max size of short neighborlist
  int *neigh3BodyCount;         // # of neighbors in short range 
                                // 3 particle forces neighbor list
  int **neigh3Body;             // neighlist for short range 3 particle forces

  void twobody_table(Param &, double, double &, int, double &);
  void setup_params();
  void create_tables();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

E: Pair style Stillinger-Weber requires atom IDs

This is a requirement to use the Stillinger-Weber potential.

E: Pair style Stillinger-Weber requires newton pair on

See the newton command.  This is a restriction to use the Stillinger-Weber
potential.

E: All pair coeffs are not set

All pair coefficients must be set in the data file or by the
pair_coeff command before running a simulation.

E: Cannot open Stillinger-Weber potential file %s

The specified Stillinger-Weber potential file cannot be opened.  Check that the path
and name are correct.

E: Incorrect format in Stillinger-Weber potential file

Incorrect number of words per line in the potential file.

E: Illegal Stillinger-Weber parameter

One or more of the coefficients defined in the potential file is
invalid.

E: Potential file has duplicate entry

The potential file has more than one entry for the same element.

E: Potential file is missing an entry

The potential file does not have a needed entry.

*/

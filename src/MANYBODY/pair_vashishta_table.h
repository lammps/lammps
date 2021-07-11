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
PairStyle(vashishta/table,PairVashishtaTable);
// clang-format on
#else

#ifndef LMP_PAIR_VASHISHITA_TABLE_H
#define LMP_PAIR_VASHISHITA_TABLE_H

#include "pair_vashishta.h"

namespace LAMMPS_NS {

class PairVashishtaTable : public PairVashishta {
 public:
  PairVashishtaTable(class LAMMPS *);
  ~PairVashishtaTable();
  void compute(int, int);
  void settings(int, char **);
  double memory_usage();

 protected:
  int ntable;
  double deltaR2;
  double oneOverDeltaR2;
  double ***forceTable;        // table of forces per element pair
  double ***potentialTable;    // table of potential energies

  void twobody_table(const Param &, double, double &, int, double &);
  void setup_params();
  void create_tables();
};

}    // namespace LAMMPS_NS

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Illegal inner cutoff for tabulation

UNDOCUMENTED

U: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

U: Pair style Vashishta requires atom IDs

This is a requirement to use the Vashishta potential.

U: Pair style Vashishta requires newton pair on

See the newton command.  This is a restriction to use the Vashishta
potential.

U: All pair coeffs are not set

All pair coefficients must be set in the data file or by the
pair_coeff command before running a simulation.

U: Cannot open Vashishta potential file %s

The specified Vashishta potential file cannot be opened.  Check that the path
and name are correct.

U: Incorrect format in Vashishta potential file

Incorrect number of words per line in the potential file.

U: Illegal Vashishta parameter

One or more of the coefficients defined in the potential file is
invalid.

U: Potential file has duplicate entry

The potential file has more than one entry for the same element.

U: Potential file is missing an entry

The potential file does not have a needed entry.

*/

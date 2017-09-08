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

PairStyle(table/rx,PairTableRX)

#else

#ifndef LMP_PAIR_TABLE_RX_H
#define LMP_PAIR_TABLE_RX_H

#include "pair_table.h"

namespace LAMMPS_NS {

class PairTableRX : public PairTable {
 public:
  PairTableRX(class LAMMPS *);
  virtual ~PairTableRX();

  virtual void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  virtual double single(int, int, int, int, double, double, double, double &);

 protected:

  int nspecies;
  char *site1, *site2;
  int isite1, isite2;
  void getMixingWeights(int, double &, double &, double &, double &);
  bool fractionalWeighting;

};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Pair distance < table inner cutoff

Two atoms are closer together than the pairwise table allows.

E: Pair distance > table outer cutoff

Two atoms are further apart than the pairwise table allows.

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Unknown table style in pair_style command

Style of table is invalid for use with pair_style table command.

E: PairTableRX requires a fix rx command

The fix rx command must come before the pair style command in the input file

E:  There are no rx species specified

There must be at least one species specified through the fix rx command

E:  Site1 name not recognized in pair coefficients

The site1 keyword does not match the species keywords specified throug the fix rx command

E: Illegal number of pair table entries

There must be at least 2 table entries.

E: Invalid pair table length

Length of read-in pair table is invalid

E: Invalid pair table cutoff

Cutoffs in pair_coeff command are not valid with read-in pair table.

E: Bitmapped table in file does not match requested table

Setting for bitmapped table in pair_coeff command must match table
in file exactly.

E: All pair coeffs are not set

All pair coefficients must be set in the data file or by the
pair_coeff command before running a simulation.

E: Cannot open file %s

The specified file cannot be opened.  Check that the path and name are
correct. If the file is a compressed file, also check that the gzip
executable can be found and run.

E: Did not find keyword in table file

Keyword used in pair_coeff command was not found in table file.

E: Bitmapped table is incorrect length in table file

Number of table entries is not a correct power of 2.

E: Invalid keyword in pair table parameters

Keyword used in list of table parameters is not recognized.

E: Pair table parameters did not set N

List of pair table parameters must include N setting.

E: Pair table cutoffs must all be equal to use with KSpace

When using pair style table with a long-range KSpace solver, the
cutoffs for all atom type pairs must all be the same, since the
long-range solver starts at that cutoff.

E:  The number of molecules in CG particle is less than 10*DBL_EPSILON

Self-explanatory.  Check the species concentrations have been properly set
and check the reaction kinetic solver parameters in fix rx to more for
sufficient accuracy.

*/

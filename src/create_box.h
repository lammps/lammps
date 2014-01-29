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

#ifdef COMMAND_CLASS

CommandStyle(create_box,CreateBox)

#else

#ifndef LMP_CREATE_BOX_H
#define LMP_CREATE_BOX_H

#include "pointers.h"

namespace LAMMPS_NS {

class CreateBox : protected Pointers {
 public:
  CreateBox(class LAMMPS *);
  void command(int, char **);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Dump file MPI-IO output not allowed with '%' in filename

UNDOCUMENTED

E: Cannot dump sort when multiple procs write the dump file

UNDOCUMENTED

E: Cannot dump sort on atom IDs with no atom IDs defined

Self-explanatory.

E: Dump sort column is invalid

Self-explanatory.

E: Too many atoms to dump sort

Cannot sort when running with more than 2^31 atoms.

E: Too much per-proc info for dump

Number of local atoms times number of columns must fit in a 32-bit
integer for dump.

E: Too much buffered per-proc info for dump

UNDOCUMENTED

E: Cannot open gzipped file

LAMMPS was compiled without support for reading and writing gzipped
files through a pipeline to the gzip program with -DLAMMPS_GZIP.

E: Cannot open dump file

The output file for the dump command cannot be opened.  Check that the
path and name are correct.

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Dump_modify buffer yes not allowed for this style

UNDOCUMENTED

E: Cannot use dump_modify fileper without % in dump file name

UNDOCUMENTED

E: Cannot use dump_modify nfile without % in dump file name

UNDOCUMENTED

*/

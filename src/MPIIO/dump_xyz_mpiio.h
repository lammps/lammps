/* ---------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef DUMP_CLASS

DumpStyle(xyz/mpiio,DumpXYZMPIIO)

#else

#ifndef LMP_DUMP_XYZ_MPIIO_H
#define LMP_DUMP_XYZ_MPII0_H

#include "dump_xyz.h"

namespace LAMMPS_NS {

class DumpXYZMPIIO : public DumpXYZ {
 public:
  DumpXYZMPIIO(class LAMMPS *, int, char **);
  virtual ~DumpXYZMPIIO();

 protected:

  bigint sumFileSize;  // size in bytes of the file up through this rank offset from the end of the header data
  char *headerBuffer; // buffer for holding header data

  MPI_File mpifh;
  MPI_Offset mpifo,offsetFromHeader,headerSize, currentFileSize;
  int performEstimate; // switch for write_data and write_header methods to use for gathering data and detemining filesize for preallocation vs actually writing the data
  char *filecurrent;  // name of file for this round (with % and * replaced)

#if defined(_OPENMP)
  int convert_string_omp(int, double *);  // multithreaded version of convert_string
#endif

  virtual void openfile();
  virtual void write_header(bigint);
  virtual void write();
  virtual void write_data(int, double *);

  virtual void init_style();

  typedef void (DumpXYZMPIIO::*FnPtrData)(int, double *);
  FnPtrData write_choice;              // ptr to write data functions
  void write_string(int, double *);

};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Cannot open dump file %s

UNDOCUMENTED

E: Too much per-proc info for dump

UNDOCUMENTED

U: No dump custom arguments specified

The dump custom command requires that atom quantities be specified to
output to dump file.

U: Invalid attribute in dump custom command

Self-explantory.

U: Dump_modify format string is too short

There are more fields to be dumped in a line of output than
your format string specifies.

U: Could not find dump custom compute ID

The compute ID needed by dump custom to compute a per-atom quantity
does not exist.

U: Could not find dump custom fix ID

Self-explanatory.

U: Dump custom and fix not computed at compatible times

The fix must produce per-atom quantities on timesteps that dump custom
needs them.

U: Could not find dump custom variable name

Self-explanatory.

U: Region ID for dump custom does not exist

Self-explanatory.

U: Threshhold for an atom property that isn't allocated

A dump threshhold has been requested on a quantity that is
not defined by the atom style used in this simulation.

U: Dumping an atom property that isn't allocated

The chosen atom style does not define the per-atom quantity being
dumped.

U: Dumping an atom quantity that isn't allocated

Only per-atom quantities that are defined for the atom style being
used are allowed.

U: Dump custom compute does not compute per-atom info

Self-explanatory.

U: Dump custom compute does not calculate per-atom vector

Self-explanatory.

U: \n

UNDOCUMENTED

U: Dump custom compute vector is accessed out-of-range

Self-explanatory.

U: Dump custom fix does not compute per-atom info

Self-explanatory.

U: Dump custom fix does not compute per-atom vector

Self-explanatory.

U: Dump custom fix does not compute per-atom array

Self-explanatory.

U: Dump custom fix vector is accessed out-of-range

Self-explanatory.

U: Dump custom variable is not atom-style variable

Only atom-style variables generate per-atom quantities, needed for
dump output.

U: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

U: Dump_modify region ID does not exist

Self-explanatory.

U: Dump modify element names do not match atom types

Number of element names must equal number of atom types.

U: Invalid attribute in dump modify command

Self-explantory.

U: Could not find dump modify compute ID

Self-explanatory.

U: Dump modify compute ID does not compute per-atom info

Self-explanatory.

U: Dump modify compute ID does not compute per-atom vector

Self-explanatory.

U: Dump modify compute ID does not compute per-atom array

Self-explanatory.

U: Dump modify compute ID vector is not large enough

Self-explanatory.

U: Could not find dump modify fix ID

Self-explanatory.

U: Dump modify fix ID does not compute per-atom info

Self-explanatory.

U: Dump modify fix ID does not compute per-atom vector

Self-explanatory.

U: Dump modify fix ID does not compute per-atom array

Self-explanatory.

U: Dump modify fix ID vector is not large enough

Self-explanatory.

U: Could not find dump modify variable name

Self-explanatory.

U: Dump modify variable is not atom-style variable

Self-explanatory.

U: Invalid dump_modify threshhold operator

Operator keyword used for threshold specification in not recognized.

U: Dump custom compute does not calculate per-atom array

Self-explanatory.

*/

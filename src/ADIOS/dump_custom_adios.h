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

#ifdef DUMP_CLASS
// clang-format off
DumpStyle(custom/adios, DumpCustomADIOS);
// clang-format on
#else

#ifndef LMP_DUMP_CUSTOM_ADIOS_H
#define LMP_DUMP_CUSTOM_ADIOS_H

#include "dump_custom.h"

namespace LAMMPS_NS {

class DumpCustomADIOSInternal;

class DumpCustomADIOS : public DumpCustom {
 public:
  DumpCustomADIOS(class LAMMPS *, int, char **);
  virtual ~DumpCustomADIOS();

 protected:
  virtual void openfile();
  virtual void write();
  virtual void init_style();

 private:
  DumpCustomADIOSInternal *internal;
};
}    // namespace LAMMPS_NS

#endif
#endif

/* ERROR/WARNING messages:

    E: Cannot open dump file %s

    The output file for the dump command cannot be opened.  Check that the
    path and name are correct.

    E: Too much per-proc info for dump

    Number of local atoms times number of columns must fit in a 32-bit
    integer for dump.

    E: Dump_modify format string is too short

    There are more fields to be dumped in a line of output than your
    format string specifies.

    E: Could not find dump custom compute ID

    Self-explanatory.

    E: Could not find dump custom fix ID

    Self-explanatory.

    E: Dump custom and fix not computed at compatible times

    The fix must produce per-atom quantities on timesteps that dump custom
    needs them.

    E: Could not find dump custom variable name

    Self-explanatory.

    E: Region ID for dump custom does not exist

    Self-explanatory.

    */

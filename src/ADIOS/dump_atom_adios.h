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
DumpStyle(atom/adios, DumpAtomADIOS);
// clang-format on
#else

#ifndef LMP_DUMP_ATOM_ADIOS_H
#define LMP_DUMP_ATOM_ADIOS_H

#include "dump_atom.h"

namespace LAMMPS_NS {

class DumpAtomADIOSInternal;

class DumpAtomADIOS : public DumpAtom {

 public:
  DumpAtomADIOS(class LAMMPS *, int, char **);
  virtual ~DumpAtomADIOS();

 protected:
  virtual void openfile();
  virtual void write();
  virtual void init_style();

 private:
  DumpAtomADIOSInternal *internal;
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

    */

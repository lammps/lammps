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

/* ----------------------------------------------------------------------
   Contributing author: Richard Berger (Temple U)
------------------------------------------------------------------------- */

#ifdef LAMMPS_ZSTD

#ifdef DUMP_CLASS
// clang-format off
DumpStyle(atom/zstd,DumpAtomZstd);
// clang-format on
#else

#ifndef LMP_DUMP_ATOM_ZSTD_H
#define LMP_DUMP_ATOM_ZSTD_H

#include "dump_atom.h"
#include "zstd_file_writer.h"

namespace LAMMPS_NS {

class DumpAtomZstd : public DumpAtom {
 public:
  DumpAtomZstd(class LAMMPS *, int, char **);
  virtual ~DumpAtomZstd();

 protected:
  ZstdFileWriter writer;

  virtual void openfile();
  virtual void write_header(bigint);
  virtual void write_data(int, double *);
  virtual void write();

  virtual int modify_param(int, char **);
};

}    // namespace LAMMPS_NS

#endif
#endif
#endif

/* ERROR/WARNING messages:

E: Dump atom/zstd only writes compressed files

The dump atom/zstd output file name must have a .zst suffix.

E: Cannot open dump file

Self-explanatory.

*/

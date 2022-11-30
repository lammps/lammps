/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef COMMAND_CLASS
// clang-format off
CommandStyle(RESET_ATOMS_MOL,ResetAtomsMol);
// clang-format on
#else

#ifndef LMP_RESET_ATOMS_MOL_H
#define LMP_RESET_ATOMS_MOL_H

#include "command.h"

namespace LAMMPS_NS {

class ResetAtomsMol : public Command {
 public:
  ResetAtomsMol(class LAMMPS *);
  ~ResetAtomsMol() override;
  void command(int, char **) override;
  void create_computes(char *, char *);
  void reset();

 private:
  std::string idfrag, idchunk;
  int nchunk;
  int groupbit;
  int compressflag;    // 1 = contiguous values for new IDs
  int singleflag;      // 0 = mol IDs of single atoms set to 0
  tagint offset;       // offset for contiguous mol ID values

  class ComputeFragmentAtom *cfa;
  class ComputeChunkAtom *cca;
};
}    // namespace LAMMPS_NS

#endif
#endif

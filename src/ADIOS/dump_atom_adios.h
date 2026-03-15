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

#ifdef DUMP_CLASS
// clang-format off
DumpStyle(atom/adios, DumpAtomADIOS);
// clang-format on
#else

#ifndef LMP_DUMP_ATOM_ADIOS_H
#define LMP_DUMP_ATOM_ADIOS_H

#include "dump_atom.h"

#include <memory>    // std::unique_ptr

namespace LAMMPS_NS {

// Forward-declare the pimpl type; defined fully in the .cpp so that
// adios2.h is not exposed in this header.
class DumpAtomADIOSInternal;

class DumpAtomADIOS : public DumpAtom {

 public:
  DumpAtomADIOS(class LAMMPS *, int, char **);

  // Defined in the .cpp where DumpAtomADIOSInternal is complete,
  // which is required for unique_ptr<incomplete type> to compile.
  ~DumpAtomADIOS() override;

  // Non-copyable: the pimpl holds ADIOS2 handle objects that must not
  // be duplicated.
  DumpAtomADIOS(const DumpAtomADIOS &) = delete;
  DumpAtomADIOS &operator=(const DumpAtomADIOS &) = delete;

 protected:
  void openfile() override;
  void write() override;
  void init_style() override;

 private:
  std::unique_ptr<DumpAtomADIOSInternal> internal;
};

}    // namespace LAMMPS_NS

#endif
#endif

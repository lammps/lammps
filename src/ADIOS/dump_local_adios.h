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
DumpStyle(local/adios, DumpLocalADIOS);
// clang-format on
#else

#ifndef LMP_DUMP_LOCAL_ADIOS_H
#define LMP_DUMP_LOCAL_ADIOS_H

#include "dump_local.h"

#include <memory>    // std::unique_ptr

namespace LAMMPS_NS {

// Forward-declare the pimpl type; defined fully in the .cpp so that
// adios2.h is not exposed in this header.
class DumpLocalADIOSInternal;

class DumpLocalADIOS : public DumpLocal {
 public:
  DumpLocalADIOS(class LAMMPS *, int, char **);

  // Defined in the .cpp where DumpLocalADIOSInternal is complete,
  // which is required for unique_ptr<incomplete type> to compile.
  ~DumpLocalADIOS() override;

  // Non-copyable: the pimpl holds ADIOS2 handle objects that must not
  // be duplicated.
  DumpLocalADIOS(const DumpLocalADIOS &) = delete;
  DumpLocalADIOS &operator=(const DumpLocalADIOS &) = delete;

 protected:
  int  modify_param(int, char **) override;
  void openfile() override;
  void write() override;
  void init_style() override;

 private:
  std::unique_ptr<DumpLocalADIOSInternal> internal;

  // dump_modify adios_precision fp32|fp64  (default: fp64)
  bool adios_use_float_{false};

  // dump_modify adios_shared_replicas yes|no  (default: no)
  // When yes and universe->nworlds > 1, all partitions write to the same
  // .bp file using universe->uworld.  Each replica writes to its own
  // per-replica ADIOS variable (e.g. local_data_replica0, ...).
  bool adios_shared_replicas_{false};
};

}    // namespace LAMMPS_NS

#endif
#endif

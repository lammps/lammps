/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef LMP_RESTARTABLE_H
#define LMP_RESTARTABLE_H

#include <cstddef>
#include <vector>

#include "pointers.h"
#include "lmptype.h"

namespace LAMMPS_NS {

class BufferReader;
class FileWriter;

class Restartable : protected Pointers {
 public:
  Restartable(class LAMMPS *lmp) : Pointers(lmp) { };
  virtual ~Restartable() = default;

  virtual void write_restart_global(FileWriter *) const;
  virtual void write_restart_local(FileWriter *) const;
  virtual void write_restart_settings(FileWriter *) const; // Not used by Fixes!

  virtual void read_restart_global(BufferReader);
  virtual void read_restart_local(BufferReader);
  virtual void read_restart_settings(BufferReader *); // Not used by Fixes!

  bool restartable_global = false; // If obj writes global data with new API
  bool restartable_local  = false; // If obj writes local data with new API

  // If obj is restartable by latest API
  [[nodiscard]] bool restartable() const {
    return restartable_global || restartable_local;
  }
  // If obj is restartable by FILE* API (a superset of restartable)
  [[nodiscard]] virtual bool fp_restartable() const {
    return restartable();
  }

  // Backwards compatibility functions
  virtual void write_restart(FILE *);
  virtual void write_restart_settings(FILE *);
  virtual void read_restart(FILE *);
  virtual void read_restart_settings(FILE *);
  virtual void restart(char *);

  // If write_restart should prefix the write with the integer size of
  // the data (backwards compatibility with Fix)
  bool write_restart_size_prefix = false;
};

} // namespace LAMMPS_NS

#endif

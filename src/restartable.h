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
class BufferReaderRootFile;
class FileWriter;

class Restartable : protected Pointers {
 public:
  Restartable(class LAMMPS *lmp) : Restartable(lmp, 0) {}
  Restartable(class LAMMPS *lmp, const int& fp_flag)
    : Pointers(lmp), restartable_file_flag(fp_flag) {}
  virtual ~Restartable() = default;

  virtual void write_restart_global(FileWriter&) const;
  virtual void write_restart_local(FileWriter&) const;
  virtual void write_restart_settings(FileWriter&) const; // Not used by Fixes!

  virtual void read_restart_global(BufferReader&);
  virtual void read_restart_local(BufferReader&);
  virtual void read_restart_settings(BufferReader&); // Not used by Fixes!

  // When restarting with fewer procs, extra proc data will be passed to some
  // existing proc(s) (typically all to rank 0). Other procs will still invoke
  // this function, but with n_extra_procs = 0
  virtual void read_restart_local_extra(int n_extra_procs, BufferReader&);

  bool restartable_global = false; // If obj writes global data with new API
  bool restartable_local  = false; // If obj writes local data with new API

  // If obj is restartable by latest API
  [[nodiscard]] bool restartable() const {
    return restartable_global || restartable_local;
  }
  // If obj is restartable by FILE* API. Since the FILE* api can be emulated by
  // the latest API, this is a superset of restartable
  [[nodiscard]] bool restartable_file() const {
    return restartable_global || restartable_local || restartable_file_flag;
  }

  // Backwards compatibility functions
  virtual void write_restart(FILE *);
  virtual void write_restart_settings(FILE *);
  virtual void read_restart(BufferReaderRootFile&);
  virtual void read_restart(FILE *);
  virtual void read_restart_settings(FILE *);
  virtual void restart(char *);

  // If write_restart should prefix the write with the integer size of
  // the data (backwards compatibility with Fix)
  bool write_restart_size_prefix = false;

  // Helpers for reading/writing to temporary variables (e.g. sub_bufs)
  void write_restart_global(FileWriter&&);
  void write_restart_local(FileWriter&&);
  void write_restart_settings(FileWriter&&);
  void read_restart_global(BufferReader&&);
  void read_restart_local(BufferReader&&);
  void read_restart_settings(BufferReader&&);
  void read_restart_local_extra(int count, BufferReader&&);

  static constexpr int ALWAYS_FILE_RESTARTABLE = 1;

 private:
  const int& restartable_file_flag; // 1 if restartable with the FILE* api
};

} // namespace LAMMPS_NS

#endif

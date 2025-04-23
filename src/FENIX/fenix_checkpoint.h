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

#ifndef LMP_FENIX_CHECKPOINT_H
#define LMP_FENIX_CHECKPOINT_H

#include "write_restart.h"
#include "error.h"

#include <string>

namespace LAMMPS_NS {

class Fenix;

class FenixCheckpoint : public WriteRestart {
 public:
  FenixCheckpoint(class LAMMPS *);
  void command(int, char**) override;

  void write(const std::string &) final;

  void checkpoint();
  void recover();

  ~FenixCheckpoint();

  using Buffer = std::pair<char*, char*>;

  template<typename T>
  static T read(Buffer&, Error*);
  template<typename T>
  static void read(T&, Buffer&, Error*);
  template<typename T>
  static void read(T*, int, Buffer&, Error*);

 protected:
  friend Fenix;

  void create_group();

  void store_meta();
  void store_global();
  void store_peratom();

  void commit();

  void restore_meta();
  void restore_global();
  void restore_peratom();

 private:
  void safe_create_member(int, void*, int, MPI_Datatype);
  void safe_store(int, int);
  void safe_storev(int, int);

  void safe_restore(int, void*, int);
  int safe_restore(int, char*&);

  void proc_grid();

  void magic_string(Buffer&);
  void header(Buffer&);
  void proc_grid(Buffer&);
  void recover_group(Buffer&);
  void type_arrays(Buffer&);
  void force_fields(Buffer&);
  int recover_modify(Buffer&);

  int group_id;
  int global_member;
  int peratom_member;
  int meta_member;

  int group_policy;
  int* group_policy_args;
};

template<typename T>
void FenixCheckpoint::read(T& var, Buffer& buf, Error* error){
  static_assert(std::is_trivially_copyable_v<T>,
    "This type needs a specialized read function"
  );
  static_assert(!std::is_pointer_v<T>,
    "Cannot (de)serialize general pointers directly"
  );

  if(buf.first + sizeof(T) > buf.second)
    error->one(FLERR, "Checkpointed data invalid");

  var = *(T*)buf.first;
  buf.first += sizeof(T);
}

template<>
inline void FenixCheckpoint::read<std::string>(
  std::string& var, Buffer& buf, Error* error
) {
  char* str_end = buf.first;
  for(; str_end < buf.second; str_end++) if(*str_end == '\0') break;

  if(str_end >= buf.second)
    error->one(FLERR, "Checkpointed data invalid");

  var = buf.first;
  buf.first = str_end+1;
}

template<>
inline void FenixCheckpoint::read<char*>(
  char*& var, Buffer& buf, Error* error
) {
  char* str_end = buf.first;
  for(; str_end < buf.second; str_end++) if(*str_end == '\0') break;

  if(str_end >= buf.second)
    error->one(FLERR, "Checkpointed data invalid");

  int len = str_end+1 - buf.first;
  var = new char[len];
  read(var, len, buf, error);
}

template<typename T>
T FenixCheckpoint::read(Buffer& buf, Error* error){
  T var;
  read(var, buf, error);
  return var;
}

template<typename T>
void FenixCheckpoint::read(T* arr, int len, Buffer& buf, Error* error){
  if(std::is_trivially_copyable_v<T> && !std::is_pointer_v<T>){
    if(buf.first + len*sizeof(T) > buf.second)
      error->one(FLERR, "Checkpointed data invalid");
    memcpy(arr, buf.first, len*sizeof(T));
    buf.first += len*sizeof(T);
  } else {
    for(int i = 0; i < len; i++){
      arr[i] = read<T>(buf, error);
    }
  }
}

}

#endif

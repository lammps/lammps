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

#ifndef LMP_LOCAL_SERIALIZER_H
#define LMP_LOCAL_SERIALIZER_H

#include "write_restart.h"
#include "error.h"

#include <string>

namespace LAMMPS_NS {

class LocalSerializer : private WriteRestart {
 public:
  LocalSerializer(class LAMMPS *);

  // <beginning, past-the-end>
  using Buffer = std::pair<char*, char*>;

  Buffer serialize();
  void deserialize(Buffer);
  
  // Read type from buffer, incrementing beginning pointer past it
  template<typename T>
  static T read(Buffer&, Error*);
  
  // As above, but write into provided location
  template<typename T>
  static void read(T&, Buffer&, Error*);

  // As above, but read N times into provided array
  template<typename T>
  static void read(T*, int, Buffer&, Error*);

  // Deserialization helpers
  void magic_string(Buffer&);
  void header(Buffer&);
  void proc_grid(Buffer&);
  void recover_group(Buffer&);
  void type_arrays(Buffer&);
  void force_fields(Buffer&);
  int recover_modify(Buffer&);
  
 private:
  // Serialization helper
  void proc_grid();
};

template<typename T>
void LocalSerializer::read(T& var, Buffer& buf, Error* error){
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
inline void LocalSerializer::read<std::string>(
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
inline void LocalSerializer::read<char*>(
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
T LocalSerializer::read(Buffer& buf, Error* error){
  T var;
  read(var, buf, error);
  return var;
}

template<typename T>
void LocalSerializer::read(T* arr, int len, Buffer& buf, Error* error){
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

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

#include "FENIX/fenix_checkpoint.h"
#include "FENIX/fenix.h"

#include "fenix.hpp"

#include <cstdio>
#include <fmt/ranges.h>

#include "memory.h"
#include "error.h"

using namespace Fenix;
using namespace Fenix::Data;

namespace LAMMPS_NS {

/* ---------------------------------------------------------------------- */

FenixCheckpoint::FenixCheckpoint(LAMMPS *lmp)
  : WriteRestart(lmp)
{
  group_id = 0;
  data_member = 0;
  metadata_member = 1;

  group_policy = FENIX_DATA_POLICY_IMR;
  group_policy_args = nullptr;
}

/* ---------------------------------------------------------------------- */

FenixCheckpoint::~FenixCheckpoint(){
  memory->destroy(group_policy_args);
}

/* ---------------------------------------------------------------------- */

void FenixCheckpoint::create_group(){
  int ret;
  group_create(group_id, world, 0, 0, group_policy, group_policy_args, &ret);
  if(ret != FENIX_SUCCESS)
    error->one(FLERR, "Failed to create data group, err={}", ret);
}


/* ---------------------------------------------------------------------- */

void FenixCheckpoint::safe_create_member(
  int member, void* ptr, int count, MPI_Datatype datatype
){
  auto members = group_members(group_id);
  if(!members){
    error->one(FLERR, "Data group does not exist!");
  }

  for(const auto& id : members.value()){
    if(id == member){
      // Member already exists, simply update pointer to latest value
      int flag;
      Fenix_Data_member_attr_set(
        group_id, member, FENIX_DATA_MEMBER_ATTRIBUTE_BUFFER, ptr, &flag
      );
      return;
    }
  }

  int ret = member_create(group_id, member, ptr, count, datatype);

  if(ret != FENIX_SUCCESS) error->one(
    FLERR, "Failed to create data member {}, err={}", member, ret
  );
}

/* ---------------------------------------------------------------------- */

void FenixCheckpoint::safe_store(int member, int length){
  DataSubset subset = length > 0 ? DataSubset({0,length-1}) : SUBSET_EMPTY;
  int ret = member_store(group_id, member, subset);
  if(ret != FENIX_SUCCESS) error->one(
    FLERR, "Failed to store data member {}, err={}", member, ret
  );
}

/* ---------------------------------------------------------------------- */

void FenixCheckpoint::safe_storev(int member, int length){
  DataSubset subset = length > 0 ? DataSubset({0,length-1}) : SUBSET_EMPTY;
  int ret = member_storev(group_id, member, subset);
  if(ret != FENIX_SUCCESS) error->one(
    FLERR, "Failed to storev data member {}, err={}", member, ret
  );
}

/* ---------------------------------------------------------------------- */

void FenixCheckpoint::safe_restore(int member, void* ptr, int length){
  int ret = member_restore(
    group_id, member, ptr, length-1, FENIX_DATA_SNAPSHOT_LATEST, SUBSET_IGNORE
  );
  if(ret != FENIX_SUCCESS) error->one(
    FLERR, "Failed to restore data member {}, err={}", member, ret
  );
}

/* ---------------------------------------------------------------------- */

// Returns length of restored data. Caller responsible to delete out_ptr
int FenixCheckpoint::safe_restore(int member, char*& out_ptr){
  DataSubset subset;
  int ret = member_restore(
    group_id, member, nullptr, 0, FENIX_DATA_SNAPSHOT_LATEST, subset
  );
  if(ret != FENIX_SUCCESS) error->one(
    FLERR, "Failed to restore resizeable data member {}, err={}", member, ret
  );

  int length = subset.max_count();
  if(length == 0){
    out_ptr = nullptr;
    return 0;
  }
  memory->create(out_ptr, length, "fenix_checkpoint:safe_restore buf");

  ret = member_lrestore(
    group_id, member, out_ptr, length, FENIX_DATA_SNAPSHOT_LATEST, subset
  );
  if(ret != FENIX_SUCCESS) error->one(
    FLERR, "Failed to lrestore resizeable data member {}, err={}", member, ret
  );
  return length;
}

/* ---------------------------------------------------------------------- */

void FenixCheckpoint::write(const std::string&){
  checkpoint();
}

/* ---------------------------------------------------------------------- */

void FenixCheckpoint::checkpoint(){
  if(!Fenix::active_controller)
    error->one(FLERR, "Fenix has not been initialized");
  store_data();
  store_metadata();
  commit();
}

/* ---------------------------------------------------------------------- */

void FenixCheckpoint::recover(){
  restore_data();
  restore_metadata();
}

/* ---------------------------------------------------------------------- */

void FenixCheckpoint::store_data(){
  auto buf = serializer.serialize();
  safe_create_member(data_member, buf.first, FENIX_RESIZEABLE, MPI_CHAR);
  try{
    safe_storev(data_member, buf.second - buf.first);
  } catch (const std::exception& e){
    memory->destroy(buf.first);
    throw;
  }
  memory->destroy(buf.first);
}

/* ---------------------------------------------------------------------- */

void FenixCheckpoint::store_metadata(){
  auto& meta = Fenix::active_controller->metadata;
  safe_create_member(metadata_member, meta, sizeof(meta), MPI_CHAR);
  safe_store(metadata_member, sizeof(meta));
}

/* ---------------------------------------------------------------------- */

void FenixCheckpoint::commit(){
  int ret = commit_barrier(group_id);
  if(ret != FENIX_SUCCESS)
    error->one(FLERR, "Unable to commit checkpointed data, code: {}", ret);
}

/* ---------------------------------------------------------------------- */

void FenixCheckpoint::restore_metadata(){
  auto& meta = Fenix::active_controller->metadata;
  safe_restore(metadata_member, &meta, sizeof(meta));
}

/* ---------------------------------------------------------------------- */

void FenixCheckpoint::restore_data(){
  char* buf_ptr;
  int buf_len = safe_restore(data_member, buf_ptr);
  try {
    serializer.deserialize({buf_ptr, buf_ptr + buf_len});
  } catch(const std::exception& e) {
    memory->destroy(buf_ptr);
    throw;
  }
  memory->destroy(buf_ptr);
}

}

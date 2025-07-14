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
#include "platform.h"
#include "comm.h"
#include "utils.h"
#include "domain.h"
#include "group.h"
#include "force.h"
#include "update.h"
#include "modify.h"
#include "atom.h"
#include "atom_vec.h"
#include "neighbor.h"
#include "label_map.h"
#include "math_extra.h"
#include "pair.h"
#include "bond.h"
#include "angle.h"
#include "dihedral.h"
#include "improper.h"
#include "fix.h"
#include "lmprestart.h"

using namespace Fenix;
using namespace Fenix::Data;

namespace LAMMPS_NS {

/* ---------------------------------------------------------------------- */

FenixCheckpoint::FenixCheckpoint(LAMMPS *lmp)
  : WriteRestart(lmp)
{
  group_id = 0;
  global_member = 0;
  peratom_member = 1;
  meta_member = 2;

  group_policy = FENIX_DATA_POLICY_IMR;
  group_policy_args = nullptr;
}

/* ---------------------------------------------------------------------- */

FenixCheckpoint::~FenixCheckpoint(){
  memory->destroy(group_policy_args);
}

/* ---------------------------------------------------------------------- */

void FenixCheckpoint::command(int, char**) {
  error->all(FLERR, "fenix_checkpoint is not an invokable command!");
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
    group_id, member, ptr, length-1, FENIX_TIME_STAMP_MAX, SUBSET_IGNORE
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
    group_id, member, nullptr, 0, FENIX_TIME_STAMP_MAX, subset
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
  out_ptr = new char[length];

  ret = member_lrestore(
    group_id, member, out_ptr, length, FENIX_TIME_STAMP_MAX, subset
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
  if(!Fenix::active_controller) error->all(FLERR,
    "Fenix has not been initialized"
  );
  for(auto& fix : modify->get_fix_list()) if(fix->restart_file) error->all(
    FLERR, "Fenix checkpoints do not yet support fixes that write their own"
    " restart files! Fix \"{}\" unsupported", fix->style
  );

  store_global();
  store_peratom();
  store_meta();

  commit();
}

/* ---------------------------------------------------------------------- */

void FenixCheckpoint::recover(){
  restore_global();
  restore_peratom();
}

/* ---------------------------------------------------------------------- */

void FenixCheckpoint::store_global(){
  // See WriteRestart::write for more explanation of these steps

  // Open a tmp file and hope it's stored in memory
  fp = std::tmpfile();
  if (fp == nullptr)
    error->one(FLERR, "Cannot open temporary file for restart: {}", utils::getsyserror());

  if(neighbor->build_once) domain->reset_box();

  try {
    // Update protected variables inside WriteRestart
    MPI_Comm_rank(world, &me);
    MPI_Comm_size(world, &nprocs);
    natoms = atom->natoms;

    WriteRestart::magic_string();
    WriteRestart::header();
    fwrite(&atom->nlocal, sizeof(atom->nlocal), 1, fp);
    proc_grid();
    group->write_restart(fp);
    WriteRestart::type_arrays();
    WriteRestart::force_fields();

    // Hacky way to get global data on to all ranks and skip communicating
    // local data
    comm->me = 0;
    MPI_Comm old_world = world; world = MPI_COMM_SELF;
    modify->write_restart(fp);
    world = old_world;
    comm->me = me;

    WriteRestart::magic_string();
  } catch(const CommException& e){
    fclose(fp);
    fp = nullptr;
    throw;
  }

  bigint globals_size = platform::ftell(fp);

  char* buf;
  memory->create(buf, globals_size, "write_restart_fenix:globals_buf");

  platform::fseek(fp, 0);
  fread(buf, sizeof(char), globals_size, fp);

  fclose(fp);
  fp = nullptr;

  safe_create_member(global_member, buf, FENIX_RESIZEABLE, MPI_CHAR);
  try{
    safe_storev(global_member, globals_size);
  } catch (const CommException& e){
    memory->destroy(buf);
    throw;
  }

  memory->destroy(buf);
}

/* ---------------------------------------------------------------------- */

void FenixCheckpoint::store_peratom(){
  // See WriteRestart::write for more explanation of these steps

  int d_count = atom->avec->size_restart();

  double* atom_data;
  memory->create(atom_data, d_count, "write_restart_fenix:peratom_buf");
  memset(atom_data, 0, d_count*sizeof(double));

  int n = 0;
  for(int i = 0; i < atom->nlocal; i++)
    n += atom->avec->pack_restart(i, &atom_data[n]);

  // It seems safe to ignore restart_pbc_any, since we are simply restoring
  // the prior state on recovery - not restoring to a different number of ranks
  // or a different set of configs etc. But if things break, fetch the code
  // for enforcing Periodic Boundary Conditions (PBC) from write_restart

  safe_create_member(peratom_member, atom_data, FENIX_RESIZEABLE, MPI_CHAR);
  try{
    safe_storev(peratom_member, d_count*sizeof(double));
  } catch (const CommException& e){
    memory->destroy(atom_data);
    throw;
  }

  memory->destroy(atom_data);
}

/* ---------------------------------------------------------------------- */

void FenixCheckpoint::store_meta(){
  auto& meta = Fenix::active_controller->metadata;
  safe_create_member(meta_member, meta, sizeof(meta), MPI_CHAR);
  safe_store(meta_member, sizeof(meta));
}

/* ---------------------------------------------------------------------- */

void FenixCheckpoint::commit(){
  int ret = commit_barrier(group_id);
  if(ret != FENIX_SUCCESS)
    error->one(FLERR, "Unable to commit checkpointed data, code: {}", ret);
}

/* ---------------------------------------------------------------------- */

void FenixCheckpoint::restore_meta(){
  auto& meta = Fenix::active_controller->metadata;
  safe_restore(meta_member, &meta, sizeof(meta));
}

/* ---------------------------------------------------------------------- */

void FenixCheckpoint::magic_string(Buffer& buf){
  if(read<std::string>(buf, error) != MAGIC_STRING)
    error->one(FLERR, "Checkpointed data is invalid");
}


/* ----------------------------------------------------------------------
  Rewritten group->read_restart to not require a FILE*
---------------------------------------------------------------------- */

void FenixCheckpoint::recover_group(Buffer& buf){
  int ngroup = read<int>(buf, error);

  group->ngroup = 0;
  for(int i = 0; i < Group::MAX_GROUP; i++){
    delete[] group->names[i];
    group->names[i] = nullptr;
  }

  for(int i = 0; i < Group::MAX_GROUP; i++){
    int len = read<int>(buf, error);
    if(len == 0) continue;

    group->names[i] = new char[len];
    read(group->names[i], len, buf, error);

    if(group->names[i][len-1] != '\0') error->all(
      FLERR, "Invalid restore file"
    );

    if(++group->ngroup == ngroup) break;
  }
}

/* ---------------------------------------------------------------------- */

void FenixCheckpoint::header(Buffer& buf){
  using std::string;
  int flag;
  while( (flag = read<int>(buf, error)) >= 0){
    if(flag == VERSION){
      read<int>(buf, error);
      lmp->restart_ver = utils::date2num(read<string>(buf, error));

    } else if(flag == SMALLINT){
      if(read<int>(buf, error) != sizeof(smallint))
        error->all(FLERR,"Smallint setting in lmptype.h is not compatible");
    } else if(flag == IMAGEINT){
      if(read<int>(buf, error) != sizeof(imageint))
        error->all(FLERR,"Imageint setting in lmptype.h is not compatible");
    } else if(flag == TAGINT){
      if(read<int>(buf, error) != sizeof(tagint))
        error->all(FLERR,"Tagint setting in lmptype.h is not compatible");
    } else if(flag == BIGINT){
      if(read<int>(buf, error) != sizeof(bigint))
        error->all(FLERR,"Bigint setting in lmptype.h is not compatible");

    } else if(flag == UNITS){
      read<int>(buf, error);
      auto style = read<char*>(buf, error);
      if(style != update->unit_style) update->set_units(style);
      delete[] style;

    } else if(flag == NTIMESTEP){
      read(update->ntimestep, buf, error);

    } else if(flag == DIMENSION){
      domain->dimension = read<int>(buf, error);
      if(domain->dimension == 2 && domain->zperiodic == 0)
        error->all(FLERR, "Cannot run 2d simulation with non-periodic Z dimension");

    } else if(flag == NPROCS){
      int nprocs = read<int>(buf, error);
      if(nprocs != comm->nprocs) error->all(
        FLERR, "Checkpointed with {} processors, but recovering with {}",
        nprocs, comm->nprocs
      );

    } else if(flag == PROCGRID){
      int dims = read<int>(buf, error);
      int procgrid[3];
      read(procgrid, 3, buf, error);

      auto& upg = comm->user_procgrid;
      for(int i = 0; i < 3; i++)
        if(upg[i] != 0 && procgrid[i] != upg[i]){
          error->warning(
            FLERR, "Checkpointed data used different 3d processor grid"
          );
          break;
        }

    } else if(flag == NEWTON_PAIR){
      // TODO: See if we can/should just recover this?
      // ReadRestart warns if different existing setting and uses existing
      read<int>(buf, error);
    } else if(flag == NEWTON_BOND){
      read(force->newton_bond, buf, error);
      force->newton = force->newton_pair || force->newton_bond ? 1 : 0;

    } else if(flag == XPERIODIC){
      read(domain->periodicity[0], buf, error);
      domain->xperiodic = domain->periodicity[0];
    } else if(flag == YPERIODIC){
      read(domain->periodicity[1], buf, error);
      domain->xperiodic = domain->periodicity[1];
    } else if(flag == ZPERIODIC){
      read(domain->periodicity[2], buf, error);
      domain->xperiodic = domain->periodicity[2];
    } else if(flag == BOUNDARY){
      read<int>(buf, error);
      int boundary[3][2];
      read(&boundary[0][0], 6, buf, error);

      domain->nonperiodic = 0;
      for(int i = 0; i < 3; i++){
        domain->boundary[i][0] = boundary[i][0];
        domain->boundary[i][1] = boundary[i][1];

        if(!domain->periodicity[i] || domain->nonperiodic){
          if(!domain->nonperiodic) domain->nonperiodic = 1;
          if(boundary[i][0] >= 2 || boundary[i][1] >= 2)
            domain->nonperiodic = 2;
        }
      }

    } else if(flag == BOUNDMIN){
      read<int>(buf, error);
      double minbound[6];
      read(minbound, 6, buf, error);
      domain->minxlo = minbound[0]; domain->minxhi = minbound[1];
      domain->minylo = minbound[2]; domain->minyhi = minbound[3];
      domain->minzlo = minbound[4]; domain->minzhi = minbound[5];

    } else if(flag == ATOM_STYLE){
      read<int>(buf, error);
      string style = read<string>(buf, error);

      std::vector<char*> arg(read<int>(buf, error));
      for(char*& c : arg){
        read<int>(buf, error);
        read(c, buf, error);
      }

      atom->create_avec(style, arg.size(), arg.data(), 1);

      for(char*& c : arg) delete[] c;

    } else if(flag == NATOMS){
      read(atom->natoms, buf, error);
    } else if(flag == NTYPES){
      read(atom->ntypes, buf, error);
    } else if(flag == NBONDS){
      read(atom->nbonds, buf, error);
    } else if (flag == NBONDTYPES) {
      read(atom->nbondtypes, buf, error);
    } else if (flag == BOND_PER_ATOM) {
      read(atom->bond_per_atom, buf, error);
    } else if (flag == NANGLES) {
      read(atom->nangles, buf, error);
    } else if (flag == NANGLETYPES) {
      read(atom->nangletypes, buf, error);
    } else if (flag == ANGLE_PER_ATOM) {
      read(atom->angle_per_atom, buf, error);
    } else if (flag == NDIHEDRALS) {
      read(atom->ndihedrals, buf, error);
    } else if (flag == NDIHEDRALTYPES) {
      read(atom->ndihedraltypes, buf, error);
    } else if (flag == DIHEDRAL_PER_ATOM) {
      read(atom->dihedral_per_atom, buf, error);
    } else if (flag == NIMPROPERS) {
      read(atom->nimpropers, buf, error);
    } else if (flag == NIMPROPERTYPES) {
      read(atom->nimpropertypes, buf, error);
    } else if (flag == IMPROPER_PER_ATOM) {
      read(atom->improper_per_atom, buf, error);

    } else if (flag == TRICLINIC) {
      read(domain->triclinic, buf, error);
    } else if (flag == BOXLO) {
      read<int>(buf, error);
      read(domain->boxlo, 3, buf, error);
    } else if (flag == BOXHI) {
      read<int>(buf, error);
      read(domain->boxhi, 3, buf, error);
    } else if (flag == XY) {
      read(domain->xy, buf, error);
    } else if (flag == XZ) {
      read(domain->xz, buf, error);
    } else if (flag == YZ) {
      read(domain->yz, buf, error);

    } else if (flag == TRICLINIC_GENERAL) {
      read(domain->triclinic_general, buf, error);
    } else if (flag == ROTATE_G2R) {
      read<int>(buf, error);
      read(&domain->rotate_g2r[0][0], 9, buf, error);
      MathExtra::transpose3(domain->rotate_g2r,domain->rotate_r2g);

    } else if (flag == SPECIAL_LJ) {
      read<int>(buf, error);
      read(&force->special_lj[1], 3, buf, error);
    } else if (flag == SPECIAL_COUL) {
      read<int>(buf, error);
      read(&force->special_coul[1], 3, buf, error);

    } else if (flag == TIMESTEP) {
      read(update->dt, buf, error);

    } else if (flag == ATOM_ID) {
      read(atom->tag_enable, buf, error);
    } else if (flag == ATOM_MAP_STYLE) {
      read(atom->map_style, buf, error);
    } else if (flag == ATOM_MAP_USER) {
      read(atom->map_user, buf, error);
    } else if (flag == ATOM_SORTFREQ) {
      read(atom->sortfreq, buf, error);
    } else if (flag == ATOM_SORTBIN) {
      read(atom->userbinsize, buf, error);

    } else if (flag == COMM_MODE) {
      read(comm->mode, buf, error);
    } else if (flag == COMM_CUTOFF) {
      read(comm->cutghostuser, buf, error);
    } else if (flag == COMM_VEL) {
      read(comm->ghost_velocity, buf, error);

    } else if (flag == EXTRA_BOND_PER_ATOM) {
      read(atom->extra_bond_per_atom, buf, error);
    } else if (flag == EXTRA_ANGLE_PER_ATOM) {
      read(atom->extra_angle_per_atom, buf, error);
    } else if (flag == EXTRA_DIHEDRAL_PER_ATOM) {
      read(atom->extra_dihedral_per_atom, buf, error);
    } else if (flag == EXTRA_IMPROPER_PER_ATOM) {
      read(atom->extra_improper_per_atom, buf, error);
    } else if (flag == ATOM_MAXSPECIAL) {
      read(atom->maxspecial, buf, error);

    } else if (flag == NELLIPSOIDS) {
      read(atom->nellipsoids, buf, error);
    } else if (flag == NLINES) {
      read(atom->nlines, buf, error);
    } else if (flag == NTRIS) {
      read(atom->ntris, buf, error);
    } else if (flag == NBODIES) {
      read(atom->nbodies, buf, error);

    } else if (flag == ATIMESTEP) {
      read(update->atimestep, buf, error);
    } else if (flag == ATIME) {
      read(update->atime, buf, error);

    } else error->all(FLERR,"Invalid flag in header section of restart file");
  }
}

/* ---------------------------------------------------------------------- */

void lmap_type(
  std::vector<std::string>& labels, int ntypes,
  std::unordered_map<std::string, int>& map,
  FenixCheckpoint::Buffer& buf, Error* error
) {
  for(int i = 0; i < ntypes; i++){
    int n = FenixCheckpoint::read<int>(buf, error);
    if(n == 0){
      labels[i] = new char[0];
      continue;
    }
    FenixCheckpoint::read(labels[i], buf, error);
    if(!labels[i].empty()) map[labels[i]] = i+1;
  }
}

/* ---------------------------------------------------------------------- */

void FenixCheckpoint::type_arrays(Buffer& buf){
  int flag;
  while( (flag = read<int>(buf, error)) >= 0){
    if(flag == MASS){
      read<int>(buf, error);
      std::vector<double> mass(atom->ntypes+1);
      read(&mass[1], atom->ntypes, buf, error);
      atom->set_mass(mass.data());

    } else if(flag == LABELMAP){
      read<int>(buf, error);
      atom->add_label_map();

      auto& l = *atom->lmap;
      lmap_type(l.typelabel, l.natomtypes, l.typelabel_map, buf, error);
      lmap_type(l.btypelabel, l.nbondtypes, l.btypelabel_map, buf, error);
      lmap_type(l.atypelabel, l.nangletypes, l.atypelabel_map, buf, error);
      lmap_type(l.dtypelabel, l.ndihedraltypes, l.dtypelabel_map, buf, error);
      lmap_type(l.itypelabel, l.nimpropertypes, l.itypelabel_map, buf, error);

    } else error->all(FLERR, "Invalid flag in checkpointed type arrays");
  }
}

/* ----------------------------------------------------------------------
  The first place we can't recover locally, since force fields take a FILE*
  directly for recovery. They also use a huge number of single-value file
  writes and MPI_Bcasts and are generally very gross for local recovery.

  On the todo list for improving, but it'll take some larger changes to LAMMPS.
  Probably manageable in a backwards compatible manner though.

  First attempt to work around the MPI communication:
    set world to MPI_COMM_SELF and comm->me to 0, so that all ranks read data
    from the input file and broadcasts are NOOPs
---------------------------------------------------------------------- */

void FenixCheckpoint::force_fields(Buffer& buf){
  int flag;
  while( (flag = read<int>(buf, error)) >= 0) {
    read<int>(buf, error);
    auto style = read<char*>(buf, error);

    FILE* fp = nullptr;
    if(flag != NO_PAIR){
      // POSIX, so should be generally portable to non-windows
      // but there's not a good windows alternative tmk
      fp = fmemopen(buf.first, buf.second - buf.first, "r");
    }

    try{

      MPI_Comm old_world = world; int old_me = comm->me;

      if(flag == PAIR){
        force->create_pair(style, 1);
        world = MPI_COMM_SELF; comm->me = 0;
        force->pair->read_restart(fp);
        world = old_world; comm->me = old_me;
      } else if(flag == NO_PAIR){
        force->create_pair("none", 0);
        force->pair_restart = style;
      } else if(flag == BOND){
        force->create_bond(style, 1);
        world = MPI_COMM_SELF; comm->me = 0;
        force->bond->read_restart(fp);
        world = old_world; comm->me = old_me;
      } else if(flag == ANGLE){
        force->create_angle(style, 1);
        world = MPI_COMM_SELF; comm->me = 0;
        force->angle->read_restart(fp);
        world = old_world; comm->me = old_me;
      } else if(flag == DIHEDRAL){
        force->create_dihedral(style, 1);
        world = MPI_COMM_SELF; comm->me = 0;
        force->dihedral->read_restart(fp);
        world = old_world; comm->me = old_me;
      } else if(flag == IMPROPER){
        force->create_improper(style, 1);
        world = MPI_COMM_SELF; comm->me = 0;
        force->improper->read_restart(fp);
        world = old_world; comm->me = old_me;
      } else {
        error->all(FLERR, "Invalid flag in checkpointed force fields");
      }

      delete[] style;
    } catch (const CommException& e){
      if(fp != nullptr) fclose(fp);
      delete[] style;
      throw;
    }

    if(flag != NO_PAIR){
      bigint offset = platform::ftell(fp);
      buf.first += offset;
      fclose(fp);
    }
  }
}

/* ---------------------------------------------------------------------- */

int FenixCheckpoint::recover_modify(Buffer& buf){
  auto& m = *modify;

  read(m.nfix_restart_global, buf, error);
  if(m.nfix_restart_global){
    m.id_restart_global = new char *[m.nfix_restart_global];
    m.style_restart_global = new char *[m.nfix_restart_global];
    m.state_restart_global = new char *[m.nfix_restart_global];

    m.used_restart_global = new int[m.nfix_restart_global];
  }

  int n;
  for(int i = 0; i < m.nfix_restart_global; i++){
    m.id_restart_global[i] = new char[n = read<int>(buf, error)];
    read(m.id_restart_global[i], n, buf, error);

    m.style_restart_global[i] = new char[n = read<int>(buf, error)];
    read(m.style_restart_global[i], n, buf, error);

    m.state_restart_global[i] = new char[n = read<int>(buf, error)];
    read(m.state_restart_global[i], n, buf, error);

    m.used_restart_global[i] = 0;
  }

  read(m.nfix_restart_peratom, buf, error);
  if(m.nfix_restart_peratom){
    m.id_restart_peratom = new char *[m.nfix_restart_peratom];
    m.style_restart_peratom = new char *[m.nfix_restart_peratom];
    m.index_restart_peratom = new int[m.nfix_restart_peratom];
    m.used_restart_peratom = new int[m.nfix_restart_peratom];
  }

  int maxsize = 0;
  for(int i = 0; i < m.nfix_restart_peratom; i++){
    m.id_restart_peratom[i] = new char[n = read<int>(buf, error)];
    read(m.id_restart_peratom[i], n, buf, error);

    m.style_restart_peratom[i] = new char[n = read<int>(buf, error)];
    read(m.style_restart_peratom[i], n, buf, error);

    maxsize += read<int>(buf, error);

    m.index_restart_peratom[i] = i;
    m.used_restart_peratom[i] = 0;
  }

  return maxsize;
}

/* ---------------------------------------------------------------------- */

void FenixCheckpoint::proc_grid(){
  auto& c = *comm;

  fwrite(c.other_procgrid, sizeof(int), 3, fp);
  fwrite(c.other_coregrid, sizeof(int), 3, fp);
  fwrite(c.procgrid, sizeof(int), 3, fp);
  fwrite(c.coregrid, sizeof(int), 3, fp);

  fwrite(c.myloc, sizeof(int), 3, fp);
  fwrite(c.procneigh, sizeof(int), 6, fp);

  for(int x = 0; x < c.procgrid[0]; x++)
    for(int y = 0; y < c.procgrid[1]; y++)
      fwrite(c.grid2proc[x][y], sizeof(int), c.procgrid[2], fp);

  fwrite(c.xsplit, sizeof(double), c.procgrid[0]+1, fp);
  fwrite(c.ysplit, sizeof(double), c.procgrid[1]+1, fp);
  fwrite(c.zsplit, sizeof(double), c.procgrid[2]+1, fp);
}

/* ---------------------------------------------------------------------- */

void FenixCheckpoint::proc_grid(Buffer& buf){
  auto& c = *comm;
  read(c.other_procgrid, 3, buf, error);
  read(c.other_coregrid, 3, buf, error);
  read(c.procgrid, 3, buf, error);
  read(c.coregrid, 3, buf, error);

  if(c.procgrid[0]*c.procgrid[1]*c.procgrid[2] != c.nprocs) error->all(
    FLERR, "Bad grid of processors during recovery"
  );
  if(domain->dimension == 2 && c.procgrid[2] != 1) error->all(
    FLERR, "Found z dimenion in procgrid for 2d simulation during recovery"
  );

  read(c.myloc, 3, buf, error);
  read(&c.procneigh[0][0], 6, buf, error);

  if(c.grid2proc) memory->destroy(c.grid2proc);
  memory->create(
    c.grid2proc, c.procgrid[0], c.procgrid[1], c.procgrid[2], "comm:grid2proc"
  );
  for(int x = 0; x < c.procgrid[0]; x++)
    for(int y = 0; y < c.procgrid[1]; y++)
      read(c.grid2proc[x][y], c.procgrid[2], buf, error);

  memory->destroy(c.xsplit);
  memory->destroy(c.ysplit);
  memory->destroy(c.zsplit);

  memory->create(c.xsplit, c.procgrid[0]+1, "comm:xsplit");
  memory->create(c.ysplit, c.procgrid[1]+1, "comm:ysplit");
  memory->create(c.zsplit, c.procgrid[2]+1, "comm:zsplit");

  read(c.xsplit, c.procgrid[0]+1, buf, error);
  read(c.ysplit, c.procgrid[1]+1, buf, error);
  read(c.zsplit, c.procgrid[2]+1, buf, error);

  if(domain->triclinic) domain->set_lamda_box();
}

/* ----------------------------------------------------------------------
  See ReadRestart::command for more detail on the steps being taken
---------------------------------------------------------------------- */

void FenixCheckpoint::restore_global(){
  // TODO: Any benefit to being able to? Less expensive restarts vs tearing everything down first?
  if(domain->box_exist) error->one(
    FLERR, "Cannot restore while simulation box is defined"
  );

  char* buf_ptr;
  int buf_len = safe_restore(global_member, buf_ptr);
  if(buf_len == 0){
    return;
  }
  Buffer buf = {buf_ptr, buf_ptr+buf_len};

  magic_string(buf);

  // TODO: We can probably ensure this once we get into localization changes
  if ((comm->me == 0) && (modify->get_fix_by_style("property/atom").size() > 0))
    error->warning(FLERR, "Fix property/atom command must be specified after read_restart "
                   "to restore its data.");

  header(buf);
  domain->box_exist = 1;
  atom->allocate_type_arrays();
  atom->deallocate_topology();

  int nlocal = read<int>(buf, error);
  atom->avec->grow(atom->avec->roundup(nlocal));

  domain->set_initial_box(0);
  domain->set_global_box();
  proc_grid(buf);
  domain->set_local_box();

  recover_group(buf);
  type_arrays(buf);
  force_fields(buf);

  int nextra = recover_modify(buf);
  atom->nextra_store = nextra;
  memory->create(atom->extra, atom->nmax, nextra, "atom:extra");

  magic_string(buf);
  memory->destroy(buf_ptr);
}

/* ---------------------------------------------------------------------- */

void FenixCheckpoint::restore_peratom(){
  char* buf_ptr;
  int buf_len = safe_restore(peratom_member, buf_ptr);

  if(buf_len % sizeof(double)) error->one(
    FLERR, "Invalid checkpoint data"
  );

  int n = buf_len / sizeof(double);
  double* buf = (double*) buf_ptr;

  int m = 0;
  while(m < n) m += atom->avec->unpack_restart(&buf[m]);

  delete[] buf_ptr;
}

}

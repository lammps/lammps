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

#include "FENIX/local_serializer.h"

#include <cstdio>

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

namespace LAMMPS_NS {
using Buffer = typename LocalSerializer::Buffer;

/* ---------------------------------------------------------------------- */

LocalSerializer::LocalSerializer(LAMMPS *lmp) : WriteRestart(lmp) {

}

/* ---------------------------------------------------------------------- */

Buffer LocalSerializer::serialize(){
  // See WriteRestart::write for more explanation of these steps
  
  for(auto& fix : modify->get_fix_list()) if(fix->restart_file) error->all(
    FLERR, "LocalSerializer does not yet support fixes that write their own"
    " restart files! Fix \"{}\" unsupported", fix->style
  );

  // Open a tmp file and hope it's stored in memory
  fp = std::tmpfile();
  if (fp == nullptr)
    error->one(FLERR, "Cannot open temporary file for restart: {}", utils::getsyserror());
  
  // Update protected variables inside WriteRestart
  MPI_Comm_rank(world, &me);
  MPI_Comm_size(world, &nprocs);
  natoms = atom->natoms;

  // Store copies of correct comm variables, in case interrupted by exception
  const MPI_Comm real_world = world;
  const int real_me = comm->me;
  
  try {
    if(neighbor->build_once) domain->reset_box();

    WriteRestart::magic_string();
    WriteRestart::header();
    fwrite(&atom->nlocal, sizeof(atom->nlocal), 1, fp);
    proc_grid();
    group->write_restart(fp);
    WriteRestart::type_arrays();
    WriteRestart::force_fields();

    // Get global data on all ranks and skip communicating local data (hacky)
    world = MPI_COMM_SELF;
    comm->me = 0;
    modify->write_restart(fp);
    comm->me = real_me;
    world = real_world;

    WriteRestart::magic_string();
  } catch(const std::exception& e){
    fclose(fp);
    fp = nullptr;
    world = real_world;
    comm->me = real_me;
    throw;
  }

  int globals_size = platform::ftell(fp);
  int locals_size = atom->avec->size_restart()*sizeof(double);
  int buf_size = globals_size + locals_size;

  // Allocate buffer
  char* buf;
  memory->create(buf, buf_size, "local_serializer:serialize_buf");

  // Write global data to buffer
  platform::fseek(fp, 0);
  fread(buf, sizeof(char), globals_size, fp);
  fclose(fp);
  fp = nullptr;
 
  // Write local data to buffer
  // It seems safe to ignore restart_pbc_any, since we are simply restoring
  // the prior state on recovery - not restoring to a different number of ranks
  // or a different set of configs etc. But if things break, fetch the code
  // for enforcing Periodic Boundary Conditions (PBC) from write_restart
  double* atom_data = (double*) (buf + globals_size);
  for(int n = 0, i = 0; i < atom->nlocal; i++)
    n += atom->avec->pack_restart(i, &atom_data[n]);

  return {buf, buf + buf_size};
}

/* ---------------------------------------------------------------------- */

void LocalSerializer::deserialize(Buffer buf){
  // TODO: Any benefit to being able to, to avoid full teardown?
  if(domain->box_exist)
    error->one(FLERR, "Cannot restore while simulation box is defined");
  // TODO: We can probably ensure this once we get into localization changes (?)
  if ((comm->me == 0) && (modify->get_fix_by_style("property/atom").size() > 0))
    error->warning(FLERR, "Fix property/atom command must be after deserialize "
                          "to restore its data.");

  // First global data, surrounded by magic strings to verify
  magic_string(buf);
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

  // Now local atom data
  double* atom_data = (double*) buf.first;

  int n = buf.second - buf.first;
  if(n % sizeof(double) != 0) error->one(FLERR, "Invalid checkpoint data");
  n = n / sizeof(double);

  int m = 0;
  while(m < n) m += atom->avec->unpack_restart(&atom_data[m]);
}

/* ---------------------------------------------------------------------- */

void LocalSerializer::proc_grid(){
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

void LocalSerializer::proc_grid(Buffer& buf){
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

/* ---------------------------------------------------------------------- */

void LocalSerializer::header(Buffer& buf){
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
    } else if (flag == ATOM_MAXEXCHANGE) {
      if(atom->avec) read(atom->avec->maxexchange, buf, error);

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

    } else error->all(FLERR,"Invalid flag {} in header section of restart file", flag);
  }
}

/* ---------------------------------------------------------------------- */

void lmap_type(
  std::vector<std::string>& labels, int ntypes,
  std::unordered_map<std::string, int>& map,
  Buffer& buf, Error* error
) {
  for(int i = 0; i < ntypes; i++){
    int n = LocalSerializer::read<int>(buf, error);
    if(n == 0){
      labels[i] = new char[0];
      continue;
    }
    LocalSerializer::read(labels[i], buf, error);
    if(!labels[i].empty()) map[labels[i]] = i+1;
  }
}

/* ---------------------------------------------------------------------- */

void LocalSerializer::type_arrays(Buffer& buf){
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

/* ---------------------------------------------------------------------- */

void LocalSerializer::force_fields(Buffer& buf){
  int flag;
  while( (flag = read<int>(buf, error)) >= 0) {
    read<int>(buf, error);
    auto style = read<char*>(buf, error);

    try {
      switch(flag) {
        case PAIR:     force->create_pair(    style, 1); break;
        case NO_PAIR:  force->create_pair(   "none", 0); break;
        case BOND:     force->create_bond(    style, 1); break;
        case ANGLE:    force->create_angle(   style, 1); break;
        case DIHEDRAL: force->create_dihedral(style, 1); break;
        case IMPROPER: force->create_improper(style, 1); break;
        default: error->all(FLERR, "Invalid flag in checkpointed force fields");
      }
    } catch (const std::exception& e){
      delete[] style;
      throw;
    }

    if(flag == NO_PAIR){
      force->pair_restart = style;
      continue;
    }
    delete[] style;

    MPI_Comm real_world = world; int real_me = comm->me;
    FILE* fp = nullptr;
    try {
      // POSIX, so should be generally portable to non-windows
      // but there's not a good windows alternative tmk
      fp = fmemopen(buf.first, buf.second - buf.first, "r");
      
      world = MPI_COMM_SELF; comm->me = 0;
      switch(flag){
        case PAIR:     force->pair->    read_restart(fp); break;
        case BOND:     force->bond->    read_restart(fp); break;
        case ANGLE:    force->angle->   read_restart(fp); break;
        case DIHEDRAL: force->dihedral->read_restart(fp); break;
        case IMPROPER: force->improper->read_restart(fp); break;
        default: error->all(FLERR, "Invalid flag in checkpointed force fields");
      }
      world = real_world; comm->me = real_me;

      buf.first += platform::ftell(fp);
      fclose(fp);
    } catch (const std::exception& e){
      world = real_world; comm->me = real_me;
      if(fp != nullptr) fclose(fp);
      throw;
    }
  }
}

/* ---------------------------------------------------------------------- */

int LocalSerializer::recover_modify(Buffer& buf){
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

/* ----------------------------------------------------------------------
  Rewritten group->read_restart to not require a FILE*
---------------------------------------------------------------------- */

void LocalSerializer::recover_group(Buffer& buf){
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

void LocalSerializer::magic_string(Buffer& buf){
  if(read<std::string>(buf, error) != MAGIC_STRING)
    error->one(FLERR, "Checkpointed data is invalid");
}

}

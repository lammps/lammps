// clang-format off
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

#include "write_restart.h"

#include "angle.h"
#include "atom.h"
#include "atom_vec.h"
#include "bond.h"
#include "comm.h"
#include "dihedral.h"
#include "domain.h"
#include "error.h"
#include "fix.h"
#include "force.h"
#include "group.h"
#include "improper.h"
#include "label_map.h"
#include "memory.h"
#include "modify.h"
#include "neighbor.h"
#include "output.h"
#include "pair.h"
#include "thermo.h"
#include "update.h"
#include "file_writer.h"
#include "file_writer_wrapper.h"
#include "file_writer_root_wrapper.h"

#include <cstring>

#include "lmprestart.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

WriteRestart::WriteRestart(LAMMPS *lmp) : Command(lmp)
{
  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);
  multiproc = 0;
  noinit = 0;
}

/* ----------------------------------------------------------------------
   called as write_restart command in input script
------------------------------------------------------------------------- */

void WriteRestart::command(int narg, char **arg)
{
  if (domain->box_exist == 0)
    error->all(FLERR,"Write_restart command before simulation box is defined" + utils::errorurl(33));
  if (narg < 1) utils::missing_cmd_args(FLERR, "write_restart", error);

  // if filename contains a "*", replace with current timestep

  std::string file = arg[0];
  std::size_t found = file.find('*');
  if (found != std::string::npos)
    file.replace(found,1,fmt::format("{}",update->ntimestep));

  // check for multiproc output and an MPI-IO filename

  if (strchr(arg[0],'%')) multiproc = nprocs;
  else multiproc = 0;
  if (utils::strmatch(arg[0],"\\.mpiio$"))
    error->all(FLERR,"MPI-IO files are no longer supported by LAMMPS");

  // setup output style and process optional args
  // also called by Output class for periodic restart files

  multiproc_options(multiproc,narg-1,&arg[1]);

  // init entire system since comm->exchange is done
  // comm::init needs neighbor::init needs pair::init needs kspace::init, etc

  if (noinit == 0) {
    if (comm->me == 0) utils::logmesg(lmp,"System init for write_restart ...\n");
    lmp->init();

    // move atoms to new processors before writing file
    // enforce PBC in case atoms are outside box
    // call borders() to rebuild atom map since exchange() destroys map
    // NOTE: removed call to setup_pre_exchange
    //   used to be needed by fixShearHistory for granular
    //   to move history info from neigh list to atoms between runs
    //   but now that is done via FIx::post_run()
    //   don't think any other fix needs this or should do it
    //   e.g. fix evaporate should not delete more atoms

    // modify->setup_pre_exchange();
    if (domain->triclinic) domain->x2lamda(atom->nlocal);
    domain->pbc();
    domain->reset_box();
    comm->setup();
    comm->exchange();
    comm->borders();
    if (domain->triclinic) domain->lamda2x(atom->nlocal+atom->nghost);
  }

  // write single restart file

  write(file);
}

/* ---------------------------------------------------------------------- */

void WriteRestart::multiproc_options(int multiproc_caller, int narg, char **arg)
{
  multiproc = multiproc_caller;

  // defaults for multiproc file writing

  nclusterprocs = nprocs;
  filewriter = 0;
  if (me == 0) filewriter = 1;
  fileproc = 0;

  if (multiproc) {
    nclusterprocs = 1;
    filewriter = 1;
    fileproc = me;
    icluster = me;
  }

  // optional args

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"fileper") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "write_restart", error);
      if (!multiproc)
        error->all(FLERR,"Cannot use write_restart fileper without % in restart file name");
      int nper = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
      if (nper <= 0) error->all(FLERR,"Invalue write_restart fileper value {}:", nper);

      multiproc = nprocs/nper;
      if (nprocs % nper) multiproc++;
      fileproc = me/nper * nper;
      int fileprocnext = MIN(fileproc+nper,nprocs);
      nclusterprocs = fileprocnext - fileproc;
      if (me == fileproc) filewriter = 1;
      else filewriter = 0;
      icluster = fileproc/nper;
      iarg += 2;

    } else if (strcmp(arg[iarg],"nfile") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "write_restart nfile", error);
      if (!multiproc)
        error->all(FLERR,"Cannot use write_restart nfile without % in restart file name");
      int nfile = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
      if (nfile <= 0) error->all(FLERR,"Invalid write_restart nfile value {}", nfile);
      nfile = MIN(nfile,nprocs);

      multiproc = nfile;
      icluster = static_cast<int> ((bigint) me * nfile/nprocs);
      fileproc = static_cast<int> ((bigint) icluster * nprocs/nfile);
      int fcluster = static_cast<int> ((bigint) fileproc * nfile/nprocs);
      if (fcluster < icluster) fileproc++;
      int fileprocnext =
        static_cast<int> ((bigint) (icluster+1) * nprocs/nfile);
      fcluster = static_cast<int> ((bigint) fileprocnext * nfile/nprocs);
      if (fcluster < icluster+1) fileprocnext++;
      nclusterprocs = fileprocnext - fileproc;
      if (me == fileproc) filewriter = 1;
      else filewriter = 0;
      iarg += 2;

    } else if (strcmp(arg[iarg],"noinit") == 0) {
      noinit = 1;
      iarg++;
    } else error->all(FLERR,"Unknown write_restart keyword: {}", arg[iarg]);
  }
}

/* ----------------------------------------------------------------------
   called from command() and directly from output within run/minimize loop
   file = final file name to write, except may contain a "%"
------------------------------------------------------------------------- */

void WriteRestart::write(const std::string &file)
{
  mode = WRITE_RESTART;

  // special case where reneighboring is not done in integrator
  //   on timestep restart file is written (due to build_once being set)
  // if box is changing, must be reset, else restart file will have
  //   wrong box size and atoms will be lost when restart file is read
  // other calls to pbc and domain and comm are not made,
  //   b/c they only make sense if reneighboring is actually performed

  if (neighbor->build_once) domain->reset_box();

  // natoms = sum of nlocal = value to write into restart file
  // if unequal and thermo lostflag is "error", don't write restart file

  bigint nblocal = atom->nlocal;
  MPI_Allreduce(&nblocal,&natoms,1,MPI_LMP_BIGINT,MPI_SUM,world);
  if (natoms != atom->natoms && output->thermo->lostflag == Thermo::ERROR)
    error->all(FLERR,"Atom count is inconsistent: {} vs {}, cannot write restart file",
               natoms, atom->natoms);

  // open single restart file or base file for multiproc case

  if (me == 0) {
    std::string base = file;
    if (multiproc) base.replace(base.find('%'),1,"base");

    fp = fopen(base.c_str(),"wb");
    if (fp == nullptr)
      error->one(FLERR, "Cannot open restart file {}: {}", base, utils::getsyserror());
  }

  // proc 0 writes magic string, global info
  {
    FileWriterRootWrapper fw_rooted(me, fp);
    global(fw_rooted);
  }

  // communication buffer for my atom info
  // max_size = largest buffer needed by any proc
  // NOTE: are assuming size_restart() returns 32-bit int
  //   for a huge one-proc problem, nlocal could be 32-bit
  //   but nlocal * doubles-peratom could overflow

  int send_size = atom->avec->size_restart();
  int max_sizes[2], sizes[2] = {send_size, modify->maxsize_restart()};
  MPI_Allreduce(sizes,max_sizes,2,MPI_INT,MPI_MAX,world);
  int max_size = max_sizes[0], max_extra = max_sizes[1];

  double *buf;
  memory->create(buf,max_size,"write_restart:buf");
  memset(buf,0,max_size*sizeof(double));

  // Proc 0 writes file layout info which may include per-proc sizes
  {
    FileWriterRootWrapper fw_rooted(me, fp);
    file_layout(fw_rooted, send_size, max_extra);
  }

  // header info is complete
  // if multiproc output:
  //   close header file, open multiname file on each writing proc,
  //   write PROCSPERFILE into new file

  int io_error = 0;
  if (multiproc) {
    if (me == 0 && fp) {
      FileWriterWrapper fw(fp);
      magic_string(fw);
      if (ferror(fp)) io_error = 1;
      fp = nullptr;             // implicitly closes file
    }

    std::string multiname = file;
    multiname.replace(multiname.find('%'),1,fmt::format("{}",icluster));

    if (filewriter) {
      fp = fopen(multiname.c_str(),"wb");
      if (fp == nullptr)
        error->one(FLERR, "Cannot open restart file {}: {}", multiname, utils::getsyserror());
    }
  }

  // check for I/O error status

  int io_all = 0;
  MPI_Allreduce(&io_error,&io_all,1,MPI_INT,MPI_MAX,world);
  if (io_all) error->all(FLERR,"I/O error while writing restart");

  // pack my atom data into buf

  AtomVec *avec = atom->avec;
  int n = 0;
  for (int i = 0; i < atom->nlocal; i++) n += avec->pack_restart(i,&buf[n]);

  // if any fix requires it, remap each atom's coords via PBC
  // is because fix changes atom coords (excepting an integrate fix)
  // just remap in buffer, not actual atoms

  if (modify->restart_pbc_any) {
    int triclinic = domain->triclinic;
    double *lo,*hi,*period;

    if (triclinic == 0) {
      lo = domain->boxlo;
      hi = domain->boxhi;
      period = domain->prd;
    } else {
      lo = domain->boxlo_lamda;
      hi = domain->boxhi_lamda;
      period = domain->prd_lamda;
    }

    int xperiodic = domain->xperiodic;
    int yperiodic = domain->yperiodic;
    int zperiodic = domain->zperiodic;

    double *x;
    int m = 0;
    for (int i = 0; i < atom->nlocal; i++) {
      x = &buf[m+1];
      if (triclinic) domain->x2lamda(x,x);

      if (xperiodic) {
        if (x[0] < lo[0]) x[0] += period[0];
        if (x[0] >= hi[0]) x[0] -= period[0];
        x[0] = MAX(x[0],lo[0]);
      }
      if (yperiodic) {
        if (x[1] < lo[1]) x[1] += period[1];
        if (x[1] >= hi[1]) x[1] -= period[1];
        x[1] = MAX(x[1],lo[1]);
      }
      if (zperiodic) {
        if (x[2] < lo[2]) x[2] += period[2];
        if (x[2] >= hi[2]) x[2] -= period[2];
        x[2] = MAX(x[2],lo[2]);
      }

      if (triclinic) domain->lamda2x(x,x);
      m += static_cast<int> (buf[m]);
    }
  }

  // output of one or more native files
  // filewriter = 1 = this proc writes to file
  // ping each proc in my cluster, receive its data, write data to file
  // else wait for ping from fileproc, send my data to fileproc

  int tmp = 0;
  int recv_size;

  if (filewriter) {
    FileWriterWrapper fw(fp);
    if (multiproc) write_val<int>(fw,PROCSPERFILE,nclusterprocs);

    MPI_Status status;
    MPI_Request request;
    for (int iproc = 0; iproc < nclusterprocs; iproc++) {
      if (iproc) {
        MPI_Irecv(buf,max_size,MPI_DOUBLE,me+iproc,0,world,&request);
        MPI_Send(&tmp,0,MPI_INT,me+iproc,0,world);
        MPI_Wait(&request,&status);
        MPI_Get_count(&status,MPI_DOUBLE,&recv_size);
      } else recv_size = send_size;

      write_vec<double>(fw,PERPROC,recv_size,buf);
    }
    magic_string(fw);
    if (ferror(fp)) io_error = 1;
    fp = nullptr;               // implicitly closes file
  } else {
    MPI_Recv(&tmp,0,MPI_INT,fileproc,0,world,MPI_STATUS_IGNORE);
    MPI_Rsend(buf,send_size,MPI_DOUBLE,fileproc,0,world);
  }

  // check again for I/O error status

  io_all = 0;
  MPI_Allreduce(&io_error,&io_all,1,MPI_INT,MPI_MAX,world);
  if (io_all) error->all(FLERR,"I/O error while writing restart");

  // clean up

  memory->destroy(buf);

  // invoke any fixes that write their own restart file

  for (const auto &fix : modify->get_fix_list())
    if (fix->restart_file)
      fix->write_restart_file(file.c_str());
}

void WriteRestart::write_restart_global(FileWriter& fw) {
  mode = WRITE_RESTART_GLOBAL;
  global(fw);
  magic_string(fw);
}

void WriteRestart::global(FileWriter& fw) {
  magic_string(fw);
  fw.writev<int>(ENDIAN);
  fw.writev<int>(FORMAT_REVISION);
  header(fw);
  write_restartable(fw, group);
  type_arrays(fw);
  force_fields(fw);
  write_restartable(fw, modify);
}

/* ----------------------------------------------------------------------
   Writes out problem description
------------------------------------------------------------------------- */

void WriteRestart::header(FileWriter& fw)
{
  write_val<std::string>(fw,VERSION,lmp->version);
  write_val<int>(fw,SMALLINT,sizeof(smallint));
  write_val<int>(fw,IMAGEINT,sizeof(imageint));
  write_val<int>(fw,TAGINT,sizeof(tagint));
  write_val<int>(fw,BIGINT,sizeof(bigint));
  write_val<std::string>(fw,UNITS,update->unit_style);
  write_val<bigint>(fw,NTIMESTEP,update->ntimestep);
  write_val<int>(fw,DIMENSION,domain->dimension);
  write_val<int>(fw,NPROCS,nprocs);
  write_vec<int>(fw,PROCGRID,3,comm->procgrid);
  write_val<int>(fw,NEWTON_PAIR,force->newton_pair);
  write_val<int>(fw,NEWTON_BOND,force->newton_bond);
  write_val<int>(fw,XPERIODIC,domain->xperiodic);
  write_val<int>(fw,YPERIODIC,domain->yperiodic);
  write_val<int>(fw,ZPERIODIC,domain->zperiodic);
  write_vec<int>(fw,BOUNDARY,6,&domain->boundary[0][0]);

  // added field for shrink-wrap boundaries with minimum - 2 Jul 2015

  double minbound[6];
  minbound[0] = domain->minxlo; minbound[1] = domain->minxhi;
  minbound[2] = domain->minylo; minbound[3] = domain->minyhi;
  minbound[4] = domain->minzlo; minbound[5] = domain->minzhi;
  write_vec<double>(fw,BOUNDMIN,6,minbound);

  // write atom_style and its args

  write_val<std::string>(fw,ATOM_STYLE,utils::strip_style_suffix(atom->atom_style,lmp));
  fw.writev<int>(atom->avec->nargcopy);
  for (int i = 0; i < atom->avec->nargcopy; i++) {
    int n = strlen(atom->avec->argcopy[i]) + 1;
    fw.writev(n);
    fw.writev(atom->avec->argcopy[i], n);
  }

  write_val<bigint>(fw,NATOMS,natoms);
  write_val<int>(fw,NTYPES,atom->ntypes);
  write_val<bigint>(fw,NBONDS,atom->nbonds);
  write_val<int>(fw,NBONDTYPES,atom->nbondtypes);
  write_val<int>(fw,BOND_PER_ATOM,atom->bond_per_atom);
  write_val<bigint>(fw,NANGLES,atom->nangles);
  write_val<int>(fw,NANGLETYPES,atom->nangletypes);
  write_val<int>(fw,ANGLE_PER_ATOM,atom->angle_per_atom);
  write_val<bigint>(fw,NDIHEDRALS,atom->ndihedrals);
  write_val<int>(fw,NDIHEDRALTYPES,atom->ndihedraltypes);
  write_val<int>(fw,DIHEDRAL_PER_ATOM,atom->dihedral_per_atom);
  write_val<bigint>(fw,NIMPROPERS,atom->nimpropers);
  write_val<int>(fw,NIMPROPERTYPES,atom->nimpropertypes);
  write_val<int>(fw,IMPROPER_PER_ATOM,atom->improper_per_atom);

  write_val<int>(fw,TRICLINIC,domain->triclinic);
  write_vec<double>(fw,BOXLO,3,domain->boxlo);
  write_vec<double>(fw,BOXHI,3,domain->boxhi);
  write_val<double>(fw,XY,domain->xy);
  write_val<double>(fw,XZ,domain->xz);
  write_val<double>(fw,YZ,domain->yz);

  write_val<int>(fw,TRICLINIC_GENERAL,domain->triclinic_general);
  if (domain->triclinic_general)
    write_vec<double>(fw,ROTATE_G2R,9,&domain->rotate_g2r[0][0]);

  write_vec<double>(fw,SPECIAL_LJ,3,&force->special_lj[1]);
  write_vec<double>(fw,SPECIAL_COUL,3,&force->special_coul[1]);

  write_val<double>(fw,TIMESTEP,update->dt);

  write_val<int>(fw,ATOM_ID,atom->tag_enable);
  write_val<int>(fw,ATOM_MAP_STYLE,atom->map_style);
  write_val<int>(fw,ATOM_MAP_USER,atom->map_user);
  write_val<int>(fw,ATOM_SORTFREQ,atom->sortfreq);
  write_val<double>(fw,ATOM_SORTBIN,atom->userbinsize);

  write_val<int>(fw,COMM_MODE,comm->mode);
  write_val<double>(fw,COMM_CUTOFF,comm->cutghostuser);
  write_val<int>(fw,COMM_VEL,comm->ghost_velocity);

  write_val<int>(fw,EXTRA_BOND_PER_ATOM,atom->extra_bond_per_atom);
  write_val<int>(fw,EXTRA_ANGLE_PER_ATOM,atom->extra_angle_per_atom);
  write_val<int>(fw,EXTRA_DIHEDRAL_PER_ATOM,atom->extra_dihedral_per_atom);
  write_val<int>(fw,EXTRA_IMPROPER_PER_ATOM,atom->extra_improper_per_atom);
  write_val<int>(fw,ATOM_MAXSPECIAL,atom->maxspecial);

  // write out AtomVec::maxexchange (extra storage for communicating
  // per-atom bond, angle, dihedral, and improper data). added 25 Oct 2025

  write_val<int>(fw,ATOM_MAXEXCHANGE,atom->avec->maxexchange);

  write_val<bigint>(fw,NELLIPSOIDS,atom->nellipsoids);
  write_val<bigint>(fw,NLINES,atom->nlines);
  write_val<bigint>(fw,NTRIS,atom->ntris);
  write_val<bigint>(fw,NBODIES,atom->nbodies);

  // write out current simulation time. added 3 May 2022

  write_val<bigint>(fw,ATIMESTEP,update->atimestep);
  write_val<double>(fw,ATIME,update->atime);

  // -1 flag signals end of header

  fw.writev<int>(-1);
}

/* ----------------------------------------------------------------------
   Writes out any type-based arrays that are defined
------------------------------------------------------------------------- */

void WriteRestart::type_arrays(FileWriter& fw)
{
  if (fw && atom->mass) write_vec<double>(fw,MASS,atom->ntypes,&atom->mass[1]);
  if (atom->labelmapflag) {
    write_val<int>(fw,LABELMAP,atom->labelmapflag);
    write_restartable(fw,atom->lmap);
  }

  // -1 flag signals end of type arrays

  fw.writev<int>(-1);
}

/* ----------------------------------------------------------------------
   Writes out force field styles and data that are defined
------------------------------------------------------------------------- */

void WriteRestart::force_fields(FileWriter& fw)
{
  if (force->pair) {
    int tag = is_restartable(force->pair) ? PAIR : NO_PAIR;
    write_val<std::string>(fw,tag,utils::strip_style_suffix(force->pair_style,lmp));
    write_restartable(fw,force->pair);
  }
  if (atom->avec->bonds_allow && force->bond) {
    write_val<std::string>(fw,BOND,utils::strip_style_suffix(force->bond_style,lmp));
    write_restartable(fw,force->bond);
  }
  if (atom->avec->angles_allow && force->angle) {
    write_val<std::string>(fw,ANGLE,utils::strip_style_suffix(force->angle_style,lmp));
    write_restartable(fw,force->angle);
  }
  if (atom->avec->dihedrals_allow && force->dihedral) {
    write_val<std::string>(fw,DIHEDRAL,utils::strip_style_suffix(force->dihedral_style,lmp));
    write_restartable(fw,force->dihedral);
  }
  if (atom->avec->impropers_allow && force->improper) {
    write_val<std::string>(fw,IMPROPER,utils::strip_style_suffix(force->improper_style,lmp));
    write_restartable(fw,force->improper);
  }

  // -1 flag signals end of force field info

  fw.writev<int>(-1);
}

/* ----------------------------------------------------------------------
   Writes out file layout info
   All procs call this method, only procs with open FileWriters write info
------------------------------------------------------------------------- */

void WriteRestart::file_layout(FileWriter& fw, int /*send_size*/, int max_extra)
{
  write_val<int>(fw,MULTIPROC,multiproc);
  write_val<int>(fw,EXTRA_FIX_PER_ATOM,max_extra);

  // -1 flag signals end of file layout info
  fw.writev<int>(-1);
}

/* ---------------------------------------------------------------------- */

void WriteRestart::write_restartable(FileWriter& fw, Restartable *r) {
  if (!is_restartable(r)) return;
  switch(mode) {
  case WRITE_RESTART:
    // "r->restartable" implies this should be called on all callers to ensure
    // all locals are aggregated as needed. Else, just call on callers with an
    // open SafeFilePtr
    if (r->restartable() || fp) r->write_restart(fp);
    break;
  case WRITE_RESTART_GLOBAL:
    r->write_restart_global(fw);
    break;
  case WRITE_RESTART_LOCAL:
    r->write_restart_local(fw);
    break;
  default: error->all(FLERR, "Unknown WriteRestart mode {}", (int)mode);
  }
}

/* ---------------------------------------------------------------------- */

bool WriteRestart::is_restartable(const Restartable *r) {
  switch (mode) {
  case WRITE_RESTART: return r->restartable_file();
  case WRITE_RESTART_GLOBAL: return r->restartable_global;
  case WRITE_RESTART_LOCAL: return r->restartable_local;
  default: error->all(FLERR, "Unknown WriteRestart mode {}", (int)mode);
  }
}

/* ---------------------------------------------------------------------- */

void WriteRestart::magic_string(FileWriter& fw)
{
  fw.writev(MAGIC_STRING, strlen(MAGIC_STRING)+1);
}

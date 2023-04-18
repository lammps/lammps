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

#include "write_data.h"

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
#include "improper.h"
#include "label_map.h"
#include "memory.h"
#include "modify.h"
#include "output.h"
#include "pair.h"
#include "thermo.h"
#include "update.h"

#include <cstring>

using namespace LAMMPS_NS;

enum{II,IJ};
enum{ELLIPSOID,LINE,TRIANGLE,BODY};   // also in AtomVecHybrid

/* ---------------------------------------------------------------------- */

WriteData::WriteData(LAMMPS *lmp) : Command(lmp)
{
  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);
}

/* ----------------------------------------------------------------------
   called as write_data command in input script
------------------------------------------------------------------------- */

void WriteData::command(int narg, char **arg)
{
  if (domain->box_exist == 0)
    error->all(FLERR,"Write_data command before simulation box is defined");

  if (narg < 1) utils::missing_cmd_args(FLERR, "write_data", error);

  // if filename contains a "*", replace with current timestep

  std::string file = arg[0];
  std::size_t found = file.find('*');
  if (found != std::string::npos)
    file.replace(found,1,fmt::format("{}",update->ntimestep));

  // read optional args
  // noinit is a hidden arg, only used by -r command-line switch

  pairflag = II;
  coeffflag = 1;
  fixflag = 1;
  lmapflag = 1;
  // store current (default) setting since we may change it.
  int types_style = atom->types_style;
  int noinit = 0;

  int iarg = 1;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"pair") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "write_data pair", error);
      if (strcmp(arg[iarg+1],"ii") == 0) pairflag = II;
      else if (strcmp(arg[iarg+1],"ij") == 0) pairflag = IJ;
      else error->all(FLERR,"Unknown write_data pair option: {}", arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"noinit") == 0) {
      noinit = 1;
      iarg++;
    } else if (strcmp(arg[iarg],"nocoeff") == 0) {
      coeffflag = 0;
      iarg++;
    } else if (strcmp(arg[iarg],"nofix") == 0) {
      fixflag = 0;
      iarg++;
    } else if (strcmp(arg[iarg],"nolabelmap") == 0) {
      lmapflag = 0;
      iarg++;
    } else if (strcmp(arg[iarg],"types") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "write_data types", error);
      if (strcmp(arg[iarg+1],"numeric") == 0) atom->types_style = Atom::NUMERIC;
      else if (strcmp(arg[iarg+1],"labels") == 0) atom->types_style = Atom::LABELS;
      else error->all(FLERR,"Unknown write_data types option: {}", arg[iarg+1]);
      iarg += 2;
    } else error->all(FLERR,"Unknown write_data keyword: {}", arg[iarg]);
  }

  // init entire system since comm->exchange is done
  // comm::init needs neighbor::init needs pair::init needs kspace::init, etc
  // exception is when called by -r command-line switch
  //   then write_data immediately follows reading of restart file
  //   assume that read_restart initialized necessary values
  //   if don't make exception:
  //     pair->init() can fail due to various unset values:
  //     e.g. pair hybrid coeffs, dpd ghost-atom velocity setting

  if (noinit == 0) {
    if (comm->me == 0) utils::logmesg(lmp,"System init for write_data ...\n");
    lmp->init();

    // move atoms to new processors before writing file
    // do setup_pre_exchange to force update of per-atom info if needed
    // enforce PBC in case atoms are outside box
    // call borders() to rebuild atom map since exchange() destroys map

    modify->setup_pre_exchange();
    if (domain->triclinic) domain->x2lamda(atom->nlocal);
    domain->pbc();
    domain->reset_box();
    comm->setup();
    comm->exchange();
    comm->borders();
    if (domain->triclinic) domain->lamda2x(atom->nlocal+atom->nghost);
  }

  write(file);
  // restore saved setting
  atom->types_style = types_style;
}

/* ----------------------------------------------------------------------
   called from command()
   might later let it be directly called within run/minimize loop
------------------------------------------------------------------------- */

void WriteData::write(const std::string &file)
{
  // special case where reneighboring is not done in integrator
  //   on timestep data file is written (due to build_once being set)
  // if box is changing, must be reset, else data file will have
  //   wrong box size and atoms will be lost when data file is read
  // other calls to pbc and domain and comm are not made,
  //   b/c they only make sense if reneighboring is actually performed

  //if (neighbor->build_once) domain->reset_box();

  // natoms = sum of nlocal = value to write into data file
  // if unequal and thermo lostflag is "error", don't write data file

  bigint nblocal = atom->nlocal;
  bigint natoms;
  MPI_Allreduce(&nblocal,&natoms,1,MPI_LMP_BIGINT,MPI_SUM,world);
  if (natoms != atom->natoms && output->thermo->lostflag == Thermo::ERROR)
    error->all(FLERR,"Atom count is inconsistent, cannot write data file");

  // sum up bond,angle,dihedral,improper counts
  // may be different than atom->nbonds,nangles, etc. if broken/turned-off

  if (atom->molecular == Atom::MOLECULAR && (atom->nbonds || atom->nbondtypes)) {
    nbonds_local = atom->avec->pack_bond(nullptr);
    MPI_Allreduce(&nbonds_local,&nbonds,1,MPI_LMP_BIGINT,MPI_SUM,world);
  }
  if (atom->molecular == Atom::MOLECULAR && (atom->nangles || atom->nangletypes)) {
    nangles_local = atom->avec->pack_angle(nullptr);
    MPI_Allreduce(&nangles_local,&nangles,1,MPI_LMP_BIGINT,MPI_SUM,world);
  }

  if (atom->molecular == Atom::MOLECULAR && (atom->ndihedrals || atom->ndihedraltypes)) {
    ndihedrals_local = atom->avec->pack_dihedral(nullptr);
    MPI_Allreduce(&ndihedrals_local,&ndihedrals,1,MPI_LMP_BIGINT,MPI_SUM,world);
  }

  if (atom->molecular == Atom::MOLECULAR && (atom->nimpropers || atom->nimpropertypes)) {
    nimpropers_local = atom->avec->pack_improper(nullptr);
    MPI_Allreduce(&nimpropers_local,&nimpropers,1,MPI_LMP_BIGINT,MPI_SUM,world);
  }

  // open data file

  if (me == 0) {
    fp = fopen(file.c_str(),"w");
    if (fp == nullptr)
      error->one(FLERR,"Cannot open data file {}: {}",
                                   file, utils::getsyserror());
  }

  // proc 0 writes header, ntype-length arrays, force fields
  // label map must come before coeffs

  if (me == 0) {
    header();
    if (lmapflag && atom->labelmapflag) atom->lmap->write_data(fp);
    type_arrays();
    if (coeffflag) force_fields();
  }

  // per atom info in Atoms and Velocities sections

  if (natoms) atoms();
  if (natoms) velocities();

  // molecular topology info if defined
  // do not write molecular topology for atom_style template

  if (atom->molecular == Atom::MOLECULAR) {
    if (atom->nbonds && nbonds) bonds();
    if (atom->nangles && nangles) angles();
    if (atom->ndihedrals) dihedrals();
    if (atom->nimpropers) impropers();
  }

  // bonus info if defined

  if (natoms && atom->ellipsoid_flag) bonus(ELLIPSOID);
  if (natoms && atom->line_flag) bonus(LINE);
  if (natoms && atom->tri_flag) bonus(TRIANGLE);
  if (natoms && atom->body_flag) bonus(BODY);

  // extra sections managed by fixes

  if (fixflag)
    for (auto &ifix : modify->get_fix_list())
      if (ifix->wd_section)
        for (int m = 0; m < ifix->wd_section; m++) fix(ifix,m);

  // close data file

  if (me == 0) fclose(fp);
}

/* ----------------------------------------------------------------------
   proc 0 writes out data file header
------------------------------------------------------------------------- */

void WriteData::header()
{
  fmt::print(fp,"LAMMPS data file via write_data, version {}, timestep = {}, units = {}\n\n",
             lmp->version, update->ntimestep, update->unit_style);

  fmt::print(fp,"{} atoms\n{} atom types\n",atom->natoms,atom->ntypes);

  // only write out number of types for atom style template

  if (atom->molecular == Atom::MOLECULAR) {
    if (atom->nbonds || atom->nbondtypes)
      fmt::print(fp,"{} bonds\n{} bond types\n",
                 nbonds,atom->nbondtypes);
    if (atom->nangles || atom->nangletypes)
      fmt::print(fp,"{} angles\n{} angle types\n",
                 nangles,atom->nangletypes);
    if (atom->ndihedrals || atom->ndihedraltypes)
      fmt::print(fp,"{} dihedrals\n{} dihedral types\n",
                 ndihedrals,atom->ndihedraltypes);
    if (atom->nimpropers || atom->nimpropertypes)
      fmt::print(fp,"{} impropers\n{} improper types\n",
                 nimpropers,atom->nimpropertypes);
  }

  if (atom->molecular == Atom::TEMPLATE) {
    if (atom->nbondtypes) fmt::print(fp,"{} bond types\n",atom->nbondtypes);
    if (atom->nangletypes) fmt::print(fp,"{} angle types\n",atom->nangletypes);
    if (atom->ndihedraltypes) fmt::print(fp,"{} dihedral types\n",atom->ndihedraltypes);
    if (atom->nimpropertypes) fmt::print(fp,"{} improper types\n",atom->nimpropertypes);
  }

  // bonus info

  if (atom->ellipsoid_flag) fmt::print(fp,"{} ellipsoids\n",atom->nellipsoids);
  if (atom->line_flag) fmt::print(fp,"{} lines\n",atom->nlines);
  if (atom->tri_flag) fmt::print(fp,"{} triangles\n",atom->ntris);
  if (atom->body_flag) fmt::print(fp,"{} bodies\n",atom->nbodies);

  // fix info

  if (fixflag)
    for (auto &ifix : modify->get_fix_list())
      if (ifix->wd_header)
        for (int m = 0; m < ifix->wd_header; m++)
          ifix->write_data_header(fp,m);

  // box info

  auto box = fmt::format("\n{} {} xlo xhi\n{} {} ylo yhi\n{} {} zlo zhi\n",
                         domain->boxlo[0],domain->boxhi[0],
                         domain->boxlo[1],domain->boxhi[1],
                         domain->boxlo[2],domain->boxhi[2]);
  if (domain->triclinic)
    box += fmt::format("{} {} {} xy xz yz\n",domain->xy,domain->xz,domain->yz);
  fputs(box.c_str(),fp);
}

/* ----------------------------------------------------------------------
   proc 0 writes out any type-based arrays that are defined
------------------------------------------------------------------------- */

void WriteData::type_arrays()
{
  if (atom->mass) {
    double *mass = atom->mass;
    fputs("\nMasses\n\n",fp);
    for (int i = 1; i <= atom->ntypes; i++)
      fmt::print(fp,"{} {:.16g}\n",i,mass[i]);
  }
}

/* ----------------------------------------------------------------------
   proc 0 writes out force field info
------------------------------------------------------------------------- */

void WriteData::force_fields()
{
  if (force->pair && force->pair->writedata) {
    if (pairflag == II) {
      if ((comm->me == 0) && (force->pair->mixed_flag == 0))
        error->warning(FLERR,"Not all mixed pair coeffs generated from mixing. "
                       "Use write_data with 'pair ij' option to store all pair coeffs.");
      fmt::print(fp,"\nPair Coeffs # {}\n\n", force->pair_style);
      force->pair->write_data(fp);
    } else if (pairflag == IJ) {
      fmt::print(fp,"\nPairIJ Coeffs # {}\n\n", force->pair_style);
      force->pair->write_data_all(fp);
    }
  }
  if (force->bond && force->bond->writedata && atom->nbondtypes) {
    fmt::print(fp,"\nBond Coeffs # {}\n\n", force->bond_style);
    force->bond->write_data(fp);
  }
  if (force->angle && force->angle->writedata && atom->nangletypes) {
    fmt::print(fp,"\nAngle Coeffs # {}\n\n", force->angle_style);
    force->angle->write_data(fp);
  }
  if (force->dihedral && force->dihedral->writedata && atom->ndihedraltypes) {
    fmt::print(fp,"\nDihedral Coeffs # {}\n\n", force->dihedral_style);
    force->dihedral->write_data(fp);
  }
  if (force->improper && force->improper->writedata && atom->nimpropertypes) {
    fmt::print(fp,"\nImproper Coeffs # {}\n\n", force->improper_style);
    force->improper->write_data(fp);
  }
}

/* ----------------------------------------------------------------------
   write out Atoms section of data file
------------------------------------------------------------------------- */

void WriteData::atoms()
{
  // communication buffer for all my Atom info
  // maxrow X ncol = largest buffer needed by any proc

  int ncol = atom->avec->size_data_atom + 3;
  int sendrow = atom->nlocal;
  int maxrow;
  MPI_Allreduce(&sendrow,&maxrow,1,MPI_INT,MPI_MAX,world);

  double **buf;
  if (me == 0) memory->create(buf,MAX(1,maxrow),ncol,"write_data:buf");
  else memory->create(buf,MAX(1,sendrow),ncol,"write_data:buf");

  // pack my atom data into buf

  atom->avec->pack_data(buf);

  // write one chunk of atoms per proc to file
  // proc 0 pings each proc, receives its chunk, writes to file
  // all other procs wait for ping, send their chunk to proc 0

  int tmp,recvrow;

  if (me == 0) {
    MPI_Status status;
    MPI_Request request;

    fmt::print(fp,"\nAtoms # {}\n\n",atom->atom_style);
    for (int iproc = 0; iproc < nprocs; iproc++) {
      if (iproc) {
        MPI_Irecv(&buf[0][0],maxrow*ncol,MPI_DOUBLE,iproc,0,world,&request);
        MPI_Send(&tmp,0,MPI_INT,iproc,0,world);
        MPI_Wait(&request,&status);
        MPI_Get_count(&status,MPI_DOUBLE,&recvrow);
        recvrow /= ncol;
      } else recvrow = sendrow;

      atom->avec->write_data(fp,recvrow,buf);
    }

  } else {
    MPI_Recv(&tmp,0,MPI_INT,0,0,world,MPI_STATUS_IGNORE);
    MPI_Rsend(&buf[0][0],sendrow*ncol,MPI_DOUBLE,0,0,world);
  }

  memory->destroy(buf);
}

/* ----------------------------------------------------------------------
   write out Velocities section of data file
------------------------------------------------------------------------- */

void WriteData::velocities()
{
  // communication buffer for all my Atom info
  // maxrow X ncol = largest buffer needed by any proc

  int ncol = atom->avec->size_velocity + 1;
  int sendrow = atom->nlocal;
  int maxrow;
  MPI_Allreduce(&sendrow,&maxrow,1,MPI_INT,MPI_MAX,world);

  double **buf;
  if (me == 0) memory->create(buf,MAX(1,maxrow),ncol,"write_data:buf");
  else memory->create(buf,MAX(1,sendrow),ncol,"write_data:buf");

  // pack my velocity data into buf

  atom->avec->pack_vel(buf);

  // write one chunk of velocities per proc to file
  // proc 0 pings each proc, receives its chunk, writes to file
  // all other procs wait for ping, send their chunk to proc 0

  int tmp,recvrow;

  if (me == 0) {
    MPI_Status status;
    MPI_Request request;

    fputs("\nVelocities\n\n",fp);
    for (int iproc = 0; iproc < nprocs; iproc++) {
      if (iproc) {
        MPI_Irecv(&buf[0][0],maxrow*ncol,MPI_DOUBLE,iproc,0,world,&request);
        MPI_Send(&tmp,0,MPI_INT,iproc,0,world);
        MPI_Wait(&request,&status);
        MPI_Get_count(&status,MPI_DOUBLE,&recvrow);
        recvrow /= ncol;
      } else recvrow = sendrow;

      atom->avec->write_vel(fp,recvrow,buf);
    }

  } else {
    MPI_Recv(&tmp,0,MPI_INT,0,0,world,MPI_STATUS_IGNORE);
    MPI_Rsend(&buf[0][0],sendrow*ncol,MPI_DOUBLE,0,0,world);
  }

  memory->destroy(buf);
}

/* ----------------------------------------------------------------------
   write out Bonds section of data file
------------------------------------------------------------------------- */

void WriteData::bonds()
{
  // communication buffer for all my Bond info
  // maxrow X ncol = largest buffer needed by any proc

  int ncol = 3;
  int sendrow = static_cast<int> (nbonds_local);
  int maxrow;
  MPI_Allreduce(&sendrow,&maxrow,1,MPI_INT,MPI_MAX,world);

  tagint **buf;
  if (me == 0) memory->create(buf,MAX(1,maxrow),ncol,"write_data:buf");
  else memory->create(buf,MAX(1,sendrow),ncol,"write_data:buf");

  // pack my bond data into buf

  atom->avec->pack_bond(buf);

  // write one chunk of info per proc to file
  // proc 0 pings each proc, receives its chunk, writes to file
  // all other procs wait for ping, send their chunk to proc 0

  int tmp,recvrow;

  int index = 1;
  if (me == 0) {
    MPI_Status status;
    MPI_Request request;

    fputs("\nBonds\n\n",fp);
    for (int iproc = 0; iproc < nprocs; iproc++) {
      if (iproc) {
        MPI_Irecv(&buf[0][0],maxrow*ncol,MPI_LMP_TAGINT,iproc,0,world,&request);
        MPI_Send(&tmp,0,MPI_INT,iproc,0,world);
        MPI_Wait(&request,&status);
        MPI_Get_count(&status,MPI_LMP_TAGINT,&recvrow);
        recvrow /= ncol;
      } else recvrow = sendrow;

      atom->avec->write_bond(fp,recvrow,buf,index);
      index += recvrow;
    }

  } else {
    MPI_Recv(&tmp,0,MPI_INT,0,0,world,MPI_STATUS_IGNORE);
    MPI_Rsend(&buf[0][0],sendrow*ncol,MPI_LMP_TAGINT,0,0,world);
  }

  memory->destroy(buf);
}

/* ----------------------------------------------------------------------
   write out Angles section of data file
------------------------------------------------------------------------- */

void WriteData::angles()
{
  // communication buffer for all my Angle info
  // maxrow X ncol = largest buffer needed by any proc

  int ncol = 4;
  int sendrow = static_cast<int> (nangles_local);
  int maxrow;
  MPI_Allreduce(&sendrow,&maxrow,1,MPI_INT,MPI_MAX,world);

  tagint **buf;
  if (me == 0) memory->create(buf,MAX(1,maxrow),ncol,"write_data:buf");
  else memory->create(buf,MAX(1,sendrow),ncol,"write_data:buf");

  // pack my angle data into buf

  atom->avec->pack_angle(buf);

  // write one chunk of info per proc to file
  // proc 0 pings each proc, receives its chunk, writes to file
  // all other procs wait for ping, send their chunk to proc 0

  int tmp,recvrow;

  int index = 1;
  if (me == 0) {
    MPI_Status status;
    MPI_Request request;

    fputs("\nAngles\n\n",fp);
    for (int iproc = 0; iproc < nprocs; iproc++) {
      if (iproc) {
        MPI_Irecv(&buf[0][0],maxrow*ncol,MPI_LMP_TAGINT,iproc,0,world,&request);
        MPI_Send(&tmp,0,MPI_INT,iproc,0,world);
        MPI_Wait(&request,&status);
        MPI_Get_count(&status,MPI_LMP_TAGINT,&recvrow);
        recvrow /= ncol;
      } else recvrow = sendrow;

      atom->avec->write_angle(fp,recvrow,buf,index);
      index += recvrow;
    }

  } else {
    MPI_Recv(&tmp,0,MPI_INT,0,0,world,MPI_STATUS_IGNORE);
    MPI_Rsend(&buf[0][0],sendrow*ncol,MPI_LMP_TAGINT,0,0,world);
  }

  memory->destroy(buf);
}

/* ----------------------------------------------------------------------
   write out Dihedrals section of data file
------------------------------------------------------------------------- */

void WriteData::dihedrals()
{
  // communication buffer for all my Dihedral info
  // maxrow X ncol = largest buffer needed by any proc

  int ncol = 5;
  int sendrow = static_cast<int> (ndihedrals_local);
  int maxrow;
  MPI_Allreduce(&sendrow,&maxrow,1,MPI_INT,MPI_MAX,world);

  tagint **buf;
  if (me == 0) memory->create(buf,MAX(1,maxrow),ncol,"write_data:buf");
  else memory->create(buf,MAX(1,sendrow),ncol,"write_data:buf");

  // pack my dihedral data into buf

  atom->avec->pack_dihedral(buf);

  // write one chunk of info per proc to file
  // proc 0 pings each proc, receives its chunk, writes to file
  // all other procs wait for ping, send their chunk to proc 0

  int tmp,recvrow;

  int index = 1;
  if (me == 0) {
    MPI_Status status;
    MPI_Request request;

    fputs("\nDihedrals\n\n",fp);
    for (int iproc = 0; iproc < nprocs; iproc++) {
      if (iproc) {
        MPI_Irecv(&buf[0][0],maxrow*ncol,MPI_LMP_TAGINT,iproc,0,world,&request);
        MPI_Send(&tmp,0,MPI_INT,iproc,0,world);
        MPI_Wait(&request,&status);
        MPI_Get_count(&status,MPI_LMP_TAGINT,&recvrow);
        recvrow /= ncol;
      } else recvrow = sendrow;

      atom->avec->write_dihedral(fp,recvrow,buf,index);
      index += recvrow;
    }

  } else {
    MPI_Recv(&tmp,0,MPI_INT,0,0,world,MPI_STATUS_IGNORE);
    MPI_Rsend(&buf[0][0],sendrow*ncol,MPI_LMP_TAGINT,0,0,world);
  }

  memory->destroy(buf);
}

/* ----------------------------------------------------------------------
   write out Impropers section of data file
------------------------------------------------------------------------- */

void WriteData::impropers()
{
  // communication buffer for all my Improper info
  // maxrow X ncol = largest buffer needed by any proc

  int ncol = 5;
  int sendrow = static_cast<int> (nimpropers_local);
  int maxrow;
  MPI_Allreduce(&sendrow,&maxrow,1,MPI_INT,MPI_MAX,world);

  tagint **buf;
  if (me == 0) memory->create(buf,MAX(1,maxrow),ncol,"write_data:buf");
  else memory->create(buf,MAX(1,sendrow),ncol,"write_data:buf");

  // pack my improper data into buf

  atom->avec->pack_improper(buf);

  // write one chunk of info per proc to file
  // proc 0 pings each proc, receives its chunk, writes to file
  // all other procs wait for ping, send their chunk to proc 0

  int tmp,recvrow;

  int index = 1;
  if (me == 0) {
    MPI_Status status;
    MPI_Request request;

    fputs("\nImpropers\n\n",fp);
    for (int iproc = 0; iproc < nprocs; iproc++) {
      if (iproc) {
        MPI_Irecv(&buf[0][0],maxrow*ncol,MPI_LMP_TAGINT,iproc,0,world,&request);
        MPI_Send(&tmp,0,MPI_INT,iproc,0,world);
        MPI_Wait(&request,&status);
        MPI_Get_count(&status,MPI_LMP_TAGINT,&recvrow);
        recvrow /= ncol;
      } else recvrow = sendrow;

      atom->avec->write_improper(fp,recvrow,buf,index);
      index += recvrow;
    }

  } else {
    MPI_Recv(&tmp,0,MPI_INT,0,0,world,MPI_STATUS_IGNORE);
    MPI_Rsend(&buf[0][0],sendrow*ncol,MPI_LMP_TAGINT,0,0,world);
  }

  memory->destroy(buf);
}

/* ----------------------------------------------------------------------
   write out Bonus sections of data file
   flag indicates which bonus section it is
------------------------------------------------------------------------- */

void WriteData::bonus(int flag)
{
  // communication buffer for all my Bonus info
  // maxvalues = largest buffer needed by any proc

  int nvalues = atom->avec->pack_data_bonus(nullptr,flag);
  int maxvalues;
  MPI_Allreduce(&nvalues,&maxvalues,1,MPI_INT,MPI_MAX,world);

  double *buf = nullptr;
  if (me == 0) memory->create(buf,MAX(1,maxvalues),"write_data:buf");
  else memory->create(buf,MAX(1,nvalues),"write_data:buf");

  // pack my bonus data into buf

  atom->avec->pack_data_bonus(buf,flag);

  // write one chunk of info per proc to file
  // proc 0 pings each proc, receives its chunk, writes to file
  // all other procs wait for ping, send their chunk to proc 0

  int tmp;

  if (me == 0) {
    MPI_Status status;
    MPI_Request request;

    if (flag == ELLIPSOID) fputs("\nEllipsoids\n\n",fp);
    if (flag == LINE)      fputs("\nLines\n\n",fp);
    if (flag == TRIANGLE)  fputs("\nTriangles\n\n",fp);
    if (flag == BODY)      fputs("\nBodies\n\n",fp);

    for (int iproc = 0; iproc < nprocs; iproc++) {
      if (iproc) {
        MPI_Irecv(buf,maxvalues,MPI_DOUBLE,iproc,0,world,&request);
        MPI_Send(&tmp,0,MPI_INT,iproc,0,world);
        MPI_Wait(&request,&status);
        MPI_Get_count(&status,MPI_DOUBLE,&nvalues);
      }

      atom->avec->write_data_bonus(fp,nvalues,buf,flag);
    }

  } else {
    MPI_Recv(&tmp,0,MPI_INT,0,0,world,MPI_STATUS_IGNORE);
    MPI_Rsend(buf,nvalues,MPI_DOUBLE,0,0,world);
  }

  memory->destroy(buf);
}

/* ----------------------------------------------------------------------
   write out Mth section of data file owned by Fix ifix
------------------------------------------------------------------------- */

void WriteData::fix(Fix *ifix, int mth)
{
  // communication buffer for Fix info
  // maxrow X ncol = largest buffer needed by any proc

  int sendrow,ncol;
  ifix->write_data_section_size(mth,sendrow,ncol);
  int maxrow;
  MPI_Allreduce(&sendrow,&maxrow,1,MPI_INT,MPI_MAX,world);

  double **buf;
  if (me == 0) memory->create(buf,MAX(1,maxrow),ncol,"write_data:buf");
  else memory->create(buf,MAX(1,sendrow),ncol,"write_data:buf");

  // pack my fix data into buf

  ifix->write_data_section_pack(mth,buf);

  // write one chunk of info per proc to file
  // proc 0 pings each proc, receives its chunk, writes to file
  // all other procs wait for ping, send their chunk to proc 0

  int tmp,recvrow;

  int index = 1;
  if (me == 0) {
    MPI_Status status;
    MPI_Request request;

    ifix->write_data_section_keyword(mth,fp);
    for (int iproc = 0; iproc < nprocs; iproc++) {
      if (iproc) {
        MPI_Irecv(&buf[0][0],maxrow*ncol,MPI_DOUBLE,iproc,0,world,&request);
        MPI_Send(&tmp,0,MPI_INT,iproc,0,world);
        MPI_Wait(&request,&status);
        MPI_Get_count(&status,MPI_DOUBLE,&recvrow);
        recvrow /= ncol;
      } else recvrow = sendrow;

      ifix->write_data_section(mth,fp,recvrow,buf,index);
      index += recvrow;
    }

  } else {
    MPI_Recv(&tmp,0,MPI_INT,0,0,world,MPI_STATUS_IGNORE);
    MPI_Rsend(&buf[0][0],sendrow*ncol,MPI_DOUBLE,0,0,world);
  }

  memory->destroy(buf);
}

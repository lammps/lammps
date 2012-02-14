/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "lmptype.h"
#include "mpi.h"
#include "string.h"
#include "write_restart.h"
#include "atom.h"
#include "atom_vec.h"
#include "atom_vec_hybrid.h"
#include "group.h"
#include "force.h"
#include "pair.h"
#include "bond.h"
#include "angle.h"
#include "dihedral.h"
#include "improper.h"
#include "update.h"
#include "neighbor.h"
#include "domain.h"
#include "modify.h"
#include "universe.h"
#include "comm.h"
#include "output.h"
#include "thermo.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

// same as read_restart.cpp and tools/restart2data.cpp

enum{VERSION,SMALLINT,TAGINT,BIGINT,
       UNITS,NTIMESTEP,DIMENSION,NPROCS,PROCGRID_0,PROCGRID_1,PROCGRID_2,
       NEWTON_PAIR,NEWTON_BOND,XPERIODIC,YPERIODIC,ZPERIODIC,
       BOUNDARY_00,BOUNDARY_01,BOUNDARY_10,BOUNDARY_11,BOUNDARY_20,BOUNDARY_21,
       ATOM_STYLE,NATOMS,NTYPES,
       NBONDS,NBONDTYPES,BOND_PER_ATOM,
       NANGLES,NANGLETYPES,ANGLE_PER_ATOM,
       NDIHEDRALS,NDIHEDRALTYPES,DIHEDRAL_PER_ATOM,
       NIMPROPERS,NIMPROPERTYPES,IMPROPER_PER_ATOM,
       BOXLO_0,BOXHI_0,BOXLO_1,BOXHI_1,BOXLO_2,BOXHI_2,
       SPECIAL_LJ_1,SPECIAL_LJ_2,SPECIAL_LJ_3,
       SPECIAL_COUL_1,SPECIAL_COUL_2,SPECIAL_COUL_3,
       XY,XZ,YZ};
enum{MASS};
enum{PAIR,BOND,ANGLE,DIHEDRAL,IMPROPER};

enum{IGNORE,WARN,ERROR};                    // same as thermo.cpp

/* ---------------------------------------------------------------------- */

WriteRestart::WriteRestart(LAMMPS *lmp) : Pointers(lmp)
{
  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);
}

/* ----------------------------------------------------------------------
   called as write_restart command in input script
------------------------------------------------------------------------- */

void WriteRestart::command(int narg, char **arg)
{
  if (domain->box_exist == 0)
    error->all(FLERR,"Write_restart command before simulation box is defined");
  if (narg != 1) error->all(FLERR,"Illegal write_restart command");

  // if filename contains a "*", replace with current timestep

  char *ptr;
  int n = strlen(arg[0]) + 16;
  char *file = new char[n];

  if (ptr = strchr(arg[0],'*')) {
    *ptr = '\0';
    sprintf(file,"%s" BIGINT_FORMAT "%s",arg[0],update->ntimestep,ptr+1);
  } else strcpy(file,arg[0]);

  // init entire system since comm->exchange is done
  // comm::init needs neighbor::init needs pair::init needs kspace::init, etc

  if (comm->me == 0 && screen)
    fprintf(screen,"System init for write_restart ...\n");
  lmp->init();

  // move atoms to new processors before writing file
  // enforce PBC before in case atoms are outside box
  // call borders() to rebuild atom map since exchange() destroys map

  if (domain->triclinic) domain->x2lamda(atom->nlocal);
  domain->pbc();
  domain->reset_box();
  comm->setup();
  comm->exchange();
  comm->borders();
  if (domain->triclinic) domain->lamda2x(atom->nlocal+atom->nghost);

  write(file);
  delete [] file;
}

/* ----------------------------------------------------------------------
   called from command() and directly from output within run/minimize loop
   file = final file name to write, except may contain a "%"
------------------------------------------------------------------------- */

void WriteRestart::write(char *file)
{
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
  if (natoms != atom->natoms && output->thermo->lostflag == ERROR) 
    error->all(FLERR,"Atom count is inconsistent, cannot write restart file");

  // check if filename contains "%"

  int multiproc;
  if (strchr(file,'%')) multiproc = 1;
  else multiproc = 0;

  // open single restart file or base file for multiproc case

  if (me == 0) {
    char *hfile;
    if (multiproc) {
      hfile = new char[strlen(file) + 16];
      char *ptr = strchr(file,'%');
      *ptr = '\0';
      sprintf(hfile,"%s%s%s",file,"base",ptr+1);
      *ptr = '%';
    } else hfile = file;
    fp = fopen(hfile,"wb");
    if (fp == NULL) {
      char str[128];
      sprintf(str,"Cannot open restart file %s",hfile);
      error->one(FLERR,str);
    }
    if (multiproc) delete [] hfile;
  }

  // proc 0 writes header, groups, ntype-length arrays, force field
  // all procs write fix info

  if (me == 0) {
    header();
    group->write_restart(fp);
    type_arrays();
    force_fields();
  }

  modify->write_restart(fp);

  // communication buffer for all my atom's info
  // max_size = largest buffer needed by any proc

  int max_size;
  int send_size = atom->avec->size_restart();
  MPI_Allreduce(&send_size,&max_size,1,MPI_INT,MPI_MAX,world);

  double *buf;
  if (me == 0) memory->create(buf,max_size,"write_restart:buf");
  else memory->create(buf,send_size,"write_restart:buf");

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

  // if single file:
  //   write one chunk of atoms per proc to file
  //   proc 0 pings each proc, receives its chunk, writes to file
  //   all other procs wait for ping, send their chunk to proc 0
  // else if one file per proc:
  //   each proc opens its own file and writes its chunk directly

  if (multiproc == 0) {
    int tmp,recv_size;
    MPI_Status status;
    MPI_Request request;

    if (me == 0) {
      for (int iproc = 0; iproc < nprocs; iproc++) {
	if (iproc) {
	  MPI_Irecv(buf,max_size,MPI_DOUBLE,iproc,0,world,&request);
	  MPI_Send(&tmp,0,MPI_INT,iproc,0,world);
	  MPI_Wait(&request,&status);
	  MPI_Get_count(&status,MPI_DOUBLE,&recv_size);
	} else recv_size = send_size;
	
	fwrite(&recv_size,sizeof(int),1,fp);
	fwrite(buf,sizeof(double),recv_size,fp);
      }
      fclose(fp);

    } else {
      MPI_Recv(&tmp,0,MPI_INT,0,0,world,&status);
      MPI_Rsend(buf,send_size,MPI_DOUBLE,0,0,world);
    }

  } else {
    if (me == 0) fclose(fp);

    char *perproc = new char[strlen(file) + 16];
    char *ptr = strchr(file,'%');
    *ptr = '\0';
    sprintf(perproc,"%s%d%s",file,me,ptr+1);
    *ptr = '%';
    fp = fopen(perproc,"wb");
    if (fp == NULL) {
      char str[128];
      sprintf(str,"Cannot open restart file %s",perproc);
      error->one(FLERR,str);
    }
    delete [] perproc;
    fwrite(&send_size,sizeof(int),1,fp);
    fwrite(buf,sizeof(double),send_size,fp);
    fclose(fp);
  }
    
  memory->destroy(buf);
}

/* ----------------------------------------------------------------------
   proc 0 writes out problem description 
------------------------------------------------------------------------- */

void WriteRestart::header()
{
  write_char(VERSION,universe->version);
  write_int(SMALLINT,sizeof(smallint));
  write_int(TAGINT,sizeof(tagint));
  write_int(BIGINT,sizeof(bigint));
  write_char(UNITS,update->unit_style);
  write_bigint(NTIMESTEP,update->ntimestep);
  write_int(DIMENSION,domain->dimension);
  write_int(NPROCS,nprocs);
  write_int(PROCGRID_0,comm->procgrid[0]);
  write_int(PROCGRID_1,comm->procgrid[1]);
  write_int(PROCGRID_2,comm->procgrid[2]);
  write_int(NEWTON_PAIR,force->newton_pair);
  write_int(NEWTON_BOND,force->newton_bond);
  write_int(XPERIODIC,domain->xperiodic);
  write_int(YPERIODIC,domain->yperiodic);
  write_int(ZPERIODIC,domain->zperiodic);
  write_int(BOUNDARY_00,domain->boundary[0][0]);
  write_int(BOUNDARY_01,domain->boundary[0][1]);
  write_int(BOUNDARY_10,domain->boundary[1][0]);
  write_int(BOUNDARY_11,domain->boundary[1][1]);
  write_int(BOUNDARY_20,domain->boundary[2][0]);
  write_int(BOUNDARY_21,domain->boundary[2][1]);

  // atom_style must be written before atom class values
  // so read_restart can create class before reading class values
  // if style = hybrid, also write sub-class styles

  write_char(ATOM_STYLE,atom->atom_style);

  if (strcmp(atom->atom_style,"hybrid") == 0) {
    AtomVecHybrid *avec_hybrid = (AtomVecHybrid *) atom->avec;
    int nstyles = avec_hybrid->nstyles;
    char **keywords = avec_hybrid->keywords;
    fwrite(&nstyles,sizeof(int),1,fp);
    for (int i = 0; i < nstyles; i++) {
      int n = strlen(keywords[i]) + 1;
      fwrite(&n,sizeof(int),1,fp);
      fwrite(keywords[i],sizeof(char),n,fp);
    }
  }

  write_bigint(NATOMS,natoms);
  write_int(NTYPES,atom->ntypes);
  write_bigint(NBONDS,atom->nbonds);
  write_int(NBONDTYPES,atom->nbondtypes);
  write_int(BOND_PER_ATOM,atom->bond_per_atom);
  write_bigint(NANGLES,atom->nangles);
  write_int(NANGLETYPES,atom->nangletypes);
  write_int(ANGLE_PER_ATOM,atom->angle_per_atom);
  write_bigint(NDIHEDRALS,atom->ndihedrals);
  write_int(NDIHEDRALTYPES,atom->ndihedraltypes);
  write_int(DIHEDRAL_PER_ATOM,atom->dihedral_per_atom);
  write_bigint(NIMPROPERS,atom->nimpropers);
  write_int(NIMPROPERTYPES,atom->nimpropertypes);
  write_int(IMPROPER_PER_ATOM,atom->improper_per_atom);

  write_double(BOXLO_0,domain->boxlo[0]);
  write_double(BOXHI_0,domain->boxhi[0]);
  write_double(BOXLO_1,domain->boxlo[1]);
  write_double(BOXHI_1,domain->boxhi[1]);
  write_double(BOXLO_2,domain->boxlo[2]);
  write_double(BOXHI_2,domain->boxhi[2]);

  write_double(SPECIAL_LJ_1,force->special_lj[1]);
  write_double(SPECIAL_LJ_2,force->special_lj[2]);
  write_double(SPECIAL_LJ_3,force->special_lj[3]);
  write_double(SPECIAL_COUL_1,force->special_coul[1]);
  write_double(SPECIAL_COUL_2,force->special_coul[2]);
  write_double(SPECIAL_COUL_3,force->special_coul[3]);

  if (domain->triclinic) {
    write_double(XY,domain->xy);
    write_double(XZ,domain->xz);
    write_double(YZ,domain->yz);
  }

  // -1 flag signals end of header

  int flag = -1;
  fwrite(&flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
   proc 0 writes out any type-based arrays that are defined
------------------------------------------------------------------------- */

void WriteRestart::type_arrays()
{
  if (atom->mass) {
    int flag = MASS;
    fwrite(&flag,sizeof(int),1,fp);
    fwrite(&atom->mass[1],sizeof(double),atom->ntypes,fp);
  }

  // -1 flag signals end of type arrays

  int flag = -1;
  fwrite(&flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
   proc 0 writes out and force field styles and data that are defined
------------------------------------------------------------------------- */

void WriteRestart::force_fields()
{
  if (force->pair) {
    int flag = PAIR;
    fwrite(&flag,sizeof(int),1,fp);
    int n = strlen(force->pair_style) + 1;
    fwrite(&n,sizeof(int),1,fp);
    fwrite(force->pair_style,sizeof(char),n,fp);
    force->pair->write_restart(fp);
  }
  if (atom->avec->bonds_allow && force->bond) {
    int flag = BOND;
    fwrite(&flag,sizeof(int),1,fp);
    int n = strlen(force->bond_style) + 1;
    fwrite(&n,sizeof(int),1,fp);
    fwrite(force->bond_style,sizeof(char),n,fp);
    force->bond->write_restart(fp);
  }
  if (atom->avec->angles_allow && force->angle) {
    int flag = ANGLE;
    fwrite(&flag,sizeof(int),1,fp);
    int n = strlen(force->angle_style) + 1;
    fwrite(&n,sizeof(int),1,fp);
    fwrite(force->angle_style,sizeof(char),n,fp);
    force->angle->write_restart(fp);
  }
  if (atom->avec->dihedrals_allow && force->dihedral) {
    int flag = DIHEDRAL;
    fwrite(&flag,sizeof(int),1,fp);
    int n = strlen(force->dihedral_style) + 1;
    fwrite(&n,sizeof(int),1,fp);
    fwrite(force->dihedral_style,sizeof(char),n,fp);
    force->dihedral->write_restart(fp);
  }
  if (atom->avec->impropers_allow && force->improper) {
    int flag = IMPROPER;
    fwrite(&flag,sizeof(int),1,fp);
    int n = strlen(force->improper_style) + 1;
    fwrite(&n,sizeof(int),1,fp);
    fwrite(force->improper_style,sizeof(char),n,fp);
    force->improper->write_restart(fp);
  }

  // -1 flag signals end of force field info

  int flag = -1;
  fwrite(&flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
   write a flag and an int into restart file 
------------------------------------------------------------------------- */

void WriteRestart::write_int(int flag, int value)
{
  fwrite(&flag,sizeof(int),1,fp);
  fwrite(&value,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
   write a flag and a double into restart file 
------------------------------------------------------------------------- */

void WriteRestart::write_double(int flag, double value)
{
  fwrite(&flag,sizeof(int),1,fp);
  fwrite(&value,sizeof(double),1,fp);
}

/* ----------------------------------------------------------------------
   write a flag and a char str into restart file 
------------------------------------------------------------------------- */

void WriteRestart::write_char(int flag, char *value)
{
  fwrite(&flag,sizeof(int),1,fp);
  int n = strlen(value) + 1;
  fwrite(&n,sizeof(int),1,fp);
  fwrite(value,sizeof(char),n,fp);
}

/* ----------------------------------------------------------------------
   write a flag and a bigint into restart file 
------------------------------------------------------------------------- */

void WriteRestart::write_bigint(int flag, bigint value)
{
  fwrite(&flag,sizeof(int),1,fp);
  fwrite(&value,sizeof(bigint),1,fp);
}

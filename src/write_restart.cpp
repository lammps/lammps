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
#include "domain.h"
#include "modify.h"
#include "universe.h"
#include "comm.h"
#include "output.h"
#include "thermo.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

enum{IGNORE,WARN,ERROR};     // same as thermo.cpp

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
    error->all("Write_restart command before simulation box is defined");
  if (narg != 1) error->all("Illegal write_restart command");

  // if filename contains a "*", replace with current timestep

  char *ptr;
  int n = strlen(arg[0]) + 16;
  char *file = new char[n];

  if (ptr = strchr(arg[0],'*')) {
    *ptr = '\0';
    sprintf(file,"%s%d%s",arg[0],update->ntimestep,ptr+1);
  } else strcpy(file,arg[0]);

  // init entire system since comm->exchange is done
  // comm::init needs neighbor::init needs pair::init needs kspace::init, etc

  if (comm->me == 0 && screen)
    fprintf(screen,"System init for write_restart ...\n");
  lmp->init();

  // move atoms to new processors before writing file
  // enforce PBC before in case atoms are outside box

  if (domain->triclinic) domain->x2lamda(atom->nlocal);
  domain->pbc();
  domain->reset_box();
  comm->setup();
  comm->exchange();
  if (domain->triclinic) domain->lamda2x(atom->nlocal);

  write(file);
  delete [] file;
}

/* ----------------------------------------------------------------------
   called from output within run/minimize loop
   file = final file name to write, except may contain a "%"
------------------------------------------------------------------------- */

void WriteRestart::write(char *file)
{
  // natoms = sum of nlocal = value to write into restart file
  // if unequal and thermo lostflag is "error", don't write restart file

  double rlocal = atom->nlocal;
  MPI_Allreduce(&rlocal,&natoms,1,MPI_DOUBLE,MPI_SUM,world);
  if (natoms != atom->natoms && output->thermo->lostflag == ERROR) 
    error->all("Atom count is inconsistent, cannot write restart file");

  // check if filename contains "%"

  int multiproc;
  if (strchr(file,'%')) multiproc = 1;
  else multiproc = 0;

  // open single restart file or base file for multiproc case

  if (me == 0) {
    char *hfile;
    if (multiproc) {
      char *hfile = new char[strlen(file) + 16];
      char *ptr = strchr(file,'%');
      *ptr = '\0';
      sprintf(hfile,"%s%s%s",file,"base",ptr+1);
      *ptr = '%';
    } else hfile = file;
    fp = fopen(hfile,"wb");
    if (fp == NULL) {
      char str[128];
      sprintf(str,"Cannot open restart file %s",hfile);
      error->one(str);
    }
    if (multiproc) delete [] hfile;
  }

  // write header, groups, ntype-length arrays, force field, fix info

  if (me == 0) {
    header();
    group->write_restart(fp);
    if (atom->mass) mass();
    if (atom->dipole) dipole();
    force_fields();
  }

  modify->write_restart(fp);

  // communication buffer for all my atom's info
  // max_size = largest buffer needed by any proc

  int max_size;
  int send_size = atom->avec->size_restart();
  MPI_Allreduce(&send_size,&max_size,1,MPI_INT,MPI_MAX,world);

  double *buf;
  if (me == 0) 
    buf = (double *) 
      memory->smalloc(max_size*sizeof(double),"write_restart:buf");
  else
    buf = (double *) 
      memory->smalloc(send_size*sizeof(double),"write_restart:buf");

  // pack my atom data into buf

  AtomVec *avec = atom->avec;
  int n = 0;
  for (int i = 0; i < atom->nlocal; i++) n += avec->pack_restart(i,&buf[n]);

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
    delete [] perproc;
    if (fp == NULL) error->one("Cannot open restart file");

    fwrite(&send_size,sizeof(int),1,fp);
    fwrite(buf,sizeof(double),send_size,fp);
    fclose(fp);
  }
    
  memory->sfree(buf);
}

/* ----------------------------------------------------------------------
   proc 0 writes out problem description 
------------------------------------------------------------------------- */

void WriteRestart::header()
{
  write_char(0,universe->version);
  write_char(1,update->unit_style);
  write_int(2,update->ntimestep);
  write_int(3,force->dimension);
  write_int(4,nprocs);
  write_int(5,comm->procgrid[0]);
  write_int(6,comm->procgrid[1]);
  write_int(7,comm->procgrid[2]);
  write_int(8,force->newton_pair);
  write_int(9,force->newton_bond);
  write_int(10,domain->xperiodic);
  write_int(11,domain->yperiodic);
  write_int(12,domain->zperiodic);
  write_int(13,domain->boundary[0][0]);
  write_int(14,domain->boundary[0][1]);
  write_int(15,domain->boundary[1][0]);
  write_int(16,domain->boundary[1][1]);
  write_int(17,domain->boundary[2][0]);
  write_int(18,domain->boundary[2][1]);

  // atom_style must be written before atom class values
  // so read_restart can create class before reading class values
  // if style = hybrid, also write sub-class styles

  write_char(19,atom->atom_style);

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

  write_double(20,natoms);
  write_int(21,atom->ntypes);
  write_int(22,atom->nbonds);
  write_int(23,atom->nbondtypes);
  write_int(24,atom->bond_per_atom);
  write_int(25,atom->nangles);
  write_int(26,atom->nangletypes);
  write_int(27,atom->angle_per_atom);
  write_int(28,atom->ndihedrals);
  write_int(29,atom->ndihedraltypes);
  write_int(30,atom->dihedral_per_atom);
  write_int(31,atom->nimpropers);
  write_int(32,atom->nimpropertypes);
  write_int(33,atom->improper_per_atom);

  write_double(34,domain->boxxlo);
  write_double(35,domain->boxxhi);
  write_double(36,domain->boxylo);
  write_double(37,domain->boxyhi);
  write_double(38,domain->boxzlo);
  write_double(39,domain->boxzhi);

  write_double(40,force->special_lj[1]);
  write_double(41,force->special_lj[2]);
  write_double(42,force->special_lj[3]);
  write_double(43,force->special_coul[1]);
  write_double(44,force->special_coul[2]);
  write_double(45,force->special_coul[3]);

  if (domain->triclinic) {
    write_double(46,domain->xy);
    write_double(47,domain->xz);
    write_double(48,domain->yz);
  }

  // -1 flag signals end of header

  int flag = -1;
  fwrite(&flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
   proc 0 writes out atom masses 
------------------------------------------------------------------------- */

void WriteRestart::mass()
{
  fwrite(&atom->mass[1],sizeof(double),atom->ntypes,fp);
}

/* ----------------------------------------------------------------------
   proc 0 writes out dipole moments 
------------------------------------------------------------------------- */

void WriteRestart::dipole()
{
  fwrite(&atom->dipole[1],sizeof(double),atom->ntypes,fp);
}

/* ----------------------------------------------------------------------
   proc 0 writes out force field styles and data 
------------------------------------------------------------------------- */

void WriteRestart::force_fields()
{
  int n;

  if (force->pair) n = strlen(force->pair_style) + 1;
  else n = 0;
  fwrite(&n,sizeof(int),1,fp);
  if (n) {
    fwrite(force->pair_style,sizeof(char),n,fp);
    force->pair->write_restart(fp);
  }

  if (atom->avec->bonds_allow) {
    if (force->bond) n = strlen(force->bond_style) + 1;
    else n = 0;
    fwrite(&n,sizeof(int),1,fp);
    if (n) {
      fwrite(force->bond_style,sizeof(char),n,fp);
      force->bond->write_restart(fp);
    }
  }

  if (atom->avec->angles_allow) {
    if (force->angle) n = strlen(force->angle_style) + 1;
    else n = 0;
    fwrite(&n,sizeof(int),1,fp);
    if (n) {
      fwrite(force->angle_style,sizeof(char),n,fp);
      force->angle->write_restart(fp);
    }
  }

  if (atom->avec->dihedrals_allow) {
    if (force->dihedral) n = strlen(force->dihedral_style) + 1;
    else n = 0;
    fwrite(&n,sizeof(int),1,fp);
    if (n) {
      fwrite(force->dihedral_style,sizeof(char),n,fp);
      force->dihedral->write_restart(fp);
    }
  }

  if (atom->avec->impropers_allow) {
    if (force->improper) n = strlen(force->improper_style) + 1;
    else n = 0;
    fwrite(&n,sizeof(int),1,fp);
    if (n) {
      fwrite(force->improper_style,sizeof(char),n,fp);
      force->improper->write_restart(fp);
    }
  }
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

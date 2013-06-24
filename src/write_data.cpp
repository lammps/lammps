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
#include "write_data.h"
#include "atom.h"
#include "atom_vec.h"
#include "group.h"
#include "force.h"
#include "pair.h"
#include "bond.h"
#include "angle.h"
#include "dihedral.h"
#include "improper.h"
#include "update.h"
#include "modify.h"
#include "domain.h"
#include "universe.h"
#include "comm.h"
#include "output.h"
#include "thermo.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

enum{IGNORE,WARN,ERROR};                    // same as thermo.cpp
enum{II,IJ};

/* ---------------------------------------------------------------------- */

WriteData::WriteData(LAMMPS *lmp) : Pointers(lmp)
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

  if (narg < 1) error->all(FLERR,"Illegal write_data command");

  // if filename contains a "*", replace with current timestep

  char *ptr;
  int n = strlen(arg[0]) + 16;
  char *file = new char[n];

  if (ptr = strchr(arg[0],'*')) {
    *ptr = '\0';
    sprintf(file,"%s" BIGINT_FORMAT "%s",arg[0],update->ntimestep,ptr+1);
  } else strcpy(file,arg[0]);

  // read optional args

  pairflag = II;

  int iarg = 1;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"pair") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal write_data command");
      if (strcmp(arg[iarg+1],"ii") == 0) pairflag = II;
      else if (strcmp(arg[iarg+1],"ij") == 0) pairflag = IJ;
      else error->all(FLERR,"Illegal write_data command");
      iarg += 2;
    } else error->all(FLERR,"Illegal write_data command");
  }

  // init entire system since comm->exchange is done
  // comm::init needs neighbor::init needs pair::init needs kspace::init, etc

  if (comm->me == 0 && screen)
    fprintf(screen,"System init for write_data ...\n");
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

  write(file);
  delete [] file;
}

/* ----------------------------------------------------------------------
   called from command()
   later might let it be directly called within run/minimize loop
------------------------------------------------------------------------- */

void WriteData::write(char *file)
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
  if (natoms != atom->natoms && output->thermo->lostflag == ERROR)
    error->all(FLERR,"Atom count is inconsistent, cannot write data file");

  // sum up bond,angle counts
  // may be different than atom->nbonds,nangles if broken/turned-off

  if (atom->nbonds || atom->nbondtypes) {
    nbonds_local = atom->avec->pack_bond(NULL);
    MPI_Allreduce(&nbonds_local,&nbonds,1,MPI_LMP_BIGINT,MPI_SUM,world);
  }
  if (atom->nangles || atom->nangletypes) {
    nangles_local = atom->avec->pack_angle(NULL);
    MPI_Allreduce(&nangles_local,&nangles,1,MPI_LMP_BIGINT,MPI_SUM,world);
  }

  // open data file

  if (me == 0) {
    fp = fopen(file,"w");
    if (fp == NULL) {
      char str[128];
      sprintf(str,"Cannot open data file %s",file);
      error->one(FLERR,str);
    }
  }

  // proc 0 writes header, ntype-length arrays, force fields

  if (me == 0) {
    header();
    type_arrays();
    force_fields();
  }

  // per atom info

  if (natoms) atoms();
  if (natoms) velocities();
  if (atom->nbonds && nbonds) bonds();
  if (atom->nangles && nangles) angles();
  if (atom->ndihedrals) dihedrals();
  if (atom->nimpropers) impropers();

  // close data file

  if (me == 0) fclose(fp);
}

/* ----------------------------------------------------------------------
   proc 0 writes out data file header
------------------------------------------------------------------------- */

void WriteData::header()
{
  fprintf(fp,"LAMMPS data file via write_data, version %s, "
          "timestep = " BIGINT_FORMAT "\n",
          universe->version,update->ntimestep);

  fprintf(fp,"\n");

  fprintf(fp,BIGINT_FORMAT " atoms\n",atom->natoms);
  fprintf(fp,"%d atom types\n",atom->ntypes);

  if (atom->nbonds || atom->nbondtypes) {
    fprintf(fp,BIGINT_FORMAT " bonds\n",nbonds);
    fprintf(fp,"%d bond types\n",atom->nbondtypes);
  }
  if (atom->nangles || atom->nangletypes) {
    fprintf(fp,BIGINT_FORMAT " angles\n",nangles);
    fprintf(fp,"%d angle types\n",atom->nangletypes);
  }
  if (atom->ndihedrals || atom->ndihedraltypes) {
    fprintf(fp,BIGINT_FORMAT " dihedrals\n",atom->ndihedrals);
    fprintf(fp,"%d dihedral types\n",atom->ndihedraltypes);
  }
  if (atom->nimpropers || atom->nimpropertypes) {
    fprintf(fp,BIGINT_FORMAT " impropers\n",atom->nimpropers);
    fprintf(fp,"%d improper types\n",atom->nimpropertypes);
  }

  fprintf(fp,"\n");

  fprintf(fp,"%-1.16e %-1.16e xlo xhi\n",domain->boxlo[0],domain->boxhi[0]);
  fprintf(fp,"%-1.16e %-1.16e ylo yhi\n",domain->boxlo[1],domain->boxhi[1]);
  fprintf(fp,"%-1.16e %-1.16e zlo zhi\n",domain->boxlo[2],domain->boxhi[2]);

  if (domain->triclinic)
    fprintf(fp,"%-1.16e %-1.16e %-1.16e xy xz yz\n",
            domain->xy,domain->xz,domain->yz);
}

/* ----------------------------------------------------------------------
   proc 0 writes out any type-based arrays that are defined
------------------------------------------------------------------------- */

void WriteData::type_arrays()
{
  if (atom->mass) {
    double *mass = atom->mass;
    fprintf(fp,"\nMasses\n\n");
    for (int i = 1; i <= atom->ntypes; i++) fprintf(fp,"%d %g\n",i,mass[i]);
  }
}

/* ----------------------------------------------------------------------
   proc 0 writes out force field info
------------------------------------------------------------------------- */

void WriteData::force_fields()
{
  if (force->pair && force->pair->writedata) {
    if (pairflag == II) {
      fprintf(fp,"\nPair Coeffs\n\n");
      force->pair->write_data(fp);
    } else if (pairflag == IJ) {
      fprintf(fp,"\nPairIJ Coeffs\n\n");
      force->pair->write_data_all(fp);
    }
  }
  if (atom->avec->bonds_allow && force->bond && force->bond->writedata) {
    fprintf(fp,"\nBond Coeffs\n\n");
    force->bond->write_data(fp);
  }
  if (atom->avec->angles_allow && force->angle && force->angle->writedata) {
    fprintf(fp,"\nAngle Coeffs\n\n");
    force->angle->write_data(fp);
  }
  if (atom->avec->dihedrals_allow && force->dihedral && 
      force->dihedral->writedata) {
    fprintf(fp,"\nDihedral Coeffs\n\n");
    force->dihedral->write_data(fp);
  }
  if (atom->avec->impropers_allow && force->improper && 
      force->improper->writedata) {
    fprintf(fp,"\nImproper Coeffs\n\n");
    force->improper->write_data(fp);
  }
}

/* ----------------------------------------------------------------------
   write out Atoms section of data file
------------------------------------------------------------------------- */

void WriteData::atoms()
{
  // communication buffer for all my Atom info
  // max_size = largest buffer needed by any proc

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
  MPI_Status status;
  MPI_Request request;

  if (me == 0) {
    fprintf(fp,"\nAtoms\n\n");
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
    MPI_Recv(&tmp,0,MPI_INT,0,0,world,&status);
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
  // max_size = largest buffer needed by any proc

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
  MPI_Status status;
  MPI_Request request;

  if (me == 0) {
    fprintf(fp,"\nVelocities\n\n");
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
    MPI_Recv(&tmp,0,MPI_INT,0,0,world,&status);
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

  int ncol = 3;
  int sendrow = static_cast<int> (nbonds_local);
  int maxrow;
  MPI_Allreduce(&sendrow,&maxrow,1,MPI_INT,MPI_MAX,world);

  int **buf;
  if (me == 0) memory->create(buf,MAX(1,maxrow),ncol,"write_data:buf");
  else memory->create(buf,MAX(1,sendrow),ncol,"write_data:buf");

  // pack my bond data into buf

  int foo = atom->avec->pack_bond(buf);

  // write one chunk of atoms per proc to file
  // proc 0 pings each proc, receives its chunk, writes to file
  // all other procs wait for ping, send their chunk to proc 0

  int tmp,recvrow;
  MPI_Status status;
  MPI_Request request;

  int index = 1;
  if (me == 0) {
    fprintf(fp,"\nBonds\n\n");
    for (int iproc = 0; iproc < nprocs; iproc++) {
      if (iproc) {
        MPI_Irecv(&buf[0][0],maxrow*ncol,MPI_INT,iproc,0,world,&request);
        MPI_Send(&tmp,0,MPI_INT,iproc,0,world);
        MPI_Wait(&request,&status);
        MPI_Get_count(&status,MPI_INT,&recvrow);
        recvrow /= ncol;
      } else recvrow = sendrow;

      atom->avec->write_bond(fp,recvrow,buf,index);
      index += recvrow;
    }
    
  } else {
    MPI_Recv(&tmp,0,MPI_INT,0,0,world,&status);
    MPI_Rsend(&buf[0][0],sendrow*ncol,MPI_INT,0,0,world);
  }

  memory->destroy(buf);
}

/* ----------------------------------------------------------------------
   write out Angles section of data file
------------------------------------------------------------------------- */

void WriteData::angles()
{
  // communication buffer for all my Angle info

  int ncol = 4;
  int sendrow = static_cast<int> (nangles_local);
  int maxrow;
  MPI_Allreduce(&sendrow,&maxrow,1,MPI_INT,MPI_MAX,world);

  int **buf;
  if (me == 0) memory->create(buf,MAX(1,maxrow),ncol,"write_data:buf");
  else memory->create(buf,MAX(1,sendrow),ncol,"write_data:buf");

  // pack my angle data into buf

  atom->avec->pack_angle(buf);

  // write one chunk of atoms per proc to file
  // proc 0 pings each proc, receives its chunk, writes to file
  // all other procs wait for ping, send their chunk to proc 0

  int tmp,recvrow;
  MPI_Status status;
  MPI_Request request;

  int index = 1;
  if (me == 0) {
    fprintf(fp,"\nAngles\n\n");
    for (int iproc = 0; iproc < nprocs; iproc++) {
      if (iproc) {
        MPI_Irecv(&buf[0][0],maxrow*ncol,MPI_INT,iproc,0,world,&request);
        MPI_Send(&tmp,0,MPI_INT,iproc,0,world);
        MPI_Wait(&request,&status);
        MPI_Get_count(&status,MPI_INT,&recvrow);
        recvrow /= ncol;
      } else recvrow = sendrow;
      
      atom->avec->write_angle(fp,recvrow,buf,index);
      index += recvrow;
    }
    
  } else {
    MPI_Recv(&tmp,0,MPI_INT,0,0,world,&status);
    MPI_Rsend(&buf[0][0],sendrow*ncol,MPI_INT,0,0,world);
  }

  memory->destroy(buf);
}

/* ----------------------------------------------------------------------
   write out Dihedrals section of data file
------------------------------------------------------------------------- */

void WriteData::dihedrals()
{
  // communication buffer for all my Dihedral info
  // max_size = largest buffer needed by any proc

  int ncol = 5;

  int *tag = atom->tag;
  int *num_dihedral = atom->num_dihedral;
  int **dihedral_atom2 = atom->dihedral_atom2;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;

  int i,j;
  int sendrow = 0;
  if (newton_bond) {
    for (i = 0; i < nlocal; i++)
      sendrow += num_dihedral[i];
  } else {
    for (i = 0; i < nlocal; i++)
      for (j = 0; j < num_dihedral[i]; j++)
        if (tag[i] == dihedral_atom2[i][j]) sendrow++;
  }

  int maxrow;
  MPI_Allreduce(&sendrow,&maxrow,1,MPI_INT,MPI_MAX,world);

  int **buf;
  if (me == 0) memory->create(buf,MAX(1,maxrow),ncol,"write_data:buf");
  else memory->create(buf,MAX(1,sendrow),ncol,"write_data:buf");

  // pack my dihedral data into buf

  atom->avec->pack_dihedral(buf);

  // write one chunk of atoms per proc to file
  // proc 0 pings each proc, receives its chunk, writes to file
  // all other procs wait for ping, send their chunk to proc 0

  int tmp,recvrow;
  MPI_Status status;
  MPI_Request request;

  int index = 1;
  if (me == 0) {
    fprintf(fp,"\nDihedrals\n\n");
    for (int iproc = 0; iproc < nprocs; iproc++) {
      if (iproc) {
        MPI_Irecv(&buf[0][0],maxrow*ncol,MPI_INT,iproc,0,world,&request);
        MPI_Send(&tmp,0,MPI_INT,iproc,0,world);
        MPI_Wait(&request,&status);
        MPI_Get_count(&status,MPI_INT,&recvrow);
        recvrow /= ncol;
      } else recvrow = sendrow;
      
      atom->avec->write_dihedral(fp,recvrow,buf,index);
      index += recvrow;
    }
    
  } else {
    MPI_Recv(&tmp,0,MPI_INT,0,0,world,&status);
    MPI_Rsend(&buf[0][0],sendrow*ncol,MPI_INT,0,0,world);
  }

  memory->destroy(buf);
}

/* ----------------------------------------------------------------------
   write out Impropers section of data file
------------------------------------------------------------------------- */

void WriteData::impropers()
{
  // communication buffer for all my Improper info
  // max_size = largest buffer needed by any proc

  int ncol = 5;

  int *tag = atom->tag;
  int *num_improper = atom->num_improper;
  int **improper_atom2 = atom->improper_atom2;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;

  int i,j;
  int sendrow = 0;
  if (newton_bond) {
    for (i = 0; i < nlocal; i++)
      sendrow += num_improper[i];
  } else {
    for (i = 0; i < nlocal; i++)
      for (j = 0; j < num_improper[i]; j++)
        if (tag[i] == improper_atom2[i][j]) sendrow++;
  }

  int maxrow;
  MPI_Allreduce(&sendrow,&maxrow,1,MPI_INT,MPI_MAX,world);

  int **buf;
  if (me == 0) memory->create(buf,MAX(1,maxrow),ncol,"write_data:buf");
  else memory->create(buf,MAX(1,sendrow),ncol,"write_data:buf");

  // pack my improper data into buf

  atom->avec->pack_improper(buf);

  // write one chunk of atoms per proc to file
  // proc 0 pings each proc, receives its chunk, writes to file
  // all other procs wait for ping, send their chunk to proc 0

  int tmp,recvrow;
  MPI_Status status;
  MPI_Request request;

  int index = 1;
  if (me == 0) {
    fprintf(fp,"\nImpropers\n\n");
    for (int iproc = 0; iproc < nprocs; iproc++) {
      if (iproc) {
        MPI_Irecv(&buf[0][0],maxrow*ncol,MPI_INT,iproc,0,world,&request);
        MPI_Send(&tmp,0,MPI_INT,iproc,0,world);
        MPI_Wait(&request,&status);
        MPI_Get_count(&status,MPI_INT,&recvrow);
        recvrow /= ncol;
      } else recvrow = sendrow;
      
      atom->avec->write_improper(fp,recvrow,buf,index);
      index += recvrow;
    }
    
  } else {
    MPI_Recv(&tmp,0,MPI_INT,0,0,world,&status);
    MPI_Rsend(&buf[0][0],sendrow*ncol,MPI_INT,0,0,world);
  }

  memory->destroy(buf);
}

// lmpqst = umbrella driver to couple LAMMPS + Quest
//          for MD using quantum forces

// Syntax: lmpqst Niter in.lammps in.quest
//         Niter = # of MD iterations
//         in.lammps = LAMMPS input script
//         in.quest = Quest input script

#include "mpi.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "stdint.h"

#include "many2one.h"
#include "one2many.h"
#include "files.h"
#include "memory.h"
#include "error.h"

#define QUOTE_(x) #x
#define QUOTE(x) QUOTE_(x)

#include "lmppath.h"
#include QUOTE(LMPPATH/src/lammps.h)
#include QUOTE(LMPPATH/src/library.h)
#include QUOTE(LMPPATH/src/input.h)
#include QUOTE(LMPPATH/src/modify.h)
#include QUOTE(LMPPATH/src/fix.h)
#include QUOTE(LMPPATH/src/fix_external.h)

#include "qstexe.h"

using namespace LAMMPS_NS;

#define ANGSTROM_per_BOHR 0.529
#define EV_per_RYDBERG 13.6056923

void quest_callback(void *, bigint, int, int *, double **, double **);

struct Info {
  int me;
  Memory *memory;
  LAMMPS *lmp;
  char *quest_input;
};

/* ---------------------------------------------------------------------- */

int main(int narg, char **arg)
{
  int n;
  char str[128];

  // setup MPI

  MPI_Init(&narg,&arg);
  MPI_Comm comm = MPI_COMM_WORLD;

  int me,nprocs;
  MPI_Comm_rank(comm,&me);
  MPI_Comm_size(comm,&nprocs);

  Memory *memory = new Memory(comm);
  Error *error = new Error(comm);

  // command-line args

  if (narg != 4) error->all("Syntax: lmpqst Niter in.lammps in.quest");

  int niter = atoi(arg[1]);
  n = strlen(arg[2]) + 1;
  char *lammps_input = new char[n];
  strcpy(lammps_input,arg[2]);
  n = strlen(arg[3]) + 1;
  char *quest_input = new char[n];
  strcpy(quest_input,arg[3]);

  // instantiate LAMMPS

  LAMMPS *lmp = new LAMMPS(0,NULL,MPI_COMM_WORLD);

  // create simulation in LAMMPS from in.lammps

  lmp->input->file(lammps_input);

  // make info avaiable to callback function

  Info info;
  info.me = me;
  info.memory = memory;
  info.lmp = lmp;
  info.quest_input = quest_input;

  // set callback to Quest inside fix external
  // this could also be done thru Python, using a ctypes callback

  int ifix = lmp->modify->find_fix("2");
  FixExternal *fix = (FixExternal *) lmp->modify->fix[ifix];
  fix->set_callback(quest_callback,&info);

  // run LAMMPS for Niter
  // each time it needs forces, it will invoke quest_callback

  sprintf(str,"run %d",niter);
  lmp->input->one(str);

  // clean up

  delete lmp;

  delete memory;
  delete error;

  delete [] lammps_input;
  delete [] quest_input;

  MPI_Finalize();
}

/* ----------------------------------------------------------------------
   callback to Quest with atom IDs and coords from each proc
   invoke Quest to compute forces, load them into f for LAMMPS to use
   f can be NULL if proc owns no atoms
------------------------------------------------------------------------- */

void quest_callback(void *ptr, bigint ntimestep,
		    int nlocal, int *id, double **x, double **f)
{
  int i,j;
  char str[128];

  Info *info = (Info *) ptr;

  // boxlines = LAMMPS box size converted into Quest lattice vectors

  char **boxlines = NULL;
  if (info->me == 0) {
    boxlines = new char*[3];
    for (i = 0; i < 3; i++) boxlines[i] = new char[128];
  }

  double boxxlo = *((double *) lammps_extract_global(info->lmp,"boxxlo"));
  double boxxhi = *((double *) lammps_extract_global(info->lmp,"boxxhi"));
  double boxylo = *((double *) lammps_extract_global(info->lmp,"boxylo"));
  double boxyhi = *((double *) lammps_extract_global(info->lmp,"boxyhi"));
  double boxzlo = *((double *) lammps_extract_global(info->lmp,"boxzlo"));
  double boxzhi = *((double *) lammps_extract_global(info->lmp,"boxzhi"));
  double boxxy = *((double *) lammps_extract_global(info->lmp,"xy"));
  double boxxz = *((double *) lammps_extract_global(info->lmp,"xz"));
  double boxyz = *((double *) lammps_extract_global(info->lmp,"yz"));

  double xprd = (boxxhi-boxxlo)/ANGSTROM_per_BOHR;
  double yprd = (boxyhi-boxylo)/ANGSTROM_per_BOHR;
  double zprd = (boxzhi-boxzlo)/ANGSTROM_per_BOHR;
  double xy = boxxy/ANGSTROM_per_BOHR;
  double xz = boxxz/ANGSTROM_per_BOHR;
  double yz = boxyz/ANGSTROM_per_BOHR;

  if (info->me == 0) {
    sprintf(boxlines[0],"%g %g %g\n",xprd,0.0,0.0);
    sprintf(boxlines[1],"%g %g %g\n",xy,yprd,0.0);
    sprintf(boxlines[2],"%g %g %g\n",xz,yz,zprd);
  }

  // xlines = x for atoms on each proc converted to text lines
  // xlines is suitable for insertion into Quest input file
  // convert LAMMPS Angstroms to Quest bohr
  
  int natoms;
  MPI_Allreduce(&nlocal,&natoms,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);

  Many2One *lmp2qst = new Many2One(MPI_COMM_WORLD);
  lmp2qst->setup(nlocal,id,natoms);

  char **xlines = NULL;
  double **xquest = NULL;
  if (info->me == 0) {
    xquest = info->memory->create_2d_double_array(natoms,3,"lmpqst:xquest");
    xlines = new char*[natoms];
    for (i = 0; i < natoms; i++) xlines[i] = new char[128];
  }

  if (info->me == 0) lmp2qst->gather(&x[0][0],3,&xquest[0][0]);
  else lmp2qst->gather(&x[0][0],3,NULL);

  if (info->me == 0) {
    for (i = 0; i < natoms; i++) {
      xquest[i][0] /= ANGSTROM_per_BOHR;
      xquest[i][1] /= ANGSTROM_per_BOHR;
      xquest[i][2] /= ANGSTROM_per_BOHR;
    }
    for (i = 0; i < natoms; i++) {
      sprintf(xlines[i],"%d %d %g %g %g\n",i+1,1,
	      xquest[i][0],xquest[i][1],xquest[i][2]);
    }
  }

  // one-processor tasks:
  // whack all lcao.* files
  // cp quest_input to lcao.in
  // replace atom coords section of lcao.in with new atom coords
  // run Quest on one proc, save screen output to file
  // flines = atom forces extracted from Quest screen file
  // fquest = atom forces
  // convert Quest Ryd/bohr to LAMMPS eV/Angstrom
  
  char **flines = NULL;
  double **fquest = NULL;
  if (info->me == 0) {
    fquest = info->memory->create_2d_double_array(natoms,3,"lmpqst:fquest");
    flines = new char*[natoms];
    for (i = 0; i < natoms; i++) flines[i] = new char[128];
  }

  if (info->me == 0) {
    system("rm lcao.*");
    sprintf(str,"cp %s lcao.in",info->quest_input);
    system(str);
    sprintf(str,"cp %s lcao.x",QUOTE(QUEST));
    system(str);
    replace("lcao.in","primitive lattice vectors",3,boxlines);
    replace("lcao.in","atom, type, position vector",natoms,xlines);
    system("lcao.x > lcao.screen");
    extract("lcao.screen","atom       x force          "
	    "y force          z force",natoms,flines);
    
    int itmp;
    for (i = 0; i < natoms; i++)
      sscanf(flines[i],"%d %lg %lg %lg",&itmp,
	     &fquest[i][0],&fquest[i][1],&fquest[i][2]);

    for (i = 0; i < natoms; i++) {
      fquest[i][0] *= EV_per_RYDBERG / ANGSTROM_per_BOHR;
      fquest[i][1] *= EV_per_RYDBERG / ANGSTROM_per_BOHR;
      fquest[i][2] *= EV_per_RYDBERG / ANGSTROM_per_BOHR;
    }
  }

  // convert fquest on one proc into f for atoms on each proc

  One2Many *qst2lmp = new One2Many(MPI_COMM_WORLD);
  qst2lmp->setup(natoms,nlocal,id);
  double *fvec = NULL;
  if (f) fvec = &f[0][0];
  if (info->me == 0) qst2lmp->scatter(&fquest[0][0],3,fvec);
  else qst2lmp->scatter(NULL,3,fvec);

  // clean up
  // some data only exists on proc 0

  delete lmp2qst;
  delete qst2lmp;

  info->memory->destroy_2d_double_array(xquest);
  info->memory->destroy_2d_double_array(fquest);

  if (boxlines) {
    for (i = 0; i < 3; i++) delete [] boxlines[i];
    delete [] boxlines;
  }
  if (xlines) {
    for (i = 0; i < natoms; i++) delete [] xlines[i];
    delete [] xlines;
  }
  if (flines) {
    for (i = 0; i < natoms; i++) delete [] flines[i];
    delete [] flines;
  }
}

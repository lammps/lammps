// lmpspk = umbrella driver to couple LAMMPS + SPPARKS
//          for a strain-induced grain growth model

// Syntax: lmpspk Niter Ndelta Sfactor in.spparks
//         Niter = # of outer iterations
//         Ndelta = time to run MC in each iteration
//         Sfactor = multiplier on strain effect
//         in.spparks = SPPARKS input script

#include "mpi.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"

#include "lammps_data_write.h"
#include "many2many.h"
#include "memory.h"
#include "error.h"

#define QUOTE_(x) #x
#define QUOTE(x) QUOTE_(x)

#include "spkpath.h"
#include QUOTE(SPKPATH/src/spparks.h)
#include QUOTE(SPKPATH/src/library.h)
#include QUOTE(SPKPATH/src/input.h)

#include "lmppath.h"
#include QUOTE(LMPPATH/src/lammps.h)
#include QUOTE(LMPPATH/src/library.h)
#include QUOTE(LMPPATH/src/input.h)
#include QUOTE(LMPPATH/src/modify.h)
#include QUOTE(LMPPATH/src/compute.h)

using namespace SPPARKS_NS;
using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

int main(int narg, char **arg)
{
  int i,n;
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

  if (narg != 5) 
    error->all("Syntax: lmpspk Niter Ndelta Sfactor in.spparks");

  int niter = atoi(arg[1]);
  double delta = atof(arg[2]);
  double sfactor = atof(arg[3]);
  n = strlen(arg[4]) + 1;
  char *spparks_input = new char[n];
  strcpy(spparks_input,arg[4]);

  // instantiate LAMMPS and SPPARKS

  SPPARKS *spk = new SPPARKS(0,NULL,MPI_COMM_WORLD);
  LAMMPS *lmp = new LAMMPS(0,NULL,MPI_COMM_WORLD);

  // create simulation in SPPARKS from in.spparks

  spk->input->file(spparks_input);

  // extract permanent info from SPPARKS

  int dimension,nglobal,nlocal_spparks,nspins;
  double boxxlo,boxxhi,boxylo,boxyhi,boxzlo,boxzhi;
  int *id_spparks,*spins;
  double **xyz;
  double *strain;

  dimension = *((int *) spparks_extract(spk,"dimension"));
  nglobal = *((int *) spparks_extract(spk,"nglobal"));
  nlocal_spparks = *((int *) spparks_extract(spk,"nlocal"));

  boxxlo = *((double *) spparks_extract(spk,"boxxlo"));
  boxxhi = *((double *) spparks_extract(spk,"boxxhi"));
  boxylo = *((double *) spparks_extract(spk,"boxylo"));
  boxyhi = *((double *) spparks_extract(spk,"boxyhi"));
  if (dimension == 3) {
    boxzlo = *((double *) spparks_extract(spk,"boxzlo"));
    boxzhi = *((double *) spparks_extract(spk,"boxzhi"));
  } else {
    boxzlo = -0.5;
    boxzhi = 0.5;
  }

  id_spparks = (int *) spparks_extract(spk,"id");
  spins = (int *) spparks_extract(spk,"site");
  xyz = (double **) spparks_extract(spk,"xyz");

  nspins = *((int *) spparks_extract(spk,"nspins"));
  strain = (double *) spparks_extract(spk,"strain");

  // write a LAMMPS input script using SPPARKS params

  if (me == 0) {
    FILE *fp = fopen("in.lammps","w");
    if (fp == NULL) error->one("Could not create LAMMPS input script");

    fprintf(fp,"units lj\n");
    sprintf(str,"dimension %d\n",dimension);
    fprintf(fp,str);
    fprintf(fp,"atom_style atomic\n\n");
    
    fprintf(fp,"read_data data.lammps\n");
    fprintf(fp,"mass * 1.0\n\n");
    
    fprintf(fp,"pair_style lj/cut 2.5\n");
    fprintf(fp,"pair_coeff * * 1.0 1.2\n");
    for (i = 0; i < nspins; i++) {
      sprintf(str,"pair_coeff %d %d 1.0 1.0\n",i+1,i+1);
      fprintf(fp,str);
    }
    fprintf(fp,"\n");

    fprintf(fp,"compute da all displace/atom\n\n");
    fprintf(fp,"dump 1 all atom 10 dump.md\n");
    fprintf(fp,"thermo 1\n");

    fclose(fp);
  }

  // write a LAMMPS data file using SPPARKS data

  LAMMPSDataWrite *lwd = new LAMMPSDataWrite(MPI_COMM_WORLD);
  lwd->file("data.lammps");
  lwd->header("%d atoms",nglobal);
  lwd->header("%d atom types",nspins);
  lwd->header("%g %g xlo xhi",boxxlo,boxxhi);
  lwd->header("%g %g ylo yhi",boxylo,boxyhi);
  lwd->header("%g %g zlo zhi",boxzlo,boxzhi);
  lwd->atoms(nlocal_spparks);
  lwd->atoms(id_spparks);
  lwd->atoms(spins);
  lwd->atoms(3,xyz);
  lwd->execute();
  delete lwd;

  // create simulation in LAMMPS from created input script

  lmp->input->file("in.lammps");

  // create transfer operators

  Many2Many *spk2lmp = new Many2Many(MPI_COMM_WORLD);
  Many2Many *lmp2spk = new Many2Many(MPI_COMM_WORLD);

  // timestep loop
  // run SPPARKS for delta time
  // use SPPARKS spins to reset LAMMPS atom types
  // perform LAMMPS minimization
  // use atom displacements to reset strain values in SPPARKS
  // realloc displace as necessary since nlocal_lammps may change
  // re-create both xfers every iteration since LAMMPS may migrate atoms

  int nmax = 0;
  double *displace = NULL;
  
  int nlocal_lammps;
  int *id_lammps,*type;
  double **displace_lammps;

  for (int iter = 0; iter < niter; iter++) {
    sprintf(str,"run %g",delta);
    spk->input->one(str);

    nlocal_lammps = *((int *) lammps_extract_global(lmp,"nlocal"));
    id_lammps = (int *) lammps_extract_atom(lmp,"id");
    type = (int *) lammps_extract_atom(lmp,"type");

    spk2lmp->setup(nlocal_spparks,id_spparks,nlocal_lammps,id_lammps);
    spk2lmp->exchange(spins,type);

    lmp->input->one("minimize 0.001 0.001 10 1000");

    nlocal_lammps = *((int *) lammps_extract_global(lmp,"nlocal"));
    id_lammps = (int *) lammps_extract_atom(lmp,"id");
    displace_lammps = (double **) lammps_extract_compute(lmp,"da",1,2);

    if (nlocal_lammps > nmax) {
      memory->sfree(displace);
      nmax = nlocal_lammps;
      displace = (double *) 
	memory->smalloc(nmax*sizeof(double),"lmpspk:displace");
    }

    for (i = 0; i < nlocal_lammps; i++)
      displace[i] = sfactor*displace_lammps[i][3];

    lmp2spk->setup(nlocal_lammps,id_lammps,nlocal_spparks,id_spparks);
    lmp2spk->exchange(displace,strain);
  }

  memory->sfree(displace);

  // clean up

  delete spk2lmp;
  delete lmp2spk;
  delete spk;
  delete lmp;

  delete [] spparks_input;
  delete memory;
  delete error;

  MPI_Finalize();
}

/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* c_driver = simple example of how an umbrella program
              can invoke LAMMPS as a library on some subset of procs
   Syntax: simpleC P in.lammps
           P = # of procs to run LAMMPS on
               must be <= # of procs the driver code itself runs on
           in.lammps = LAMMPS input script
   See README for compilation instructions */

#include <mpi.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* define so interface to lammps_open() is available,
   since we will run on split communicator */
#define LAMMPS_LIB_MPI 1

#include "liblammpsplugin.h"  /* this is the include for the plugin loader */

int main(int narg, char **arg)
{
  int n,me,nprocs;
  int nprocs_lammps, lammps;
  char line[1024];
  MPI_Comm comm_lammps;
  FILE *fp = NULL;
  liblammpsplugin_t *plugin = NULL;
  void *lmp = NULL;
  double *x = NULL;
  double *v = NULL;
  char *strtwo;
  char *cmds[2];
  int *type = NULL;

  /* setup MPI and various communicators
     driver runs on all procs in MPI_COMM_WORLD
     comm_lammps only has 1st P procs (could be all or any subset) */

  MPI_Init(&narg,&arg);

  if (narg != 4) {
    printf("Syntax: %s P in.lammps /path/to/liblammps.so\n", arg[0]);
    exit(1);
  }

  MPI_Comm_rank(MPI_COMM_WORLD,&me);
  MPI_Comm_size(MPI_COMM_WORLD,&nprocs);

  nprocs_lammps = atoi(arg[1]);
  if (nprocs_lammps > nprocs) {
    if (me == 0)
      printf("ERROR: LAMMPS cannot use more procs than available\n");
    MPI_Abort(MPI_COMM_WORLD,1);
  }

  if (me < nprocs_lammps) lammps = 1;
  else lammps = MPI_UNDEFINED;
  MPI_Comm_split(MPI_COMM_WORLD,lammps,0,&comm_lammps);

  /* open LAMMPS input script */

  if (me == 0) {
    fp = fopen(arg[2],"r");
    if (fp == NULL) {
      printf("ERROR: Could not open LAMMPS input script\n");
      MPI_Abort(MPI_COMM_WORLD,1);
    }
  }

  /* run the input script thru LAMMPS one line at a time until end-of-file
     driver proc 0 reads a line, Bcasts it to all procs
     (could just send it to proc 0 of comm_lammps and let it Bcast)
     all LAMMPS procs call lammps_command() on the line */

  if (lammps == 1) {
    plugin = liblammpsplugin_load(arg[3]);
    if (plugin == NULL) {
      if (me == 0) printf("ERROR: Could not load shared LAMMPS library\n");
      MPI_Abort(MPI_COMM_WORLD,1);
    }
    /* must match the plugin ABI version */
    if (plugin->abiversion != LAMMPSPLUGIN_ABI_VERSION) {
      if (me == 0) printf("ERROR: Plugin abi version does not match: %d vs %d\n",
                          plugin->abiversion, LAMMPSPLUGIN_ABI_VERSION);
      MPI_Abort(MPI_COMM_WORLD,1);
    }
  }
  if (lammps == 1) {
    if (plugin->open == NULL) {
      printf("ERROR: liblammpsmpi.c must be compiled with -DLAMMPS_LIB_MPI=1 for this program\n");
      MPI_Abort(MPI_COMM_WORLD,2);
    }
    lmp = plugin->open(0,NULL,comm_lammps,NULL);
  }

  while (1) {
    if (me == 0) {
      if (fgets(line,1024,fp) == NULL) n = 0;
      else n = strlen(line) + 1;
      if (n == 0) fclose(fp);
    }
    MPI_Bcast(&n,1,MPI_INT,0,MPI_COMM_WORLD);
    if (n == 0) break;
    MPI_Bcast(line,n,MPI_CHAR,0,MPI_COMM_WORLD);
    if (lammps == 1) plugin->command(lmp,line);
  }

  /* run 10 more steps
     get coords from LAMMPS
     change coords of 1st atom
     put coords back into LAMMPS
     run a single step with changed coords */

  if (lammps == 1) {
    plugin->command(lmp,"run 10");

    int natoms = plugin->get_natoms(lmp);
    x = (double *) malloc(3*natoms*sizeof(double));
    plugin->gather_atoms(lmp,(char *)"x",1,3,x);
    v = (double *) malloc(3*natoms*sizeof(double));
    plugin->gather_atoms(lmp,(char *)"v",1,3,v);
    double epsilon = 0.1;
    x[0] += epsilon;
    plugin->scatter_atoms(lmp,(char *)"x",1,3,x);

    plugin->command(lmp,"run 1");
  }

  /* extract force on single atom two different ways */

  if (lammps == 1) {
    double **f = (double **) plugin->extract_atom(lmp,"f");
    printf("Force on 1 atom via extract_atom: %g\n",f[0][0]);

    double *fx = (double *) plugin->extract_variable(lmp,"fx",(char *)"all");
    printf("Force on 1 atom via extract_variable: %g\n",fx[0]);
    plugin->free(fx);
  }

  /* use commands_string() and commands_list() to invoke more commands */

  strtwo = (char *)"run 10\nrun 20";
  if (lammps == 1) plugin->commands_string(lmp,strtwo);

  cmds[0] = (char *)"run 10";
  cmds[1] = (char *)"run 20";
  if (lammps == 1) plugin->commands_list(lmp,2,(const char **)cmds);

  /* delete all atoms
     create_atoms() to create new ones with old coords, vels
     initial thermo should be same as step 20 */

  if (lammps == 1) {
    int natoms = plugin->get_natoms(lmp);
    type = (int *) malloc(natoms*sizeof(double));
    int i;
    for (i = 0; i < natoms; i++) type[i] = 1;

    plugin->command(lmp,"delete_atoms group all");
    plugin->create_atoms(lmp,natoms,NULL,type,x,v,NULL,0);
    plugin->command(lmp,"run 10");
  }

  if (x) free(x);
  if (v) free(v);
  if (type) free(type);

  /* close down LAMMPS */

  if (lammps == 1) {
    plugin->close(lmp);
    MPI_Barrier(comm_lammps);
    MPI_Comm_free(&comm_lammps);
    liblammpsplugin_release(plugin);
  }

  /* close down MPI */

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
}

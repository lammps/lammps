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
#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "neb.h"
#include "universe.h"
#include "atom.h"
#include "update.h"
#include "domain.h"
#include "min.h"
#include "modify.h"
#include "fix.h"
#include "fix_neb.h"
#include "output.h"
#include "thermo.h"
#include "finish.h"
#include "timer.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

#define CHUNK 1000
#define MAXLINE 256

/* ---------------------------------------------------------------------- */

NEB::NEB(LAMMPS *lmp) : Pointers(lmp) {}

/* ---------------------------------------------------------------------- */

NEB::~NEB()
{
  MPI_Comm_free(&roots);
  memory->destroy_2d_double_array(all);
  delete [] rdist;
}

/* ----------------------------------------------------------------------
   perform NEB on multiple replicas
------------------------------------------------------------------------- */

void NEB::command(int narg, char **arg)
{
  if (domain->box_exist == 0) 
    error->all("NEB command before simulation box is defined");

  if (narg != 5) error->universe_all("Illegal NEB command");
  
  double ftol = atof(arg[0]);
  int n1steps = atoi(arg[1]);
  int n2steps = atoi(arg[2]);
  int nevery = atoi(arg[3]);
  char *infile = arg[4];

  // error checks

  if (ftol < 0.0) error->all("Illegal NEB command");
  if (nevery == 0) error->universe_all("Illegal NEB command");
  if (n1steps % nevery || n2steps % nevery)
    error->universe_all("Illegal NEB command");

  // replica info

  nreplica = universe->nworlds;
  ireplica = universe->iworld;
  me_universe = universe->me;
  uworld = universe->uworld;
  MPI_Comm_rank(world,&me);

  // create MPI communicator for root proc from each world

  int color;
  if (me == 0) color = 0;
  else color = 1;
  MPI_Comm_split(uworld,color,0,&roots);

  // error checks

  if (nreplica == 1) error->all("Cannot use NEB with a single replica");
  if (nreplica != universe->nprocs)
    error->all("Can only use NEB with 1-processor replicas");
  if (atom->sortfreq > 0)
    error->all("Cannot use NEB with atom_modify sort enabled");
  if (atom->map_style == 0) 
    error->all("Cannot use NEB unless atom map exists");

  int ineb,idamp;
  for (ineb = 0; ineb < modify->nfix; ineb++)
    if (strcmp(modify->fix[ineb]->style,"neb") == 0) break;
  if (ineb == modify->nfix) error->all("NEB requires use of fix neb");

  fneb = (FixNEB *) modify->fix[ineb];
  all = memory->create_2d_double_array(nreplica,3,"neb:all");
  rdist = new double[nreplica];

  // initialize LAMMPS

  update->whichflag = 2;
  update->etol = 0.0;
  update->ftol = ftol;
  update->multireplica = 1;

  lmp->init();

  if (update->minimize->searchflag)
    error->all("NEB requires damped dynamics minimizer");

  // read in file of final state atom coords and reset my coords

  readfile(infile);

  // setup regular NEB minimizaiton

  if (me_universe == 0 && universe->uscreen)
    fprintf(universe->uscreen,"Setting up regular NEB ...\n");
  
  update->beginstep = update->firststep = update->ntimestep;
  update->endstep = update->laststep = update->firststep + n1steps;
  update->nsteps = n1steps;
  update->max_eval = n1steps;

  update->minimize->setup();
  
  if (me_universe == 0) {
    if (universe->uscreen)
      fprintf(universe->uscreen,"Step RD1 PE1 RD2 PE2 ... RDN PEN\n");
    if (universe->ulogfile)
      fprintf(universe->ulogfile,"Step RD1 PE1 RD2 PE2 ... RDN PEN\n");
  }
  print_status();
  
  // perform regular NEB for n1steps or until replicas converge
  // retrieve PE values from fix NEB and print every nevery iterations
  // break induced if converged
  // damped dynamic min styles insure all replicas converge together

  int flag,flagall;

  timer->barrier_start(TIME_LOOP);
  
  while (update->minimize->niter < n1steps) {
    update->minimize->run(nevery);
    print_status();
    if (update->minimize->stop_condition) break;
  }
	
  timer->barrier_stop(TIME_LOOP);

  update->minimize->cleanup();

  Finish finish(lmp);
  finish.end(1);

  // switch fix NEB to climbing mode
  // top = replica that becomes hill climber

  double vmax = all[0][0];
  int top = 0;
  for (int m = 1; m < nreplica; m++)
    if (vmax < all[m][0]) {
      vmax = all[m][0];
      top = m;
    }

  // setup climbing NEB minimization
  // must reinitialize minimizer so it re-creates its fix MINIMIZE

  if (me_universe == 0 && universe->uscreen)
    fprintf(universe->uscreen,"Setting up climbing NEB ...\n");
  
  update->beginstep = update->firststep = update->ntimestep;
  update->endstep = update->laststep = update->firststep + n2steps;
  update->nsteps = n2steps;
  update->max_eval = n2steps;

  update->minimize->init();
  fneb->rclimber = top;
  update->minimize->setup();

  if (me_universe == 0) {
    if (universe->uscreen)
      fprintf(universe->uscreen,"Step RD1 PE1 RD2 PE2 ... RDN PEN\n");
    if (universe->ulogfile)
      fprintf(universe->ulogfile,"Step RD1 PE1 RD2 PE2 ... RDN PEN\n");
  }
  print_status();
  
  // perform climbing NEB for n2steps or until replicas converge
  // retrieve PE values from fix NEB and print every nevery iterations
  // break induced if converged
  // damped dynamic min styles insure all replicas converge together
  
  timer->barrier_start(TIME_LOOP);
  
  while (update->minimize->niter < n2steps) {
    update->minimize->run(nevery);
    print_status();
    if (update->minimize->stop_condition) break;
  }
	
  timer->barrier_stop(TIME_LOOP);

  update->minimize->cleanup();

  finish.end(1);

  update->whichflag = 0;
  update->multireplica = 0;
  update->firststep = update->laststep = 0;
  update->beginstep = update->endstep = 0;
}

/* ----------------------------------------------------------------------
   read target coordinates from file, store with appropriate atom
   adjust coords of each atom based on ireplica
   new coord = replica fraction between current and final state
------------------------------------------------------------------------- */

void NEB::readfile(char *file)
{
  if (me_universe == 0) {
    if (screen) fprintf(screen,"Reading NEB coordinate file %s ...\n",file);
    open(file);
  }

  double fraction = ireplica/(nreplica-1.0);

  double **x = atom->x;
  int *image = atom->image;
  int nlocal = atom->nlocal;

  char *buffer = new char[CHUNK*MAXLINE];
  char *ptr,*next,*bufptr;
  int i,m,nlines,tag;
  double xx,yy,zz,delx,dely,delz;

  int firstline = 1;
  int ncount = 0;
  int eof = 0;

  while (!eof) {
    if (me_universe == 0) {
      m = 0;
      for (nlines = 0; nlines < CHUNK; nlines++) {
	ptr = fgets(&buffer[m],MAXLINE,fp);
	if (ptr == NULL) break;
	m += strlen(&buffer[m]);
      }
      if (ptr == NULL) eof = 1;
      buffer[m++] = '\n';
    }

    MPI_Bcast(&eof,1,MPI_INT,0,uworld);
    MPI_Bcast(&nlines,1,MPI_INT,0,uworld);
    MPI_Bcast(&m,1,MPI_INT,0,uworld);
    MPI_Bcast(buffer,m,MPI_CHAR,0,uworld);

    bufptr = buffer;
    for (i = 0; i < nlines; i++) {
      next = strchr(bufptr,'\n');
      *next = '\0';

      if (firstline) {
	if (atom->count_words(bufptr) == 4) firstline = 0;
	else error->all("Incorrect format in NEB coordinate file");
      }

      sscanf(bufptr,"%d %lg %lg %lg",&tag,&xx,&yy,&zz);

      // adjust atom coord based on replica fraction
      // ignore image flags of final x
      // new x is displacement from old x
      // if final x is across periodic boundary:
      //   new x may be outside box
      //   will be remapped back into box when simulation starts
      //   its image flags will be adjusted appropriately

      m = atom->map(tag);
      if (m >= 0 && m < nlocal) {
	delx = xx - x[m][0];
	dely = yy - x[m][1];
	delz = zz - x[m][2];
	domain->minimum_image(delx,dely,delz);
	x[m][0] += fraction*delx;
	x[m][1] += fraction*dely;
	x[m][2] += fraction*delz;
	ncount++;
      }

      bufptr = next + 1;
    }
  }

  // clean up

  delete [] buffer;

  if (me_universe == 0) {
    if (compressed) pclose(fp);
    else fclose(fp);
  }
}

/* ----------------------------------------------------------------------
   universe proc 0 opens NEB data file
   test if gzipped
------------------------------------------------------------------------- */

void NEB::open(char *file)
{
  compressed = 0;
  char *suffix = file + strlen(file) - 3;
  if (suffix > file && strcmp(suffix,".gz") == 0) compressed = 1;
  if (!compressed) fp = fopen(file,"r");
  else {
#ifdef LAMMPS_GZIP
    char gunzip[128];
    sprintf(gunzip,"gunzip -c %s",file);
    fp = popen(gunzip,"r");
#else
    error->one("Cannot open gzipped file");
#endif
  }

  if (fp == NULL) {
    char str[128];
    sprintf(str,"Cannot open file %s",file);
    error->one(str);
  }
}

/* ----------------------------------------------------------------------
   query fix NEB for PE of each replica
   proc 0 prints current NEB status
------------------------------------------------------------------------- */

void NEB::print_status()
{
  double one[3];
  one[0] = fneb->veng;
  one[1] = fneb->plen;
  one[2] = fneb->nlen;
  if (output->thermo->normflag) one[0] /= atom->natoms;
  if (me == 0) MPI_Allgather(one,3,MPI_DOUBLE,&all[0][0],3,MPI_DOUBLE,roots);

  rdist[0] = 0.0;
  for (int i = 1; i < nreplica; i++)
    rdist[i] = rdist[i-1] + all[i][1];
  double endpt = rdist[nreplica-1] = rdist[nreplica-2] + all[nreplica-2][2];
  for (int i = 1; i < nreplica; i++)
    rdist[i] /= endpt;
  
  if (me_universe == 0) {
    if (universe->uscreen) {
      fprintf(universe->uscreen,"%d ",update->ntimestep);
      for (int i = 0; i < nreplica; i++) 
	fprintf(universe->uscreen,"%g %g ",rdist[i],all[i][0]);
      fprintf(universe->uscreen,"\n");
    }
    if (universe->ulogfile) {
      fprintf(universe->ulogfile,"%d ",update->ntimestep);
      for (int i = 0; i < nreplica; i++)
	fprintf(universe->ulogfile,"%g %g ",rdist[i],all[i][0]);
      fprintf(universe->ulogfile,"\n");
      fflush(universe->ulogfile);
    }
  }
}

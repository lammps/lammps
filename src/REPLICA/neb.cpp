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

// lmptype.h must be first b/c this file uses MAXBIGINT and includes mpi.h
// due to OpenMPI bug which sets INT64_MAX via its mpi.h
//   before lmptype.h can set flags to insure it is done correctly

#include "neb.h"
#include "lmptype.h"
#include <mpi.h>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include "universe.h"
#include "atom.h"
#include "update.h"
#include "domain.h"
#include "comm.h"
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
#include "force.h"
#include "math_const.h"

using namespace LAMMPS_NS;
using namespace MathConst;

#define MAXLINE 256
#define CHUNK 1024
#define ATTRIBUTE_PERLINE 4

/* ---------------------------------------------------------------------- */

NEB::NEB(LAMMPS *lmp) : Pointers(lmp) {}

/* ----------------------------------------------------------------------
   internal NEB constructor, called from TAD
------------------------------------------------------------------------- */

NEB::NEB(LAMMPS *lmp, double etol_in, double ftol_in, int n1steps_in,
         int n2steps_in, int nevery_in, double *buf_init, double *buf_final)
  : Pointers(lmp)
{
  double delx,dely,delz;

  etol = etol_in;
  ftol = ftol_in;
  n1steps = n1steps_in;
  n2steps = n2steps_in;
  nevery = nevery_in;

  // replica info

  nreplica = universe->nworlds;
  ireplica = universe->iworld;
  me_universe = universe->me;
  uworld = universe->uworld;
  MPI_Comm_rank(world,&me);

  // generate linear interpolate replica

  double fraction = ireplica/(nreplica-1.0);

  double **x = atom->x;
  int nlocal = atom->nlocal;

  int ii = 0;
  for (int i = 0; i < nlocal; i++) {
    delx = buf_final[ii] - buf_init[ii];
    dely = buf_final[ii+1] - buf_init[ii+1];
    delz = buf_final[ii+2] - buf_init[ii+2];
    domain->minimum_image(delx,dely,delz);
    x[i][0] = buf_init[ii] + fraction*delx;
    x[i][1] = buf_init[ii+1] + fraction*dely;
    x[i][2] = buf_init[ii+2] + fraction*delz;
    ii += 3;
  }
}

/* ---------------------------------------------------------------------- */

NEB::~NEB()
{
  MPI_Comm_free(&roots);
  memory->destroy(all);
  delete [] rdist;
}

/* ----------------------------------------------------------------------
   perform NEB on multiple replicas
------------------------------------------------------------------------- */

void NEB::command(int narg, char **arg)
{
  if (domain->box_exist == 0)
    error->all(FLERR,"NEB command before simulation box is defined");

  if (narg < 6) error->universe_all(FLERR,"Illegal NEB command");

  etol = force->numeric(FLERR,arg[0]);
  ftol = force->numeric(FLERR,arg[1]);
  n1steps = force->inumeric(FLERR,arg[2]);
  n2steps = force->inumeric(FLERR,arg[3]);
  nevery = force->inumeric(FLERR,arg[4]);

  // error checks

  if (etol < 0.0) error->all(FLERR,"Illegal NEB command");
  if (ftol < 0.0) error->all(FLERR,"Illegal NEB command");
  if (nevery <= 0) error->universe_all(FLERR,"Illegal NEB command");
  if (n1steps % nevery || n2steps % nevery)
    error->universe_all(FLERR,"Illegal NEB command");

  // replica info

  nreplica = universe->nworlds;
  ireplica = universe->iworld;
  me_universe = universe->me;
  uworld = universe->uworld;
  MPI_Comm_rank(world,&me);

  // error checks

  if (nreplica == 1) error->all(FLERR,"Cannot use NEB with a single replica");
  if (atom->map_style == 0)
    error->all(FLERR,"Cannot use NEB unless atom map exists");

  // process file-style setting to setup initial configs for all replicas

  if (strcmp(arg[5],"final") == 0) {
    if (narg != 7 && narg !=8) error->universe_all(FLERR,"Illegal NEB command");
    inpfile = arg[6];
    readfile(inpfile,0);
  } else if (strcmp(arg[5],"each") == 0) {
    if (narg != 7 && narg !=8) error->universe_all(FLERR,"Illegal NEB command");
    inpfile = arg[6];
    readfile(inpfile,1);
  } else if (strcmp(arg[5],"none") == 0) {
    if (narg != 6 && narg !=7) error->universe_all(FLERR,"Illegal NEB command");
  } else error->universe_all(FLERR,"Illegal NEB command");

  verbose=false;
  if (strcmp(arg[narg-1],"verbose") == 0) verbose=true;
  // run the NEB calculation

  run();
}

/* ----------------------------------------------------------------------
   run NEB on multiple replicas
------------------------------------------------------------------------- */

void NEB::run()
{
  // create MPI communicator for root proc from each world

  int color;
  if (me == 0) color = 0;
  else color = 1;
  MPI_Comm_split(uworld,color,0,&roots);

  int ineb;
  for (ineb = 0; ineb < modify->nfix; ineb++)
    if (strcmp(modify->fix[ineb]->style,"neb") == 0) break;
  if (ineb == modify->nfix) error->all(FLERR,"NEB requires use of fix neb");

  fneb = (FixNEB *) modify->fix[ineb];
  if (verbose) numall =7;
  else  numall = 4;
  memory->create(all,nreplica,numall,"neb:all");
  rdist = new double[nreplica];

  // initialize LAMMPS

  update->whichflag = 2;
  update->etol = etol;
  update->ftol = ftol;
  update->multireplica = 1;

  lmp->init();

  if (update->minimize->searchflag)
    error->all(FLERR,"NEB requires damped dynamics minimizer");

  // setup regular NEB minimization
  FILE *uscreen = universe->uscreen;
  FILE *ulogfile = universe->ulogfile;

  if (me_universe == 0 && uscreen)
    fprintf(uscreen,"Setting up regular NEB ...\n");

  update->beginstep = update->firststep = update->ntimestep;
  update->endstep = update->laststep = update->firststep + n1steps;
  update->nsteps = n1steps;
  update->max_eval = n1steps;
  if (update->laststep < 0)
    error->all(FLERR,"Too many timesteps for NEB");

  update->minimize->setup();

  if (me_universe == 0) {
    if (uscreen) {
      if (verbose) {
        fprintf(uscreen,"Step MaxReplicaForce MaxAtomForce "
                "GradV0 GradV1 GradVc EBF EBR RDT RD1 PE1 RD2 PE2 ... "
                "RDN PEN pathangle1 angletangrad1 anglegrad1 gradV1 "
                "ReplicaForce1 MaxAtomForce1 pathangle2 angletangrad2 "
                "... ReplicaForceN MaxAtomForceN\n");
      } else {
        fprintf(uscreen,"Step MaxReplicaForce MaxAtomForce "
                "GradV0 GradV1 GradVc EBF EBR RDT RD1 PE1 RD2 PE2 ... "
                "RDN PEN\n");
      }
    }

    if (ulogfile) {
      if (verbose) {
        fprintf(ulogfile,"Step MaxReplicaForce MaxAtomForce "
                "GradV0 GradV1 GradVc EBF EBR RDT RD1 PE1 RD2 PE2 ... "
                "RDN PEN pathangle1 angletangrad1 anglegrad1 gradV1 "
                "ReplicaForce1 MaxAtomForce1 pathangle2 angletangrad2 "
                "... ReplicaForceN MaxAtomForceN\n");
      } else {
        fprintf(ulogfile,"Step MaxReplicaForce MaxAtomForce "
                "GradV0 GradV1 GradVc EBF EBR RDT RD1 PE1 RD2 PE2 ... "
                "RDN PEN\n");
      }
    }
  }
  print_status();

  // perform regular NEB for n1steps or until replicas converge
  // retrieve PE values from fix NEB and print every nevery iterations
  // break out of while loop early if converged
  // damped dynamic min styles insure all replicas converge together

  timer->init();
  timer->barrier_start();

  while (update->minimize->niter < n1steps) {
    update->minimize->run(nevery);
    print_status();
    if (update->minimize->stop_condition) break;
  }

  timer->barrier_stop();

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

  if (me_universe == 0 && uscreen)
    fprintf(uscreen,"Setting up climbing ...\n");

  if (me_universe == 0) {
    if (uscreen)
      fprintf(uscreen,"Climbing replica = %d\n",top+1);
    if (ulogfile)
      fprintf(ulogfile,"Climbing replica = %d\n",top+1);
  }

  update->beginstep = update->firststep = update->ntimestep;
  update->endstep = update->laststep = update->firststep + n2steps;
  update->nsteps = n2steps;
  update->max_eval = n2steps;
  if (update->laststep < 0)
    error->all(FLERR,"Too many timesteps");

  update->minimize->init();
  fneb->rclimber = top;
  update->minimize->setup();

  if (me_universe == 0) {
    if (uscreen) {
      if (verbose) {
        fprintf(uscreen,"Step MaxReplicaForce MaxAtomForce "
                "GradV0 GradV1 GradVc EBF EBR RDT "
                "RD1 PE1 RD2 PE2 ... RDN PEN "
                "pathangle1 angletangrad1 anglegrad1 gradV1 "
                "ReplicaForce1 MaxAtomForce1 pathangle2 angletangrad2 "
                "... ReplicaForceN MaxAtomForceN\n");
      } else {
        fprintf(uscreen,"Step MaxReplicaForce MaxAtomForce "
                "GradV0 GradV1 GradVc "
                "EBF EBR RDT "
                "RD1 PE1 RD2 PE2 ... RDN PEN\n");
      }
    }
    if (ulogfile) {
      if (verbose) {
        fprintf(ulogfile,"Step MaxReplicaForce MaxAtomForce "
                "GradV0 GradV1 GradVc EBF EBR RDT "
                "RD1 PE1 RD2 PE2 ... RDN PEN "
                "pathangle1 angletangrad1 anglegrad1 gradV1 "
                "ReplicaForce1 MaxAtomForce1 pathangle2 angletangrad2 "
                "... ReplicaForceN MaxAtomForceN\n");
      } else {
        fprintf(ulogfile,"Step MaxReplicaForce MaxAtomForce "
                "GradV0 GradV1 GradVc "
                "EBF EBR RDT "
                "RD1 PE1 RD2 PE2 ... RDN PEN\n");
      }
    }
  }
  print_status();

  // perform climbing NEB for n2steps or until replicas converge
  // retrieve PE values from fix NEB and print every nevery iterations
  // break induced if converged
  // damped dynamic min styles insure all replicas converge together

  timer->init();
  timer->barrier_start();

  while (update->minimize->niter < n2steps) {
    update->minimize->run(nevery);
    print_status();
    if (update->minimize->stop_condition) break;
  }

  timer->barrier_stop();

  update->minimize->cleanup();

  finish.end(1);

  update->whichflag = 0;
  update->multireplica = 0;
  update->firststep = update->laststep = 0;
  update->beginstep = update->endstep = 0;
}

/* ----------------------------------------------------------------------
   read initial config atom coords from file
   flag = 0
   only first replica opens file and reads it
   first replica bcasts lines to all replicas
   final replica stores coords
   intermediate replicas interpolate from coords
   new coord = replica fraction between current and final state
   initial replica does nothing
   flag = 1
   each replica (except first) opens file and reads it
   each replica stores coords
   initial replica does nothing
------------------------------------------------------------------------- */

void NEB::readfile(char *file, int flag)
{
  int i,j,m,nchunk,eofflag,nlines;
  tagint tag;
  char *eof,*start,*next,*buf;
  char line[MAXLINE];
  double xx,yy,zz,delx,dely,delz;

  if (me_universe == 0 && screen)
    fprintf(screen,"Reading NEB coordinate file(s) ...\n");

  // flag = 0, universe root reads header of file, bcast to universe
  // flag = 1, each replica's root reads header of file, bcast to world
  //   but explicitly skip first replica

  if (flag == 0) {
    if (me_universe == 0) {
      open(file);
      while (1) {
        eof = fgets(line,MAXLINE,fp);
        if (eof == NULL) error->one(FLERR,"Unexpected end of neb file");
        start = &line[strspn(line," \t\n\v\f\r")];
        if (*start != '\0' && *start != '#') break;
      }
      sscanf(line,"%d",&nlines);
    }
    MPI_Bcast(&nlines,1,MPI_INT,0,uworld);

  } else {
    if (me == 0) {
      if (ireplica) {
        open(file);
        while (1) {
          eof = fgets(line,MAXLINE,fp);
          if (eof == NULL) error->one(FLERR,"Unexpected end of neb file");
          start = &line[strspn(line," \t\n\v\f\r")];
          if (*start != '\0' && *start != '#') break;
        }
        sscanf(line,"%d",&nlines);
      } else nlines = 0;
    }
    MPI_Bcast(&nlines,1,MPI_INT,0,world);
  }

  char *buffer = new char[CHUNK*MAXLINE];
  char **values = new char*[ATTRIBUTE_PERLINE];

  double fraction = ireplica/(nreplica-1.0);

  double **x = atom->x;
  int nlocal = atom->nlocal;

  // loop over chunks of lines read from file
  // two versions of read_lines_from_file() for world vs universe bcast
  // count # of atom coords changed so can check for invalid atom IDs in file

  int ncount = 0;

  int nread = 0;
  while (nread < nlines) {
    nchunk = MIN(nlines-nread,CHUNK);
    if (flag == 0)
      eofflag = comm->read_lines_from_file_universe(fp,nchunk,MAXLINE,buffer);
    else
      eofflag = comm->read_lines_from_file(fp,nchunk,MAXLINE,buffer);
    if (eofflag) error->all(FLERR,"Unexpected end of neb file");

    buf = buffer;
    next = strchr(buf,'\n');
    *next = '\0';
    int nwords = atom->count_words(buf);
    *next = '\n';

    if (nwords != ATTRIBUTE_PERLINE)
      error->all(FLERR,"Incorrect atom format in neb file");

    // loop over lines of atom coords
    // tokenize the line into values

    for (i = 0; i < nchunk; i++) {
      next = strchr(buf,'\n');

      values[0] = strtok(buf," \t\n\r\f");
      for (j = 1; j < nwords; j++)
        values[j] = strtok(NULL," \t\n\r\f");

      // adjust atom coord based on replica fraction
      // for flag = 0, interpolate for intermediate and final replicas
      // for flag = 1, replace existing coord with new coord
      // ignore image flags of final x
      // for interpolation:
      //   new x is displacement from old x via minimum image convention
      //   if final x is across periodic boundary:
      //     new x may be outside box
      //     will be remapped back into box when simulation starts
      //     its image flags will then be adjusted

      tag = ATOTAGINT(values[0]);
      m = atom->map(tag);
      if (m >= 0 && m < nlocal) {
        ncount++;
        xx = atof(values[1]);
        yy = atof(values[2]);
        zz = atof(values[3]);

        if (flag == 0) {
          delx = xx - x[m][0];
          dely = yy - x[m][1];
          delz = zz - x[m][2];
          domain->minimum_image(delx,dely,delz);
          x[m][0] += fraction*delx;
          x[m][1] += fraction*dely;
          x[m][2] += fraction*delz;
        } else {
          x[m][0] = xx;
          x[m][1] = yy;
          x[m][2] = zz;
        }
      }

      buf = next + 1;
    }

    nread += nchunk;
  }

  // check that all atom IDs in file were found by a proc

  if (flag == 0) {
    int ntotal;
    MPI_Allreduce(&ncount,&ntotal,1,MPI_INT,MPI_SUM,uworld);
    if (ntotal != nreplica*nlines)
      error->universe_all(FLERR,"Invalid atom IDs in neb file");
  } else {
    int ntotal;
    MPI_Allreduce(&ncount,&ntotal,1,MPI_INT,MPI_SUM,world);
    if (ntotal != nlines)
      error->all(FLERR,"Invalid atom IDs in neb file");
  }

  // clean up

  delete [] buffer;
  delete [] values;

  if (flag == 0) {
    if (me_universe == 0) {
      if (compressed) pclose(fp);
      else fclose(fp);
    }
  } else {
    if (me == 0 && ireplica) {
      if (compressed) pclose(fp);
      else fclose(fp);
    }
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
    snprintf(gunzip,128,"gzip -c -d %s",file);

#ifdef _WIN32
    fp = _popen(gunzip,"rb");
#else
    fp = popen(gunzip,"r");
#endif

#else
    error->one(FLERR,"Cannot open gzipped file");
#endif
  }

  if (fp == NULL) {
    char str[128];
    snprintf(str,128,"Cannot open file %s",file);
    error->one(FLERR,str);
  }
}

/* ----------------------------------------------------------------------
   query fix NEB for info on each replica
   universe proc 0 prints current NEB status
------------------------------------------------------------------------- */

void NEB::print_status()
{
  double fnorm2 = sqrt(update->minimize->fnorm_sqr());
  double fmaxreplica;
  MPI_Allreduce(&fnorm2,&fmaxreplica,1,MPI_DOUBLE,MPI_MAX,roots);
  double fnorminf = update->minimize->fnorm_inf();
  double fmaxatom;
  MPI_Allreduce(&fnorminf,&fmaxatom,1,MPI_DOUBLE,MPI_MAX,roots);

  if (verbose) {
    freplica = new double[nreplica];
    MPI_Allgather(&fnorm2,1,MPI_DOUBLE,&freplica[0],1,MPI_DOUBLE,roots);
    fmaxatomInRepl = new double[nreplica];
    MPI_Allgather(&fnorminf,1,MPI_DOUBLE,&fmaxatomInRepl[0],1,MPI_DOUBLE,roots);
  }

  double one[7];
  one[0] = fneb->veng;
  one[1] = fneb->plen;
  one[2] = fneb->nlen;
  one[3] = fneb->gradlen;

  if (verbose) {
    one[4] = fneb->dotpath;
    one[5] = fneb->dottangrad;
    one[6] = fneb->dotgrad;
  }

  if (output->thermo->normflag) one[0] /= atom->natoms;
  if (me == 0)
    MPI_Allgather(one,numall,MPI_DOUBLE,&all[0][0],numall,MPI_DOUBLE,roots);
  MPI_Bcast(&all[0][0],numall*nreplica,MPI_DOUBLE,0,world);

  rdist[0] = 0.0;
  for (int i = 1; i < nreplica; i++)
    rdist[i] = rdist[i-1] + all[i][1];
  double endpt = rdist[nreplica-1] = rdist[nreplica-2] + all[nreplica-2][2];
  for (int i = 1; i < nreplica; i++)
    rdist[i] /= endpt;

  // look up GradV for the initial, final, and climbing replicas
  // these are identical to fnorm2, but to be safe we
  // take them straight from fix_neb

  double gradvnorm0, gradvnorm1, gradvnormc;

  int irep;
  irep = 0;
  gradvnorm0 = all[irep][3];
  irep = nreplica-1;
  gradvnorm1 = all[irep][3];
  irep = fneb->rclimber;
  if (irep > -1) {
    gradvnormc = all[irep][3];
    ebf = all[irep][0]-all[0][0];
    ebr = all[irep][0]-all[nreplica-1][0];
  } else {
    double vmax = all[0][0];
    int top = 0;
    for (int m = 1; m < nreplica; m++)
      if (vmax < all[m][0]) {
        vmax = all[m][0];
        top = m;
      }
    irep = top;
    gradvnormc = all[irep][3];
    ebf = all[irep][0]-all[0][0];
    ebr = all[irep][0]-all[nreplica-1][0];
  }

  if (me_universe == 0) {
    const double todeg=180.0/MY_PI;
    FILE *uscreen = universe->uscreen;
    FILE *ulogfile = universe->ulogfile;
    if (uscreen) {
      fprintf(uscreen,BIGINT_FORMAT " %12.8g %12.8g ",
              update->ntimestep,fmaxreplica,fmaxatom);
      fprintf(uscreen,"%12.8g %12.8g %12.8g ",
              gradvnorm0,gradvnorm1,gradvnormc);
      fprintf(uscreen,"%12.8g %12.8g %12.8g ",ebf,ebr,endpt);
      for (int i = 0; i < nreplica; i++)
        fprintf(uscreen,"%12.8g %12.8g ",rdist[i],all[i][0]);
      if (verbose) {
        fprintf(uscreen,"%12.5g %12.5g %12.5g %12.5g %12.5g %12.5g",
                NAN,180-acos(all[0][5])*todeg,180-acos(all[0][6])*todeg,
                all[0][3],freplica[0],fmaxatomInRepl[0]);
        for (int i = 1; i < nreplica-1; i++)
          fprintf(uscreen,"%12.5g %12.5g %12.5g %12.5g %12.5g %12.5g",
                  180-acos(all[i][4])*todeg,180-acos(all[i][5])*todeg,
                  180-acos(all[i][6])*todeg,all[i][3],freplica[i],
                  fmaxatomInRepl[i]);
        fprintf(uscreen,"%12.5g %12.5g %12.5g %12.5g %12.5g %12.5g",
                NAN,180-acos(all[nreplica-1][5])*todeg,NAN,all[nreplica-1][3],
                freplica[nreplica-1],fmaxatomInRepl[nreplica-1]);
      }
      fprintf(uscreen,"\n");
    }

    if (ulogfile) {
      fprintf(ulogfile,BIGINT_FORMAT " %12.8g %12.8g ",
              update->ntimestep,fmaxreplica,fmaxatom);
      fprintf(ulogfile,"%12.8g %12.8g %12.8g ",
              gradvnorm0,gradvnorm1,gradvnormc);
      fprintf(ulogfile,"%12.8g %12.8g %12.8g ",ebf,ebr,endpt);
      for (int i = 0; i < nreplica; i++)
        fprintf(ulogfile,"%12.8g %12.8g ",rdist[i],all[i][0]);
      if (verbose) {
        fprintf(ulogfile,"%12.5g %12.5g %12.5g %12.5g %12.5g %12.5g",
                NAN,180-acos(all[0][5])*todeg,180-acos(all[0][6])*todeg,
                all[0][3],freplica[0],fmaxatomInRepl[0]);
        for (int i = 1; i < nreplica-1; i++)
          fprintf(ulogfile,"%12.5g %12.5g %12.5g %12.5g %12.5g %12.5g",
                  180-acos(all[i][4])*todeg,180-acos(all[i][5])*todeg,
                  180-acos(all[i][6])*todeg,all[i][3],freplica[i],
                  fmaxatomInRepl[i]);
        fprintf(ulogfile,"%12.5g %12.5g %12.5g %12.5g %12.5g %12.5g",
                NAN,180-acos(all[nreplica-1][5])*todeg,NAN,all[nreplica-1][3],
                freplica[nreplica-1],fmaxatomInRepl[nreplica-1]);
      }
      fprintf(ulogfile,"\n");
      fflush(ulogfile);
    }
  }
}

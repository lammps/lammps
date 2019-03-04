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
// before lmptype.h can set flags to insure it is done correctly

#include "lmptype.h"
#include <mpi.h>
#include <cmath>
#include <cstdlib>
#include <cstring>
//#include "neb.h"
// test spin
#include "neb_spin.h"
#include "compute.h"
#include "force.h"

#include "universe.h"
#include "atom.h"
#include "update.h"
#include "domain.h"
#include "comm.h"
#include "min.h"
#include "modify.h"
#include "fix.h"
//#include "fix_neb.h"
// test spin
#include "fix_neb_spin.h"

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
//#define ATTRIBUTE_PERLINE 4
// 8 attributes: tag, spin norm, position (3), spin direction (3)
#define ATTRIBUTE_PERLINE 8

/* ---------------------------------------------------------------------- */

NEB_spin::NEB_spin(LAMMPS *lmp) : Pointers(lmp) {}

/* ----------------------------------------------------------------------
   internal NEB_spin constructor, called from TAD
------------------------------------------------------------------------- */

NEB_spin::NEB_spin(LAMMPS *lmp, double etol_in, double ftol_in, int n1steps_in,
         int n2steps_in, int nevery_in, double *buf_init, double *buf_final)
  : Pointers(lmp)
{
  //double delx,dely,delz;
  double delspx,delspy,delspz;

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
  // spin quantitites
  double **sp = atom->sp;
  int nlocal = atom->nlocal;

  int ii = 0;
  double spinit[3],spfinal[3];
  for (int i = 0; i < nlocal; i++) {

    spinit[0] = buf_init[ii];
    spinit[1] = buf_init[ii+1];
    spinit[2] = buf_init[ii+2];
    spfinal[0] = buf_final[ii];
    spfinal[1] = buf_final[ii+1];
    spfinal[2] = buf_final[ii+2];

    initial_rotation(spinit,spfinal,fraction);

    sp[i][0] = spfinal[0];
    sp[i][1] = spfinal[1];
    sp[i][2] = spfinal[2];

    //delx = buf_final[ii] - buf_init[ii];
    //dely = buf_final[ii+1] - buf_init[ii+1];
    //delz = buf_final[ii+2] - buf_init[ii+2];
    
    // adjust distance if pbc
    // not implemented yet
    //domain->minimum_image(delx,dely,delz);
   
    // need to define a procedure for circular initialization

    //x[i][0] = buf_init[ii] + fraction*delx;
    //x[i][1] = buf_init[ii+1] + fraction*dely;
    //x[i][2] = buf_init[ii+2] + fraction*delz;
    ii += 3;
  }
}

/* ---------------------------------------------------------------------- */

NEB_spin::~NEB_spin()
{
  MPI_Comm_free(&roots);
  memory->destroy(all);
  delete [] rdist;
}

/* ----------------------------------------------------------------------
   perform NEB_spin on multiple replicas
------------------------------------------------------------------------- */

void NEB_spin::command(int narg, char **arg)
{

  //printf("test 1 \n");

  // test 1
  double **sp1 = atom->sp;
  //printf("test 1 atom: i=%d,%g,%g,%g \n",1,sp1[1][0],sp1[1][1],sp1[1][2]);
  //error->all(FLERR,"end neb_spin test");


  if (domain->box_exist == 0)
    error->all(FLERR,"NEB_spin command before simulation box is defined");

  if (narg < 6) error->universe_all(FLERR,"Illegal NEB_spin command");

  etol = force->numeric(FLERR,arg[0]);
  ftol = force->numeric(FLERR,arg[1]);
  n1steps = force->inumeric(FLERR,arg[2]);
  n2steps = force->inumeric(FLERR,arg[3]);
  nevery = force->inumeric(FLERR,arg[4]);

  // error checks

  if (etol < 0.0) error->all(FLERR,"Illegal NEB_spin command");
  if (ftol < 0.0) error->all(FLERR,"Illegal NEB_spin command");
  if (nevery <= 0) error->universe_all(FLERR,"Illegal NEB_spin command");
  if (n1steps % nevery || n2steps % nevery)
    error->universe_all(FLERR,"Illegal NEB_spin command");

  // replica info

  nreplica = universe->nworlds;
  ireplica = universe->iworld;
  me_universe = universe->me;
  uworld = universe->uworld;
  MPI_Comm_rank(world,&me);

  // check metal units and spin atom/style 

  if (!atom->sp_flag)
    error->all(FLERR,"neb/spin requires atom/spin style");
  if (strcmp(update->unit_style,"metal") != 0)
    error->all(FLERR,"neb/spin simulation requires metal unit style");

  // error checks

  if (nreplica == 1) error->all(FLERR,"Cannot use NEB_spin with a single replica");
  if (atom->map_style == 0)
    error->all(FLERR,"Cannot use NEB_spin unless atom map exists");

  // process file-style setting to setup initial configs for all replicas

  // check what options are available
  if (strcmp(arg[5],"final") == 0) {
    if (narg != 7 && narg !=8) error->universe_all(FLERR,"Illegal NEB_spin command");
    infile = arg[6];
    readfile(infile,0);
  } else if (strcmp(arg[5],"each") == 0) {
    if (narg != 7 && narg !=8) error->universe_all(FLERR,"Illegal NEB_spin command");
    infile = arg[6];
    readfile(infile,1);
  } else if (strcmp(arg[5],"none") == 0) {
    if (narg != 6 && narg !=7) error->universe_all(FLERR,"Illegal NEB_spin command");
  } else error->universe_all(FLERR,"Illegal NEB_spin command");

  verbose=false;
  if (strcmp(arg[narg-1],"verbose") == 0) verbose=true;
  

  // test 1
  double **sp = atom->sp;
  //printf("test 2 atom: i=%d,%g,%g,%g \n",1,sp[1][0],sp[1][1],sp[1][2]);
  //error->all(FLERR,"end neb_spin test");

  // run the NEB_spin calculation

  run();
}

/* ----------------------------------------------------------------------
   run NEB_spin on multiple replicas
------------------------------------------------------------------------- */

void NEB_spin::run()
{
  // create MPI communicator for root proc from each world

  int color;
  if (me == 0) color = 0;
  else color = 1;
  MPI_Comm_split(uworld,color,0,&roots);

  // search for neb_spin fix, allocate it
  int ineb;
  for (ineb = 0; ineb < modify->nfix; ineb++)
    if (strcmp(modify->fix[ineb]->style,"neb/spin") == 0) break;
  if (ineb == modify->nfix) error->all(FLERR,"NEB_spin requires use of fix neb/spin");

  fneb = (FixNEB_spin *) modify->fix[ineb];
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

  // put flag to check gilbert damping procedure is set
  
  //if (update->minimize->searchflag)
  //  error->all(FLERR,"NEB_spin requires damped dynamics minimizer");

  // setup regular NEB_spin minimization
  FILE *uscreen = universe->uscreen;
  FILE *ulogfile = universe->ulogfile;
  
  //printf("test before run 1 \n");

  if (me_universe == 0 && uscreen)
    fprintf(uscreen,"Setting up regular NEB_spin ...\n");

  update->beginstep = update->firststep = update->ntimestep;
  update->endstep = update->laststep = update->firststep + n1steps;
  update->nsteps = n1steps;
  update->max_eval = n1steps;
  if (update->laststep < 0)
    error->all(FLERR,"Too many timesteps for NEB_spin");
  
  //printf("test before run 2 \n");

  update->minimize->setup();
  
  //printf("test before run 3 \n");

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
  //printf("test before run 4 \n");
  print_status();
  //printf("test before run 5 \n");

  // perform regular NEB_spin for n1steps or until replicas converge
  // retrieve PE values from fix NEB_spin and print every nevery iterations
  // break out of while loop early if converged
  // damped dynamic min styles insure all replicas converge together

  timer->init();
  timer->barrier_start();

  // test import fix_nve scheme

  //printf("test 2 atom: i=%d,%g,%g,%g \n",1,sp[1][0],sp[1][1],sp[1][2]);
  //error->all(FLERR,"end neb_spin test");
  double dts;
  while (update->minimize->niter < n1steps) {
    //dts = evaluate_dt();
    //advance_spins(dts);
    //fneb->
    
    //dts = fneb->evaluate_dt();
    //fneb->advance_spins(dts);
   
    // no minimizer for spins
    update->minimize->run(nevery);


    // no minimizer for spins
    //update->minimize->run(nevery);
    //
    // evaluate dts
    // loop on spins, damped advance
    //
    print_status();
    if (update->minimize->stop_condition) break;
  }

  // test neb end
  //error->all(FLERR,"end neb_spin test");
  
  timer->barrier_stop();

  update->minimize->cleanup();

  Finish finish(lmp);
  finish.end(1);

  // switch fix NEB_spin to climbing mode
  // top = replica that becomes hill climber

  double vmax = all[0][0];
  int top = 0;
  for (int m = 1; m < nreplica; m++)
    if (vmax < all[m][0]) {
      vmax = all[m][0];
      top = m;
    }

  // setup climbing NEB_spin minimization
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
  printf("test print 6.2 \n");
  update->minimize->setup();

  printf("test print 6.3 \n");
  
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

  // perform climbing NEB_spin for n2steps or until replicas converge
  // retrieve PE values from fix NEB_spin and print every nevery iterations
  // break induced if converged
  // damped dynamic min styles insure all replicas converge together

  timer->init();
  timer->barrier_start();

  while (update->minimize->niter < n2steps) {
    //dts = evaluate_dt();
    //advance_spins(dts);
    //dts = fneb->evaluate_dt();
    //fneb->advance_spins(dts);
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
   geodesic distance calculation (Vincenty's formula)
------------------------------------------------------------------------- */

//double NEB_spin::geodesic_distance2(double spi[3], double spj[3])
//{
//  double dist;
//  double crossx,crossy,crossz;
//  double dotx,doty,dotz;
//  double crosslen,dots;
//
//  crossx = spi[1]*spj[2]-spi[2]*spj[1];
//  crossy = spi[2]*spj[0]-spi[0]*spj[2];
//  crossz = spi[0]*spj[1]-spi[1]*spj[0];
//  crosslen = sqrt(crossx*crossx + crossy*crossy + crossz*crossz);
//  dotx = spi[0]*spj[0];
//  doty = spi[1]*spj[1];
//  dotz = spi[2]*spj[2];
//  dots = dotx+doty+dotz;
//
//  dist = atan2(crosslen,dots);
//
//  return dist;
//}

/* ----------------------------------------------------------------------
   evaluate max timestep
---------------------------------------------------------------------- */

//double NEB_spin::evaluate_dt()
//{
//  double dtmax;
//  double fmsq;
//  double fmaxsqone,fmaxsqloc,fmaxsqall;
//  int nlocal = atom->nlocal;
//  int *mask = atom->mask;
//  double **fm = atom->fm;
//
//  // finding max fm on this proc. 
//  
//  fmsq = fmaxsqone = fmaxsqloc = fmaxsqall = 0.0;
//  for (int i = 0; i < nlocal; i++)
//    if (mask[i] & groupbit) {
//      fmsq = fm[i][0]*fm[i][0]+fm[i][1]*fm[i][1]+fm[i][2]*fm[i][2];
//      fmaxsqone = MAX(fmaxsqone,fmsq);
//    }
//
//  // finding max fm on this replica
// 
//  fmaxsqloc = fmaxsqone; 
//  MPI_Allreduce(&fmaxsqone,&fmaxsqloc,1,MPI_DOUBLE,MPI_MAX,world); 
//
//  // finding max fm over all replicas, if necessary
//  // this communicator would be invalid for multiprocess replicas
//
//  if (update->multireplica == 1) {
//    fmaxsqall = fmaxsqloc;
//    MPI_Allreduce(&fmaxsqloc,&fmaxsqall,1,MPI_DOUBLE,MPI_MAX,universe->uworld);
//  }
//   
//  if (fmaxsqall < fmaxsqloc)
//    error->all(FLERR,"Incorrect fmaxall calc.");
//
//  // define max timestep
//  // dividing by 10 the inverse of max frequency
//
//  dtmax = MY_2PI/(10.0*sqrt(fmaxsqall));
//
//  return dtmax;
//}

/* ----------------------------------------------------------------------
   geometric damped advance os spins
---------------------------------------------------------------------- */

//void NEB_spin::advance_spins(double dts)
//{
//  //int j=0;
//  //int *sametag = atom->sametag;
//  int nlocal = atom->nlocal;
//  int *mask = atom->mask;
//  double **sp = atom->sp;
//  double **fm = atom->fm;
//  double tdampx,tdampy,tdampz;
//  double msq,scale,fm2,energy,dts2;
//  double alpha;
//  double spi[3],fmi[3];
//  double cp[3],g[3];
//
//  //cp[0] = cp[1] = cp[2] = 0.0;
//  //g[0] = g[1] = g[2] = 0.0;
//  dts2 = dts*dts;		
//
//  // fictitious Gilbert damping of 1
//  alpha = 1.0;
//
//  // loop on all spins on proc.
//
//  if (ireplica != nreplica-1 && ireplica != 0)
//    for (int i = 0; i < nlocal; i++)
//      if (mask[i] & groupbit) {
//
//	spi[0] = sp[i][0];
//	spi[1] = sp[i][1];
//	spi[2] = sp[i][2];
//
//	fmi[0] = fm[i][0];
//	fmi[1] = fm[i][1];
//	fmi[2] = fm[i][2];
//
//        // calc. damping torque
//      
//        tdampx = -alpha*(fmi[1]*spi[2] - fmi[2]*spi[1]);
//        tdampy = -alpha*(fmi[2]*spi[0] - fmi[0]*spi[2]);
//        tdampz = -alpha*(fmi[0]*spi[1] - fmi[1]*spi[0]);
//      
//        // apply advance algorithm (geometric, norm preserving)
//        
//        fm2 = (tdampx*tdampx+tdampy*tdampy+tdampz*tdampz);
//        energy = (sp[i][0]*tdampx)+(sp[i][1]*tdampy)+(sp[i][2]*tdampz);
//        
//        cp[0] = tdampy*sp[i][2]-tdampz*sp[i][1];
//        cp[1] = tdampz*sp[i][0]-tdampx*sp[i][2];
//        cp[2] = tdampx*sp[i][1]-tdampy*sp[i][0];
//      
//        g[0] = sp[i][0]+cp[0]*dts;
//        g[1] = sp[i][1]+cp[1]*dts;
//        g[2] = sp[i][2]+cp[2]*dts;
//      			
//        g[0] += (fm[i][0]*energy-0.5*sp[i][0]*fm2)*0.5*dts2;
//        g[1] += (fm[i][1]*energy-0.5*sp[i][1]*fm2)*0.5*dts2;
//        g[2] += (fm[i][2]*energy-0.5*sp[i][2]*fm2)*0.5*dts2;
//      			
//        g[0] /= (1+0.25*fm2*dts2);
//        g[1] /= (1+0.25*fm2*dts2);
//        g[2] /= (1+0.25*fm2*dts2);
//      
//        sp[i][0] = g[0];
//        sp[i][1] = g[1];
//        sp[i][2] = g[2];			
//      
//        // renormalization (check if necessary)
//      
//        msq = g[0]*g[0] + g[1]*g[1] + g[2]*g[2];
//        scale = 1.0/sqrt(msq);
//        sp[i][0] *= scale;
//        sp[i][1] *= scale;
//        sp[i][2] *= scale;
//      
//        // comm. sp[i] to atoms with same tag (for serial algo)
//      
//        // no need for simplecticity
//        //if (sector_flag == 0) {
//        //  if (sametag[i] >= 0) {
//        //    j = sametag[i];
//        //    while (j >= 0) {
//        //      sp[j][0] = sp[i][0];
//        //      sp[j][1] = sp[i][1];
//        //      sp[j][2] = sp[i][2];
//        //      j = sametag[j];
//        //    }
//        //  }
//        //}
//        //
//
//      }
//}

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

void NEB_spin::readfile(char *file, int flag)
{
  int i,j,m,nchunk,eofflag,nlines;
  tagint tag;
  char *eof,*start,*next,*buf;
  char line[MAXLINE];
  double xx,yy,zz,delx,dely,delz;
  // spin quantities
  double musp,spx,spy,spz;
    //,delx,dely,delz;

  if (me_universe == 0 && screen)
    fprintf(screen,"Reading NEB_spin coordinate file(s) ...\n");

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
  // spin quantities
  double **sp = atom->sp;
  double spinit[3],spfinal[3];
  int nlocal = atom->nlocal;

  // test 1.2
  //double **sp = atom->sp;
  //printf("test 1.2 atom: i=%d,%g,%g,%g \n",1,sp[1][0],sp[1][1],sp[1][2]);
  //error->all(FLERR,"end neb_spin test");


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
	musp = atof(values[1]);
	xx = atof(values[2]);
        yy = atof(values[3]);
        zz = atof(values[4]);
        spx = atof(values[5]);
        spy = atof(values[6]);
        spz = atof(values[7]);
	//xx = atof(values[1]);
        //yy = atof(values[2]);
        //zz = atof(values[3]);

        if (flag == 0) {
          
	  // here, function interp. spin states

	  //spinit[0] = x[m][0];
	  //spinit[1] = x[m][1];
	  //spinit[2] = x[m][2];
	  spinit[0] = sp[m][0];
	  spinit[1] = sp[m][1];
	  spinit[2] = sp[m][2];
	  spfinal[0] = spx;
	  spfinal[1] = spy;
	  spfinal[2] = spz;
          //domain->minimum_image(delx,dely,delz);

	  // test 
	  //printf("spinit: %g %g %g \n",spinit[0],spinit[1],spinit[2]);
	  //printf("spfinal bef: %g %g %g \n",spfinal[0],spfinal[1],spfinal[2]);

	  initial_rotation(spinit,spfinal,fraction);
	  
	  // test 
	  //printf("spfinal aft: %g %g %g \n",spfinal[0],spfinal[1],spfinal[2]);

	  sp[m][0] = spfinal[0];
	  sp[m][1] = spfinal[1];
	  sp[m][2] = spfinal[2];
	  sp[m][3] = musp;
	  //delx = xx - x[m][0];
          //dely = yy - x[m][1];
          //delz = zz - x[m][2];
          //x[m][0] += fraction*delx;
          //x[m][1] += fraction*dely;
          //x[m][2] += fraction*delz;
        } else {
          sp[m][3] = musp;
	  x[m][0] = xx;
          x[m][1] = yy;
          x[m][2] = zz;
	  sp[m][0] = spx;
	  sp[m][1] = spy; 
	  sp[m][2] = spz; 
        }
      }

      buf = next + 1;
    }

    nread += nchunk;
  }

  // test 1.3
  //double **sp = atom->sp;
  //printf("test 1.3 atom: i=%d,%g,%g,%g \n",1,sp[1][0],sp[1][1],sp[1][2]);
  //error->all(FLERR,"end neb_spin test");


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
   initial configuration of spin sploc using Rodrigues' formula
   interpolates between initial (spi) and final (stored in sploc) 
------------------------------------------------------------------------- */

void NEB_spin::initial_rotation(double *spi, double *sploc, double fraction)
{
  
  // implementing initial rotation using atan2
  //atan2(crosslen,dots);

  // this is not a sufficient routine, 
  // we need more accurate verifications


  // initial and final and inter ang. values 
  double itheta,iphi,ftheta,fphi,ktheta,kphi;
  double spix,spiy,spiz,spfx,spfy,spfz;
  double spkx,spky,spkz,iknorm;

  spix = spi[0];
  spiy = spi[1];
  spiz = spi[2];

  spfx = sploc[0];
  spfy = sploc[1];
  spfz = sploc[2];

  iphi = acos(spiz);
  itheta = acos(spix/sin(iphi));

  fphi = acos(spfz);
  ftheta = acos(spfx/sin(fphi));
 
  kphi = iphi + fraction*(fphi-iphi);
  ktheta = itheta + fraction*(ftheta-itheta);
 
  spkx = cos(ktheta)*sin(kphi);
  spky = sin(ktheta)*sin(kphi);
  spkz = cos(kphi);

  iknorm = spkx*spkx+spky*spky+spkz*spkz;

  spkx *= iknorm;
  spky *= iknorm;
  spkz *= iknorm;

  sploc[0] = spkx;
  sploc[1] = spky;
  sploc[2] = spkz;
  
  //double theta,spdot;
  //double inormdot,ispinorm;
  //double kix,kiy,kiz;
  //double kinorm, ikinorm;
  //double crossx,crossy,crossz;

  ////printf("inside rot, spi %g, spf %g \n",spi[0],sploc[0]);

  //spdot = spi[0]*sploc[0]+spi[1]*sploc[1]+spi[2]*sploc[2];
  //theta = fraction*acos(spdot);
 
  //printf("inside rot, theta %g \n",theta);

  //kix = spi[1]*sploc[2]-spi[2]*sploc[1];
  //kiy = spi[2]*sploc[0]-spi[0]*sploc[2];
  //kiz = spi[0]*sploc[1]-spi[1]*sploc[0];
  //
  ////printf("inside rot1.1, ki %g %g %g \n",kix,kiy,kiz);

  //inormdot = 1.0/sqrt(spdot);
  //kinorm = kix*kix+kiy*kiy+kiz*kiz;
  //if (kinorm == 0.0) {
  //  kix = 0.0; 
  //  kiy = 0.0; 
  //  kiz = 0.0;
  //} else {
  //  ikinorm = 1.0/kinorm; 
  //  kix *= ikinorm;
  //  kiy *= ikinorm;
  //  kiz *= ikinorm;
  //}
 
  ////printf("inside rot1.2, kin %g %g %g \n",kix,kiy,kiz);

  //crossx = kiy*spi[2]-kiz*spi[1];
  //crossy = kiz*spi[0]-kix*spi[2];
  //crossz = kix*spi[1]-kiy*spi[0];
  //
  ////printf("inside rot1.3, cross %g %g %g \n",crossx,crossy,crossz);

  //sploc[0] = spi[0]*cos(theta)+crossx*sin(theta);
  //sploc[1] = spi[1]*cos(theta)+crossy*sin(theta);
  //sploc[2] = spi[2]*cos(theta)+crossz*sin(theta);
  //
  ////printf("inside rot2, spf %g %g %g \n",sploc[0],sploc[1],sploc[2]);

  //ispinorm = 1.0/sqrt(sploc[0]*sploc[0]+sploc[1]*sploc[1]+sploc[2]*sploc[2]);

  //sploc[0] *= ispinorm;
  //sploc[1] *= ispinorm;
  //sploc[2] *= ispinorm;
  //printf("inside rot2, spf %g %g %g \n",sploc[0],sploc[1],sploc[2]);
}

/* ----------------------------------------------------------------------
   universe proc 0 opens NEB_spin data file
   test if gzipped
------------------------------------------------------------------------- */

void NEB_spin::open(char *file)
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
   query fix NEB_spin for info on each replica
   universe proc 0 prints current NEB_spin status
------------------------------------------------------------------------- */

void NEB_spin::print_status()
{
  
  //double fnorm2 = sqrt(update->minimize->fnorm_sqr());
  
  // test fmax spin
  int nlocal = atom->nlocal;
  double **fm = atom->fm;
  double fnorm2;
  for (int i = 0; i < nlocal; i++)
    fnorm2 += (fm[i][0]*fm[i][0]+fm[i][1]*fm[i][1]+fm[i][2]*fm[i][2]);
  //for (int i = 0; i < nlocal; i++)
  //  if (mask[i] & groupbit) {
  //    fnorm2 += (fm[i][0]*fm[i][0]+fm[i][1]*fm[i][1]+fm[i][2]*fm[i][2]);
  //  }
  
  double fmaxreplica;
  MPI_Allreduce(&fnorm2,&fmaxreplica,1,MPI_DOUBLE,MPI_MAX,roots);
  
  // no minimize->fnorm_inf for spins
  //double fnorminf = update->minimize->fnorm_inf();
  //double fmaxatom;
  //MPI_Allreduce(&fnorminf,&fmaxatom,1,MPI_DOUBLE,MPI_MAX,roots);
  double fnorminf = 0.0;
  double fmaxatom = 0.0;

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

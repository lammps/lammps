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

/* ----------------------------------------------------------------------
   Sources: Numerical Recipes frprmn routine
            "Conjugate Gradient Method Without the Agonizing Pain" by
            JR Shewchuk, http://www-2.cs.cmu.edu/~jrs/jrspapers.html#cg
------------------------------------------------------------------------- */

#include "math.h"
#include "string.h"
#include "mpi.h"
#include "min_cg.h"
#include "neighbor.h"
#include "domain.h"
#include "comm.h"
#include "atom.h"
#include "force.h"
#include "pair.h"
#include "bond.h"
#include "angle.h"
#include "dihedral.h"
#include "improper.h"
#include "kspace.h"
#include "output.h"
#include "thermo.h"
#include "update.h"
#include "modify.h"
#include "compute.h"
#include "fix_minimize.h"
#include "thermo.h"
#include "timer.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

#define MIN(A,B) ((A) < (B)) ? (A) : (B)
#define MAX(A,B) ((A) > (B)) ? (A) : (B)

#define EPS_ENERGY 1.0e-8
#define BACKTRACK_SLOPE 0.5
#define ALPHA_REDUCE 0.5

enum{FAIL,MAXITER,MAXEVAL,ETOL,FTOL};   // same as in other min classes

char *stopstrings[] = {"failed linesearch","max iterations",
		       "max force evaluations","energy tolerance",
		       "force tolerance"};

/* ---------------------------------------------------------------------- */

MinCG::MinCG(LAMMPS *lmp) : Min(lmp) {}

/* ---------------------------------------------------------------------- */

void MinCG::init()
{
  // create fix needed for storing atom-based gradient vectors
  // will delete it at end of run

  char **fixarg = new char*[3];
  fixarg[0] = (char *) "MINIMIZE";
  fixarg[1] = (char *) "all";
  fixarg[2] = (char *) "MINIMIZE";
  modify->add_fix(3,fixarg);
  delete [] fixarg;
  fix_minimize = (FixMinimize *) modify->fix[modify->nfix-1];

  // zero gradient vectors before first atom exchange

  setup_vectors();
  for (int i = 0; i < ndof; i++) h[i] = g[i] = 0.0;

  // virial_style:
  // 1 if computed explicitly by pair->compute via sum over pair interactions
  // 2 if computed implicitly by pair->virial_compute via sum over ghost atoms

  if (force->newton_pair) virial_style = 2;
  else virial_style = 1;

  // setup lists of computes for global and per-atom PE and pressure

  ev_setup();

  // set flags for what arrays to clear in force_clear()
  // need to clear torques if array exists

  torqueflag = 0;
  if (atom->torque) torqueflag = 1;

  // orthogonal vs triclinic simulation box

  triclinic = domain->triclinic;

  // reset reneighboring criteria if necessary

  neigh_every = neighbor->every;
  neigh_delay = neighbor->delay;
  neigh_dist_check = neighbor->dist_check;

  if (neigh_every != 1 || neigh_delay != 0 || neigh_dist_check != 1) {
    if (comm->me == 0) 
      error->warning("Resetting reneighboring criteria during minimization");
  }

  neighbor->every = 1;
  neighbor->delay = 0;
  neighbor->dist_check = 1;

  // set ptr to linemin function

  linemin = &MinCG::linemin_backtrack;
}

/* ----------------------------------------------------------------------
   perform minimization, with setup first
------------------------------------------------------------------------- */

void MinCG::run()
{
  double tmp,*f;

  // set initial force & energy
  // normalize energy if thermo PE does

  setup();
  setup_vectors();

  int id = modify->find_compute("thermo_pe");
  if (id < 0) error->all("Minimization could not find thermo_pe compute");
  pe_compute = modify->compute[id];

  ecurrent = pe_compute->compute_scalar();
  if (output->thermo->normflag) ecurrent /= atom->natoms;

  // stats for Finish to print
	
  einitial = ecurrent;

  f = NULL;
  if (ndof) f = atom->f[0];
  tmp = 0.0;
  for (int i = 0; i < ndof; i++) tmp += f[i]*f[i];
  MPI_Allreduce(&tmp,&fnorm2_init,1,MPI_DOUBLE,MPI_SUM,world);
  fnorm2_init = sqrt(fnorm2_init);

  tmp = 0.0;
  for (int i = 0; i < ndof; i++) tmp = MAX(fabs(f[i]),tmp);
  MPI_Allreduce(&tmp,&fnorminf_init,1,MPI_DOUBLE,MPI_MAX,world);

  // minimizer iterations

  timer->barrier_start(TIME_LOOP);
  int stop_condition = iterate(update->nsteps);
  stopstr = stopstrings[stop_condition];

  // account for early exit from iterate loop due to convergence
  // set niter/nsteps for Finish stats to print
  // set output->next values to this timestep
  // call eng_force to insure vflag is set when forces computed
  // output->write does final output for thermo, dump, restart files
  // add ntimestep to ALL computes that store invocation times
  //   since just hardwired call to thermo/dumps and they may not be ready

  if (niter < update->nsteps) {
    niter++;
    update->nsteps = niter;

    for (int idump = 0; idump < output->ndump; idump++)
      output->next_dump[idump] = update->ntimestep;
    output->next_dump_any = update->ntimestep;
    if (output->restart_every) output->next_restart = update->ntimestep;
    output->next_thermo = update->ntimestep;

    modify->addstep_compute_all(update->ntimestep);

    int ntmp;
    double *xtmp,*htmp,etmp;
    eng_force(&ntmp,&xtmp,&htmp,&etmp);
    output->write(update->ntimestep);
  }

  timer->barrier_stop(TIME_LOOP);

  // delete fix at end of run, so its atom arrays won't persist

  modify->delete_fix("MINIMIZE");

  // reset reneighboring criteria

  neighbor->every = neigh_every;
  neighbor->delay = neigh_delay;
  neighbor->dist_check = neigh_dist_check;

  // stats for Finish to print
	
  efinal = ecurrent;

  f = NULL;
  if (ndof) f = atom->f[0];
  tmp = 0.0;
  for (int i = 0; i < ndof; i++) tmp += f[i]*f[i];
  MPI_Allreduce(&tmp,&fnorm2_final,1,MPI_DOUBLE,MPI_SUM,world);
  fnorm2_final = sqrt(fnorm2_final);

  tmp = 0.0;
  for (int i = 0; i < ndof; i++) tmp = MAX(fabs(f[i]),tmp);
  MPI_Allreduce(&tmp,&fnorminf_final,1,MPI_DOUBLE,MPI_MAX,world);
}

/* ----------------------------------------------------------------------
   setup before run
------------------------------------------------------------------------- */

void MinCG::setup()
{
  if (comm->me == 0 && screen) fprintf(screen,"Setting up minimization ...\n");

  // setup domain, communication and neighboring
  // acquire ghosts
  // build neighbor lists
  // reset gradient vector ptrs

  if (triclinic) domain->x2lamda(atom->nlocal);
  domain->pbc();
  domain->reset_box();
  comm->setup();
  if (neighbor->style) neighbor->setup_bins();
  comm->exchange();
  comm->borders();
  if (triclinic) domain->lamda2x(atom->nlocal+atom->nghost);
  neighbor->build();
  neighbor->ncalls = 0;
  setup_vectors();

  // compute all forces

  ev_set(update->ntimestep);
  force_clear();
  
  if (force->pair) force->pair->compute(eflag,vflag);

  if (atom->molecular) {
    if (force->bond) force->bond->compute(eflag,vflag);
    if (force->angle) force->angle->compute(eflag,vflag);
    if (force->dihedral) force->dihedral->compute(eflag,vflag);
    if (force->improper) force->improper->compute(eflag,vflag);
  }

  if (force->kspace) {
    force->kspace->setup();
    force->kspace->compute(eflag,vflag);
  }

  if (force->newton) comm->reverse_communicate();

  modify->setup(vflag);
  output->setup(1);
}

/* ----------------------------------------------------------------------
   minimization via conjugate gradient iterations
   Polak-Ribiere formulation
------------------------------------------------------------------------- */

int MinCG::iterate(int n)
{
  int i,fail,ntimestep;
  double beta,gg,dot[2],dotall[2];

  double *x = NULL;
  double *f = NULL;

  if (ndof) f = atom->f[0];
  for (i = 0; i < ndof; i++) h[i] = g[i] = f[i];

  dot[0] = 0.0;
  for (i = 0; i < ndof; i++) dot[0] += f[i]*f[i];
  MPI_Allreduce(dot,&gg,1,MPI_DOUBLE,MPI_SUM,world);

  neval = 0;

  for (niter = 0; niter < n; niter++) {

    ntimestep = ++update->ntimestep;

    // line minimization along direction h from current atom->x

    eprevious = ecurrent;
    if (ndof) x = atom->x[0];
    fail = (this->*linemin)(ndof,x,h,ecurrent,dmax,alpha_final,neval);
    if (fail) return FAIL;

    // function evaluation criterion

    if (neval >= update->max_eval) return MAXEVAL;

    // energy tolerance criterion

    if (fabs(ecurrent-eprevious) < 
	update->etol * 0.5*(fabs(ecurrent) + fabs(eprevious) + EPS_ENERGY))
      return ETOL;

    // force tolerance criterion

    if (ndof) f = atom->f[0];
    dot[0] = dot[1] = 0.0;
    for (i = 0; i < ndof; i++) {
      dot[0] += f[i]*f[i];
      dot[1] += f[i]*g[i];
    }
    MPI_Allreduce(dot,dotall,2,MPI_DOUBLE,MPI_SUM,world);

    if (dotall[0] < update->ftol * update->ftol) return FTOL;

    // update h from new f = -Grad(x) and old g
    // beta = dotall[0]/gg would be Fletcher-Reeves CG
    
    beta = MAX(0.0,(dotall[0] - dotall[1])/gg);
    gg = dotall[0];

    for (i = 0; i < ndof; i++) {
      g[i] = f[i];
      h[i] = g[i] + beta*h[i];
    }

    // output for thermo, dump, restart files

    if (output->next == ntimestep) {
      timer->stamp();
      output->write(ntimestep);
      timer->stamp(TIME_OUTPUT);
    }
  }

  return MAXITER;
}

/* ----------------------------------------------------------------------
   set ndof and vector pointers after atoms have migrated
------------------------------------------------------------------------- */

void MinCG::setup_vectors()
{
  ndof = 3 * atom->nlocal;
  if (ndof) g = fix_minimize->gradient[0];
  else g = NULL;
  if (ndof) h = fix_minimize->searchdir[0];
  else h = NULL;
}

/* ----------------------------------------------------------------------
   evaluate potential energy and forces
   may migrate atoms
   new energy stored in ecurrent and returned (in case caller not in class)
   negative gradient will be stored in atom->f
------------------------------------------------------------------------- */

void MinCG::eng_force(int *pndof, double **px, double **ph, double *peng)
{
  // check for reneighboring
  // always communicate since minimizer moved atoms
  // if reneighbor, have to setup_vectors() since atoms migrated
  
  int nflag = neighbor->decide();

  if (nflag == 0) {
    timer->stamp();
    comm->communicate();
    timer->stamp(TIME_COMM);
  } else {
    if (triclinic) domain->x2lamda(atom->nlocal);
    domain->pbc();
    if (domain->box_change) {
      domain->reset_box();
      comm->setup();
      if (neighbor->style) neighbor->setup_bins();
    }
    timer->stamp();
    comm->exchange();
    comm->borders();
    if (triclinic) domain->lamda2x(atom->nlocal+atom->nghost);
    timer->stamp(TIME_COMM);
    neighbor->build();
    timer->stamp(TIME_NEIGHBOR);
    setup_vectors();
  }

  ev_set(update->ntimestep);
  force_clear();

  timer->stamp();

  if (force->pair) {
    force->pair->compute(eflag,vflag);
    timer->stamp(TIME_PAIR);
  }

  if (atom->molecular) {
    if (force->bond) force->bond->compute(eflag,vflag);
    if (force->angle) force->angle->compute(eflag,vflag);
    if (force->dihedral) force->dihedral->compute(eflag,vflag);
    if (force->improper) force->improper->compute(eflag,vflag);
    timer->stamp(TIME_BOND);
  }

  if (force->kspace) {
    force->kspace->compute(eflag,vflag);
    timer->stamp(TIME_KSPACE);
  }

  if (force->newton) {
    comm->reverse_communicate();
    timer->stamp(TIME_COMM);
  }

  // fixes that affect minimization

  if (modify->n_min_post_force) modify->min_post_force(vflag);

  // compute potential energy of system
  // normalize if thermo PE does

  ecurrent = pe_compute->compute_scalar();
  if (output->thermo->normflag) ecurrent /= atom->natoms;

  // return updated ptrs to caller since atoms may have migrated

  *pndof = ndof;
  if (ndof) *px = atom->x[0];
  else *px = NULL;
  *ph = h;
  *peng = ecurrent;
}

/* ----------------------------------------------------------------------
   clear force on own & ghost atoms
   setup and clear other arrays as needed
------------------------------------------------------------------------- */

void MinCG::force_clear()
{
  int i;

  // clear global force array
  // nall includes ghosts only if either newton flag is set

  int nall;
  if (force->newton) nall = atom->nlocal + atom->nghost;
  else nall = atom->nlocal;

  double **f = atom->f;
  for (i = 0; i < nall; i++) {
    f[i][0] = 0.0;
    f[i][1] = 0.0;
    f[i][2] = 0.0;
  }

  if (torqueflag) {
    double **torque = atom->torque;
    for (i = 0; i < nall; i++) {
      torque[i][0] = 0.0;
      torque[i][1] = 0.0;
      torque[i][2] = 0.0;
    }
  }
}

/* ----------------------------------------------------------------------
   line minimization methods
   find minimum-energy starting at x along dir direction
   input: n = # of degrees of freedom on this proc
          x = ptr to atom->x[0] as vector
	  dir = search direction as vector
	  eng = current energy at initial x
	  maxdist = max distance to move any atom coord
   output: return 0 if successful move, non-zero alpha alpha
           return 1 if failed, alpha = 0.0
           alpha = distance moved along dir to set x to minimun eng config
           caller has several quantities set via last call to eng_force()
	     must insure last call to eng_force() is consistent with returns
	       if fail, eng_force() of original x
	       if succeed, eng_force() at x + alpha*dir
             atom->x = coords at new configuration
	     atom->f = force (-Grad) is evaulated at new configuration
	     ecurrent = energy of new configuration
   NOTE: when call eng_force: n,x,dir,eng may change due to atom migration
	 updated values are returned by eng_force()
	 b/c of migration, linemin routines CANNOT store atom-based quantities
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   linemin: backtracking line search (Proc 3.1, p 41 in Nocedal and Wright)
   uses no gradient info, but should be very robust
   start at maxdist, backtrack until energy decrease is sufficient
------------------------------------------------------------------------- */

int MinCG::linemin_backtrack(int n, double *x, double *dir, double eng,
			     double maxdist, double &alpha, int &nfunc)
{
  int i;

  // stopping criterion, must be scaled by normflag

  double *f = NULL;
  if (n) f = atom->f[0];
  double fdotdirme = 0.0;
  for (i = 0; i < n; i++) fdotdirme += f[i]*dir[i];
  double fdotdirall;
  MPI_Allreduce(&fdotdirme,&fdotdirall,1,MPI_DOUBLE,MPI_SUM,world);
  if (output->thermo->normflag) fdotdirall /= atom->natoms;

  // alphamax = step that moves some atom coord by maxdist

  double fme = 0.0;
  for (i = 0; i < n; i++) fme = MAX(fme,fabs(dir[i]));
  double fmax;
  MPI_Allreduce(&fme,&fmax,1,MPI_DOUBLE,MPI_MAX,world);
  if (fmax == 0.0) return 1;

  double alphamax = maxdist/fmax;

  // reduce alpha by ALPHA_REDUCE until energy decrease is sufficient

  double eoriginal = eng;
  alpha = alphamax;
  double alpha_previous = 0.0;
  double delta;

  while (1) {
    delta = alpha - alpha_previous;
    for (i = 0; i < n; i++) x[i] += delta*dir[i];
    eng_force(&n,&x,&dir,&eng);
    nfunc++;

    if (eng <= eoriginal - BACKTRACK_SLOPE*alpha*fdotdirall) return 0;

    alpha_previous = alpha;
    alpha *= ALPHA_REDUCE;
    if (alpha == 0.0) {
      for (i = 0; i < n; i++) x[i] -= alpha_previous*dir[i];
      eng_force(&n,&x,&dir,&eng);
      return 1;
    }
  }    
}

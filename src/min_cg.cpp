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

// ALPHA_REDUCE = reduction ratio, should be in range [0.5,1)
// BACKTRACK_SLOPE, should be in range (0,0.5]
// QUADRATIC_TOL = tolerance on alpha0, should be in range [0.1,1)
// QUADRATIC_ENERGY = threshhold for switch to quadratic, small value > 0
// EPS_ENERGY = minimum normalization for energy tolerance
// EPS_QUAD = tolerance for quadratic projection

#define BACKTRACK_SLOPE 0.01
#define ALPHA_REDUCE 0.5
#define QUADRATIC_TOL 0.1
#define QUADRATIC_ENERGY 1.0e-4
#define EPS_ENERGY 1.0e-8
#define EPS_QUAD 1.0e-28

enum{FAIL,MAXITER,MAXEVAL,ETOL,FTOL};   // same as in other min classes

/* ---------------------------------------------------------------------- */

MinCG::MinCG(LAMMPS *lmp) : Min(lmp)
{
  fextra = gextra = hextra = NULL;
}

/* ---------------------------------------------------------------------- */

MinCG::~MinCG()
{
  delete [] fextra;
  delete [] gextra;
  delete [] hextra;
}

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
  // clear torques if array exists

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

  if (linestyle == 0) linemin = &MinCG::linemin_backtrack;
  else if (linestyle == 1) linemin = &MinCG::linemin_quadratic;
}

/* ----------------------------------------------------------------------
   perform minimization, with setup first
------------------------------------------------------------------------- */

void MinCG::run()
{
  int i;
  double tmp,*f;

  // possible stop conditions

  char *stopstrings[] = {"failed linesearch","max iterations",
			 "max force evaluations","energy tolerance",
			 "force tolerance"};

  // set initial force & energy
  // normalize energy if thermo PE does

  setup();
  setup_vectors();

  // setup any extra dof due to fixes
  // can't be done until now b/c update init() comes before modify init()

  delete [] fextra;
  delete [] gextra;
  delete [] hextra;

  nextra = modify->min_dof();
  if (nextra) {
    fextra = new double[nextra];
    gextra = new double[nextra];
    hextra = new double[nextra];
  }

  // compute potential energy of system
  // normalize if thermo PE does

  int id = modify->find_compute("thermo_pe");
  if (id < 0) error->all("Minimization could not find thermo_pe compute");
  pe_compute = modify->compute[id];

  ecurrent = pe_compute->compute_scalar();
  if (nextra) ecurrent += modify->min_energy(fextra);
  if (output->thermo->normflag) ecurrent /= atom->natoms;
  
  // stats for Finish to print
	
  einitial = ecurrent;

  f = NULL;
  if (ndof) f = atom->f[0];
  tmp = 0.0;
  for (i = 0; i < ndof; i++) tmp += f[i]*f[i];
  MPI_Allreduce(&tmp,&fnorm2_init,1,MPI_DOUBLE,MPI_SUM,world);
  if (nextra)
    for (i = 0; i < nextra; i++) fnorm2_init += fextra[i]*fextra[i];
  fnorm2_init = sqrt(fnorm2_init);

  tmp = 0.0;
  for (i = 0; i < ndof; i++) tmp = MAX(fabs(f[i]),tmp);
  MPI_Allreduce(&tmp,&fnorminf_init,1,MPI_DOUBLE,MPI_MAX,world);
  if (nextra)
    for (i = 0; i < nextra; i++)
      fnorminf_init = MAX(fabs(fextra[i]),fnorminf_init);

  // minimizer iterations

  timer->barrier_start(TIME_LOOP);
  int stop_condition = iterate(update->nsteps);
  stopstr = stopstrings[stop_condition];

  // account for early exit from iterate loop due to convergence
  // set niter/nsteps for Finish stats to print
  // set output->next values to this timestep
  // call eng_force to insure vflag is set when forces computed
  // output->write does final output for thermo, dump, restart files
  // add ntimestep to all computes that store invocation times
  //   since are hardwireing call to thermo/dumps and computes may not be ready

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
    double *xtmp,*htmp,*x0tmp,etmp;
    eng_force(&ntmp,&xtmp,&htmp,&x0tmp,&etmp);
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
  for (i = 0; i < ndof; i++) tmp += f[i]*f[i];
  MPI_Allreduce(&tmp,&fnorm2_final,1,MPI_DOUBLE,MPI_SUM,world);
  if (nextra)
    for (i = 0; i < nextra; i++)
      fnorm2_final += fextra[i]*fextra[i];
  fnorm2_final = sqrt(fnorm2_final);

  tmp = 0.0;
  for (i = 0; i < ndof; i++) tmp = MAX(fabs(f[i]),tmp);
  MPI_Allreduce(&tmp,&fnorminf_final,1,MPI_DOUBLE,MPI_MAX,world);
  if (nextra)
    for (i = 0; i < nextra; i++)
      fnorminf_final = MAX(fabs(fextra[i]),fnorminf_final);
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
------------------------------------------------------------------------- */

int MinCG::iterate(int n)
{
  int i,fail,ntimestep;
  double beta,gg,dot[2],dotall[2];

  double *x = NULL;
  double *f = NULL;

  if (ndof) f = atom->f[0];
  for (i = 0; i < ndof; i++) h[i] = g[i] = f[i];
  if (nextra)
    for (i = 0; i < nextra; i++)
      hextra[i] = gextra[i] = fextra[i];
  
  dot[0] = 0.0;
  for (i = 0; i < ndof; i++) dot[0] += f[i]*f[i];
  MPI_Allreduce(dot,&gg,1,MPI_DOUBLE,MPI_SUM,world);
  if (nextra)
    for (i = 0; i < nextra; i++) gg += fextra[i]*fextra[i];

  neval = 0;

  for (niter = 0; niter < n; niter++) {

    ntimestep = ++update->ntimestep;

    // line minimization along direction h from current atom->x

    eprevious = ecurrent;
    if (ndof) x = atom->x[0];
    fail = (this->*linemin)(ndof,x,h,x0,ecurrent,dmax,alpha_final,neval);
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
    if (nextra)
      for (i = 0; i < nextra; i++) {
	dotall[0] += fextra[i]*fextra[i];
	dotall[1] += fextra[i]*gextra[i];
      }

    if (dotall[0] < update->ftol * update->ftol) return FTOL;

    // update new search direction h from new f = -Grad(x) and old g
    // this is Polak-Ribieri formulation
    // beta = dotall[0]/gg would be Fletcher-Reeves

    beta = MAX(0.0,(dotall[0] - dotall[1])/gg);
    gg = dotall[0];

    for (i = 0; i < ndof; i++) {
      g[i] = f[i];
      h[i] = g[i] + beta*h[i];
    }
    if (nextra)
      for (i = 0; i < nextra; i++) {
	gextra[i] = fextra[i];
	hextra[i] = gextra[i] + beta*hextra[i];
      }

    // reinitialize CG if new search direction h is not downhill
    // also reinitialize CG every ndof iterations

    dot[0] = 0.0;
    for (i = 0; i < ndof; i++) dot[0] += g[i]*h[i];
    MPI_Allreduce(dot,dotall,1,MPI_DOUBLE,MPI_SUM,world);
    if (nextra)
      for (i = 0; i < nextra; i++)
	dotall[0] += gextra[i]*hextra[i];

    if (dotall[0] <= 0.0) {
      for (i = 0; i < ndof; i++) h[i] = g[i];
      if (nextra)
	for (i = 0; i < nextra; i++)
	  hextra[i] = gextra[i];
    }

    if (niter % (ndof+nextra) == 0) beta = 0.0;

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
  if (ndof) x0 = fix_minimize->x0[0];
  else x0 = NULL;
}

/* ----------------------------------------------------------------------
   evaluate potential energy and forces
   may migrate atoms
   new energy stored in ecurrent and returned (in case caller not in class)
   negative gradient will be stored in atom->f
------------------------------------------------------------------------- */

void MinCG::eng_force(int *pndof, double **px, double **ph, double **px0,
		      double *peng)
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
  if (nextra) ecurrent += modify->min_energy(fextra);
  if (output->thermo->normflag) ecurrent /= atom->natoms;

  // return updated ptrs to caller since atoms may have migrated

  *pndof = ndof;
  if (ndof) *px = atom->x[0];
  else *px = NULL;
  *ph = h;
  *px0 = x0;
  *peng = ecurrent;
}

/* ----------------------------------------------------------------------
   clear force on own & ghost atoms
   setup and clear other arrays as needed
------------------------------------------------------------------------- */

void MinCG::force_clear()
{
  // clear global force array
  // nall includes ghosts only if either newton flag is set

  int nall;
  if (force->newton) nall = atom->nlocal + atom->nghost;
  else nall = atom->nlocal;

  double **f = atom->f;
  for (int i = 0; i < nall; i++) {
    f[i][0] = 0.0;
    f[i][1] = 0.0;
    f[i][2] = 0.0;
  }

  if (torqueflag) {
    double **torque = atom->torque;
    for (int i = 0; i < nall; i++) {
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
          x0 = ptr to fix->x0[0] as vector, for storing initial coords
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
   NOTE: when call eng_force: n,x,dir,x0,eng may change due to atom migration
	 updated values are returned by eng_force()
	 b/c of migration, linemin routines CANNOT store atom-based quantities
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   linemin: backtracking line search (Proc 3.1, p 41 in Nocedal and Wright)
   uses no gradient info, but should be very robust
   start at maxdist, backtrack until energy decrease is sufficient
------------------------------------------------------------------------- */

int MinCG::linemin_backtrack(int n, double *x, double *dir,
			     double *x0, double eng, double maxdist,
			     double &alpha, int &nfunc)
{
  int i,m;
  double fdotdirall,fdotdirme,hmax,hme,eoriginal;
  double de_ideal,de;

  double *f = NULL;
  if (n) f = atom->f[0];

  // fdotdirall = projection of search dir along downhill gradient
  // if search direction is not downhill, exit with error

  fdotdirme = 0.0;
  for (i = 0; i < n; i++) fdotdirme += f[i]*dir[i];
  MPI_Allreduce(&fdotdirme,&fdotdirall,1,MPI_DOUBLE,MPI_SUM,world);
  if (nextra)
    for (i = 0; i < nextra; i++) fdotdirall += fextra[i]*hextra[i];
  if (output->thermo->normflag) fdotdirall /= atom->natoms;
  if (fdotdirall <= 0.0) return 1;

  // initial alpha = stepsize to change any atom coord by maxdist
  // limit alpha <= 1.0 else backtrack from huge value when forces are tiny
  // if all search dir components are already 0.0, exit with error

  hme = 0.0;
  for (i = 0; i < n; i++) hme = MAX(hme,fabs(dir[i]));
  MPI_Allreduce(&hme,&hmax,1,MPI_DOUBLE,MPI_MAX,world);
  if (nextra)
    for (i = 0; i < nextra; i++)
      hmax = MAX(hmax,fabs(hextra[i]));
  if (hmax == 0.0) return 1;
  alpha = MIN(1.0,maxdist/hmax);

  // store coords and other dof at start of linesearch

  for (i = 0; i < n; i++) x0[i] = x[i];
  if (nextra) modify->min_store();

  // eoriginal = energy at start of linesearch

  eng_force(&n,&x,&dir,&x0,&eng);
  nfunc++;
  eoriginal = eng;

  // backtrack with alpha until energy decrease is sufficient

  while (1) {
    for (i = 0; i < n; i++) x[i] = x0[i];
    if (nextra) modify->min_step(0.0,hextra);
    for (i = 0; i < n; i++) x[i] += alpha*dir[i];
    if (nextra) modify->min_step(alpha,hextra);
    eng_force(&n,&x,&dir,&x0,&eng);
    nfunc++;

    // if energy change is better than ideal, exit with success

    de_ideal = -BACKTRACK_SLOPE*alpha*fdotdirall;
    de = eng - eoriginal;
    if (de <= de_ideal) return 0;

    // reduce alpha

    alpha *= ALPHA_REDUCE;

    // backtracked all the way to 0.0
    // reset to starting point, exit with error

    if (alpha <= 0.0 || de_ideal >= -1.0e-20) {
      for (i = 0; i < n; i++) x[i] = x0[i];
      if (nextra) modify->min_step(0.0,hextra);
      eng_force(&n,&x,&dir,&x0,&eng);
      nfunc++;
      return 1;
    }
  }
}

/* ----------------------------------------------------------------------
   linemin: quadratic line search (adapted from Dennis and Schnabel)
   basic idea is to backtrack until change in energy is sufficiently small
     based on ENERGY_QUADRATIC, then use a quadratic approximation
     using forces at two alpha values to project to minimum
   use forces rather than energy change to do projection
   this is b/c the forces are going to zero and can become very small
     unlike energy differences which are the difference of two finite
     values and are thus limited by machine precision
   two changes that were critical to making this method work:
     a) limit maximum step to alpha <= 1
     b) ignore energy criterion if delE <= ENERGY_QUADRATIC
   several other ideas also seemed to help:
     c) making each step from starting point (alpha = 0), not previous alpha
     d) quadratic model based on forces, not energy
     e) exiting immediately if f.dir <= 0 (search direction not downhill)
        so that CG can restart
   a,c,e were also adopted for the backtracking linemin function
------------------------------------------------------------------------- */

int MinCG::linemin_quadratic(int n, double *x, double *dir,
			     double *x0, double eng, double maxdist,
			     double &alpha, int &nfunc)
{
  int i,m;
  double fdotdirall,fdotdirme,hmax,hme,alphamax,eoriginal;
  double de_ideal,de;
  double delfh,engprev,relerr,alphaprev,fhprev,ff,fh,alpha0;
  double dot[2],dotall[2];	
  double *f = atom->f[0];

  // fdotdirall = projection of search dir along downhill gradient
  // if search direction is not downhill, exit with error

  fdotdirme = 0.0;
  for (i = 0; i < n; i++) fdotdirme += f[i]*dir[i];
  MPI_Allreduce(&fdotdirme,&fdotdirall,1,MPI_DOUBLE,MPI_SUM,world);
  if (nextra)
    for (i = 0; i < nextra; i++) fdotdirall += fextra[i]*hextra[i];
  if (output->thermo->normflag) fdotdirall /= atom->natoms;
  if (fdotdirall <= 0.0) return 1;

  // initial alpha = stepsize to change any atom coord by maxdist
  // limit alpha <= 1.0 else backtrack from huge value when forces are tiny
  // if all search dir components are already 0.0, exit with error

  hme = 0.0;
  for (i = 0; i < n; i++) hme = MAX(hme,fabs(dir[i]));
  MPI_Allreduce(&hme,&hmax,1,MPI_DOUBLE,MPI_MAX,world);
  if (nextra)
    for (i = 0; i < nextra; i++)
      hmax = MAX(hmax,fabs(hextra[i]));
  if (hmax == 0.0) return 1;
  alpha = MIN(1.0,maxdist/hmax);

  // store coords and other dof at start of linesearch

  for (i = 0; i < n; i++) x0[i] = x[i];
  if (nextra) modify->min_store();

  // eoriginal = energy at start of linesearch

  eng_force(&n,&x,&dir,&x0,&eng);
  nfunc++;
  eoriginal = eng;

  // backtrack with alpha until energy decrease is sufficient
  // or until get to small energy change, then perform quadratic projection

  fhprev = fdotdirall;
  engprev = eoriginal;
  alphaprev = 0.0;

  while (1) {
    for (i = 0; i < n; i++) x[i] = x0[i];
    if (nextra) modify->min_step(0.0,hextra);
    for (i = 0; i < n; i++) x[i] += alpha*dir[i];
    if (nextra) modify->min_step(alpha,hextra);
    eng_force(&n,&x,&dir,&x0,&eng);
    nfunc++;

    // compute new fh, alpha, delfh

    dot[0] = dot[1] = 0.0;
    for (i = 0; i < ndof; i++) {
      dot[0] += f[i]*f[i];
      dot[1] += f[i]*dir[i];
    }
    MPI_Allreduce(dot,dotall,2,MPI_DOUBLE,MPI_SUM,world);
    if (nextra) {
      for (i = 0; i < nextra; i++) {
	dotall[0] += fextra[i]*fextra[i];
	dotall[1] += fextra[i]*hextra[i];
      }
    }
    ff = dotall[0];
    fh = dotall[1];
    if (output->thermo->normflag) {
      ff /= atom->natoms;
      fh /= atom->natoms;
    }

    delfh = fh - fhprev;

    // if fh or delfh is epsilon, reset to starting point, exit with error

    if (fabs(fh) < EPS_QUAD || fabs(delfh) < EPS_QUAD) {
      for (i = 0; i < n; i++) x[i] = x0[i];
      if (nextra) modify->min_step(0.0,hextra);
      eng_force(&n,&x,&dir,&x0,&eng);
      nfunc++;
      return 1;
    }

    // check if ready for quadratic projection, equivalent to secant method
    // alpha0 = projected alpha, perform step, exit with success

    de_ideal = -BACKTRACK_SLOPE*alpha*fdotdirall;
    de = eng - eoriginal;
    relerr = fabs(1.0+(0.5*alpha*(alpha-alphaprev)*(fh+fhprev)-eng)/engprev);
    alpha0 = alpha-(alpha-alphaprev)*fh/delfh;

    if (de_ideal >= -QUADRATIC_ENERGY && relerr <= QUADRATIC_TOL && 
	alpha0 > 0.0) {
      for (i = 0; i < n; i++) x[i] = x0[i];
      if (nextra) modify->min_step(0.0,hextra);

      alpha = alpha0;
      for (i = 0; i < n; i++) x[i] += alpha*dir[i];
      if (nextra) modify->min_step(alpha,hextra);
      eng_force(&n,&x,&dir,&x0,&eng);
      nfunc++;
      
      return 0;
    }

    // if backtracking energy change is better than ideal, exit with success
    
    if (de <= de_ideal) return 0;

    // save previous state

    fhprev = fh;
    engprev = eng;
    alphaprev = alpha;

    // reduce alpha

    alpha *= ALPHA_REDUCE;

    // backtracked all the way to 0.0
    // reset to starting point, exit with error

    if (alpha <= 0.0 || de_ideal >= -1.0e-20) {
      for (i = 0; i < n; i++) x[i] = x0[i];
      if (nextra) modify->min_step(0.0,hextra);
      eng_force(&n,&x,&dir,&x0,&eng);
      nfunc++;
      return 1;
    }
  }
}

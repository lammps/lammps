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
#include "fix_minimize.h"
#include "timer.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

#define MIN(A,B) ((A) < (B)) ? (A) : (B)
#define MAX(A,B) ((A) > (B)) ? (A) : (B)

#define EPS         1.0e-6
#define SCAN_FACTOR 2.0
#define SECANT_EPS  1.0e-3

#define SCAN   0   // same as in min.cpp
#define SECANT 1

/* ---------------------------------------------------------------------- */

MinCG::MinCG(LAMMPS *lmp) : Min(lmp) {}

/* ---------------------------------------------------------------------- */

void MinCG::init()
{
  // create fix needed for storing atom-based gradient vectors
  // will delete it at end of run

  char **fixarg = new char*[3];
  fixarg[0] = "MINIMIZE";
  fixarg[1] = "all";
  fixarg[2] = "MINIMIZE";
  modify->add_fix(3,fixarg);
  delete [] fixarg;
  fix_minimize = (FixMinimize *) modify->fix[modify->nfix-1];

  // zero gradient vectors before first atom exchange

  setup_vectors();
  for (int i = 0; i < ndof; i++) h[i] = g[i] = 0.0;

  // virial_thermo is how virial should be computed on thermo timesteps
  // 1 = computed explicity by pair, 2 = computed implicitly by pair

  if (force->newton_pair) virial_thermo = 2;
  else virial_thermo = 1;

  // set flags for what arrays to clear in force_clear()
  // need to clear torques if atom_style is dipole
  // need to clear phia if atom_style is granular
  // don't need to clear f_pair if atom_style is only granular (no virial)

  torqueflag = 0;
  if (atom->check_style("dipole")) torqueflag = 1;
  granflag = 0;
  if (atom->check_style("granular")) granflag = 1;
  pairflag = 1;
  if (strcmp(atom->atom_style,"granular") == 0) pairflag = 0;

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

  if (linestyle == SCAN) linemin = &MinCG::linemin_scan;
  else if (linestyle == SECANT) linemin = &MinCG::linemin_secant;

  // local versions of Update quantities

  maxpair = update->maxpair;
  f_pair = update->f_pair;
}

/* ----------------------------------------------------------------------
   perform minimization, with setup first
------------------------------------------------------------------------- */

void MinCG::run()
{
  double tmp,*f;

  // set initial force & energy

  setup();
  setup_vectors();
  output->thermo->compute_pe();
  ecurrent = output->thermo->potential_energy;

  // stats for Finish to print
	
  einitial = ecurrent;

  f = atom->f[0];
  tmp = 0.0;
  for (int i = 0; i < ndof; i++) tmp += f[i]*f[i];
  MPI_Allreduce(&tmp,&gnorm2_init,1,MPI_DOUBLE,MPI_SUM,world);
  gnorm2_init = sqrt(gnorm2_init);

  tmp = 0.0;
  for (int i = 0; i < ndof; i++) tmp = MAX(fabs(f[i]),tmp);
  MPI_Allreduce(&tmp,&gnorminf_init,1,MPI_DOUBLE,MPI_MAX,world);

  // minimizer iterations

  timer->barrier_start(TIME_LOOP);
  iterate(update->nsteps);

  // account for early exit from iterate loop due to convergence
  // set niter/nsteps for Finish stats to print
  // set output->next values to this timestep
  // call eng_force to insure vflag is set when forces computed
  // output->write does final output for thermo, dump, restart files

  if (niter < update->nsteps) {
    niter++;
    update->nsteps = niter;
    for (int idump = 0; idump < output->ndump; idump++)
      output->next_dump[idump] = update->ntimestep;
    output->next_dump_any = update->ntimestep;
    if (output->restart_every) output->next_restart = update->ntimestep;
    output->next_thermo = update->ntimestep;
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

  f = atom->f[0];
  tmp = 0.0;
  for (int i = 0; i < ndof; i++) tmp += f[i]*f[i];
  MPI_Allreduce(&tmp,&gnorm2_final,1,MPI_DOUBLE,MPI_SUM,world);
  gnorm2_final = sqrt(gnorm2_final);

  tmp = 0.0;
  for (int i = 0; i < ndof; i++) tmp = MAX(fabs(f[i]),tmp);
  MPI_Allreduce(&tmp,&gnorminf_final,1,MPI_DOUBLE,MPI_MAX,world);
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

  domain->pbc();
  domain->reset_box();
  comm->setup();
  if (neighbor->style) neighbor->setup_bins();
  comm->exchange();
  comm->borders();
  neighbor->build();
  neighbor->ncalls = 0;
  setup_vectors();

  // compute all forces

  int eflag = 1;
  int vflag = virial_thermo;
  force_clear(vflag);
  
  if (atom->molecular) {
    if (force->bond) force->bond->compute(eflag,vflag);
    if (force->angle) force->angle->compute(eflag,vflag);
    if (force->dihedral) force->dihedral->compute(eflag,vflag);
    if (force->improper) force->improper->compute(eflag,vflag);
  }

  if (force->pair) force->pair->compute(eflag,vflag);

  if (force->kspace) {
    force->kspace->setup();
    force->kspace->compute(eflag,vflag);
  }

  if (force->newton) comm->reverse_communicate();

  modify->setup();
  output->setup(1);
}

/* ----------------------------------------------------------------------
   minimization via conjugate gradient iterations
   Polak-Ribiere formulation
------------------------------------------------------------------------- */

void MinCG::iterate(int n)
{
  int i,gradsearch,fail;
  double alpha,beta,gg,dot[2],dotall[2];
  double *f;

  f = atom->f[0];
  for (int i = 0; i < ndof; i++) h[i] = g[i] = f[i];

  dot[0] = 0.0;
  for (i = 0; i < ndof; i++) dot[0] += f[i]*f[i];
  MPI_Allreduce(dot,&gg,1,MPI_DOUBLE,MPI_SUM,world);

  neval = 0;
  gradsearch = 1;

  for (niter = 0; niter < n; niter++) {

    update->ntimestep++;

    // line minimization along direction h from current atom->x

    eprevious = ecurrent;
    fail = (this->*linemin)(ndof,atom->x[0],h,ecurrent,dmin,dmax,alpha,neval);

    // if max_eval exceeded, all done
    // if linemin failed or energy did not decrease sufficiently:
    //   all done if searched in grad direction
    //   else force next search to be in grad direction (CG restart)

    if (neval >= update->max_eval) break;

    if (fail || fabs(ecurrent-eprevious) <= 
    	update->tolerance * 0.5*(fabs(ecurrent) + fabs(eprevious) + EPS)) {
      if (gradsearch == 1) break;
      gradsearch = -1;
    }

    // update h from new f = -Grad(x) and old g
    // old g,h must have migrated with atoms to do this correctly
    // done if size sq of grad vector < EPS
    // force new search dir to be grad dir if need to restart CG
    // set gradsearch to 1 if will search in grad dir on next iteration

    f = atom->f[0];
    dot[0] = dot[1] = 0.0;
    for (i = 0; i < ndof; i++) {
      dot[0] += f[i]*f[i];
      dot[1] += f[i]*g[i];
    }
    MPI_Allreduce(dot,dotall,2,MPI_DOUBLE,MPI_SUM,world);

    beta = MAX(0.0,(dotall[0] - dotall[1])/gg);
    gg = dotall[0];
    if (gg < EPS) break;

    if (gradsearch == -1) beta = 0.0;
    if (beta == 0.0) gradsearch = 1;
    else gradsearch = 0;

    for (i = 0; i < ndof; i++) {
      g[i] = f[i];
      h[i] = g[i] + beta*h[i];
    }

    // output for thermo, dump, restart files

    if (output->next == update->ntimestep) {
      timer->stamp();
      output->write(update->ntimestep);
      timer->stamp(TIME_OUTPUT);
    }
  }
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
    domain->pbc();
    if (domain->box_change) {
      domain->reset_box();
      comm->setup();
      if (neighbor->style) neighbor->setup_bins();
    }
    timer->stamp();
    comm->exchange();
    comm->borders();
    timer->stamp(TIME_COMM);
    neighbor->build();
    timer->stamp(TIME_NEIGHBOR);
    setup_vectors();
  }

  // eflag is always set, since minimizer needs potential energy

  int eflag = 1;
  int vflag = 0;
  if (output->next_thermo == update->ntimestep) vflag = virial_thermo;
  force_clear(vflag);

  timer->stamp();
  if (atom->molecular) {
    if (force->bond) force->bond->compute(eflag,vflag);
    if (force->angle) force->angle->compute(eflag,vflag);
    if (force->dihedral) force->dihedral->compute(eflag,vflag);
    if (force->improper) force->improper->compute(eflag,vflag);
    timer->stamp(TIME_BOND);
  }

  if (force->pair) {
    force->pair->compute(eflag,vflag);
    timer->stamp(TIME_PAIR);
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

  // compute potential energy of system via Thermo

  output->thermo->compute_pe();
  ecurrent = output->thermo->potential_energy;

  // return updated ptrs to caller since atoms may have migrated

  *pndof = ndof;
  *px = atom->x[0];
  *ph = h;
  *peng = ecurrent;
}

/* ----------------------------------------------------------------------
   clear force on own & ghost atoms
   setup and clear other arrays as needed
------------------------------------------------------------------------- */

void MinCG::force_clear(int vflag)
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

  if (granflag) {
    double **phia = atom->phia;
    for (i = 0; i < nall; i++) {
      phia[i][0] = 0.0;
      phia[i][1] = 0.0;
      phia[i][2] = 0.0;
    }
  }

  // clear f_pair array if using it this timestep to compute virial

  if (vflag == 2 && pairflag) {
    if (atom->nmax > maxpair) {
      maxpair = atom->nmax;
      memory->destroy_2d_double_array(f_pair);
      f_pair = memory->create_2d_double_array(maxpair,3,"min:f_pair");
      update->maxpair = maxpair;
      update->f_pair = f_pair;
    }
    for (i = 0; i < nall; i++) {
      f_pair[i][0] = 0.0;
      f_pair[i][1] = 0.0;
      f_pair[i][2] = 0.0;
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
	  min/max dist = min/max distance to move any atom coord
   output: return 0 if successful move, set alpha
           return 1 if failed, no move, no need to set alpha
           alpha = distance moved along dir to set x to min-eng config
           caller has several quantities set via last call to eng_force()
	     INSURE last call to eng_force() is consistent with returns
	       if fail, eng_force() of original x
	       if succeed, eng_force() at x + alpha*dir
             atom->x = coords at new configuration
	     atom->f = force (-Grad) is evaulated at new configuration
	     ecurrent = energy of new configuration
   NOTE: when call eng_force: n,x,dir,eng may change due to atom migration
	 updated values are returned by eng_force()
	 this routine CANNOT store atom-based quantities b/c of migration
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   linemin: scan forward by larger and larger steps (SCAN_FACTOR)
   uses no gradient info, but should be very robust
   start at mindist, continue until maxdist
   quit as soon as energy starts to rise
------------------------------------------------------------------------- */

int MinCG::linemin_scan(int n, double *x, double *dir, double eng,
			double mindist, double maxdist,
			double &alpha, int &nfunc)
{
  int i;
  double fmax,fme,elowest,alphamin,alphamax,alphalast;

  // alphamin = step that moves some atom coord by mindist
  // alphamax = step that moves some atom coord by maxdist

  fme = 0.0;
  for (i = 0; i < n; i++) fme = MAX(fme,fabs(dir[i]));
  MPI_Allreduce(&fme,&fmax,1,MPI_DOUBLE,MPI_MAX,world);
  if (fmax == 0.0) return 1;

  alphamin = mindist/fmax;
  alphamax = maxdist/fmax;

  // if minstep is already uphill, fail
  // if eng increases, stop and return previous alpha
  // if alphamax, stop and return alphamax

  elowest = eng;
  alpha = alphamin;

  while (1) {
    for (i = 0; i < n; i++) x[i] += alpha*dir[i];
    eng_force(&n,&x,&dir,&eng);
    nfunc++;

    if (alpha == alphamin && eng >= elowest) {
      for (i = 0; i < n; i++) x[i] -= alpha*dir[i];
      eng_force(&n,&x,&dir,&eng);
      nfunc++;
      return 1;
    }
    if (eng > elowest) {
      for (i = 0; i < n; i++) x[i] += (alphalast-alpha)*dir[i];
      eng_force(&n,&x,&dir,&eng);
      nfunc++;
      alpha = alphalast;
      return 0;
    }      
    if (alpha == alphamax) return 0;

    elowest = eng;
    alphalast = alpha;
    alpha *= SCAN_FACTOR;
    if (alpha > alphamax) alpha = alphamax;
  }
}

/* ----------------------------------------------------------------------
   linemin: use secant approximation to estimate parabola minimum at each step
   should converge more quickly/accurately than "scan", but may be less robust
   how to prevent func evals at new x bigger than maxdist?
   initial secant from two points: 0 and sigma0 = mindist
------------------------------------------------------------------------- */

int MinCG::linemin_secant(int n, double *x, double *dir, double eng,
			  double mindist, double maxdist,
			  double &alpha, int &nfunc)
{
  int i,iter;
  double eta,eta_prev,sigma0,alphadelta,fme,fmax,dsq,e0,tmp;
  double *f;
  double epssq = SECANT_EPS * SECANT_EPS;

  // stopping criterion for secant iterations
  // change this

  fme = 0.0;
  for (i = 0; i < n; i++) fme += dir[i]*dir[i];
  MPI_Allreduce(&fme,&dsq,1,MPI_DOUBLE,MPI_SUM,world);

  // sigma0 = smallest allowed step of mindist
  // eval func at sigma0
  // test if minstep is already uphill

  fme = 0.0;
  for (i = 0; i < n; i++) fme = MAX(fme,fabs(dir[i]));
  MPI_Allreduce(&fme,&fmax,1,MPI_DOUBLE,MPI_MAX,world);
  if (fmax == 0.0) return 1;

  sigma0 = mindist/fmax;

  e0 = eng;
  for (i = 0; i < n; i++) x[i] += sigma0*dir[i];
  eng_force(&n,&x,&dir,&eng);
  nfunc++;

  if (eng >= e0) {
    for (i = 0; i < n; i++) x[i] -= sigma0*dir[i];
    eng_force(&n,&x,&dir,&eng);
    nfunc++;
    return 1;
  }

  f = atom->f[0];
  tmp = 0.0;
  for (i = 0; i < n; i++) tmp -= f[i]*dir[i];
  MPI_Allreduce(&tmp,&eta_prev,1,MPI_DOUBLE,MPI_SUM,world);

  // secant iterations
  // alphadelta = new increment to move, alpha = accumulated move

  alpha = sigma0;
  alphadelta = -sigma0;

  for (iter = 0; iter < lineiter; iter++) {
    alpha += alphadelta;
    for (i = 0; i < n; i++) x[i] += alphadelta*dir[i];
    eng_force(&n,&x,&dir,&eng);
    nfunc++;

    f = atom->f[0];
    tmp = 0.0;
    for (i = 0; i < n; i++) tmp -= f[i]*dir[i];
    MPI_Allreduce(&tmp,&eta,1,MPI_DOUBLE,MPI_SUM,world);

    alphadelta *= eta / (eta_prev - eta);
    eta_prev = eta;
    if (alphadelta*alphadelta*dsq <= epssq) break;
  }

  // if exited loop on first iteration, func eval was at alpha = 0.0
  // else successful line search

  if (iter == 0) return 1;
  return 0;
}

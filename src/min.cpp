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
   Contributing author: Aidan Thompson (SNL)
                        improved CG and backtrack ls, added quadratic ls
   Sources: Numerical Recipes frprmn routine
            "Conjugate Gradient Method Without the Agonizing Pain" by
            JR Shewchuk, http://www-2.cs.cmu.edu/~jrs/jrspapers.html#cg
------------------------------------------------------------------------- */

#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "min.h"
#include "atom.h"
#include "domain.h"
#include "comm.h"
#include "update.h"
#include "modify.h"
#include "fix_minimize.h"
#include "compute.h"
#include "neighbor.h"
#include "force.h"
#include "pair.h"
#include "bond.h"
#include "angle.h"
#include "dihedral.h"
#include "improper.h"
#include "kspace.h"
#include "output.h"
#include "thermo.h"
#include "timer.h"
#include "error.h"

using namespace LAMMPS_NS;

// ALPHA_MAX = max alpha allowed to avoid long backtracks
// ALPHA_REDUCE = reduction ratio, should be in range [0.5,1)
// BACKTRACK_SLOPE, should be in range (0,0.5]
// QUADRATIC_TOL = tolerance on alpha0, should be in range [0.1,1)
// IDEAL_TOL = ideal energy tolerance for backtracking
// EPS_QUAD = tolerance for quadratic projection

#define ALPHA_MAX 1.0
#define ALPHA_REDUCE 0.5
#define BACKTRACK_SLOPE 0.4
#define QUADRATIC_TOL 0.1
#define IDEAL_TOL 1.0e-8
#define EPS_QUAD 1.0e-28

// same as in other min classes

enum{MAXITER,MAXEVAL,ETOL,FTOL,DOWNHILL,ZEROALPHA,ZEROFORCE,ZEROQUAD};

#define MIN(A,B) ((A) < (B)) ? (A) : (B)
#define MAX(A,B) ((A) > (B)) ? (A) : (B)

/* ---------------------------------------------------------------------- */

Min::Min(LAMMPS *lmp) : Pointers(lmp)
{
  dmax = 0.1;
  linestyle = 0;

  elist_atom = NULL;
  vlist_global = vlist_atom = NULL;

  fextra = gextra = hextra = NULL;
}

/* ---------------------------------------------------------------------- */

Min::~Min()
{
  delete [] elist_atom;
  delete [] vlist_global;
  delete [] vlist_atom;

  delete [] fextra;
  delete [] gextra;
  delete [] hextra;
}

/* ---------------------------------------------------------------------- */

void Min::init()
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

  if (linestyle == 0) linemin = &Min::linemin_backtrack;
  else if (linestyle == 1) linemin = &Min::linemin_quadratic;
}

/* ----------------------------------------------------------------------
   setup before run
------------------------------------------------------------------------- */

void Min::setup()
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
   perform minimization, with setup first
------------------------------------------------------------------------- */

void Min::run()
{
  int i;
  double tmp,*f;

  // possible stop conditions

  char *stopstrings[] = {"max iterations","max force evaluations",
                         "energy tolerance","force tolerance",
			 "search direction is not downhill",
			 "linesearch alpha is zero",
			 "forces are zero","quadratic factors are zero"};

  // set initial force & energy

  setup();

  // setup any extra dof due to fixes
  // can't be done until now b/c update init() comes before modify init()

  delete [] fextra;
  delete [] gextra;
  delete [] hextra;
  fextra = NULL;
  gextra = NULL;
  hextra = NULL;

  nextra = modify->min_dof();
  if (nextra) {
    fextra = new double[nextra];
    gextra = new double[nextra];
    hextra = new double[nextra];
  }

  // compute potential energy of system
  // normalize energy if thermo PE does

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
    eng_force(&ntmp,&xtmp,&htmp,&x0tmp,&etmp,0);
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
   evaluate potential energy and forces
   may migrate atoms
   if resetflag = 1, update x0 by PBC for atoms that migrate
   new energy stored in ecurrent and returned (in case caller not in class)
   negative gradient will be stored in atom->f
------------------------------------------------------------------------- */

void Min::eng_force(int *pndof, double **px, double **ph, double **px0,
		    double *peng, int resetflag)
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

    // update x0 for atoms that migrated
    // must do minimum_image on box size when x0 was stored
    // domain->set_global_box() changes to x0 box, then restores current box

    if (resetflag) {
      box_swap();
      domain->set_global_box();

      double **x = atom->x;
      double **x0 = fix_minimize->x0;
      int nlocal = atom->nlocal;

      double dx,dy,dz,dx0,dy0,dz0;
      for (int i = 0; i < nlocal; i++) {
	dx = dx0 = x[i][0] - x0[i][0];
	dy = dy0 = x[i][1] - x0[i][1];
	dz = dz0 = x[i][2] - x0[i][2];
	domain->minimum_image(dx,dy,dz);
	if (dx != dx0) x0[i][0] = x[i][0] - dx;
	if (dy != dy0) x0[i][1] = x[i][1] - dy;
	if (dz != dz0) x0[i][2] = x[i][2] - dz;
      }

      box_swap();
      domain->set_global_box();
    }
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
   set ndof and vector pointers after atoms have migrated
------------------------------------------------------------------------- */

void Min::setup_vectors()
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
   clear force on own & ghost atoms
   setup and clear other arrays as needed
------------------------------------------------------------------------- */

void Min::force_clear()
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
	  eoriginal = energy at initial x
	  maxdist = max distance to move any atom coord
   output: return 0 if successful move, non-zero alpha
           return non-zero if failed
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

int Min::linemin_backtrack(int n, double *x, double *dir,
			   double *x0, double eoriginal, double maxdist,
			   double &alpha, int &nfunc)
{
  int i,m;
  double fdotdirall,fdotdirme,hmax,hme,alpha_extra;
  double eng,de_ideal,de;

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
  if (fdotdirall <= 0.0) return DOWNHILL;

  // initial alpha = stepsize to change any atom coord by maxdist
  // alpha <= ALPHA_MAX, else backtrack from huge value when forces are tiny
  // if all search dir components are already 0.0, exit with error

  hme = 0.0;
  for (i = 0; i < n; i++) hme = MAX(hme,fabs(dir[i]));
  MPI_Allreduce(&hme,&hmax,1,MPI_DOUBLE,MPI_MAX,world);
  alpha = MIN(ALPHA_MAX,maxdist/hmax);
  if (nextra) {
    double alpha_extra = modify->max_alpha(hextra);
    alpha = MIN(alpha,alpha_extra);
    for (i = 0; i < nextra; i++)
      hmax = MAX(hmax,fabs(hextra[i]));
  }
  if (hmax == 0.0) return ZEROFORCE;

  // store coords and other dof at start of linesearch

  box_store();
  for (i = 0; i < n; i++) x0[i] = x[i];
  if (nextra) modify->min_store();

  // backtrack with alpha until energy decrease is sufficient

  while (1) {
    if (nextra) modify->min_step(0.0,hextra);
    for (i = 0; i < n; i++) x[i] = x0[i];
    if (nextra) modify->min_step(alpha,hextra);
    for (i = 0; i < n; i++) x[i] += alpha*dir[i];
    eng_force(&n,&x,&dir,&x0,&eng,1);
    nfunc++;

    // if energy change is better than ideal, exit with success

    de_ideal = -BACKTRACK_SLOPE*alpha*fdotdirall;
    de = eng - eoriginal;
    if (de <= de_ideal) return 0;

    // reduce alpha

    alpha *= ALPHA_REDUCE;

    // backtracked all the way to 0.0
    // reset to starting point, exit with error

    if (alpha <= 0.0 || de_ideal >= -IDEAL_TOL) {
      if (nextra) modify->min_step(0.0,hextra);
      for (i = 0; i < n; i++) x[i] = x0[i];
      eng_force(&n,&x,&dir,&x0,&eng,0);
      nfunc++;
      return ZEROALPHA;
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

int Min::linemin_quadratic(int n, double *x, double *dir,
			   double *x0, double eoriginal, double maxdist,
			   double &alpha, int &nfunc)
{
  int i,m;
  double fdotdirall,fdotdirme,hmax,hme,alphamax,alpha_extra;
  double eng,de_ideal,de;
  double delfh,engprev,relerr,alphaprev,fhprev,ff,fh,alpha0,fh0,ff0;
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
  if (fdotdirall <= 0.0) return DOWNHILL;

  // initial alpha = stepsize to change any atom coord by maxdist
  // alpha <= ALPHA_MAX, else backtrack from huge value when forces are tiny
  // if all search dir components are already 0.0, exit with error

  hme = 0.0;
  for (i = 0; i < n; i++) hme = MAX(hme,fabs(dir[i]));
  MPI_Allreduce(&hme,&hmax,1,MPI_DOUBLE,MPI_MAX,world);
  alpha = MIN(ALPHA_MAX,maxdist/hmax);
  if (nextra) {
    double alpha_extra = modify->max_alpha(hextra);
    alpha = MIN(alpha,alpha_extra);
    for (i = 0; i < nextra; i++)
      hmax = MAX(hmax,fabs(hextra[i]));
  }
  if (hmax == 0.0) return ZEROFORCE;

  // store coords and other dof at start of linesearch

  box_store();
  for (i = 0; i < n; i++) x0[i] = x[i];
  if (nextra) modify->min_store();

  // backtrack with alpha until energy decrease is sufficient
  // or until get to small energy change, then perform quadratic projection

  fhprev = fdotdirall;
  engprev = eoriginal;
  alphaprev = 0.0;

  while (1) {
    if (nextra) modify->min_step(0.0,hextra);
    for (i = 0; i < n; i++) x[i] = x0[i];
    if (nextra) modify->min_step(alpha,hextra);
    for (i = 0; i < n; i++) x[i] += alpha*dir[i];
    eng_force(&n,&x,&dir,&x0,&eng,1);
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
      if (nextra) modify->min_step(0.0,hextra);
      for (i = 0; i < n; i++) x[i] = x0[i];
      eng_force(&n,&x,&dir,&x0,&eng,0);
      nfunc++;
      return ZEROQUAD;
    }

    // check if ready for quadratic projection, equivalent to secant method
    // alpha0 = projected alpha

    relerr = fabs(1.0+(0.5*alpha*(alpha-alphaprev)*(fh+fhprev)-eng)/engprev);
    alpha0 = alpha - (alpha-alphaprev)*fh/delfh;

    if (relerr <= QUADRATIC_TOL && alpha0 > 0.0) {
      if (nextra) modify->min_step(0.0,hextra);
      for (i = 0; i < n; i++) x[i] = x0[i];

      if (nextra) modify->min_step(alpha0,hextra);
      for (i = 0; i < n; i++) x[i] += alpha0*dir[i];
      eng_force(&n,&x,&dir,&x0,&eng,1);
      nfunc++;

      // if backtracking energy change is better than ideal, exit with success

      de_ideal = -BACKTRACK_SLOPE*alpha0*fdotdirall;
      de = eng - eoriginal;
      if (de <= de_ideal || de_ideal >= -IDEAL_TOL) return 0;

      // drop back from alpha0 to alpha

      if (nextra) modify->min_step(0.0,hextra);
      for (i = 0; i < n; i++) x[i] = x0[i];
      if (nextra) modify->min_step(alpha,hextra);
      for (i = 0; i < n; i++) x[i] += alpha*dir[i];
      eng_force(&n,&x,&dir,&x0,&eng,1);
      nfunc++;
    }

    // if backtracking energy change is better than ideal, exit with success

    de_ideal = -BACKTRACK_SLOPE*alpha*fdotdirall;
    de = eng - eoriginal;
    if (de <= de_ideal) return 0;

    // save previous state

    fhprev = fh;
    engprev = eng;
    alphaprev = alpha;

    // reduce alpha

    alpha *= ALPHA_REDUCE;

    // backtracked all the way to 0.0
    // reset to starting point, exit with error

    if (alpha <= 0.0 || de_ideal >= -IDEAL_TOL) {
      if (nextra) modify->min_step(0.0,hextra);
      for (i = 0; i < n; i++) x[i] = x0[i];
      eng_force(&n,&x,&dir,&x0,&eng,0);
      nfunc++;
      return ZEROALPHA;
    }
  }
}

/* ---------------------------------------------------------------------- */

void Min::modify_params(int narg, char **arg)
{
  if (narg == 0) error->all("Illegal min_modify command");

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"dmax") == 0) {
      if (iarg+2 > narg) error->all("Illegal min_modify command");
      dmax = atof(arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"line") == 0) {
      if (iarg+2 > narg) error->all("Illegal min_modify command");
      if (strcmp(arg[iarg+1],"backtrack") == 0) linestyle = 0;
      else if (strcmp(arg[iarg+1],"quadratic") == 0) linestyle = 1;
      else error->all("Illegal min_modify command");
      iarg += 2;
    } else error->all("Illegal min_modify command");
  }
}

/* ----------------------------------------------------------------------
   setup lists of computes for global and per-atom PE and pressure
------------------------------------------------------------------------- */

void Min::ev_setup()
{
  delete [] elist_atom;
  delete [] vlist_global;
  delete [] vlist_atom;
  elist_atom = NULL;
  vlist_global = vlist_atom = NULL;

  nelist_atom = 0;
  nvlist_global = nvlist_atom = 0;
  for (int i = 0; i < modify->ncompute; i++) {
    if (modify->compute[i]->peatomflag) nelist_atom++;
    if (modify->compute[i]->pressflag) nvlist_global++;
    if (modify->compute[i]->pressatomflag) nvlist_atom++;
  }

  if (nelist_atom) elist_atom = new Compute*[nelist_atom];
  if (nvlist_global) vlist_global = new Compute*[nvlist_global];
  if (nvlist_atom) vlist_atom = new Compute*[nvlist_atom];

  nelist_atom = 0;
  nvlist_global = nvlist_atom = 0;
  for (int i = 0; i < modify->ncompute; i++) {
    if (modify->compute[i]->peatomflag)
      elist_atom[nelist_atom++] = modify->compute[i];
    if (modify->compute[i]->pressflag)
      vlist_global[nvlist_global++] = modify->compute[i];
    if (modify->compute[i]->pressatomflag)
      vlist_atom[nvlist_atom++] = modify->compute[i];
  }
}

/* ----------------------------------------------------------------------
   set eflag,vflag for current iteration with ntimestep
   always set eflag_global = 1, since need energy every iteration
   eflag = 0 = no energy computation
   eflag = 1 = global energy only
   eflag = 2 = per-atom energy only
   eflag = 3 = both global and per-atom energy
   vflag = 0 = no virial computation (pressure)
   vflag = 1 = global virial with pair portion via sum of pairwise interactions
   vflag = 2 = global virial with pair portion via F dot r including ghosts
   vflag = 4 = per-atom virial only
   vflag = 5 or 6 = both global and per-atom virial
------------------------------------------------------------------------- */

void Min::ev_set(int ntimestep)
{
  int i;

  int eflag_global = 1;

  int eflag_atom = 0;
  for (i = 0; i < nelist_atom; i++)
    if (elist_atom[i]->matchstep(ntimestep)) break;
  if (i < nelist_atom) eflag_atom = 2;

  if (eflag_global) update->eflag_global = update->ntimestep;
  if (eflag_atom) update->eflag_atom = update->ntimestep;
  eflag = eflag_global + eflag_atom;

  int vflag_global = 0;
  for (i = 0; i < nvlist_global; i++)
    if (vlist_global[i]->matchstep(ntimestep)) break;
  if (i < nvlist_global) vflag_global = virial_style;

  int vflag_atom = 0;
  for (i = 0; i < nvlist_atom; i++)
    if (vlist_atom[i]->matchstep(ntimestep)) break;
  if (i < nvlist_atom) vflag_atom = 4;

  if (vflag_global) update->vflag_global = update->ntimestep;
  if (vflag_atom) update->vflag_atom = update->ntimestep;
  vflag = vflag_global + vflag_atom;
}

/* ----------------------------------------------------------------------
   store box size at beginning of line search
------------------------------------------------------------------------- */

void Min::box_store()
{
  boxlo0[0] = domain->boxlo[0];
  boxlo0[1] = domain->boxlo[1];
  boxlo0[2] = domain->boxlo[2];

  boxhi0[0] = domain->boxhi[0];
  boxhi0[1] = domain->boxhi[1];
  boxhi0[2] = domain->boxhi[2];
}

/* ----------------------------------------------------------------------
   swap current box size with stored box size
------------------------------------------------------------------------- */

void Min::box_swap()
{
  double tmp;

  tmp = boxlo0[0];
  boxlo0[0] = domain->boxlo[0];
  domain->boxlo[0] = tmp;
  tmp = boxlo0[1];
  boxlo0[1] = domain->boxlo[1];
  domain->boxlo[1] = tmp;
  tmp = boxlo0[2];
  boxlo0[2] = domain->boxlo[2];
  domain->boxlo[2] = tmp;

  tmp = boxhi0[0];
  boxhi0[0] = domain->boxhi[0];
  domain->boxhi[0] = tmp;
  tmp = boxhi0[1];
  boxhi0[1] = domain->boxhi[1];
  domain->boxhi[1] = tmp;
  tmp = boxhi0[2];
  boxhi0[2] = domain->boxhi[2];
  domain->boxhi[2] = tmp;
}

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

  // zero local vectors before first atom exchange

  set_local_vectors();
  for (int i = 0; i < ndof; i++) h[i] = g[i] = 0.0;

  // virial_thermo = how virial computed on thermo timesteps
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
  int i;
  double tmp;

  // set initial force & energy

  setup();
  output->thermo->compute_pe();
  energy = output->thermo->potential_energy;
  energy += energy_extra;

  // stats for Finish to print
	
  einitial = energy;

  tmp = 0.0;
  for (i = 0; i < ndof; i++) tmp += f[i]*f[i];
  MPI_Allreduce(&tmp,&gnorm2_init,1,MPI_DOUBLE,MPI_SUM,world);
  for (i = 0; i < ndof_extra; i++) gnorm2_init += fextra[i]*fextra[i];
  gnorm2_init = sqrt(gnorm2_init);

  tmp = 0.0;
  for (i = 0; i < ndof; i++) tmp = MAX(fabs(f[i]),tmp);
  MPI_Allreduce(&tmp,&gnorminf_init,1,MPI_DOUBLE,MPI_MAX,world);
  for (i = 0; i < ndof_extra; i++)
    gnorminf_init = MAX(gnorminf_init,fabs(fextra[i]));

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
    eng_force();
    output->write(update->ntimestep);
  }
  timer->barrier_stop(TIME_LOOP);

  // delete fix at end of run, so its atom arrays won't persist
  // delete extra arrays

  modify->delete_fix("MINIMIZE");
  delete [] xextra;
  delete [] fextra;
  delete [] gextra;
  delete [] hextra;

  // reset reneighboring criteria

  neighbor->every = neigh_every;
  neighbor->delay = neigh_delay;
  neighbor->dist_check = neigh_dist_check;

  // stats for Finish to print
	
  efinal = energy;

  f = atom->f[0];
  tmp = 0.0;
  for (i = 0; i < ndof; i++) tmp += f[i]*f[i];
  MPI_Allreduce(&tmp,&gnorm2_final,1,MPI_DOUBLE,MPI_SUM,world);
  for (i = 0; i < ndof_extra; i++) gnorm2_final += fextra[i]*fextra[i];
  gnorm2_final = sqrt(gnorm2_final);

  tmp = 0.0;
  for (i = 0; i < ndof; i++) tmp = MAX(fabs(f[i]),tmp);
  MPI_Allreduce(&tmp,&gnorminf_final,1,MPI_DOUBLE,MPI_MAX,world);
  for (i = 0; i < ndof_extra; i++)
    gnorminf_final = MAX(gnorminf_final,fabs(fextra[i]));
}

/* ----------------------------------------------------------------------
   setup before run
------------------------------------------------------------------------- */

void MinCG::setup()
{
  // allocate extra arrays
  // couldn't do in init(), b/c modify and fixes weren't yet init()
  // set initial xextra values via fixes

  ndof_extra = modify->min_dof();
  xextra = new double[ndof_extra];
  fextra = new double[ndof_extra];
  gextra = new double[ndof_extra];
  hextra = new double[ndof_extra];

  modify->min_xinitial(xextra);

  // perform usual setup

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
  set_local_vectors();

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
  energy_extra = modify->min_energy(xextra,fextra);
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

  for (i = 0; i < ndof; i++) h[i] = g[i] = f[i];
  for (i = 0; i < ndof_extra; i++) hextra[i] = gextra[i] = fextra[i];

  dot[0] = 0.0;
  for (i = 0; i < ndof; i++) dot[0] += f[i]*f[i];
  MPI_Allreduce(dot,&gg,1,MPI_DOUBLE,MPI_SUM,world);
  for (i = 0; i < ndof_extra; i++) gg += fextra[i]*fextra[i];

  neval = 0;
  gradsearch = 1;

  for (niter = 0; niter < n; niter++) {

    update->ntimestep++;

    // line minimization along direction h from current atom->x

    eprevious = energy;
    fail = (this->*linemin)(neval);

    // if max_eval exceeded, all done
    // if linemin failed or energy did not decrease sufficiently:
    //   all done if searched in grad direction
    //   else force next search to be in grad direction (CG restart)

    if (neval >= update->max_eval) break;

    if (fail || fabs(energy-eprevious) <= 
    	update->tolerance * 0.5*(fabs(energy) + fabs(eprevious) + EPS)) {
      if (gradsearch == 1) break;
      gradsearch = -1;
    }

    // update h from new f = -Grad(x) and old g
    // old g,h must have migrated with atoms to do this correctly
    // done if size sq of grad vector < EPS
    // force new search dir to be grad dir if need to restart CG
    // set gradsearch to 1 if will search in grad dir on next iteration

    dot[0] = dot[1] = 0.0;
    for (i = 0; i < ndof; i++) {
      dot[0] += f[i]*f[i];
      dot[1] += f[i]*g[i];
    }
    MPI_Allreduce(dot,dotall,2,MPI_DOUBLE,MPI_SUM,world);
    for (i = 0; i < ndof_extra; i++) {
      dotall[0] += fextra[i]*fextra[i];
      dotall[1] += fextra[i]*gextra[i];
    }
    
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
    for (i = 0; i < ndof_extra; i++) {
      gextra[i] = fextra[i];
      hextra[i] = gextra[i] + beta*hextra[i];
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

void MinCG::set_local_vectors()
{
  ndof = 3 * atom->nlocal;
  x = atom->x[0];
  f = atom->f[0];
  if (ndof) g = fix_minimize->gradient[0];
  else g = NULL;
  if (ndof) h = fix_minimize->searchdir[0];
  else h = NULL;
}

/* ----------------------------------------------------------------------
   evaluate potential energy and forces
   may migrate atoms
   energy = new objective function energy = poteng of atoms + eng_extra
   atom->f, fextra = negative gradient of objective function
------------------------------------------------------------------------- */

void MinCG::eng_force()
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
    set_local_vectors();
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

  // min_post_force = forces on atoms that affect minimization
  // min_energy = energy, forces on extra degrees of freedom

  if (modify->n_min_post_force) modify->min_post_force(vflag);
  if (modify->n_min_energy) energy_extra = modify->min_energy(xextra,fextra);

  // compute potential energy of system via Thermo

  output->thermo->compute_pe();
  energy = output->thermo->potential_energy;
  energy += energy_extra;
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
   find minimum-energy starting at x along h direction
   update atom->x by alpha, call eng_force() for result
   alpha = distance moved along h to set x to minimun-energy configuration
   return 0 if successful move, 1 if failed (no move)
   insure last call to eng_force() is consistent with return
     if fail, eng_force() of original x
     if succeed, eng_force() at x + alpha*h
   eng_force() may migrate atoms due to neighbor list build
     therefore linemin routines CANNOT store atom-based quantities
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   linemin: scan forward by larger and larger steps (SCAN_FACTOR)
   uses no gradient info, but should be very robust
   start at mindist, continue until maxdist
   quit as soon as energy starts to rise
------------------------------------------------------------------------- */

int MinCG::linemin_scan(int &nfunc)
{
  int i;
  double fmax,fme,elowest,alpha,alphamin,alphamax,alphalast;

  // alphamin = step that moves some atom coord by mindist
  // alphamax = step that moves some atom coord by maxdist

  fme = 0.0;
  for (i = 0; i < ndof; i++) fme = MAX(fme,fabs(h[i]));
  MPI_Allreduce(&fme,&fmax,1,MPI_DOUBLE,MPI_MAX,world);
  for (i = 0; i < ndof_extra; i++) fmax = MAX(fmax,fabs(hextra[i]));

  if (fmax == 0.0) return 1;
  alphamin = dmin/fmax;
  alphamax = dmax/fmax;

  // if minstep is already uphill, fail
  // if eng increases, stop and return previous alpha
  // if alphamax, stop and return alphamax

  elowest = energy;
  alpha = alphamin;

  while (1) {
    for (i = 0; i < ndof; i++) x[i] += alpha*h[i];
    for (i = 0; i < ndof_extra; i++) xextra[i] += alpha*hextra[i];
    eng_force();
    nfunc++;

    if (alpha == alphamin && energy >= elowest) {
      for (i = 0; i < ndof; i++) x[i] -= alpha*h[i];
      for (i = 0; i < ndof_extra; i++) xextra[i] -= alpha*hextra[i];
      eng_force();
      nfunc++;
      return 1;
    }
    if (energy > elowest) {
      for (i = 0; i < ndof; i++) x[i] += (alphalast-alpha)*h[i];
      for (i = 0; i < ndof_extra; i++)
	xextra[i] += (alphalast-alpha)*hextra[i];
      eng_force();
      nfunc++;
      alpha = alphalast;
      return 0;
    }      
    if (alpha == alphamax) return 0;

    elowest = energy;
    alphalast = alpha;
    alpha *= SCAN_FACTOR;
    if (alpha > alphamax) alpha = alphamax;
  }
}

/* ----------------------------------------------------------------------
   linemin: use secant approximation to estimate parabola minimum at each step
   should converge more quickly/accurately than "scan", but may be less robust
   initial secant from two points: 0 and sigma0 = mindist
   prevents successvive func evals further apart in x than maxdist
------------------------------------------------------------------------- */

int MinCG::linemin_secant(int &nfunc)
{
  int i,iter;
  double eta,eta_prev,sigma0,sigmamax,alpha,alphadelta,fme,fmax,dsq,e0,tmp;
  double epssq = SECANT_EPS * SECANT_EPS;

  // stopping criterion for secant iterations

  fme = 0.0;
  for (i = 0; i < ndof; i++) fme += h[i]*h[i];
  MPI_Allreduce(&fme,&dsq,1,MPI_DOUBLE,MPI_SUM,world);
  for (i = 0; i < ndof_extra; i++) dsq += hextra[i]*hextra[i];

  // sigma0 = smallest allowed step of mindist
  // sigmamax = largest allowed step (in single iteration) of maxdist

  fme = 0.0;
  for (i = 0; i < ndof; i++) fme = MAX(fme,fabs(h[i]));
  MPI_Allreduce(&fme,&fmax,1,MPI_DOUBLE,MPI_MAX,world);
  for (i = 0; i < ndof_extra; i++) fmax = MAX(fmax,fabs(hextra[i]));

  if (fmax == 0.0) return 1;
  sigma0 = dmin/fmax;
  sigmamax = dmax/fmax;

  // eval func at sigma0
  // test if minstep is already uphill

  e0 = energy;
  for (i = 0; i < ndof; i++) x[i] += sigma0*h[i];
  for (i = 0; i < ndof_extra; i++) xextra[i] += sigma0*hextra[i];
  eng_force();
  nfunc++;

  if (energy >= e0) {
    for (i = 0; i < ndof; i++) x[i] -= sigma0*h[i];
    for (i = 0; i < ndof_extra; i++) xextra[i] -= sigma0*hextra[i];
    eng_force();
    nfunc++;
    return 1;
  }

  // secant iterations
  // alphadelta = new increment to move, alpha = accumulated move

  tmp = 0.0;
  for (i = 0; i < ndof; i++) tmp -= f[i]*h[i];
  MPI_Allreduce(&tmp,&eta_prev,1,MPI_DOUBLE,MPI_SUM,world);
  for (i = 0; i < ndof_extra; i++) eta_prev -= fextra[i]*hextra[i];

  alpha = sigma0;
  alphadelta = -sigma0;

  for (iter = 0; iter < lineiter; iter++) {
    alpha += alphadelta;
    for (i = 0; i < ndof; i++) x[i] += alphadelta*h[i];
    for (i = 0; i < ndof_extra; i++) xextra[i] += alphadelta*hextra[i];
    eng_force();
    nfunc++;

    tmp = 0.0;
    for (i = 0; i < ndof; i++) tmp -= f[i]*h[i];
    MPI_Allreduce(&tmp,&eta,1,MPI_DOUBLE,MPI_SUM,world);
    for (i = 0; i < ndof_extra; i++) eta -= fextra[i]*hextra[i];

    alphadelta *= eta / (eta_prev - eta);
    eta_prev = eta;
    if (alphadelta*alphadelta*dsq <= epssq) break;
    if (fabs(alphadelta) > sigmamax) {
      if (alphadelta > 0.0) alphadelta = sigmamax;
      else alphadelta = -sigmamax;
    }
  }

  // if exited loop on first iteration, func eval was at alpha = 0.0
  // else successful line search

  if (iter == 0) return 1;
  return 0;
}

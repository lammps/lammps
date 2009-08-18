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
#include "min_linesearch.h"
#include "atom.h"
#include "update.h"
#include "neighbor.h"
#include "domain.h"
#include "modify.h"
#include "fix_minimize.h"
#include "pair.h"
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

MinLineSearch::MinLineSearch(LAMMPS *lmp) : Min(lmp)
{
  gextra = hextra = NULL;
  x0extra_atom = gextra_atom = hextra_atom = NULL;
}

/* ---------------------------------------------------------------------- */

MinLineSearch::~MinLineSearch()
{
  delete [] gextra;
  delete [] hextra;
  delete [] x0extra_atom;
  delete [] gextra_atom;
  delete [] hextra_atom;
}

/* ---------------------------------------------------------------------- */

void MinLineSearch::init_style()
{
  if (linestyle == 0) linemin = &MinLineSearch::linemin_backtrack;
  else if (linestyle == 1) linemin = &MinLineSearch::linemin_quadratic;

  delete [] gextra;
  delete [] hextra;
  gextra = hextra = NULL;

  delete [] x0extra_atom;
  delete [] gextra_atom;
  delete [] hextra_atom;
  x0extra_atom = gextra_atom = hextra_atom = NULL;
}

/* ---------------------------------------------------------------------- */

void MinLineSearch::setup_style()
{
  // memory for x0,g,h for atomic dof

  fix_minimize->add_vector(3);
  fix_minimize->add_vector(3);
  fix_minimize->add_vector(3);

  // memory for g,h for extra global dof, fix stores x0

  if (nextra_global) {
    gextra = new double[nextra_global];
    hextra = new double[nextra_global];
  }

  // memory for x0,g,h for extra per-atom dof

  if (nextra_atom) {
    x0extra_atom = new double*[nextra_atom];
    gextra_atom = new double*[nextra_atom];
    hextra_atom = new double*[nextra_atom];

    for (int m = 0; m < nextra_atom; m++) {
      fix_minimize->add_vector(extra_peratom[m]);
      fix_minimize->add_vector(extra_peratom[m]);
      fix_minimize->add_vector(extra_peratom[m]);
    }
  }
}

/* ----------------------------------------------------------------------
   set current vector lengths and pointers
   called after atoms have migrated
------------------------------------------------------------------------- */

void MinLineSearch::reset_vectors()
{
  // atomic dof

  n3 = 3 * atom->nlocal;
  if (n3) x = atom->x[0];
  if (n3) f = atom->f[0];
  x0 = fix_minimize->request_vector(0);
  g = fix_minimize->request_vector(1);
  h = fix_minimize->request_vector(2);

  // extra per-atom dof

  if (nextra_atom) {
    int n = 3;
    for (int m = 0; m < nextra_atom; m++) {
      extra_nlen[m] = extra_peratom[m] * atom->nlocal;
      requestor[m]->min_pointers(&xextra_atom[m],&fextra_atom[m]);
      x0extra_atom[m] = fix_minimize->request_vector(n++);
      gextra_atom[m] = fix_minimize->request_vector(n++);
      hextra_atom[m] = fix_minimize->request_vector(n++);
    }
  }
}

/* ----------------------------------------------------------------------
   line minimization methods
   find minimum-energy starting at x along h direction
   input args:   eoriginal = energy at initial x
   input extra:  n,x,x0,f,h for atomic, extra global, extra per-atom dof
   output args:  return 0 if successful move, non-zero alpha
                 return non-zero if failed
                 alpha = distance moved along h for x at min eng config
		 nfunc = updated counter of eng/force function evals
   output extra: if fail, energy_force() of original x
	         if succeed, energy_force() at x + alpha*h
                 atom->x = coords at new configuration
	         atom->f = force at new configuration
	         ecurrent = energy of new configuration
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   linemin: backtracking line search (Proc 3.1, p 41 in Nocedal and Wright)
   uses no gradient info, but should be very robust
   start at maxdist, backtrack until energy decrease is sufficient
------------------------------------------------------------------------- */

int MinLineSearch::linemin_backtrack(double eoriginal, double &alpha, 
				     int &nfunc)
{
  int i,m,n;
  double fdothall,fdothme,hme,hmax,hmaxall;
  double de_ideal,de;
  double *xatom,*x0atom,*fatom,*hatom;

  // fdothall = projection of search dir along downhill gradient
  // if search direction is not downhill, exit with error

  fdothme = 0.0;
  for (i = 0; i < n3; i++) fdothme += f[i]*h[i];
  if (nextra_atom)
    for (m = 0; m < nextra_atom; m++) {
      fatom = fextra_atom[m];
      hatom = hextra_atom[m];
      n = extra_nlen[m];
      for (i = 0; i < n; i++) fdothme += fatom[i]*hatom[i];
    }
  MPI_Allreduce(&fdothme,&fdothall,1,MPI_DOUBLE,MPI_SUM,world);
  if (nextra_global)
    for (i = 0; i < nextra_global; i++) fdothall += fextra[i]*hextra[i];
  if (output->thermo->normflag) fdothall /= atom->natoms;
  if (fdothall <= 0.0) return DOWNHILL;

  // set alpha so no dof is changed by more than max allowed amount
  // for atom coords, max amount = dmax
  // for extra per-atom dof, max amount = extra_max[]
  // for extra global dof, max amount is set by fix
  // also insure alpha <= ALPHA_MAX
  // else will have to backtrack from huge value when forces are tiny
  // if all search dir components are already 0.0, exit with error

  hme = 0.0;
  for (i = 0; i < n3; i++) hme = MAX(hme,fabs(h[i]));
  MPI_Allreduce(&hme,&hmaxall,1,MPI_DOUBLE,MPI_MAX,world);
  alpha = MIN(ALPHA_MAX,dmax/hmaxall);
  if (nextra_atom)
    for (m = 0; m < nextra_atom; m++) {
      hme = 0.0;
      fatom = fextra_atom[m];
      n = extra_nlen[m];
      for (i = 0; i < n; i++) hme = MAX(hme,fabs(hatom[i]));
      MPI_Allreduce(&hme,&hmax,1,MPI_DOUBLE,MPI_MAX,world);
      alpha = MIN(alpha,extra_max[m]/hmax);
      hmaxall = MAX(hmaxall,hmax);
    }
  if (nextra_global) {
    double alpha_extra = modify->max_alpha(hextra);
    alpha = MIN(alpha,alpha_extra);
    for (i = 0; i < nextra_global; i++)
      hmaxall = MAX(hmaxall,fabs(hextra[i]));
  }
  if (hmaxall == 0.0) return ZEROFORCE;

  // store box and values of all dof at start of linesearch

  fix_minimize->store_box();
  for (i = 0; i < n3; i++) x0[i] = x[i];
  if (nextra_atom)
    for (m = 0; m < nextra_atom; m++) {
      xatom = xextra_atom[m];
      x0atom = x0extra_atom[m];
      n = extra_nlen[m];
      for (i = 0; i < n; i++) x0atom[i] = xatom[i];
    }
  if (nextra_global) modify->min_store();

  // backtrack with alpha until energy decrease is sufficient

  while (1) {
    ecurrent = alpha_step(alpha,1,nfunc);

    // if energy change is better than ideal, exit with success

    de_ideal = -BACKTRACK_SLOPE*alpha*fdothall;
    de = ecurrent - eoriginal;
    if (de <= de_ideal) return 0;

    // reduce alpha

    alpha *= ALPHA_REDUCE;

    // backtracked all the way to 0.0
    // reset to starting point, exit with error

    if (alpha <= 0.0 || de_ideal >= -IDEAL_TOL) {
      ecurrent = alpha_step(0.0,0,nfunc);
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

int MinLineSearch::linemin_quadratic(double eoriginal, double &alpha, 
				     int &nfunc)
{
  int i,m,n;
  double fdothall,fdothme,hme,hmax,hmaxall;
  double de_ideal,de;
  double delfh,engprev,relerr,alphaprev,fhprev,ff,fh,alpha0,fh0,ff0;
  double dot[2],dotall[2];	
  double *xatom,*x0atom,*fatom,*hatom;

  // fdothall = projection of search dir along downhill gradient
  // if search direction is not downhill, exit with error

  fdothme = 0.0;
  for (i = 0; i < n3; i++) fdothme += f[i]*h[i];
  if (nextra_atom)
    for (m = 0; m < nextra_atom; m++) {
      fatom = fextra_atom[m];
      hatom = hextra_atom[m];
      n = extra_nlen[m];
      for (i = 0; i < n; i++) fdothme += fatom[i]*hatom[i];
    }
  MPI_Allreduce(&fdothme,&fdothall,1,MPI_DOUBLE,MPI_SUM,world);
  if (nextra_global)
    for (i = 0; i < nextra_global; i++) fdothall += fextra[i]*hextra[i];
  if (output->thermo->normflag) fdothall /= atom->natoms;
  if (fdothall <= 0.0) return DOWNHILL;

  // set alpha so no dof is changed by more than max allowed amount
  // for atom coords, max amount = dmax
  // for extra per-atom dof, max amount = extra_max[]
  // for extra global dof, max amount is set by fix
  // also insure alpha <= ALPHA_MAX
  // else will have to backtrack from huge value when forces are tiny
  // if all search dir components are already 0.0, exit with error

  hme = 0.0;
  for (i = 0; i < n3; i++) hme = MAX(hme,fabs(h[i]));
  MPI_Allreduce(&hme,&hmaxall,1,MPI_DOUBLE,MPI_MAX,world);
  alpha = MIN(ALPHA_MAX,dmax/hmaxall);
  if (nextra_atom)
    for (m = 0; m < nextra_atom; m++) {
      hme = 0.0;
      fatom = fextra_atom[m];
      n = extra_nlen[m];
      for (i = 0; i < n; i++) hme = MAX(hme,fabs(hatom[i]));
      MPI_Allreduce(&hme,&hmax,1,MPI_DOUBLE,MPI_MAX,world);
      alpha = MIN(alpha,extra_max[m]/hmax);
      hmaxall = MAX(hmaxall,hmax);
    }
  if (nextra_global) {
    double alpha_extra = modify->max_alpha(hextra);
    alpha = MIN(alpha,alpha_extra);
    for (i = 0; i < nextra_global; i++)
      hmaxall = MAX(hmaxall,fabs(hextra[i]));
  }
  if (hmaxall == 0.0) return ZEROFORCE;

  // store box and values of all dof at start of linesearch

  fix_minimize->store_box();
  for (i = 0; i < n3; i++) x0[i] = x[i];
  if (nextra_atom)
    for (m = 0; m < nextra_atom; m++) {
      xatom = xextra_atom[m];
      x0atom = x0extra_atom[m];
      n = extra_nlen[m];
      for (i = 0; i < n; i++) x0atom[i] = xatom[i];
    }
  if (nextra_global) modify->min_store();

  // backtrack with alpha until energy decrease is sufficient
  // or until get to small energy change, then perform quadratic projection

  fhprev = fdothall;
  engprev = eoriginal;
  alphaprev = 0.0;

  while (1) {
    ecurrent = alpha_step(alpha,1,nfunc);

    // compute new fh, alpha, delfh

    dot[0] = dot[1] = 0.0;
    for (i = 0; i < n3; i++) {
      dot[0] += f[i]*f[i];
      dot[1] += f[i]*h[i];
    }
    if (nextra_atom)
      for (m = 0; m < nextra_atom; m++) {
	xatom = xextra_atom[m];
	hatom = hextra_atom[m];
	n = extra_nlen[m];
	for (i = 0; i < n; i++) {
	  dot[0] += fatom[i]*fatom[i];
	  dot[1] += fatom[i]*hatom[i];
	}
      }
    MPI_Allreduce(dot,dotall,2,MPI_DOUBLE,MPI_SUM,world);
    if (nextra_global) {
      for (i = 0; i < nextra_global; i++) {
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
      ecurrent = alpha_step(0.0,0,nfunc);
      return ZEROQUAD;
    }

    // check if ready for quadratic projection, equivalent to secant method
    // alpha0 = projected alpha

    relerr = fabs(1.0+(0.5*alpha*(alpha-alphaprev)*
		       (fh+fhprev)-ecurrent)/engprev);
    alpha0 = alpha - (alpha-alphaprev)*fh/delfh;

    if (relerr <= QUADRATIC_TOL && alpha0 > 0.0) {
      ecurrent = alpha_step(alpha0,1,nfunc);
      return 0;
    }

    // if backtracking energy change is better than ideal, exit with success

    de_ideal = -BACKTRACK_SLOPE*alpha*fdothall;
    de = ecurrent - eoriginal;
    if (de <= de_ideal) return 0;

    // save previous state

    fhprev = fh;
    engprev = ecurrent;
    alphaprev = alpha;

    // reduce alpha

    alpha *= ALPHA_REDUCE;

    // backtracked all the way to 0.0
    // reset to starting point, exit with error

    if (alpha <= 0.0 || de_ideal >= -IDEAL_TOL) {
      ecurrent = alpha_step(0.0,0,nfunc);
      return ZEROALPHA;
    }
  }
}

/* ---------------------------------------------------------------------- */

double MinLineSearch::alpha_step(double alpha, int resetflag, int &nfunc)
{
  int i,n,m;
  double *xatom,*x0atom,*hatom;

  // reset to starting point

  if (nextra_global) modify->min_step(0.0,hextra);
  for (i = 0; i < n3; i++) x[i] = x0[i];
  if (nextra_atom)
      for (m = 0; m < nextra_atom; m++) {
	xatom = xextra_atom[m];
	x0atom = x0extra_atom[m];
	n = extra_nlen[m];
	for (i = 0; i < n; i++) xatom[i] = x0atom[i];
      }
  
  // step forward along h

  if (alpha > 0.0) {
    if (nextra_global) modify->min_step(alpha,hextra);
    for (i = 0; i < n3; i++) x[i] += alpha*h[i];
    if (nextra_atom)
      for (m = 0; m < nextra_atom; m++) {
	xatom = xextra_atom[m];
	hatom = hextra_atom[m];
	n = extra_nlen[m];
	for (i = 0; i < n; i++) xatom[i] += alpha*hatom[i];
      }
  }

  // compute and return new energy

  nfunc++;
  return energy_force(resetflag);
}

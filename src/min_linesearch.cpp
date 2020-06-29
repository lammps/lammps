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
   Contributing author: Asad Hasan (CMU)
                        added forcezero ls
   Sources: Numerical Recipes frprmn routine
            "Conjugate Gradient Method Without the Agonizing Pain" by
            JR Shewchuk, http://www-2.cs.cmu.edu/~jrs/jrspapers.html#cg
------------------------------------------------------------------------- */

#include "min_linesearch.h"
#include <mpi.h>
#include <cmath>
#include "atom.h"
#include "modify.h"
#include "fix_minimize.h"
#include "pair.h"
#include "output.h"
#include "thermo.h"

using namespace LAMMPS_NS;

// ALPHA_MAX = max alpha allowed to avoid long backtracks
// ALPHA_REDUCE = reduction ratio, should be in range [0.5,1)
// BACKTRACK_SLOPE, should be in range (0,0.5]
// QUADRATIC_TOL = tolerance on alpha0, should be in range [0.1,1)
// EMACH = machine accuracy limit of energy changes (1.0e-8)
// EPS_QUAD = tolerance for quadratic projection

#define ALPHA_MAX 1.0
#define ALPHA_REDUCE 0.5
#define BACKTRACK_SLOPE 0.4
#define QUADRATIC_TOL 0.1
//#define EMACH 1.0e-8
#define EMACH 1.0e-8
#define EPS_QUAD 1.0e-28

/* ---------------------------------------------------------------------- */

MinLineSearch::MinLineSearch(LAMMPS *lmp) : Min(lmp)
{
  searchflag = 1;
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

void MinLineSearch::init()
{
  Min::init();

  if (linestyle == 0) linemin = &MinLineSearch::linemin_backtrack;
  else if (linestyle == 1) linemin = &MinLineSearch::linemin_quadratic;
  else if (linestyle == 2) linemin = &MinLineSearch::linemin_forcezero;

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

  nvec = 3 * atom->nlocal;
  if (nvec) xvec = atom->x[0];
  if (nvec) fvec = atom->f[0];
  x0 = fix_minimize->request_vector(0);
  g = fix_minimize->request_vector(1);
  h = fix_minimize->request_vector(2);

  // extra per-atom dof

  if (nextra_atom) {
    int n = 3;
    for (int m = 0; m < nextra_atom; m++) {
      extra_nlen[m] = extra_peratom[m] * atom->nlocal;
      requestor[m]->min_xf_pointers(m,&xextra_atom[m],&fextra_atom[m]);
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
                 update neval counter of eng/force function evaluations
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

int MinLineSearch::linemin_backtrack(double eoriginal, double &alpha)
{
  int i,m,n;
  double fdothall,fdothme,hme,hmax,hmaxall;
  double de_ideal,de;
  double *xatom,*x0atom,*fatom,*hatom;

  // fdothall = projection of search dir along downhill gradient
  // if search direction is not downhill, exit with error

  fdothme = 0.0;
  for (i = 0; i < nvec; i++) fdothme += fvec[i]*h[i];
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
  for (i = 0; i < nvec; i++) hme = MAX(hme,fabs(h[i]));
  MPI_Allreduce(&hme,&hmaxall,1,MPI_DOUBLE,MPI_MAX,world);
  alpha = MIN(ALPHA_MAX,dmax/hmaxall);
  if (nextra_atom)
    for (m = 0; m < nextra_atom; m++) {
      hatom = hextra_atom[m];
      n = extra_nlen[m];
      hme = 0.0;
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
  for (i = 0; i < nvec; i++) x0[i] = xvec[i];
  if (nextra_atom)
    for (m = 0; m < nextra_atom; m++) {
      xatom = xextra_atom[m];
      x0atom = x0extra_atom[m];
      n = extra_nlen[m];
      for (i = 0; i < n; i++) x0atom[i] = xatom[i];
    }
  if (nextra_global) modify->min_store();

  // // important diagnostic: test the gradient against energy
  // double etmp;
  // double alphatmp = alpha*1.0e-4;
  // etmp = alpha_step(alphatmp,1);
  // printf("alpha = %g dele = %g dele_force = %g err = %g\n",
  //        alphatmp,etmp-eoriginal,-alphatmp*fdothall,
  //        etmp-eoriginal+alphatmp*fdothall);
  // alpha_step(0.0,1);

  // backtrack with alpha until energy decrease is sufficient

  while (1) {
    ecurrent = alpha_step(alpha,1);

    // if energy change is better than ideal, exit with success

    de_ideal = -BACKTRACK_SLOPE*alpha*fdothall;
    de = ecurrent - eoriginal;
    if (de <= de_ideal) {
      if (nextra_global) {
        int itmp = modify->min_reset_ref();
        if (itmp) ecurrent = energy_force(1);
      }
      return 0;
    }

    // reduce alpha

    alpha *= ALPHA_REDUCE;

    // backtracked too much
    // reset to starting point
    // if de is positive, exit with error
    // if de is negative, exit with ETOL

    if (alpha <= 0.0 || de_ideal >= -EMACH) {
      ecurrent = alpha_step(0.0,0);
      if (de < 0.0) return ETOL;
      else return ZEROALPHA;
    }
  }
}

/* ----------------------------------------------------------------------
    // linemin: quadratic line search (adapted from Dennis and Schnabel)
    // The objective function is approximated by a quadratic
    // function in alpha, for sufficiently small alpha.
    // This idea is the same as that used in the well-known secant
    // method. However, since the change in the objective function
    // (difference of two finite numbers) is not known as accurately
    // as the gradient (which is close to zero), all the expressions
    // are written in terms of gradients. In this way, we can converge
    // the LAMMPS forces much closer to zero.
    //
    // We know E,Eprev,fh,fhprev. The Taylor series about alpha_prev
    // truncated at the quadratic term is:
    //
    //     E = Eprev - del_alpha*fhprev + (1/2)del_alpha^2*Hprev
    //
    // and
    //
    //     fh = fhprev - del_alpha*Hprev
    //
    // where del_alpha = alpha-alpha_prev
    //
    // We solve these two equations for Hprev and E=Esolve, giving:
    //
    //     Esolve = Eprev - del_alpha*(f+fprev)/2
    //
    // We define relerr to be:
    //
    //      relerr = |(Esolve-E)/Eprev|
    //             = |1.0 - (0.5*del_alpha*(f+fprev)+E)/Eprev|
    //
    // If this is accurate to within a reasonable tolerance, then
    // we go ahead and use a secant step to fh = 0:
    //
    //      alpha0 = alpha - (alpha-alphaprev)*fh/delfh;
    //
------------------------------------------------------------------------- */

int MinLineSearch::linemin_quadratic(double eoriginal, double &alpha)
{
  int i,m,n;
  double fdothall,fdothme,hme,hmax,hmaxall;
  double de_ideal,de;
  double delfh,engprev,relerr,alphaprev,fhprev,ff,fh,alpha0;
  double dot[2],dotall[2];
  double *xatom,*x0atom,*fatom,*hatom;
  double alphamax;

  // fdothall = projection of search dir along downhill gradient
  // if search direction is not downhill, exit with error

  fdothme = 0.0;
  for (i = 0; i < nvec; i++) fdothme += fvec[i]*h[i];
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

  // set alphamax so no dof is changed by more than max allowed amount
  // for atom coords, max amount = dmax
  // for extra per-atom dof, max amount = extra_max[]
  // for extra global dof, max amount is set by fix
  // also insure alphamax <= ALPHA_MAX
  // else will have to backtrack from huge value when forces are tiny
  // if all search dir components are already 0.0, exit with error

  hme = 0.0;
  for (i = 0; i < nvec; i++) hme = MAX(hme,fabs(h[i]));
  MPI_Allreduce(&hme,&hmaxall,1,MPI_DOUBLE,MPI_MAX,world);
  alphamax = MIN(ALPHA_MAX,dmax/hmaxall);
  if (nextra_atom)
    for (m = 0; m < nextra_atom; m++) {
      hatom = hextra_atom[m];
      n = extra_nlen[m];
      hme = 0.0;
      for (i = 0; i < n; i++) hme = MAX(hme,fabs(hatom[i]));
      MPI_Allreduce(&hme,&hmax,1,MPI_DOUBLE,MPI_MAX,world);
      alphamax = MIN(alphamax,extra_max[m]/hmax);
      hmaxall = MAX(hmaxall,hmax);
    }
  if (nextra_global) {
    double alpha_extra = modify->max_alpha(hextra);
    alphamax = MIN(alphamax,alpha_extra);
    for (i = 0; i < nextra_global; i++)
      hmaxall = MAX(hmaxall,fabs(hextra[i]));
  }

  if (hmaxall == 0.0) return ZEROFORCE;

  // store box and values of all dof at start of linesearch

  fix_minimize->store_box();
  for (i = 0; i < nvec; i++) x0[i] = xvec[i];
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

  alpha = alphamax;
  fhprev = fdothall;
  engprev = eoriginal;
  alphaprev = 0.0;

  // // important diagnostic: test the gradient against energy
  // double etmp;
  // double alphatmp = alphamax*1.0e-4;
  // etmp = alpha_step(alphatmp,1);
  // printf("alpha = %g dele = %g dele_force = %g err = %g\n",
  //        alphatmp,etmp-eoriginal,-alphatmp*fdothall,
  //        etmp-eoriginal+alphatmp*fdothall);
  // alpha_step(0.0,1);

  while (1) {
    ecurrent = alpha_step(alpha,1);

    // compute new fh, alpha, delfh

    dot[0] = dot[1] = 0.0;
    for (i = 0; i < nvec; i++) {
      dot[0] += fvec[i]*fvec[i];
      dot[1] += fvec[i]*h[i];
    }
    if (nextra_atom)
      for (m = 0; m < nextra_atom; m++) {
        fatom = fextra_atom[m];
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
      ecurrent = alpha_step(0.0,0);
      return ZEROQUAD;
    }

    // Check if ready for quadratic projection, equivalent to secant method
    // alpha0 = projected alpha

    relerr = fabs(1.0-(0.5*(alpha-alphaprev)*(fh+fhprev)+ecurrent)/engprev);
    alpha0 = alpha - (alpha-alphaprev)*fh/delfh;

    if (relerr <= QUADRATIC_TOL && alpha0 > 0.0 && alpha0 < alphamax) {
      ecurrent = alpha_step(alpha0,1);
      if (ecurrent - eoriginal < EMACH) {
        if (nextra_global) {
          int itmp = modify->min_reset_ref();
          if (itmp) ecurrent = energy_force(1);
        }
        return 0;
      }
    }

    // if backtracking energy change is better than ideal, exit with success

    de_ideal = -BACKTRACK_SLOPE*alpha*fdothall;
    de = ecurrent - eoriginal;

    if (de <= de_ideal) {
      if (nextra_global) {
        int itmp = modify->min_reset_ref();
        if (itmp) ecurrent = energy_force(1);
      }
      return 0;
    }

    // save previous state

    fhprev = fh;
    engprev = ecurrent;
    alphaprev = alpha;

    // reduce alpha

    alpha *= ALPHA_REDUCE;

    // backtracked all the way to 0.0
    // reset to starting point, exit with error

    if (alpha <= 0.0 || de_ideal >= -EMACH) {
      ecurrent = alpha_step(0.0,0);
      return ZEROALPHA;
    }
  }
}

/* ----------------------------------------------------------------------

forcezero linesearch method - seeks a zero of force in a robust manner.
    (motivated by a line minimization routine of f77 DYNAMO code)

central idea:
  In each linesearch we attempt to converge to a zero of force
  (usual case) or reduces forces (worst case).
  Energy does not play any role in the search procedure,
  except we ensure that it doesn't increase.

pseudo code:
  i)  Fix an alpha max:
        // also account for nextra atom & global
        alpha_max <= dmax/hmaxall

  ii) Initialize:
        fhCurr = current_force.dot.search_direction
        fhoriginal = fhCurr
        // try decreasing the energy to 1/10 of initial
        alpha_init = 0.1*fabs(eoriginal)/fhCurr;
        // initial alpha is smaller than alpha_max
        alpha_del = MIN(alpha_init, 0.5*alpha_max);
        alpha = 0.0
  iii) Loop:
        backtrack = false
        alpha += alpha_del
        if (alpha > alpha_max):
           // we have done enough in the search space
           EXIT with success

        Step with the new alpha
        Compute:
           current energy and 'fhCurr'
           de = ecurrent - eprev

           // ZERO_ENERGY = 1e-12, is max allowed energy increase
           if (de > ZERO_ENERGY):
              bactrack = true

           // GRAD_TOL = 0.1
           if ( (not backtrack) && (fabs(fhCurr/fh0) <= GRAD_TOL) ):
              // forces sufficiently reduced without energy increase
              EXIT with success

           // projected force changed sign but didn't become small enough
           if ( fhCurr < 0):
              backtrack = true

           if (bactrack):
              // forces along search direction changed sign
              if (fhCurr < 0):
                 Get alpha_del by solving for zero
                    of force (1D Newton's Method)
              else:
                 // force didn't change sign but only energy increased,
                 // we overshot a minimum which is very close to a
                 // maximum (or there is an inflection point)

                 // New alpha_del should be much smaller
                 // ALPHA_FACT = 0.1
                 alpha_del *= ALPHA_FACT

                 // Check to see if new 'alpha_del' isn't too small
                 if (alpha_del < MIN_ALPHA):
                    EXIT with failure("linesearch alpha is zero")

               Undo the step of alpha.

           // continue the loop with a new alpha_del
           else:
              Get new alpha_del by linearizing force and solving for its zero

 ---------------------------------------------------------------------- */

int MinLineSearch::linemin_forcezero(double eoriginal, double &alpha)
{
  int i,m,n;
  double fdothall,fdothme,hme,hmax,hmaxall;
  double de;
  double *xatom,*x0atom,*fatom,*hatom;

  double alpha_max, alpha_init, alpha_del;
  // projection of: force on itself, current force on search direction,
  double ffCurr, fhCurr;
  // previous force on search direction, initial force on search direction
  double fhPrev, fhoriginal;
  // current energy, previous energy
  double engCurr, engPrev;
  bool backtrack;

  // hardcoded constants

  // factor by which alpha is reduced when backtracking
  double ALPHA_FACT = 0.1;
  // maximum amount by which we'll permit energy increase
  double ZERO_ENERGY = 1e-12;
  // fraction to which we want to reduce the directional derivative
  double GRAD_TOL = 0.1;
  // largest alpha increment which will trigger a failed_linesearch
  double MIN_ALPHA_FAC = 1e-14;
  double LIMIT_BOOST = 4.0;

  // fdothall = projection of search dir along downhill gradient
  // if search direction is not downhill, exit with error

  fdothme = 0.0;
  for (i = 0; i < nvec; i++) fdothme += fvec[i]*h[i];
  if (nextra_atom)
    for (m = 0; m < nextra_atom; m++) {
      fatom = fextra_atom[m];
      hatom = hextra_atom[m];
      n = extra_nlen[m];
      for (i = 0; i < n; i++)
        fdothme += fatom[i]*hatom[i];
    }

  MPI_Allreduce(&fdothme,&fdothall,1,MPI_DOUBLE,MPI_SUM,world);
  if (nextra_global)
    for (i = 0; i < nextra_global; i++)
      fdothall += fextra[i]*hextra[i];
  if (output->thermo->normflag) fdothall /= atom->natoms;
  if (fdothall <= 0.0) return DOWNHILL;

  // set alpha so no dof is changed by more than max allowed amount
  // for atom coords, max amount = dmax
  // for extra per-atom dof, max amount = extra_max[]
  // for extra global dof, max amount is set by fix

  // also insure alpha <= ALPHA_MAX else will have
  // to backtrack from huge value when forces are tiny

  // if all search dir components are already 0.0, exit with error

  hme = 0.0;
  for (i = 0; i < nvec; i++)
    hme = MAX(hme,fabs(h[i]));

  MPI_Allreduce(&hme,&hmaxall,1,MPI_DOUBLE,MPI_MAX,world);
  alpha_max = dmax/hmaxall;
  if (nextra_atom)
    for (m = 0; m < nextra_atom; m++) {
      hatom = hextra_atom[m];
      n = extra_nlen[m];
      hme = 0.0;
      for (i = 0; i < n; i++) hme = MAX(hme,fabs(hatom[i]));
      MPI_Allreduce(&hme,&hmax,1,MPI_DOUBLE,MPI_MAX,world);
      alpha_max = MIN(alpha_max,extra_max[m]/hmax);
      hmaxall = MAX(hmaxall,hmax);
    }

  if (nextra_global) {
    double alpha_extra = modify->max_alpha(hextra);
    alpha_max = MIN(alpha_max, alpha_extra);
    for (i = 0; i < nextra_global; i++)
      hmaxall = MAX(hmaxall,fabs(hextra[i]));
  }

  if (hmaxall == 0.0) return ZEROFORCE;

  // store box and values of all dof at start of linesearch

  fix_minimize->store_box();

  for (i = 0; i < nvec; i++) x0[i] = xvec[i];
  if (nextra_atom)
    for (m = 0; m < nextra_atom; m++) {
      xatom = xextra_atom[m];
      x0atom = x0extra_atom[m];
      n = extra_nlen[m];
      for (i = 0; i < n; i++) x0atom[i] = xatom[i];
    }

  if (nextra_global) modify->min_store();

  // initialize important variables before main linesearch loop

  ffCurr = 0.0;
  fhCurr = fdothall;
  fhoriginal = fhCurr;
  engCurr = eoriginal;

  // stores energy difference due to the current move

  de = 0.0;

  // choosing the initial alpha that we'll use
  // rough estimate that'll decrease energy to 1/10

  alpha_init = 0.1*fabs(eoriginal)/fdothall;

  // initialize aplha to 0.0

  alpha = 0.0;

  // compute increment to alpha, ensure that we
  // don't take the largest allowed alpha
  // first alpha that will actually apply

  alpha_del = MIN(alpha_init,0.5*alpha_max);

  // main linesearch loop

  while (1) {
    backtrack = false;
    fhPrev = fhCurr;
    engPrev = engCurr;

    // apply the increment to alpha, but first
    // check whether we are still in allowed search space

    alpha += alpha_del;
    if (alpha > alpha_max) {

      // undo the increment

      alpha -= alpha_del;
      if (nextra_global) {
        int itmp = modify->min_reset_ref();
        if (itmp) ecurrent = energy_force(1);
      }

      // exit linesearch with success: have done
      // enough in allowed search space

      return 0;
    }

    // move the system

    // '1' updates coordinates of atoms which cross PBC

    engCurr = alpha_step(alpha,1);
    ecurrent = engCurr;

    // compute the new directional derivative and also f_dot_f

    fhCurr = compute_dir_deriv(ffCurr);

    // energy change

    de = engCurr - engPrev;

    // if the function value increases measurably,
    // then we have to reduce alpha

    if (de >= ZERO_ENERGY)
      backtrack = true;

    // check if the directional derivative has sufficiently decreased
    // NOTE: the fabs is essential here

    if ((!backtrack) && (fabs(fhCurr/fhoriginal) <= GRAD_TOL)) {
      if (nextra_global) {
        int itmp = modify->min_reset_ref();
        if (itmp) ecurrent = energy_force(1);
      }

      // we are done

      return 0;
    }

    // check if the directional derivative changed sign
    // but it's not small: we overshot the minima -- BACKTRACK

    if (fhCurr < 0.0)
      backtrack = true;

    // backtrack by undoing step and choosing a new alpha

    if (backtrack) {

      // move back

      alpha -= alpha_del;

      // choose new alpha
      // if the force changed sign, linearize force and
      // solve for new alpha_del

      if (fhCurr < 0.0)
        alpha_del *= fhPrev/(fhPrev - fhCurr);
      else

    // force didn't change sign but only energy increased,
    // we overshot a minimum which is very close to a maxima
    // (or there is an inflection point)
    // new alpha_del should be much smaller

        alpha_del *= ALPHA_FACT;

      // since we moved back ...

      engCurr = engPrev;
      ecurrent = engCurr;
      fhCurr = fhPrev;

      // if new move is too small then we have failed;
      // exit with 'failed_linesearch'

      if (hmaxall*alpha_del <= MIN_ALPHA_FAC) {

        // undo all line minization moves

        engCurr = alpha_step(0.0,1);
        ecurrent= engCurr;
        return ZEROALPHA;
      }

    } else {

      // get a new alpha by linearizing force and start over

      double boostFactor = LIMIT_BOOST;

      // avoids problems near an energy inflection point

      if (fhPrev > fhCurr)
        boostFactor = fhCurr/(fhPrev - fhCurr);

      // don't want to boost too much

      boostFactor = MIN(boostFactor, LIMIT_BOOST);
      alpha_del *= boostFactor;
    }
  }
}

/* ---------------------------------------------------------------------- */

double MinLineSearch::alpha_step(double alpha, int resetflag)
{
  int i,n,m;
  double *xatom,*x0atom,*hatom;

  // reset to starting point

  if (nextra_global) modify->min_step(0.0,hextra);
  for (i = 0; i < nvec; i++) xvec[i] = x0[i];
  if (nextra_atom)
    for (m = 0; m < nextra_atom; m++) {
      xatom = xextra_atom[m];
      x0atom = x0extra_atom[m];
      n = extra_nlen[m];
      for (i = 0; i < n; i++) xatom[i] = x0atom[i];
      requestor[m]->min_x_set(m);
    }

  // step forward along h

  if (alpha > 0.0) {
    if (nextra_global) modify->min_step(alpha,hextra);
    for (i = 0; i < nvec; i++) xvec[i] += alpha*h[i];
    if (nextra_atom)
      for (m = 0; m < nextra_atom; m++) {
        xatom = xextra_atom[m];
        hatom = hextra_atom[m];
        n = extra_nlen[m];
        for (i = 0; i < n; i++) xatom[i] += alpha*hatom[i];
        requestor[m]->min_x_set(m);
      }
  }

  // compute and return new energy

  neval++;
  return energy_force(resetflag);
}

/* ---------------------------------------------------------------------- */

// compute projection of force on: itself and the search direction

double MinLineSearch::compute_dir_deriv(double &ff)
{
   int i,m,n;
   double *hatom, *fatom;
   double dot[2],dotall[2];
   double fh;

   // compute new fh, alpha, delfh

    dot[0] = dot[1] = 0.0;
    for (i = 0; i < nvec; i++) {
      dot[0] += fvec[i]*fvec[i];
      dot[1] += fvec[i]*h[i];
    }

    if (nextra_atom)
      for (m = 0; m < nextra_atom; m++) {
        fatom = fextra_atom[m];
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

    return fh;
}

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
   Contributing author: Stan Moore (SNL)
------------------------------------------------------------------------- */

#include "min_linesearch_kokkos.h"
#include <mpi.h>
#include <cmath>
#include "atom_kokkos.h"
#include "modify.h"
#include "fix_minimize_kokkos.h"
#include "pair.h"
#include "output.h"
#include "thermo.h"
#include "error.h"
#include "atom_masks.h"

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

MinLineSearchKokkos::MinLineSearchKokkos(LAMMPS *lmp) : MinKokkos(lmp)
{
  searchflag = 1;
  atomKK = (AtomKokkos *) atom;
}

/* ---------------------------------------------------------------------- */

MinLineSearchKokkos::~MinLineSearchKokkos()
{

}

/* ---------------------------------------------------------------------- */

void MinLineSearchKokkos::init()
{
  MinKokkos::init();

  if (linestyle == 1) linemin = &MinLineSearchKokkos::linemin_quadratic;
  else error->all(FLERR,"Kokkos minimize only supports the 'min_modify line "
   "quadratic' option");
}

/* ---------------------------------------------------------------------- */

void MinLineSearchKokkos::setup_style()
{
  // memory for x0,g,h for atomic dof

  fix_minimize_kk->add_vector_kokkos();
  fix_minimize_kk->add_vector_kokkos();
  fix_minimize_kk->add_vector_kokkos();
}

/* ----------------------------------------------------------------------
   set current vector lengths and pointers
   called after atoms have migrated
------------------------------------------------------------------------- */

void MinLineSearchKokkos::reset_vectors()
{
  // atomic dof

  nvec = 3 * atom->nlocal;
  atomKK->sync(Device,F_MASK|X_MASK);
  auto d_x = atomKK->k_x.d_view;
  auto d_f = atomKK->k_f.d_view;

  if (nvec) xvec = DAT::t_ffloat_1d(d_x.data(),d_x.size());
  if (nvec) fvec = DAT::t_ffloat_1d(d_f.data(),d_f.size());
  x0 = fix_minimize_kk->request_vector_kokkos(0);
  g = fix_minimize_kk->request_vector_kokkos(1);
  h = fix_minimize_kk->request_vector_kokkos(2);

  auto h_fvec = Kokkos::create_mirror_view(fvec);
  Kokkos::deep_copy(h_fvec,fvec);
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

int MinLineSearchKokkos::linemin_quadratic(double eoriginal, double &alpha)
{
  double fdothall,fdothme,hme,hmaxall;
  double de_ideal,de;
  double delfh,engprev,relerr,alphaprev,fhprev,ff,fh,alpha0;
  double dot[2],dotall[2];
  double alphamax;

  fix_minimize_kk->k_vectors.sync<LMPDeviceType>();
  fix_minimize_kk->k_vectors.modify<LMPDeviceType>();

  // fdothall = projection of search dir along downhill gradient
  // if search direction is not downhill, exit with error

  fdothme = 0.0;
  {
    // local variables for lambda capture

    auto l_fvec = fvec;
    auto l_h = h;

    Kokkos::parallel_reduce(nvec, LAMMPS_LAMBDA(const int& i, double& fdothme) {
      fdothme += l_fvec[i]*l_h[i];
    },fdothme);
  }
  MPI_Allreduce(&fdothme,&fdothall,1,MPI_DOUBLE,MPI_SUM,world);
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
  {
    // local variables for lambda capture

    auto l_h = h;

    Kokkos::parallel_reduce(nvec, LAMMPS_LAMBDA(const int& i, double& hme) {
      hme = MAX(hme,fabs(l_h[i]));
    },Kokkos::Max<double>(hme));
  }
  MPI_Allreduce(&hme,&hmaxall,1,MPI_DOUBLE,MPI_MAX,world);
  alphamax = MIN(ALPHA_MAX,dmax/hmaxall);

  if (hmaxall == 0.0) return ZEROFORCE;

  // store box and values of all dof at start of linesearch

  {
    // local variables for lambda capture

    auto l_xvec = xvec;
    auto l_x0 = x0;

    fix_minimize_kk->store_box();
    Kokkos::parallel_for(nvec, LAMMPS_LAMBDA(const int& i) {
      l_x0[i] = l_xvec[i];
    });
  }

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

    s_double2 sdot;
    {
      // local variables for lambda capture

      auto l_fvec = fvec;
      auto l_h = h;

      Kokkos::parallel_reduce(nvec, LAMMPS_LAMBDA(const int& i, s_double2& sdot) {
        sdot.d0 += l_fvec[i]*l_fvec[i];
        sdot.d1 += l_fvec[i]*l_h[i];
      },sdot);
    }
    dot[0] = sdot.d0;
    dot[1] = sdot.d1;

    MPI_Allreduce(dot,dotall,2,MPI_DOUBLE,MPI_SUM,world);
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
        return 0;
      }
    }

    // if backtracking energy change is better than ideal, exit with success

    de_ideal = -BACKTRACK_SLOPE*alpha*fdothall;
    de = ecurrent - eoriginal;

    if (de <= de_ideal)
      return 0;

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

/* ---------------------------------------------------------------------- */

double MinLineSearchKokkos::alpha_step(double alpha, int resetflag)
{
  // reset to starting point

  atomKK->k_x.clear_sync_state(); // ignore if host positions since device
                                  //  positions will be reset below
  {
    // local variables for lambda capture

    auto l_xvec = xvec;
    auto l_x0 = x0;

    Kokkos::parallel_for(nvec, LAMMPS_LAMBDA(const int& i) {
      l_xvec[i] = l_x0[i];
    });
  }

  // step forward along h

  if (alpha > 0.0) {
    // local variables for lambda capture

    auto l_xvec = xvec;
    auto l_h = h;

    Kokkos::parallel_for(nvec, LAMMPS_LAMBDA(const int& i) {
      l_xvec[i] += alpha*l_h[i];
    });
  }

  atomKK->modified(Device,X_MASK);

  // compute and return new energy

  neval++;
  return energy_force(resetflag);
}

/* ---------------------------------------------------------------------- */

// compute projection of force on: itself and the search direction

double MinLineSearchKokkos::compute_dir_deriv(double &ff)
{
  double dot[2],dotall[2];
  double fh;

  // compute new fh, alpha, delfh

  s_double2 sdot;
  {
    // local variables for lambda capture

    auto l_fvec = fvec;
    auto l_h = h;

    Kokkos::parallel_reduce(nvec, LAMMPS_LAMBDA(const int& i, s_double2& sdot) {
      sdot.d0 += l_fvec[i]*l_fvec[i];
      sdot.d1 += l_fvec[i]*l_h[i];
    },sdot);
  }
  dot[0] = sdot.d0;
  dot[1] = sdot.d1;

  MPI_Allreduce(dot,dotall,2,MPI_DOUBLE,MPI_SUM,world);

  ff = dotall[0];
  fh = dotall[1];
  if (output->thermo->normflag) {
    ff /= atom->natoms;
    fh /= atom->natoms;
  }

  return fh;
}

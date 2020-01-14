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

/* ------------------------------------------------------------------------
   Contributing authors: Aleksei Ivanov (University of Iceland)
                         Julien Tranchida (SNL)

   Please cite the related publication:
   Ivanov, A. V., Uzdin, V. M., & Jónsson, H. (2019). Fast and Robust
   Algorithm for the Minimisation of the Energy of Spin Systems. arXiv
   preprint arXiv:1904.02669.
------------------------------------------------------------------------- */

#include <mpi.h>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include "min_spin_cg.h"
#include "universe.h"
#include "atom.h"
#include "citeme.h"
#include "comm.h"
#include "force.h"
#include "update.h"
#include "output.h"
#include "timer.h"
#include "error.h"
#include "memory.h"
#include "modify.h"
#include "math_special.h"
#include "math_const.h"
#include "universe.h"

using namespace LAMMPS_NS;
using namespace MathConst;

static const char cite_minstyle_spin_cg[] =
  "min_style spin/cg command:\n\n"
  "@article{ivanov2019fast,\n"
  "title={Fast and Robust Algorithm for the Minimisation of the Energy of "
  "Spin Systems},\n"
  "author={Ivanov, A. V and Uzdin, V. M. and J{\'o}nsson, H.},\n"
  "journal={arXiv preprint arXiv:1904.02669},\n"
  "year={2019}\n"
  "}\n\n";

// EPS_ENERGY = minimum normalization for energy tolerance

#define EPS_ENERGY 1.0e-8

#define DELAYSTEP 5

/* ---------------------------------------------------------------------- */

MinSpinCG::MinSpinCG(LAMMPS *lmp) :
  Min(lmp), g_old(NULL), g_cur(NULL), p_s(NULL), sp_copy(NULL)
{
  if (lmp->citeme) lmp->citeme->add(cite_minstyle_spin_cg);
  nlocal_max = 0;

  // nreplica = number of partitions
  // ireplica = which world I am in universe

  nreplica = universe->nworlds;
  ireplica = universe->iworld;
  use_line_search = 0;  // no line search as default option for CG

  discrete_factor = 10.0;
}

/* ---------------------------------------------------------------------- */

MinSpinCG::~MinSpinCG()
{
  memory->destroy(g_old);
  memory->destroy(g_cur);
  memory->destroy(p_s);
  if (use_line_search)
    memory->destroy(sp_copy);
}

/* ---------------------------------------------------------------------- */

void MinSpinCG::init()
{
  local_iter = 0;
  der_e_cur = 0.0;
  der_e_pr = 0.0;

  Min::init();

  // warning if line_search combined to gneb

  if ((nreplica >= 1) && (linestyle != 4) && (comm->me == 0))
    error->warning(FLERR,"Line search incompatible gneb");

  // set back use_line_search to 0 if more than one replica

  if (linestyle == 3 && nreplica == 1){
    use_line_search = 1;
  }
  else{
    use_line_search = 0;
  }

  dts = dt = update->dt;
  last_negative = update->ntimestep;

  // allocate tables

  nlocal_max = atom->nlocal;
  memory->grow(g_old,3*nlocal_max,"min/spin/cg:g_old");
  memory->grow(g_cur,3*nlocal_max,"min/spin/cg:g_cur");
  memory->grow(p_s,3*nlocal_max,"min/spin/cg:p_s");
  if (use_line_search)
    memory->grow(sp_copy,nlocal_max,3,"min/spin/cg:sp_copy");
}

/* ---------------------------------------------------------------------- */

void MinSpinCG::setup_style()
{
  double **v = atom->v;
  int nlocal = atom->nlocal;

  // check if the atom/spin style is defined

  if (!atom->sp_flag)
    error->all(FLERR,"min spin/cg requires atom/spin style");

  for (int i = 0; i < nlocal; i++)
    v[i][0] = v[i][1] = v[i][2] = 0.0;
}

/* ---------------------------------------------------------------------- */

int MinSpinCG::modify_param(int narg, char **arg)
{
  if (strcmp(arg[0],"discrete_factor") == 0) {
    if (narg < 2) error->all(FLERR,"Illegal fix_modify command");
    discrete_factor = force->numeric(FLERR,arg[1]);
    return 2;
  }
  return 0;
}

/* ----------------------------------------------------------------------
   set current vector lengths and pointers
   called after atoms have migrated
------------------------------------------------------------------------- */

void MinSpinCG::reset_vectors()
{
  // atomic dof

  // size sp is 4N vector
  nvec = 4 * atom->nlocal;
  if (nvec) spvec = atom->sp[0];

  nvec = 3 * atom->nlocal;
  if (nvec) fmvec = atom->fm[0];

  if (nvec) xvec = atom->x[0];
  if (nvec) fvec = atom->f[0];
}

/* ----------------------------------------------------------------------
   minimization via orthogonal spin optimisation
------------------------------------------------------------------------- */

int MinSpinCG::iterate(int maxiter)
{
  int nlocal = atom->nlocal;
  bigint ntimestep;
  double fmdotfm,fmsq;
  int flag, flagall;
  double **sp = atom->sp;
  double der_e_cur_tmp = 0.0;

  if (nlocal_max < nlocal) {
    local_iter = 0;
    nlocal_max = nlocal;
    memory->grow(g_old,3*nlocal_max,"min/spin/cg:g_old");
    memory->grow(g_cur,3*nlocal_max,"min/spin/cg:g_cur");
    memory->grow(p_s,3*nlocal_max,"min/spin/cg:p_s");
    if (use_line_search)
      memory->grow(sp_copy,nlocal_max,3,"min/spin/cg:sp_copy");
  }

  for (int iter = 0; iter < maxiter; iter++) {

    if (timer->check_timeout(niter))
      return TIMEOUT;

    ntimestep = ++update->ntimestep;
    niter++;

    // optimize timestep accross processes / replicas
    // need a force calculation for timestep optimization

    if (use_line_search) {

      // here we need to do line search
      if (local_iter == 0){
        calc_gradient();
      }

      calc_search_direction();
      der_e_cur = 0.0;
      for (int i = 0; i < 3 * nlocal; i++)
        der_e_cur += g_cur[i] * p_s[i];
      MPI_Allreduce(&der_e_cur,&der_e_cur_tmp,1,MPI_DOUBLE,MPI_SUM,world);
      der_e_cur = der_e_cur_tmp;
      if (update->multireplica == 1) {
        MPI_Allreduce(&der_e_cur_tmp,&der_e_cur,1,MPI_DOUBLE,MPI_SUM,universe->uworld);
      }
      for (int i = 0; i < nlocal; i++)
        for (int j = 0; j < 3; j++)
          sp_copy[i][j] = sp[i][j];

      eprevious = ecurrent;
      der_e_pr = der_e_cur;
      calc_and_make_step(0.0, 1.0, 0);
    }
    else{

      // here we don't do line search
      // if gneb calc., nreplica > 1
      // then calculate gradients and advance spins
      // of intermediate replicas only
      calc_gradient();
      calc_search_direction();
      advance_spins();
      neval++;
      eprevious = ecurrent;
      ecurrent = energy_force(0);
    }

    // energy tolerance criterion
    // only check after DELAYSTEP elapsed since velocties reset to 0
    // sync across replicas if running multi-replica minimization

    if (update->etol > 0.0 && ntimestep-last_negative > DELAYSTEP) {
      if (update->multireplica == 0) {
        if (fabs(ecurrent-eprevious) <
            update->etol * 0.5*(fabs(ecurrent) + fabs(eprevious) + EPS_ENERGY))
          return ETOL;
      } else {
        if (fabs(ecurrent-eprevious) <
            update->etol * 0.5*(fabs(ecurrent) + fabs(eprevious) + EPS_ENERGY))
          flag = 0;
        else flag = 1;
        MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_SUM,universe->uworld);
        if (flagall == 0) return ETOL;
      }
    }

    // magnetic torque tolerance criterion
    // sync across replicas if running multi-replica minimization

    fmdotfm = fmsq = 0.0;
    if (update->ftol > 0.0) {
      if (normstyle == MAX) fmsq = max_torque();        // max torque norm
      else if (normstyle == INF) fmsq = inf_torque();   // inf torque norm
      else if (normstyle == TWO) fmsq = total_torque(); // Euclidean torque 2-norm
      else error->all(FLERR,"Illegal min_modify command");
      fmdotfm = fmsq*fmsq;
      if (update->multireplica == 0) {
        if (fmdotfm < update->ftol*update->ftol) return FTOL;
      } else {
        if (fmdotfm < update->ftol*update->ftol) flag = 0;
        else flag = 1;
        MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_SUM,universe->uworld);
        if (flagall == 0) return FTOL;
      }
    }

    // output for thermo, dump, restart files

    if (output->next == ntimestep) {
      timer->stamp();
      output->write(ntimestep);
      timer->stamp(Timer::OUTPUT);
    }
  }

  return MAXITER;
}

/* ----------------------------------------------------------------------
   calculate gradients
---------------------------------------------------------------------- */

void MinSpinCG::calc_gradient()
{
  int nlocal = atom->nlocal;
  double **sp = atom->sp;
  double **fm = atom->fm;
  double hbar = force->hplanck/MY_2PI;
  double factor;

  if (use_line_search)
    factor = hbar;
  else factor = evaluate_dt();

  // loop on all spins on proc.

  for (int i = 0; i < nlocal; i++) {
    g_cur[3 * i + 0] = (fm[i][0]*sp[i][1] - fm[i][1]*sp[i][0]) * factor;
    g_cur[3 * i + 1] = -(fm[i][2]*sp[i][0] - fm[i][0]*sp[i][2]) * factor;
    g_cur[3 * i + 2] = (fm[i][1]*sp[i][2] - fm[i][2]*sp[i][1]) * factor;
  }
}

/* ----------------------------------------------------------------------
   search direction:
   The Fletcher-Reeves conj. grad. method
   See Jorge Nocedal and Stephen J. Wright 'Numerical
   Optimization' Second Edition, 2006 (p. 121)
---------------------------------------------------------------------- */

void MinSpinCG::calc_search_direction()
{
  int nlocal = atom->nlocal;
  double g2old = 0.0;
  double g2 = 0.0;
  double beta = 0.0;

  double g2_global = 0.0;
  double g2old_global = 0.0;

  double factor = 1.0;

  // for multiple replica do not move end points
  if (nreplica > 1)
    if (ireplica == 0 || ireplica == nreplica - 1)
      factor = 0.0;


  if (local_iter == 0 || local_iter % 5 == 0){  // steepest descent direction
    for (int i = 0; i < 3 * nlocal; i++) {
      p_s[i] = -g_cur[i] * factor;
      g_old[i] = g_cur[i] * factor;
    }
  } else {                              // conjugate direction
    for (int i = 0; i < 3 * nlocal; i++) {
      g2old += g_old[i] * g_old[i];
      g2 += g_cur[i] * g_cur[i];
    }

    // now we need to collect/broadcast beta on this replica
    // need to check what is beta for GNEB

    MPI_Allreduce(&g2,&g2_global,1,MPI_DOUBLE,MPI_SUM,world);
    MPI_Allreduce(&g2old,&g2old_global,1,MPI_DOUBLE,MPI_SUM,world);

    // Sum over all replicas. Good for GNEB.

    if (nreplica > 1) {
      g2 = g2_global * factor;
      g2old = g2old_global * factor;
      MPI_Allreduce(&g2,&g2_global,1,MPI_DOUBLE,MPI_SUM,universe->uworld);
      MPI_Allreduce(&g2old,&g2old_global,1,MPI_DOUBLE,MPI_SUM,universe->uworld);
    }
    if (fabs(g2_global) < 1.0e-60) beta = 0.0;
    else beta = g2_global / g2old_global;

    // calculate conjugate direction

    for (int i = 0; i < 3 * nlocal; i++) {
      p_s[i] = (beta * p_s[i] - g_cur[i]) * factor;
      g_old[i] = g_cur[i] * factor;
    }
  }

  local_iter++;
}

/* ----------------------------------------------------------------------
   rotation of spins along the search direction
---------------------------------------------------------------------- */

void MinSpinCG::advance_spins()
{
  int nlocal = atom->nlocal;
  double **sp = atom->sp;
  double rot_mat[9];    // exponential of matrix made of search direction
  double s_new[3];

  // loop on all spins on proc.

  for (int i = 0; i < nlocal; i++) {
    rodrigues_rotation(p_s + 3 * i, rot_mat);

    // rotate spins

    vm3(rot_mat, sp[i], s_new);
    for (int j = 0; j < 3; j++) sp[i][j] = s_new[j];
  }
}

/* ----------------------------------------------------------------------
  calculate 3x3 matrix exponential using Rodrigues' formula
  (R. Murray, Z. Li, and S. Shankar Sastry,
  A Mathematical Introduction to
  Robotic Manipulation (1994), p. 28 and 30).

  upp_tr - vector x, y, z so that one calculate
  U = exp(A) with A= [[0, x, y],
                      [-x, 0, z],
                      [-y, -z, 0]]
------------------------------------------------------------------------- */

void MinSpinCG::rodrigues_rotation(const double *upp_tr, double *out)
{
  double theta,A,B,D,x,y,z;
  double s1,s2,s3,a1,a2,a3;

  if (fabs(upp_tr[0]) < 1.0e-40 &&
      fabs(upp_tr[1]) < 1.0e-40 &&
      fabs(upp_tr[2]) < 1.0e-40){

    // if upp_tr is zero, return unity matrix

    for(int k = 0; k < 3; k++){
      for(int m = 0; m < 3; m++){
        if (m == k) out[3 * k + m] = 1.0;
        else out[3 * k + m] = 0.0;
      }
    }
    return;
  }

  theta = sqrt(upp_tr[0] * upp_tr[0] +
               upp_tr[1] * upp_tr[1] +
               upp_tr[2] * upp_tr[2]);

  A = cos(theta);
  B = sin(theta);
  D = 1.0 - A;
  x = upp_tr[0]/theta;
  y = upp_tr[1]/theta;
  z = upp_tr[2]/theta;

  // diagonal elements of U

  out[0] = A + z * z * D;
  out[4] = A + y * y * D;
  out[8] = A + x * x * D;

  // off diagonal of U

  s1 = -y * z *D;
  s2 = x * z * D;
  s3 = -x * y * D;

  a1 = x * B;
  a2 = y * B;
  a3 = z * B;

  out[1] = s1 + a1;
  out[3] = s1 - a1;
  out[2] = s2 + a2;
  out[6] = s2 - a2;
  out[5] = s3 + a3;
  out[7] = s3 - a3;

}

/* ----------------------------------------------------------------------
  out = vector^T x m,
  m -- 3x3 matrix , v -- 3-d vector
------------------------------------------------------------------------- */

void MinSpinCG::vm3(const double *m, const double *v, double *out)
{
  for(int i = 0; i < 3; i++){
    out[i] = 0.0;
    for(int j = 0; j < 3; j++) out[i] += *(m + 3 * j + i) * v[j];
  }
}

/* ----------------------------------------------------------------------
  advance spins
------------------------------------------------------------------------- */

void MinSpinCG::make_step(double c, double *energy_and_der)
{
  double p_scaled[3];
  int nlocal = atom->nlocal;
  double rot_mat[9]; // exponential of matrix made of search direction
  double s_new[3];
  double **sp = atom->sp;
  double der_e_cur_tmp = 0.0;

  for (int i = 0; i < nlocal; i++) {

    // scale the search direction

    for (int j = 0; j < 3; j++) p_scaled[j] = c * p_s[3 * i + j];

    // calculate rotation matrix

    rodrigues_rotation(p_scaled, rot_mat);

    // rotate spins

    vm3(rot_mat, sp[i], s_new);
    for (int j = 0; j < 3; j++) sp[i][j] = s_new[j];
  }

  ecurrent = energy_force(0);
  calc_gradient();
  neval++;
  der_e_cur = 0.0;
  for (int i = 0; i < 3 * nlocal; i++) {
    der_e_cur += g_cur[i] * p_s[i];
  }
  MPI_Allreduce(&der_e_cur,&der_e_cur_tmp,1,MPI_DOUBLE,MPI_SUM,world);
  der_e_cur = der_e_cur_tmp;
  if (update->multireplica == 1) {
    MPI_Allreduce(&der_e_cur_tmp,&der_e_cur,1,MPI_DOUBLE,MPI_SUM,universe->uworld);
  }

  energy_and_der[0] = ecurrent;
  energy_and_der[1] = der_e_cur;
}

/* ----------------------------------------------------------------------
  Calculate step length which satisfies approximate Wolfe conditions
  using the cubic interpolation
------------------------------------------------------------------------- */

int MinSpinCG::calc_and_make_step(double a, double b, int index)
{
  double e_and_d[2] = {0.0,0.0};
  double alpha,c1,c2,c3;
  double **sp = atom->sp;
  int nlocal = atom->nlocal;

  make_step(b,e_and_d);
  ecurrent = e_and_d[0];
  der_e_cur = e_and_d[1];
  index++;

  if (adescent(eprevious,e_and_d[0]) || index == 5){
    MPI_Bcast(&b,1,MPI_DOUBLE,0,world);
    for (int i = 0; i < 3 * nlocal; i++) {
      p_s[i] = b * p_s[i];
    }
    return 1;
  }
  else {
    double r,f0,f1,df0,df1;
    r = b - a;
    f0 = eprevious;
    f1 = ecurrent;
    df0 = der_e_pr;
    df1 = der_e_cur;

    c1 = -2.0*(f1-f0)/(r*r*r)+(df1+df0)/(r*r);
    c2 = 3.0*(f1-f0)/(r*r)-(df1+2.0*df0)/(r);
    c3 = df0;

    // f(x) = c1 x^3 + c2 x^2 + c3 x^1 + c4
    // has minimum at alpha below. We do not check boundaries.

    alpha = (-c2 + sqrt(c2*c2 - 3.0*c1*c3))/(3.0*c1);
    MPI_Bcast(&alpha,1,MPI_DOUBLE,0,world);

    if (alpha < 0.0) alpha = r/2.0;

    for (int i = 0; i < nlocal; i++) {
      for (int j = 0; j < 3; j++) sp[i][j] = sp_copy[i][j];
    }
    calc_and_make_step(0.0, alpha, index);
   }

  return 0;
}

/* ----------------------------------------------------------------------
  Approximate descent
------------------------------------------------------------------------- */

int MinSpinCG::adescent(double phi_0, double phi_j){

  double eps = 1.0e-6;

  if (phi_j<=phi_0+eps*fabs(phi_0))
    return 1;
  else
    return 0;
}

/* ----------------------------------------------------------------------
   evaluate max timestep
---------------------------------------------------------------------- */

double MinSpinCG::evaluate_dt()
{
  double dtmax;
  double fmsq;
  double fmaxsqone,fmaxsqloc,fmaxsqall;
  int nlocal = atom->nlocal;
  double **fm = atom->fm;

  // finding max fm on this proc.

  fmsq = fmaxsqone = fmaxsqloc = fmaxsqall = 0.0;
  for (int i = 0; i < nlocal; i++) {
    fmsq = fm[i][0]*fm[i][0]+fm[i][1]*fm[i][1]+fm[i][2]*fm[i][2];
    fmaxsqone = MAX(fmaxsqone,fmsq);
  }

  // finding max fm on this replica

  fmaxsqloc = fmaxsqone;
  MPI_Allreduce(&fmaxsqone,&fmaxsqloc,1,MPI_DOUBLE,MPI_MAX,world);

  // finding max fm over all replicas, if necessary
  // this communicator would be invalid for multiprocess replicas

  fmaxsqall = fmaxsqloc;
  if (update->multireplica == 1) {
    fmaxsqall = fmaxsqloc;
    MPI_Allreduce(&fmaxsqloc,&fmaxsqall,1,MPI_DOUBLE,MPI_MAX,universe->uworld);
  }

  if (fmaxsqall == 0.0)
    error->all(FLERR,"Incorrect fmaxsqall calculation");

  // define max timestep by dividing by the
  // inverse of max frequency by discrete_factor

  dtmax = MY_2PI/(discrete_factor*sqrt(fmaxsqall));

  return dtmax;
}

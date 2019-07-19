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
#include "min_spin_oso_lbfgs_ls.h"
#include "universe.h"
#include "atom.h"
#include "citeme.h"
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
#include <iostream>


using namespace LAMMPS_NS;
using namespace MathConst;

static const char cite_minstyle_spin_oso_lbfgs_ls[] =
  "min_style spin/oso_lbfgs_ls command:\n\n"
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

MinSpinOSO_LBFGS_LS::MinSpinOSO_LBFGS_LS(LAMMPS *lmp) :
  Min(lmp), g_old(NULL), g_cur(NULL), p_s(NULL), ds(NULL), dy(NULL), rho(NULL), sp_copy(NULL)
{
  if (lmp->citeme) lmp->citeme->add(cite_minstyle_spin_oso_lbfgs_ls);
  nlocal_max = 0;

  // nreplica = number of partitions
  // ireplica = which world I am in universe

  nreplica = universe->nworlds;
  ireplica = universe->iworld;
  use_line_search = 1;

}

/* ---------------------------------------------------------------------- */

MinSpinOSO_LBFGS_LS::~MinSpinOSO_LBFGS_LS()
{
    memory->destroy(g_old);
    memory->destroy(g_cur);
    memory->destroy(p_s);
    memory->destroy(ds);
    memory->destroy(dy);
    memory->destroy(rho);
    if (use_line_search)
      memory->destroy(sp_copy);
}

/* ---------------------------------------------------------------------- */

void MinSpinOSO_LBFGS_LS::init()
{
  num_mem = 3;
  local_iter = 0;
  der_e_cur = 0.0;
  der_e_pr = 0.0;

  Min::init();

  last_negative = update->ntimestep;

  // allocate tables

  nlocal_max = atom->nlocal;
  memory->grow(g_old,3*nlocal_max,"min/spin/oso/lbfgs_ls:g_old");
  memory->grow(g_cur,3*nlocal_max,"min/spin/oso/lbfgs_ls:g_cur");
  memory->grow(p_s,3*nlocal_max,"min/spin/oso/lbfgs_ls:p_s");
  memory->grow(rho,num_mem,"min/spin/oso/lbfgs_ls:rho");
  memory->grow(ds,num_mem,3*nlocal_max,"min/spin/oso/lbfgs_ls:ds");
  memory->grow(dy,num_mem,3*nlocal_max,"min/spin/oso/lbfgs_ls:dy");
  if (use_line_search)
    memory->grow(sp_copy,nlocal_max,3,"min/spin/oso/lbfgs_ls:sp_copy");

}

/* ---------------------------------------------------------------------- */

void MinSpinOSO_LBFGS_LS::setup_style()
{
  double **v = atom->v;
  int nlocal = atom->nlocal;

  // check if the atom/spin style is defined

  if (!atom->sp_flag)
    error->all(FLERR,"min/spin_oso_lbfgs_ls requires atom/spin style");

  for (int i = 0; i < nlocal; i++)
    v[i][0] = v[i][1] = v[i][2] = 0.0;
}

/* ---------------------------------------------------------------------- */

int MinSpinOSO_LBFGS_LS::modify_param(int narg, char **arg)
{

  if (strcmp(arg[0],"line_search") == 0) {
    if (narg < 2) error->all(FLERR,"Illegal fix_modify command");
    use_line_search = force->numeric(FLERR,arg[1]);
    return 2;
  }
  return 0;
}

/* ----------------------------------------------------------------------
   set current vector lengths and pointers
   called after atoms have migrated
------------------------------------------------------------------------- */

void MinSpinOSO_LBFGS_LS::reset_vectors()
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
   minimization via damped spin dynamics
------------------------------------------------------------------------- */

int MinSpinOSO_LBFGS_LS::iterate(int maxiter)
{
  int nlocal = atom->nlocal;
  bigint ntimestep;
  double fmdotfm;
  int flag, flagall;
  double **sp = atom->sp;
  double der_e_cur_tmp = 0.0;

  if (nlocal_max < nlocal) {
    nlocal_max = nlocal;
    memory->grow(g_old,3*nlocal_max,"min/spin/oso/lbfgs_ls:g_old");
    memory->grow(g_cur,3*nlocal_max,"min/spin/oso/lbfgs_ls:g_cur");
    memory->grow(p_s,3*nlocal_max,"min/spin/oso/lbfgs_ls:p_s");
    memory->grow(rho,num_mem,"min/spin/oso/lbfgs_ls:rho");
    memory->grow(ds,num_mem,3*nlocal_max,"min/spin/oso/lbfgs_ls:ds");
    memory->grow(dy,num_mem,3*nlocal_max,"min/spin/oso/lbfgs_ls:dy");
    if (use_line_search)
      memory->grow(sp_copy,nlocal_max,3,"min/spin/oso/lbfgs_ls:sp_copy");
  }

  for (int iter = 0; iter < maxiter; iter++) {
    
    if (timer->check_timeout(niter))
      return TIMEOUT;
  
    ntimestep = ++update->ntimestep;
    niter++;
  
    // optimize timestep accross processes / replicas
    // need a force calculation for timestep optimization

    if (local_iter == 0){
        ecurrent = energy_force(0);
        calc_gradient();
        neval++;
    }
    calc_search_direction();

    if (use_line_search) {

      // here we need to do line search

      der_e_cur = 0.0;
      for (int i = 0; i < 3 * nlocal; i++) {
        der_e_cur += g_cur[i] * p_s[i];
      }
      MPI_Allreduce(&der_e_cur,&der_e_cur_tmp,1,MPI_DOUBLE,MPI_SUM,world);
      der_e_cur = der_e_cur_tmp;
      if (update->multireplica == 1) {
        MPI_Allreduce(&der_e_cur_tmp,&der_e_cur,1,MPI_DOUBLE,MPI_SUM,universe->uworld);
      }
      for (int i = 0; i < nlocal; i++) {
        for (int j = 0; j < 3; j++) sp_copy[i][j] = sp[i][j];
      }
      eprevious = ecurrent;
      der_e_pr = der_e_cur;
      calc_and_make_step(0.0, 1.0, 0);
    }
    else{

      // here we don't do line search

      advance_spins();
      eprevious = ecurrent;
      ecurrent = energy_force(0);
      neval++;
    }

    //// energy tolerance criterion
    //// only check after DELAYSTEP elapsed since velocties reset to 0
    //// sync across replicas if running multi-replica minimization
  
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

    if (update->ftol > 0.0) {
      fmdotfm = fmnorm_sqr();
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

void MinSpinOSO_LBFGS_LS::calc_gradient()
{
  int nlocal = atom->nlocal;
  double **sp = atom->sp;
  double **fm = atom->fm;
  double tdampx, tdampy, tdampz;
  
  // loop on all spins on proc.

  for (int i = 0; i < nlocal; i++) {
    
    // calculate gradients
    
    g_cur[3 * i + 0] = (fm[i][0]*sp[i][1] - fm[i][1]*sp[i][0]);
    g_cur[3 * i + 1] = -(fm[i][2]*sp[i][0] - fm[i][0]*sp[i][2]);
    g_cur[3 * i + 2] = (fm[i][1]*sp[i][2] - fm[i][2]*sp[i][1]);
  }
}

/* ----------------------------------------------------------------------
   search direction:
   Limited-memory BFGS.
   See Jorge Nocedal and Stephen J. Wright 'Numerical
   Optimization' Second Edition, 2006 (p. 177)
---------------------------------------------------------------------- */

void MinSpinOSO_LBFGS_LS::calc_search_direction()
{
  int nlocal = atom->nlocal;

  double dyds = 0.0;
  double sq = 0.0;
  double yy = 0.0;
  double yr = 0.0;
  double beta = 0.0;

  double dyds_global = 0.0;
  double sq_global = 0.0;
  double yy_global = 0.0;
  double yr_global = 0.0;
  double beta_global = 0.0;

  int m_index = local_iter % num_mem; // memory index
  int c_ind = 0;
  double *q;
  double *alpha;

  double factor;

  if (nreplica > 1) {
    if (ireplica == 0 || ireplica == nreplica - 1) {
      factor = 0.0;
    }
    else factor = 1.0;
  }else{
    factor = 1.0;
  }

  q = (double *) calloc(3*nlocal, sizeof(double));
  alpha = (double *) calloc(num_mem, sizeof(double));

  if (local_iter == 0){ 	// steepest descent direction
    for (int i = 0; i < 3 * nlocal; i++) {
      p_s[i] = -g_cur[i];
      g_old[i] = g_cur[i];
      for (int k = 0; k < num_mem; k++){
        ds[k][i] = 0.0;
        dy[k][i] = 0.0;
        rho[k] = 0.0;
      }
    }
  } else {
    dyds = 0.0;
    for (int i = 0; i < 3 * nlocal; i++) {
      ds[m_index][i] = p_s[i];
      dy[m_index][i] = g_cur[i] - g_old[i];
      dyds += ds[m_index][i] * dy[m_index][i];
    }
    MPI_Allreduce(&dyds, &dyds_global, 1, MPI_DOUBLE, MPI_SUM, world);

    if (update->multireplica == 1) {
      dyds_global *= factor;
      dyds = dyds_global;
      MPI_Allreduce(&dyds, &dyds_global, 1,MPI_DOUBLE,MPI_SUM,universe->uworld);
    }

    if (fabs(dyds_global) > 1.0e-60) rho[m_index] = 1.0 / dyds_global;
    else rho[m_index] = 1.0e60;

    if (rho[m_index] < 0.0){
      local_iter = 0;
      for (int k = 0; k < num_mem; k++){
	      for (int i = 0; i < nlocal; i ++){
      ds[k][i] = 0.0;
      dy[k][i] = 0.0;
        }
      }
      return calc_search_direction();
    }

    // set the q vector

    for (int i = 0; i < 3 * nlocal; i++) {
      q[i] = g_cur[i];
    }

    // loop over last m indecies
    for(int k = num_mem - 1; k > -1; k--) {
      // this loop should run from the newest memory to the oldest one.

      c_ind = (k + m_index + 1) % num_mem;

      // dot product between dg and q

      sq = 0.0;
      for (int i = 0; i < 3 * nlocal; i++) {
        sq += ds[c_ind][i] * q[i];
      }
      MPI_Allreduce(&sq,&sq_global,1,MPI_DOUBLE,MPI_SUM,world);
      if (update->multireplica == 1) {
        sq_global *= factor;
        sq = sq_global;
        MPI_Allreduce(&sq,&sq_global,1,MPI_DOUBLE,MPI_SUM,universe->uworld);
      }

      // update alpha

      alpha[c_ind] = rho[c_ind] * sq_global;

      // update q

      for (int i = 0; i < 3 * nlocal; i++) {
        q[i] -= alpha[c_ind] * dy[c_ind][i];
      }
    }

    // dot product between dg with itself
    yy = 0.0;
    for (int i = 0; i < 3 * nlocal; i++) {
      yy += dy[m_index][i] * dy[m_index][i];
    }
    MPI_Allreduce(&yy,&yy_global,1,MPI_DOUBLE,MPI_SUM,world);
    if (update->multireplica == 1) {
      yy_global *= factor;
      yy = yy_global;
      MPI_Allreduce(&yy,&yy_global,1,MPI_DOUBLE,MPI_SUM,universe->uworld);
    }

    // calculate now search direction

    double devis = rho[m_index] * yy_global;

    if (fabs(devis) > 1.0e-60) {
      for (int i = 0; i < 3 * nlocal; i++) {
        p_s[i] = factor * q[i] / devis;
      }
    }else{
      for (int i = 0; i < 3 * nlocal; i++) {
        p_s[i] = factor * q[i] * 1.0e60;
      }
    }

    for (int k = 0; k < num_mem; k++){
      // this loop should run from the oldest memory to the newest one.

      if (local_iter < num_mem) c_ind = k;
      else c_ind = (k + m_index + 1) % num_mem;

      // dot product between p and da
      yr = 0.0;
      for (int i = 0; i < 3 * nlocal; i++) {
        yr += dy[c_ind][i] * p_s[i];
      }

      MPI_Allreduce(&yr,&yr_global,1,MPI_DOUBLE,MPI_SUM,world);
      if (update->multireplica == 1) {
        yr_global *= factor;
        yr = yr_global;
        MPI_Allreduce(&yr,&yr_global,1,MPI_DOUBLE,MPI_SUM,universe->uworld);
      }

      beta = rho[c_ind] * yr_global;
      for (int i = 0; i < 3 * nlocal; i++) {
        p_s[i] += ds[c_ind][i] * (alpha[c_ind] - beta);
      }
    }
    for (int i = 0; i < 3 * nlocal; i++) {
      p_s[i] = - factor * p_s[i];
      g_old[i] = g_cur[i];
    }
  }

  local_iter++;
  free(q);
  free(alpha);

}

/* ----------------------------------------------------------------------
   rotation of spins along the search direction
---------------------------------------------------------------------- */

void MinSpinOSO_LBFGS_LS::advance_spins()
{
  int nlocal = atom->nlocal;
  double **sp = atom->sp;
  double **fm = atom->fm;
  double tdampx, tdampy, tdampz;
  double rot_mat[9]; // exponential of matrix made of search direction
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
   compute and return ||mag. torque||_2^2
------------------------------------------------------------------------- */

double MinSpinOSO_LBFGS_LS::fmnorm_sqr()
{
  int nlocal = atom->nlocal;
  double tx,ty,tz;
  double **sp = atom->sp;
  double **fm = atom->fm;

  // calc. magnetic torques norm

  double local_norm2_sqr = 0.0;
  for (int i = 0; i < 3 * nlocal; i++) {
    local_norm2_sqr += g_cur[i]*g_cur[i];
  }

  // no extra atom calc. for spins

  if (nextra_atom)
    error->all(FLERR,"extra atom option not available yet");

  double norm2_sqr = 0.0;
  MPI_Allreduce(&local_norm2_sqr,&norm2_sqr,1,MPI_DOUBLE,MPI_SUM,world);

  return norm2_sqr;
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

void MinSpinOSO_LBFGS_LS::rodrigues_rotation(const double *upp_tr, double *out)
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
  D = 1 - A;
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

void MinSpinOSO_LBFGS_LS::vm3(const double *m, const double *v, double *out)
{
  for(int i = 0; i < 3; i++){
    out[i] *= 0.0;
    for(int j = 0; j < 3; j++)
    out[i] += *(m + 3 * j + i) * v[j];
  }
}


void MinSpinOSO_LBFGS_LS::make_step(double c, double *energy_and_der)
{
  double p_scaled[3];
  int nlocal = atom->nlocal;
  double rot_mat[9]; // exponential of matrix made of search direction
  double s_new[3];
  double **sp = atom->sp;
  double der_e_cur_tmp = 0.0;;

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
  MPI_Allreduce(&der_e_cur,&der_e_cur_tmp, 1, MPI_DOUBLE, MPI_SUM, world);
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

int MinSpinOSO_LBFGS_LS::calc_and_make_step(double a, double b, int index)
{
  double e_and_d[2] = {0.0,0.0};
  double alpha,c1,c2,c3;
  double **sp = atom->sp;
  int nlocal = atom->nlocal;

  make_step(b,e_and_d);
  ecurrent = e_and_d[0];
  der_e_cur = e_and_d[1];
  index++;

  if (awc(der_e_pr,eprevious,e_and_d[1],e_and_d[0]) || index == 5){
    MPI_Bcast(&b,1,MPI_DOUBLE,0,world);
    for (int i = 0; i < 3 * nlocal; i++) {
      p_s[i] = b * p_s[i];
    }
    return 1;
  }
  else{
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
  Approximate Wolfe conditions:
  William W. Hager and Hongchao Zhang
  SIAM J. optim., 16(1), 170-192. (23 pages)
------------------------------------------------------------------------- */

int MinSpinOSO_LBFGS_LS::awc(double der_phi_0, double phi_0, double der_phi_j, double phi_j){

  double eps = 1.0e-6;
  double delta = 0.1;
  double sigma = 0.9;

  if ((phi_j<=phi_0+eps*fabs(phi_0)) && ((2.0*delta-1.0) * der_phi_0>=der_phi_j>=sigma*der_phi_0))
    return 1;
  else
    return 0;
}

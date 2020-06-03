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

#include "fix_rx_kokkos.h"
#include <cstring>
#include "atom_masks.h"
#include "atom_kokkos.h"
#include "force.h"
#include "memory_kokkos.h"
#include "update.h"
#include "modify.h"
#include "neighbor.h"
#include "neigh_list_kokkos.h"
#include "neigh_request.h"
#include "error.h"
#include "math_special_kokkos.h"
#include "comm.h"
#include "domain.h"
#include "kokkos.h"
#include "utils.h"

#include <cfloat> // DBL_EPSILON

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathSpecialKokkos;

#ifdef DBL_EPSILON
  #define MY_EPSILON (10.0*DBL_EPSILON)
#else
  #define MY_EPSILON (10.0*2.220446049250313e-16)
#endif

#define SparseKinetics_enableIntegralReactions (true)
#define SparseKinetics_invalidIndex (-1)

// From fix_rx.cpp ... this should be lifted into fix_rx.h or fix_rx_kokkos.h?
enum{NONE,HARMONIC};
enum{LUCY};

namespace /* anonymous */
{

typedef KK_FLOAT TimerType;
TimerType getTimeStamp(void) { return MPI_Wtime(); }
KK_FLOAT getElapsedTime(const TimerType &t0, const TimerType &t1) { return t1-t0; }

} // end namespace

/* ---------------------------------------------------------------------- */

template <ExecutionSpace Space>
FixRxKokkos<Space>::FixRxKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixRX(lmp, narg, arg),
  pairDPDEKK(NULL),
  update_kinetics_data(true)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = Space;
  datamask_read = EMPTY_MASK;
  datamask_modify = EMPTY_MASK;

  k_error_flag = DAT::tdual_int_scalar("FixRxKokkos::k_error_flag");

  //printf("Inside FixRxKokkos::FixRxKokkos\n");
}

template <ExecutionSpace Space>
FixRxKokkos<Space>::~FixRxKokkos()
{
  //printf("Inside FixRxKokkos::~FixRxKokkos copymode= %d\n", copymode);
  if (copymode) return;

  if (localTempFlag)
    memoryKK->destroy_kokkos(k_dpdThetaLocal, dpdThetaLocal);

  memoryKK->destroy_kokkos(k_sumWeights, sumWeights);
  //memoryKK->destroy_kokkos(k_sumWeights);

  //delete [] scratchSpace;
  memoryKK->destroy_kokkos(d_scratchSpace);

  memoryKK->destroy_kokkos(k_cutsq);
}

/* ---------------------------------------------------------------------- */

template <ExecutionSpace Space>
void FixRxKokkos<Space>::post_constructor()
{
  // Run the parents and then reset one value.
  FixRX::post_constructor();

  // Need a copy of this
  this->my_restartFlag = modify->fix[modify->nfix-1]->restart_reset;
}

/* ---------------------------------------------------------------------- */

template <ExecutionSpace Space>
void FixRxKokkos<Space>::init()
{
  //printf("Inside FixRxKokkos::init\n");

  // Call the parent's version.
  //FixRX::init();

  pairDPDE = (PairDPDfdtEnergy *) force->pair_match("dpd/fdt/energy",1);
  if (pairDPDE == NULL)
    pairDPDE = (PairDPDfdtEnergy *) force->pair_match("dpd/fdt/energy/kk",1);

  if (pairDPDE == NULL)
    error->all(FLERR,"Must use pair_style dpd/fdt/energy with fix rx");

  pairDPDEKK = dynamic_cast<decltype(pairDPDEKK)>(pairDPDE);
  if (pairDPDEKK == NULL)
    error->all(FLERR,"Must use pair_style dpd/fdt/energy/kk with fix rx/kk");

  bool eos_flag = false;
  for (int i = 0; i < modify->nfix; i++)
    if (utils::strmatch(modify->fix[i]->style,"^eos/table/rx")) eos_flag = true;
  if(!eos_flag) error->all(FLERR,"fix rx requires fix eos/table/rx to be specified");

  if (update_kinetics_data)
    create_kinetics_data();

  // From FixRX::init()
  // need a half neighbor list
  // built whenever re-neighboring occurs

  int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->fix = 1;

  // Update the neighbor data for Kokkos.
  int neighflag = lmp->kokkos->neighflag;

  neighbor->requests[irequest]->
    kokkos_host = (Space == Host) &&
    !(Space == Device);
  neighbor->requests[irequest]->
    kokkos_device = (Space == Device);

  if (neighflag == FULL) {
    neighbor->requests[irequest]->full = 1;
    neighbor->requests[irequest]->half = 0;
  } else { //if (neighflag == HALF || neighflag == HALFTHREAD)
    neighbor->requests[irequest]->full = 0;
    neighbor->requests[irequest]->half = 1;
  }
}

/* ---------------------------------------------------------------------- */

template <ExecutionSpace Space>
void FixRxKokkos<Space>::init_list(int, class NeighList* ptr)
{
  //printf("Inside FixRxKokkos::init_list\n");
  this->list = ptr;
}

/* ---------------------------------------------------------------------- */

template <ExecutionSpace Space>
void FixRxKokkos<Space>::rk4(const KK_FLOAT t_stop, KK_FLOAT *y, KK_FLOAT *rwork, void* v_params) const
{
  KK_FLOAT *k1 = rwork;
  KK_FLOAT *k2 = k1 + nspecies;
  KK_FLOAT *k3 = k2 + nspecies;
  KK_FLOAT *k4 = k3 + nspecies;
  KK_FLOAT *yp = k4 + nspecies;

  const int numSteps = minSteps;

  const KK_FLOAT h = t_stop / KK_FLOAT(numSteps);

  // Run the requested steps with h.
  for (int step = 0; step < numSteps; step++)
  {
    // k1
    rhs(0.0,y,k1,v_params);

    // k2
    for (int ispecies = 0; ispecies < nspecies; ispecies++)
      yp[ispecies] = y[ispecies] + 0.5*h*k1[ispecies];

    rhs(0.0,yp,k2,v_params);

    // k3
    for (int ispecies = 0; ispecies < nspecies; ispecies++)
      yp[ispecies] = y[ispecies] + 0.5*h*k2[ispecies];

    rhs(0.0,yp,k3,v_params);

    // k4
    for (int ispecies = 0; ispecies < nspecies; ispecies++)
      yp[ispecies] = y[ispecies] + h*k3[ispecies];

    rhs(0.0,yp,k4,v_params);

    for (int ispecies = 0; ispecies < nspecies; ispecies++)
      y[ispecies] += h*(k1[ispecies]/6.0 + k2[ispecies]/3.0 + k3[ispecies]/3.0 + k4[ispecies]/6.0);

  } // end for (int step...

}

/* ---------------------------------------------------------------------- */

template <ExecutionSpace Space>
  template <typename VectorType, typename UserDataType>
void FixRxKokkos<Space>::k_rk4(const KK_FLOAT t_stop, VectorType& y, VectorType& rwork, UserDataType& userData) const
{
  VectorType k1(rwork);
  VectorType k2(&k1[nspecies]);
  VectorType k3(&k2[nspecies]);
  VectorType k4(&k3[nspecies]);
  VectorType yp(&k4[nspecies]);

  const int numSteps = minSteps;

  const KK_FLOAT h = t_stop / KK_FLOAT(numSteps);

  // Run the requested steps with h.
  for (int step = 0; step < numSteps; step++)
  {
    // k1
    k_rhs(0.0,y,k1, userData);

    // k2
    for (int ispecies = 0; ispecies < nspecies; ispecies++)
      yp[ispecies] = y[ispecies] + 0.5*h*k1[ispecies];

    k_rhs(0.0,yp,k2, userData);

    // k3
    for (int ispecies = 0; ispecies < nspecies; ispecies++)
      yp[ispecies] = y[ispecies] + 0.5*h*k2[ispecies];

    k_rhs(0.0,yp,k3, userData);

    // k4
    for (int ispecies = 0; ispecies < nspecies; ispecies++)
      yp[ispecies] = y[ispecies] + h*k3[ispecies];

    k_rhs(0.0,yp,k4, userData);

    for (int ispecies = 0; ispecies < nspecies; ispecies++)
      y[ispecies] += h*(k1[ispecies]/6.0 + k2[ispecies]/3.0 + k3[ispecies]/3.0 + k4[ispecies]/6.0);

  } // end for (int step...

}

/* ---------------------------------------------------------------------- */

//     f1 = dt*f(t,x)
//     f2 = dt*f(t+ c20*dt,x + c21*f1)
//     f3 = dt*f(t+ c30*dt,x + c31*f1 + c32*f2)
//     f4 = dt*f(t+ c40*dt,x + c41*f1 + c42*f2 + c43*f3)
//     f5 = dt*f(t+dt,x + c51*f1 + c52*f2 + c53*f3 + c54*f4)
//     f6 = dt*f(t+ c60*dt,x + c61*f1 + c62*f2 + c63*f3 + c64*f4 + c65*f5)
//
//     fifth-order runge-kutta integration
//        x5 = x + b1*f1 + b3*f3 + b4*f4 + b5*f5 + b6*f6
//     fourth-order runge-kutta integration
//        x  = x + a1*f1 + a3*f3 + a4*f4 + a5*f5

template <ExecutionSpace Space>
  template <typename VectorType, typename UserDataType>
void FixRxKokkos<Space>::k_rkf45_step (const int neq, const KK_FLOAT h, VectorType& y, VectorType& y_out, VectorType& rwk, UserDataType& userData) const
{
   const KK_FLOAT c21=0.25;
   const KK_FLOAT c31=0.09375;
   const KK_FLOAT c32=0.28125;
   const KK_FLOAT c41=0.87938097405553;
   const KK_FLOAT c42=-3.2771961766045;
   const KK_FLOAT c43=3.3208921256258;
   const KK_FLOAT c51=2.0324074074074;
   const KK_FLOAT c52=-8.0;
   const KK_FLOAT c53=7.1734892787524;
   const KK_FLOAT c54=-0.20589668615984;
   const KK_FLOAT c61=-0.2962962962963;
   const KK_FLOAT c62=2.0;
   const KK_FLOAT c63=-1.3816764132554;
   const KK_FLOAT c64=0.45297270955166;
   const KK_FLOAT c65=-0.275;
   const KK_FLOAT a1=0.11574074074074;
   const KK_FLOAT a3=0.54892787524366;
   const KK_FLOAT a4=0.5353313840156;
   const KK_FLOAT a5=-0.2;
   const KK_FLOAT b1=0.11851851851852;
   const KK_FLOAT b3=0.51898635477583;
   const KK_FLOAT b4=0.50613149034201;
   const KK_FLOAT b5=-0.18;
   const KK_FLOAT b6=0.036363636363636;

   // local dependent variables (5 total)
   VectorType& f1 = rwk;
   VectorType  f2(&rwk[  neq]);
   VectorType  f3(&rwk[2*neq]);
   VectorType  f4(&rwk[3*neq]);
   VectorType  f5(&rwk[4*neq]);
   VectorType  f6(&rwk[5*neq]);

   // scratch for the intermediate solution.
   VectorType& ytmp = y_out;

   // 1)
   k_rhs (0.0, y, f1, userData);

   for (int k = 0; k < neq; k++){
      f1[k] *= h;
      ytmp[k] = y[k] + c21 * f1[k];
   }

   // 2)
   k_rhs(0.0, ytmp, f2, userData);

   for (int k = 0; k < neq; k++){
      f2[k] *= h;
      ytmp[k] = y[k] + c31 * f1[k] + c32 * f2[k];
   }

   // 3)
   k_rhs(0.0, ytmp, f3, userData);

   for (int k = 0; k < neq; k++) {
      f3[k] *= h;
      ytmp[k] = y[k] + c41 * f1[k] + c42 * f2[k] + c43 * f3[k];
   }

   // 4)
   k_rhs(0.0, ytmp, f4, userData);

   for (int k = 0; k < neq; k++) {
      f4[k] *= h;
      ytmp[k] = y[k] + c51 * f1[k] + c52 * f2[k] + c53 * f3[k] + c54 * f4[k];
   }

   // 5)
   k_rhs(0.0, ytmp, f5, userData);

   for (int k = 0; k < neq; k++) {
      f5[k] *= h;
      ytmp[k] = y[k] + c61*f1[k] + c62*f2[k] + c63*f3[k] + c64*f4[k] + c65*f5[k];
   }

   // 6)
   k_rhs(0.0, ytmp, f6, userData);

   for (int k = 0; k < neq; k++)
   {
      //const KK_FLOAT f6 = h * ydot[k];
      f6[k] *= h;

      // 5th-order solution.
      const KK_FLOAT r5 = b1*f1[k] + b3*f3[k] + b4*f4[k] + b5*f5[k] + b6*f6[k];

      // 4th-order solution.
      const KK_FLOAT r4 = a1*f1[k] + a3*f3[k] + a4*f4[k] + a5*f5[k];

      // Truncation error: difference between 4th and 5th-order solutions.
      rwk[k] = fabs(r5 - r4);

      // Update solution.
    //y_out[k] = y[k] + r5; // Local extrapolation
      y_out[k] = y[k] + r4;
   }

   return;
}

template <ExecutionSpace Space>
  template <typename VectorType, typename UserDataType>
int FixRxKokkos<Space>::k_rkf45_h0
                    (const int neq, const KK_FLOAT t, const KK_FLOAT t_stop,
                     const KK_FLOAT hmin, const KK_FLOAT hmax,
                     KK_FLOAT& h0, VectorType& y, VectorType& rwk, UserDataType& userData) const
{
   // Set lower and upper bounds on h0, and take geometric mean as first trial value.
   // Exit with this value if the bounds cross each other.

   // Adjust upper bound based on ydot ...
   KK_FLOAT hg = sqrt(hmin*hmax);

   //if (hmax < hmin)
   //{
   //   h0 = hg;
   //   return;
   //}

   // Start iteration to find solution to ... {WRMS norm of (h0^2 y'' / 2)} = 1

   VectorType& ydot  = rwk;
   VectorType  y1    (&ydot[  neq]);
   VectorType  ydot1 (&ydot[2*neq]);

   const int max_iters = 10;
   bool hnew_is_ok = false;
   KK_FLOAT hnew = hg;
   int iter = 0;

   // compute ydot at t=t0
   k_rhs (t, y, ydot, userData);

   while(1)
   {
      // Estimate y'' with finite-difference ...

      for (int k = 0; k < neq; k++)
         y1[k] = y[k] + hg * ydot[k];

      // compute y' at t1
      k_rhs (t + hg, y1, ydot1, userData);

      // Compute WRMS norm of y''
      KK_FLOAT yddnrm = 0.0;
      for (int k = 0; k < neq; k++){
         KK_FLOAT ydd = (ydot1[k] - ydot[k]) / hg;
         KK_FLOAT wterr = ydd / (relTol * fabs(y[k]) + absTol);
         yddnrm += wterr * wterr;
      }

      yddnrm = sqrt(yddnrm / KK_FLOAT(neq));

      //std::cout << "iter " << _iter << " hg " << hg << " y'' " << yddnrm << std::endl;
      //std::cout << "ydot " << ydot[neq-1] << std::endl;

      // should we accept this?
      if (hnew_is_ok || iter == max_iters){
         hnew = hg;
         //if (iter == max_iters)
         //   fprintf(stderr, "ERROR_HIN_MAX_ITERS\n");
         break;
      }

      // Get the new value of h ...
      hnew = (yddnrm*hmax*hmax > 2.0) ? sqrt(2.0 / yddnrm) : sqrt(hg * hmax);

      // test the stopping conditions.
      KK_FLOAT hrat = hnew / hg;

      // Accept this value ... the bias factor should bring it within range.
      if ((hrat > 0.5) && (hrat < 2.0))
         hnew_is_ok = true;

      // If y'' is still bad after a few iterations, just accept h and give up.
      if ((iter > 1) && hrat > 2.0) {
         hnew = hg;
         hnew_is_ok = true;
      }

      //printf("iter=%d, yddnrw=%e, hnew=%e, hmin=%e, hmax=%e\n", iter, yddnrm, hnew, hmin, hmax);

      hg = hnew;
      iter ++;
   }

   // bound and bias estimate
   h0 = hnew * 0.5;
   h0 = fmax(h0, hmin);
   h0 = fmin(h0, hmax);
   //printf("h0=%e, hmin=%e, hmax=%e\n", h0, hmin, hmax);

   return (iter + 1);
}

template <ExecutionSpace Space>
  template <typename VectorType, typename UserDataType>
void FixRxKokkos<Space>::k_rkf45(const int neq, const KK_FLOAT t_stop, VectorType& y, VectorType& rwork, UserDataType& userData, CounterType& counter) const
{
  // Rounding coefficient.
  const KK_FLOAT uround = DBL_EPSILON;

  // Adaption limit (shrink or grow)
  const KK_FLOAT adaption_limit = 4.0;

  // Safety factor on the adaption. very specific but not necessary .. 0.9 is common.
  const KK_FLOAT hsafe = 0.840896415;

  // Time rounding factor.
  const KK_FLOAT tround = t_stop * uround;

  // Counters for diagnostics.
  int nst = 0; // # of steps (accepted)
  int nit = 0; // # of iterations total
  int nfe = 0; // # of RHS evaluations

  // Min/Max step-size limits.
  const KK_FLOAT h_min = 100.0 * tround;
  const KK_FLOAT h_max = (minSteps > 0) ? t_stop / KK_FLOAT(minSteps) : t_stop;

  // Set the initial step-size. 0 forces an internal estimate ... stable Euler step size.
  KK_FLOAT h = (minSteps > 0) ? t_stop / KK_FLOAT(minSteps) : 0.0;

  KK_FLOAT t = 0.0;

  if (h < h_min){
    //fprintf(stderr,"hin not implemented yet\n");
    //exit(-1);
    nfe = k_rkf45_h0 (neq, t, t_stop, h_min, h_max, h, y, rwork, userData);
  }

  //printf("t= %e t_stop= %e h= %e\n", t, t_stop, h);

  // Integrate until we reach the end time.
  while (fabs(t - t_stop) > tround)
  {
    VectorType& yout = rwork;
    VectorType  eout (&yout[neq]);

    // Take a trial step.
    k_rkf45_step (neq, h, y, yout, eout, userData);

    // Estimate the solution error.
      // ... weighted 2-norm of the error.
      KK_FLOAT err2 = 0.0;
      for (int k = 0; k < neq; k++){
        const KK_FLOAT wterr = eout[k] / (relTol * fabs(y[k]) + absTol);
        err2 += wterr * wterr;
      }

    KK_FLOAT err = fmax(uround, sqrt(err2 / KK_FLOAT(nspecies)));

    // Accept the solution?
    if (err <= 1.0 || h <= h_min){
      t += h;
      nst++;

      for (int k = 0; k < neq; k++)
        y[k] = yout[k];
    }

    // Adjust h for the next step.
    KK_FLOAT hfac = hsafe * sqrt(sqrt(1.0 / err));

    // Limit the adaption.
    hfac = fmax(hfac, 1.0 / adaption_limit);
    hfac = fmin(hfac,       adaption_limit);

    // Apply the adaption factor...
    h *= hfac;

    // Limit h.
    h = fmin(h, h_max);
    h = fmax(h, h_min);

    // Stretch h if we're within 5% ... and we didn't just fail.
    if (err <= 1.0 && (t + 1.05*h) > t_stop)
      h = t_stop - t;

    // And don't overshoot the end.
    if (t + h > t_stop)
      h = t_stop - t;

    nit++;
    nfe += 6;

    if (maxIters && nit > maxIters){
      //fprintf(stderr,"atom[%d] took too many iterations in rkf45 %d %e %e\n", id, nit, t, t_stop);
      counter.nFails ++;
      break;
      // We should set an error here so that the solution is not used!
    }

  } // end while

  counter.nSteps += nst;
  counter.nIters += nit;
  counter.nFuncs += nfe;

  //printf("id= %d nst= %d nit= %d\n", id, nst, nit);
}
/* ---------------------------------------------------------------------- */

//     f1 = dt*f(t,x)
//     f2 = dt*f(t+ c20*dt,x + c21*f1)
//     f3 = dt*f(t+ c30*dt,x + c31*f1 + c32*f2)
//     f4 = dt*f(t+ c40*dt,x + c41*f1 + c42*f2 + c43*f3)
//     f5 = dt*f(t+dt,x + c51*f1 + c52*f2 + c53*f3 + c54*f4)
//     f6 = dt*f(t+ c60*dt,x + c61*f1 + c62*f2 + c63*f3 + c64*f4 + c65*f5)
//
//     fifth-order runge-kutta integration
//        x5 = x + b1*f1 + b3*f3 + b4*f4 + b5*f5 + b6*f6
//     fourth-order runge-kutta integration
//        x  = x + a1*f1 + a3*f3 + a4*f4 + a5*f5

template <ExecutionSpace Space>
void FixRxKokkos<Space>::rkf45_step (const int neq, const KK_FLOAT h, KK_FLOAT y[], KK_FLOAT y_out[], KK_FLOAT rwk[], void* v_param) const
{
   const KK_FLOAT c21=0.25;
   const KK_FLOAT c31=0.09375;
   const KK_FLOAT c32=0.28125;
   const KK_FLOAT c41=0.87938097405553;
   const KK_FLOAT c42=-3.2771961766045;
   const KK_FLOAT c43=3.3208921256258;
   const KK_FLOAT c51=2.0324074074074;
   const KK_FLOAT c52=-8.0;
   const KK_FLOAT c53=7.1734892787524;
   const KK_FLOAT c54=-0.20589668615984;
   const KK_FLOAT c61=-0.2962962962963;
   const KK_FLOAT c62=2.0;
   const KK_FLOAT c63=-1.3816764132554;
   const KK_FLOAT c64=0.45297270955166;
   const KK_FLOAT c65=-0.275;
   const KK_FLOAT a1=0.11574074074074;
   const KK_FLOAT a3=0.54892787524366;
   const KK_FLOAT a4=0.5353313840156;
   const KK_FLOAT a5=-0.2;
   const KK_FLOAT b1=0.11851851851852;
   const KK_FLOAT b3=0.51898635477583;
   const KK_FLOAT b4=0.50613149034201;
   const KK_FLOAT b5=-0.18;
   const KK_FLOAT b6=0.036363636363636;

   // local dependent variables (5 total)
   KK_FLOAT* f1 = &rwk[    0];
   KK_FLOAT* f2 = &rwk[  neq];
   KK_FLOAT* f3 = &rwk[2*neq];
   KK_FLOAT* f4 = &rwk[3*neq];
   KK_FLOAT* f5 = &rwk[4*neq];
   KK_FLOAT* f6 = &rwk[5*neq];

   // scratch for the intermediate solution.
   //KK_FLOAT* ytmp = &rwk[6*neq];
   KK_FLOAT* ytmp = y_out;

   // 1)
   rhs (0.0, y, f1, v_param);

   for (int k = 0; k < neq; k++){
      f1[k] *= h;
      ytmp[k] = y[k] + c21 * f1[k];
   }

   // 2)
   rhs(0.0, ytmp, f2, v_param);

   for (int k = 0; k < neq; k++){
      f2[k] *= h;
      ytmp[k] = y[k] + c31 * f1[k] + c32 * f2[k];
   }

   // 3)
   rhs(0.0, ytmp, f3, v_param);

   for (int k = 0; k < neq; k++) {
      f3[k] *= h;
      ytmp[k] = y[k] + c41 * f1[k] + c42 * f2[k] + c43 * f3[k];
   }

   // 4)
   rhs(0.0, ytmp, f4, v_param);

   for (int k = 0; k < neq; k++) {
      f4[k] *= h;
      ytmp[k] = y[k] + c51 * f1[k] + c52 * f2[k] + c53 * f3[k] + c54 * f4[k];
   }

   // 5)
   rhs(0.0, ytmp, f5, v_param);

   for (int k = 0; k < neq; k++) {
      f5[k] *= h;
      ytmp[k] = y[k] + c61*f1[k] + c62*f2[k] + c63*f3[k] + c64*f4[k] + c65*f5[k];
   }

   // 6)
   rhs(0.0, ytmp, f6, v_param);

   for (int k = 0; k < neq; k++)
   {
      //const KK_FLOAT f6 = h * ydot[k];
      f6[k] *= h;

      // 5th-order solution.
      const KK_FLOAT r5 = b1*f1[k] + b3*f3[k] + b4*f4[k] + b5*f5[k] + b6*f6[k];

      // 4th-order solution.
      const KK_FLOAT r4 = a1*f1[k] + a3*f3[k] + a4*f4[k] + a5*f5[k];

      // Truncation error: difference between 4th and 5th-order solutions.
      rwk[k] = fabs(r5 - r4);

      // Update solution.
    //y_out[k] = y[k] + r5; // Local extrapolation
      y_out[k] = y[k] + r4;
   }

   return;
}

template <ExecutionSpace Space>
int FixRxKokkos<Space>::rkf45_h0
                    (const int neq, const KK_FLOAT t, const KK_FLOAT t_stop,
                     const KK_FLOAT hmin, const KK_FLOAT hmax,
                     KK_FLOAT& h0, KK_FLOAT y[], KK_FLOAT rwk[], void* v_params) const
{
   // Set lower and upper bounds on h0, and take geometric mean as first trial value.
   // Exit with this value if the bounds cross each other.

   // Adjust upper bound based on ydot ...
   KK_FLOAT hg = sqrt(hmin*hmax);

   //if (hmax < hmin)
   //{
   //   h0 = hg;
   //   return;
   //}

   // Start iteration to find solution to ... {WRMS norm of (h0^2 y'' / 2)} = 1

   KK_FLOAT *ydot  = rwk;
   KK_FLOAT *y1    = ydot + neq;
   KK_FLOAT *ydot1 = y1 + neq;

   const int max_iters = 10;
   bool hnew_is_ok = false;
   KK_FLOAT hnew = hg;
   int iter = 0;

   // compute ydot at t=t0
   rhs (t, y, ydot, v_params);

   while(1)
   {
      // Estimate y'' with finite-difference ...

      for (int k = 0; k < neq; k++)
         y1[k] = y[k] + hg * ydot[k];

      // compute y' at t1
      rhs (t + hg, y1, ydot1, v_params);

      // Compute WRMS norm of y''
      KK_FLOAT yddnrm = 0.0;
      for (int k = 0; k < neq; k++){
         KK_FLOAT ydd = (ydot1[k] - ydot[k]) / hg;
         KK_FLOAT wterr = ydd / (relTol * fabs(y[k]) + absTol);
         yddnrm += wterr * wterr;
      }

      yddnrm = sqrt(yddnrm / KK_FLOAT(neq));

      //std::cout << "iter " << _iter << " hg " << hg << " y'' " << yddnrm << std::endl;
      //std::cout << "ydot " << ydot[neq-1] << std::endl;

      // should we accept this?
      if (hnew_is_ok || iter == max_iters){
         hnew = hg;
         if (iter == max_iters)
            fprintf(stderr, "ERROR_HIN_MAX_ITERS\n");
         break;
      }

      // Get the new value of h ...
      hnew = (yddnrm*hmax*hmax > 2.0) ? sqrt(2.0 / yddnrm) : sqrt(hg * hmax);

      // test the stopping conditions.
      KK_FLOAT hrat = hnew / hg;

      // Accept this value ... the bias factor should bring it within range.
      if ((hrat > 0.5) && (hrat < 2.0))
         hnew_is_ok = true;

      // If y'' is still bad after a few iterations, just accept h and give up.
      if ((iter > 1) && hrat > 2.0) {
         hnew = hg;
         hnew_is_ok = true;
      }

      //printf("iter=%d, yddnrw=%e, hnew=%e, hmin=%e, hmax=%e\n", iter, yddnrm, hnew, hmin, hmax);

      hg = hnew;
      iter ++;
   }

   // bound and bias estimate
   h0 = hnew * 0.5;
   h0 = fmax(h0, hmin);
   h0 = fmin(h0, hmax);
   //printf("h0=%e, hmin=%e, hmax=%e\n", h0, hmin, hmax);

   return (iter + 1);
}

template <ExecutionSpace Space>
void FixRxKokkos<Space>::rkf45(const int neq, const KK_FLOAT t_stop, KK_FLOAT *y, KK_FLOAT *rwork, void *v_param, CounterType& counter) const
{
  // Rounding coefficient.
  const KK_FLOAT uround = DBL_EPSILON;

  // Adaption limit (shrink or grow)
  const KK_FLOAT adaption_limit = 4.0;

  // Safety factor on the adaption. very specific but not necessary .. 0.9 is common.
  const KK_FLOAT hsafe = 0.840896415;

  // Time rounding factor.
  const KK_FLOAT tround = t_stop * uround;

  // Counters for diagnostics.
  int nst = 0; // # of steps (accepted)
  int nit = 0; // # of iterations total
  int nfe = 0; // # of RHS evaluations

  // Min/Max step-size limits.
  const KK_FLOAT h_min = 100.0 * tround;
  const KK_FLOAT h_max = (minSteps > 0) ? t_stop / KK_FLOAT(minSteps) : t_stop;

  // Set the initial step-size. 0 forces an internal estimate ... stable Euler step size.
  KK_FLOAT h = (minSteps > 0) ? t_stop / KK_FLOAT(minSteps) : 0.0;

  KK_FLOAT t = 0.0;

  if (h < h_min){
    //fprintf(stderr,"hin not implemented yet\n");
    //exit(-1);
    nfe = rkf45_h0 (neq, t, t_stop, h_min, h_max, h, y, rwork, v_param);
  }

  //printf("t= %e t_stop= %e h= %e\n", t, t_stop, h);

  // Integrate until we reach the end time.
  while (fabs(t - t_stop) > tround){
    KK_FLOAT *yout = rwork;
    KK_FLOAT *eout = yout + neq;

    // Take a trial step.
    rkf45_step (neq, h, y, yout, eout, v_param);

    // Estimate the solution error.
      // ... weighted 2-norm of the error.
      KK_FLOAT err2 = 0.0;
      for (int k = 0; k < neq; k++){
        const KK_FLOAT wterr = eout[k] / (relTol * fabs(y[k]) + absTol);
        err2 += wterr * wterr;
      }

    KK_FLOAT err = fmax(uround, sqrt(err2 / KK_FLOAT(nspecies)));

    // Accept the solution?
    if (err <= 1.0 || h <= h_min){
      t += h;
      nst++;

      for (int k = 0; k < neq; k++)
        y[k] = yout[k];
    }

    // Adjust h for the next step.
    KK_FLOAT hfac = hsafe * sqrt(sqrt(1.0 / err));

    // Limit the adaption.
    hfac = fmax(hfac, 1.0 / adaption_limit);
    hfac = fmin(hfac,       adaption_limit);

    // Apply the adaption factor...
    h *= hfac;

    // Limit h.
    h = fmin(h, h_max);
    h = fmax(h, h_min);

    // Stretch h if we're within 5% ... and we didn't just fail.
    if (err <= 1.0 && (t + 1.05*h) > t_stop)
      h = t_stop - t;

    // And don't overshoot the end.
    if (t + h > t_stop)
      h = t_stop - t;

    nit++;
    nfe += 6;

    if (maxIters && nit > maxIters){
      //fprintf(stderr,"atom[%d] took too many iterations in rkf45 %d %e %e\n", id, nit, t, t_stop);
      counter.nFails ++;
      break;
      // We should set an error here so that the solution is not used!
    }

  } // end while

  counter.nSteps += nst;
  counter.nIters += nit;
  counter.nFuncs += nfe;

  //printf("id= %d nst= %d nit= %d\n", id, nst, nit);
}

/* ---------------------------------------------------------------------- */

template <ExecutionSpace Space>
int FixRxKokkos<Space>::rhs(KK_FLOAT t, const KK_FLOAT *y, KK_FLOAT *dydt, void *params) const
{
  // Use the sparse format instead.
  if (useSparseKinetics)
    return this->rhs_sparse(t, y, dydt, params);
  else
    return this->rhs_dense (t, y, dydt, params);
}

/* ---------------------------------------------------------------------- */

template <ExecutionSpace Space>
int FixRxKokkos<Space>::rhs_dense(KK_FLOAT t, const KK_FLOAT *y, KK_FLOAT *dydt, void *params) const
{
  UserRHSData *userData = (UserRHSData *) params;

  double *rxnRateLaw = userData->rxnRateLaw;
  double *kFor       = userData->kFor;

  //const KK_FLOAT VDPD = domain->xprd * domain->yprd * domain->zprd / atom->natoms;
  //const int nspecies = atom->nspecies_dpd;

  for(int ispecies=0; ispecies<nspecies; ispecies++)
    dydt[ispecies] = 0.0;

  // Construct the reaction rate laws
  for(int jrxn=0; jrxn<nreactions; jrxn++){
    KK_FLOAT rxnRateLawForward = kFor[jrxn];

    for(int ispecies=0; ispecies<nspecies; ispecies++){
      const KK_FLOAT concentration = y[ispecies]/VDPD;
      rxnRateLawForward *= pow(concentration, d_kineticsData.stoichReactants(jrxn,ispecies));
      //rxnRateLawForward *= pow(concentration,stoichReactants[jrxn][ispecies]);
    }
    rxnRateLaw[jrxn] = rxnRateLawForward;
  }

  // Construct the reaction rates for each species
  for(int ispecies=0; ispecies<nspecies; ispecies++)
    for(int jrxn=0; jrxn<nreactions; jrxn++)
    {
      dydt[ispecies] += d_kineticsData.stoich(jrxn,ispecies) *VDPD*rxnRateLaw[jrxn];
      //dydt[ispecies] += stoich[jrxn][ispecies]*VDPD*rxnRateLaw[jrxn];
    }

  return 0;
}

/* ---------------------------------------------------------------------- */

template <ExecutionSpace Space>
int FixRxKokkos<Space>::rhs_sparse(KK_FLOAT t, const KK_FLOAT *y, KK_FLOAT *dydt, void *v_params) const
{
   UserRHSData *userData = (UserRHSData *) v_params;

   //const KK_FLOAT VDPD = domain->xprd * domain->yprd * domain->zprd / atom->natoms;

   #define kFor         (userData->kFor)
   #define kRev         (NULL)
   #define rxnRateLaw   (userData->rxnRateLaw)
   #define conc         (dydt)
   #define maxReactants (this->sparseKinetics_maxReactants)
   #define maxSpecies   (this->sparseKinetics_maxSpecies)
   #define nuk          (this->d_kineticsData.nuk)
   #define nu           (this->d_kineticsData.nu)
   #define inu          (this->d_kineticsData.inu)
   #define isIntegral(idx) (SparseKinetics_enableIntegralReactions \
                             && this->d_kineticsData.isIntegral(idx))

   for (int k = 0; k < nspecies; ++k)
      conc[k] = y[k] / VDPD;

   // Construct the reaction rate laws
   for (int i = 0; i < nreactions; ++i)
   {
      KK_FLOAT rxnRateLawForward;
      if (isIntegral(i)){
         rxnRateLawForward = kFor[i] * powint(conc[ nuk(i,0) ], inu(i,0));
         for (int kk = 1; kk < maxReactants; ++kk){
            const int k = nuk(i,kk);
            if (k == SparseKinetics_invalidIndex) break;
            //if (k != SparseKinetics_invalidIndex)
               rxnRateLawForward *= powint(conc[k], inu(i,kk));
         }
      } else {
         rxnRateLawForward = kFor[i] * pow(conc[ nuk(i,0) ], nu(i,0));
         for (int kk = 1; kk < maxReactants; ++kk){
            const int k = nuk(i,kk);
            if (k == SparseKinetics_invalidIndex) break;
            //if (k != SparseKinetics_invalidIndex)
               rxnRateLawForward *= pow(conc[k], nu(i,kk));
         }
      }

      rxnRateLaw[i] = rxnRateLawForward;
   }

   // Construct the reaction rates for each species from the
   // Stoichiometric matrix and ROP vector.
   for (int k = 0; k < nspecies; ++k)
      dydt[k] = 0.0;

   for (int i = 0; i < nreactions; ++i){
      // Reactants ...
      dydt[ nuk(i,0) ] -= nu(i,0) * rxnRateLaw[i];
      for (int kk = 1; kk < maxReactants; ++kk){
         const int k = nuk(i,kk);
         if (k == SparseKinetics_invalidIndex) break;
         //if (k != SparseKinetics_invalidIndex)
            dydt[k] -= nu(i,kk) * rxnRateLaw[i];
      }

      // Products ...
      dydt[ nuk(i,maxReactants) ] += nu(i,maxReactants) * rxnRateLaw[i];
      for (int kk = maxReactants+1; kk < maxSpecies; ++kk){
         const int k = nuk(i,kk);
         if (k == SparseKinetics_invalidIndex) break;
         //if (k != SparseKinetics_invalidIndex)
            dydt[k] += nu(i,kk) * rxnRateLaw[i];
      }
   }

   // Add in the volume factor to convert to the proper units.
   for (int k = 0; k < nspecies; ++k)
      dydt[k] *= VDPD;

   #undef kFor
   #undef kRev
   #undef rxnRateLaw
   #undef conc
   #undef maxReactants
   #undef maxSpecies
   #undef nuk
   #undef nu
   #undef inu
   #undef isIntegral
   //#undef invalidIndex

   return 0;
}

/* ---------------------------------------------------------------------- */

template <ExecutionSpace Space>
  template <typename VectorType, typename UserDataType>
int FixRxKokkos<Space>::k_rhs(KK_FLOAT t, const VectorType& y, VectorType& dydt, UserDataType& userData) const
{
  // Use the sparse format instead.
  if (useSparseKinetics)
    return this->k_rhs_sparse(t, y, dydt, userData);
  else
    return this->k_rhs_dense (t, y, dydt, userData);
}

/* ---------------------------------------------------------------------- */

template <ExecutionSpace Space>
  template <typename VectorType, typename UserDataType>
int FixRxKokkos<Space>::k_rhs_dense(KK_FLOAT t, const VectorType& y, VectorType& dydt, UserDataType& userData) const
{
  #define rxnRateLaw (userData.rxnRateLaw)
  #define kFor       (userData.kFor     )

  //const KK_FLOAT VDPD = domain->xprd * domain->yprd * domain->zprd / atom->natoms;
  //const int nspecies = atom->nspecies_dpd;

  for(int ispecies=0; ispecies<nspecies; ispecies++)
    dydt[ispecies] = 0.0;

  // Construct the reaction rate laws
  for(int jrxn=0; jrxn<nreactions; jrxn++){
    KK_FLOAT rxnRateLawForward = kFor[jrxn];

    for(int ispecies=0; ispecies<nspecies; ispecies++){
      const KK_FLOAT concentration = y[ispecies]/VDPD;
      rxnRateLawForward *= pow(concentration, d_kineticsData.stoichReactants(jrxn,ispecies));
      //rxnRateLawForward *= pow(concentration,stoichReactants[jrxn][ispecies]);
    }
    rxnRateLaw[jrxn] = rxnRateLawForward;
  }

  // Construct the reaction rates for each species
  for(int ispecies=0; ispecies<nspecies; ispecies++)
    for(int jrxn=0; jrxn<nreactions; jrxn++)
    {
      dydt[ispecies] += d_kineticsData.stoich(jrxn,ispecies) *VDPD*rxnRateLaw[jrxn];
      //dydt[ispecies] += stoich[jrxn][ispecies]*VDPD*rxnRateLaw[jrxn];
    }

  #undef rxnRateLaw
  #undef kFor

  return 0;
}

/* ---------------------------------------------------------------------- */

template <ExecutionSpace Space>
  template <typename VectorType, typename UserDataType>
int FixRxKokkos<Space>::k_rhs_sparse(KK_FLOAT t, const VectorType& y, VectorType& dydt, UserDataType& userData) const
{
   #define kFor         (userData.kFor)
   #define kRev         (NULL)
   #define rxnRateLaw   (userData.rxnRateLaw)
   #define conc         (dydt)
   #define maxReactants (this->sparseKinetics_maxReactants)
   #define maxSpecies   (this->sparseKinetics_maxSpecies)
   #define nuk          (this->d_kineticsData.nuk)
   #define nu           (this->d_kineticsData.nu)
   #define inu          (this->d_kineticsData.inu)
   #define isIntegral(idx) (SparseKinetics_enableIntegralReactions \
                             && this->d_kineticsData.isIntegral(idx))

   for (int k = 0; k < nspecies; ++k)
      conc[k] = y[k] / VDPD;

   // Construct the reaction rate laws
   for (int i = 0; i < nreactions; ++i)
   {
      KK_FLOAT rxnRateLawForward;
      if (isIntegral(i)){
         rxnRateLawForward = kFor[i] * powint(conc[ nuk(i,0) ], inu(i,0));
         for (int kk = 1; kk < maxReactants; ++kk){
            const int k = nuk(i,kk);
            if (k == SparseKinetics_invalidIndex) break;
            //if (k != SparseKinetics_invalidIndex)
               rxnRateLawForward *= powint(conc[k], inu(i,kk));
         }
      } else {
         rxnRateLawForward = kFor[i] * pow(conc[ nuk(i,0) ], nu(i,0));
         for (int kk = 1; kk < maxReactants; ++kk){
            const int k = nuk(i,kk);
            if (k == SparseKinetics_invalidIndex) break;
            //if (k != SparseKinetics_invalidIndex)
               rxnRateLawForward *= pow(conc[k], nu(i,kk));
         }
      }

      rxnRateLaw[i] = rxnRateLawForward;
   }

   // Construct the reaction rates for each species from the
   // Stoichiometric matrix and ROP vector.
   for (int k = 0; k < nspecies; ++k)
      dydt[k] = 0.0;

   for (int i = 0; i < nreactions; ++i){
      // Reactants ...
      dydt[ nuk(i,0) ] -= nu(i,0) * rxnRateLaw[i];
      for (int kk = 1; kk < maxReactants; ++kk){
         const int k = nuk(i,kk);
         if (k == SparseKinetics_invalidIndex) break;
         //if (k != SparseKinetics_invalidIndex)
            dydt[k] -= nu(i,kk) * rxnRateLaw[i];
      }

      // Products ...
      dydt[ nuk(i,maxReactants) ] += nu(i,maxReactants) * rxnRateLaw[i];
      for (int kk = maxReactants+1; kk < maxSpecies; ++kk){
         const int k = nuk(i,kk);
         if (k == SparseKinetics_invalidIndex) break;
         //if (k != SparseKinetics_invalidIndex)
            dydt[k] += nu(i,kk) * rxnRateLaw[i];
      }
   }

   // Add in the volume factor to convert to the proper units.
   for (int k = 0; k < nspecies; ++k)
      dydt[k] *= VDPD;

   #undef kFor
   #undef kRev
   #undef rxnRateLaw
   #undef conc
   #undef maxReactants
   #undef maxSpecies
   #undef nuk
   #undef nu
   #undef inu
   #undef isIntegral
   //#undef invalidIndex

   return 0;
}

/* ---------------------------------------------------------------------- */

/*template <ExecutionSpace Space>
  template <typename SolverType>
    KOKKOS_INLINE_FUNCTION
void FixRxKokkos<Space>::operator()(SolverType, const int &i) const
{
  if (atom->mask[i] & groupbit)
  {
    KK_FLOAT *rwork = new KK_FLOAT[8*nspecies];

    UserRHSData userData;
    userData.kFor = new KK_FLOAT[nreactions];
    userData.rxnRateLaw = new KK_FLOAT[nreactions];

    int ode_counter[4] = { 0 };

    const KK_FLOAT theta = (localTempFlag) ? dpdThetaLocal[i] : atom->dpdTheta[i];

    //Compute the reaction rate constants
    for (int irxn = 0; irxn < nreactions; irxn++)
    {
      if (SolverType::setToZero)
        userData.kFor[irxn] = 0.0;
      else
        userData.kFor[irxn] = Arr[irxn]*pow(theta,nArr[irxn])*exp(-Ea[irxn]/force->boltz/theta);
    }

    if (odeIntegrationFlag == ODE_LAMMPS_RK4)
      rk4(i, rwork, &userData);
    else if (odeIntegrationFlag == ODE_LAMMPS_RKF45)
      rkf45(i, rwork, &userData, ode_counter);

    delete [] rwork;
    delete [] userData.kFor;
    delete [] userData.rxnRateLaw;
  }
} */

/* ---------------------------------------------------------------------- */

template <ExecutionSpace Space>
void FixRxKokkos<Space>::create_kinetics_data(void)
{
  //printf("Inside FixRxKokkos::create_kinetics_data\n");

  memoryKK->create_kokkos(k_kineticsData.k_Arr, nreactions, "Kinetics::Arr");
  memoryKK->create_kokkos(k_kineticsData.k_nArr, nreactions, "Kinetics::nArr");
  memoryKK->create_kokkos(k_kineticsData.k_Ea, nreactions, "Kinetics::Ea");
  h_kineticsData = KineticsHost(k_kineticsData);

  for (int i = 0; i < nreactions; ++i)
  {
    h_kineticsData.Arr[i]  = Arr[i];
    h_kineticsData.nArr[i] = nArr[i];
    h_kineticsData.Ea[i]   = Ea[i];
  }

  k_kineticsData.k_Arr.modify_host();
  k_kineticsData.k_Arr.sync_device();
  k_kineticsData.k_nArr.modify_host();
  k_kineticsData.k_nArr.sync_device();
  k_kineticsData.k_Ea.modify_host();
  k_kineticsData.k_Ea.sync_device();

  if (useSparseKinetics)
  {

    memoryKK->create_kokkos(k_kineticsData.k_nu, nreactions, sparseKinetics_maxSpecies, "Kinetics::nu");
    memoryKK->create_kokkos(k_kineticsData.k_nuk, nreactions, sparseKinetics_maxSpecies, "Kinetics::nuk");
    h_kineticsData = KineticsHost(k_kineticsData);

    for (int i = 0; i < nreactions; ++i)
      for (int k = 0; k < sparseKinetics_maxSpecies; ++k)
      {
        h_kineticsData.nu (i,k) = sparseKinetics_nu [i][k];
        h_kineticsData.nuk(i,k) = sparseKinetics_nuk[i][k];
      }

    k_kineticsData.k_nu.modify_host();
    k_kineticsData.k_nu.sync_device();
    k_kineticsData.k_nuk.modify_host();
    k_kineticsData.k_nuk.sync_device();

    if (SparseKinetics_enableIntegralReactions)
    {
      memoryKK->create_kokkos(k_kineticsData.k_inu, nreactions, sparseKinetics_maxSpecies, "Kinetics::inu");
      memoryKK->create_kokkos(k_kineticsData.k_isIntegral, nreactions, "Kinetics::isIntegral");
      h_kineticsData = KineticsHost(k_kineticsData);

      for (int i = 0; i < nreactions; ++i)
      {
        h_kineticsData.isIntegral(i) = sparseKinetics_isIntegralReaction[i];

        for (int k = 0; k < sparseKinetics_maxSpecies; ++k)
          h_kineticsData.inu(i,k) = sparseKinetics_inu[i][k];
      }

      k_kineticsData.k_inu.modify_host();
      k_kineticsData.k_inu.sync_device();
      k_kineticsData.k_isIntegral.modify_host();
      k_kineticsData.k_isIntegral.sync_device();
    }
  }

  //else
  //{

    // Dense option
    memoryKK->create_kokkos(k_kineticsData.k_stoich, nreactions, nspecies, "Kinetics::stoich");
    memoryKK->create_kokkos(k_kineticsData.k_stoichReactants, nreactions, nspecies, "Kinetics::stoichReactants");
    memoryKK->create_kokkos(k_kineticsData.k_stoichProducts, nreactions, nspecies, "Kinetics::stoichProducts");
    h_kineticsData = KineticsHost(k_kineticsData);

    for (int i = 0; i < nreactions; ++i)
      for (int k = 0; k < nspecies; ++k)
      {
        h_kineticsData.stoich(i,k) = stoich[i][k];
        h_kineticsData.stoichReactants(i,k) = stoichReactants[i][k];
        h_kineticsData.stoichProducts(i,k) = stoichProducts[i][k];
      }

    k_kineticsData.k_stoich.modify_host();
    k_kineticsData.k_stoich.sync_device();
    k_kineticsData.k_stoichReactants.modify_host();
    k_kineticsData.k_stoichReactants.sync_device();
    k_kineticsData.k_stoichProducts.modify_host();
    k_kineticsData.k_stoichProducts.sync_device();

  //}

  d_kineticsData = KineticsDevice(k_kineticsData);
  update_kinetics_data = false;
}

/* ---------------------------------------------------------------------- */

template <ExecutionSpace Space>
void FixRxKokkos<Space>::setup_pre_force(int vflag)
{
  //printf("Inside FixRxKokkos<Space>::setup_pre_force restartFlag= %d\n", my_restartFlag);

  if (my_restartFlag)
    my_restartFlag = 0;
  else
    this->solve_reactions(vflag, false);
}

/* ---------------------------------------------------------------------- */

template <ExecutionSpace Space>
void FixRxKokkos<Space>::pre_force(int vflag)
{
  //printf("Inside FixRxKokkos<Space>::pre_force localTempFlag= %d\n", localTempFlag);

  this->solve_reactions(vflag, true);
}

/* ---------------------------------------------------------------------- */

template <ExecutionSpace Space>
  KOKKOS_INLINE_FUNCTION
void FixRxKokkos<Space>::operator()(Tag_FixRxKokkos_zeroCounterViews, const int& i) const
{
  d_diagnosticCounterPerODEnSteps(i) = 0;
  d_diagnosticCounterPerODEnFuncs(i) = 0;
}

/* ---------------------------------------------------------------------- */

template <ExecutionSpace Space>
  template <bool ZERO_RATES>
  KOKKOS_INLINE_FUNCTION
void FixRxKokkos<Space>::operator()(Tag_FixRxKokkos_solveSystems<ZERO_RATES>, const int& i, CounterType& counter) const
{
  if (d_mask(i) & groupbit)
  {
    StridedArrayType<SPACE_FLOAT,1> y(d_scratchSpace.data() + scratchSpaceSize * i);
    StridedArrayType<SPACE_FLOAT,1> rwork(&y[nspecies]);

    UserRHSDataKokkos<1> userData;
    userData.kFor.m_data = &(rwork[7*nspecies]);
    userData.rxnRateLaw.m_data = &(userData.kFor[ nreactions ]);

    CounterType counter_i;

    const KK_FLOAT theta = (localTempFlag) ? d_dpdThetaLocal(i) : d_dpdTheta(i);

    //Compute the reaction rate constants
    for (int irxn = 0; irxn < nreactions; irxn++)
    {
      if (ZERO_RATES)
        userData.kFor[irxn] = 0.0;
      else
      {
        userData.kFor[irxn] = d_kineticsData.Arr(irxn) *
                               pow(theta, d_kineticsData.nArr(irxn)) *
                               exp(-d_kineticsData.Ea(irxn) / boltz / theta);
      }
    }

    // Update ConcOld and initialize the ODE solution vector y[].
    for (int ispecies = 0; ispecies < nspecies; ispecies++)
    {
      const KK_FLOAT tmp = d_dvector(ispecies, i);
      d_dvector(ispecies+nspecies, i) = tmp;
      y[ispecies] = tmp;
    }

    // Solver the ODE system.
    if (odeIntegrationFlag == ODE_LAMMPS_RK4)
    {
      k_rk4(t_stop, y, rwork, userData);
    }
    else if (odeIntegrationFlag == ODE_LAMMPS_RKF45)
    {
      k_rkf45(nspecies, t_stop, y, rwork, userData, counter_i);

      if (diagnosticFrequency == 1)
      {
        d_diagnosticCounterPerODEnSteps(i) = counter_i.nSteps;
        d_diagnosticCounterPerODEnFuncs(i) = counter_i.nFuncs;
      }
    }

    // Store the solution back in dvector.
    for (int ispecies = 0; ispecies < nspecies; ispecies++)
    {
      if (y[ispecies] < -1.0e-10)
      {
        //error->one(FLERR,"Computed concentration in RK solver is < -1.0e-10");
        DualViewHelper<Space>::view(k_error_flag)() = 2;
        // This should be an atomic update.
      }
      else if (y[ispecies] < MY_EPSILON)
        y[ispecies] = 0.0;

      d_dvector(ispecies,i) = y[ispecies];
    }

    // Update the iteration statistics counter. Is this unique for each iteration?
    counter += counter_i;

  } // if
}

/* ---------------------------------------------------------------------- */

template <ExecutionSpace Space>
void FixRxKokkos<Space>::solve_reactions(const int vflag, const bool isPreForce)
{
  //printf("Inside FixRxKokkos<Space>::solve_reactions localTempFlag= %d isPreForce= %s\n", localTempFlag, isPreForce ? "True" : "false");

  copymode = 1;

  if (update_kinetics_data)
    create_kinetics_data();

  TimerType timer_start = getTimeStamp();

  //const int nlocal = atom->nlocal;
  this->nlocal = atom->nlocal;
  const int nghost = atom->nghost;
  const int newton_pair = force->newton_pair;

  // Set the forward rates to zero if acting as setup_pre_force.
  const bool setRatesToZero = (isPreForce == false);

  if (localTempFlag)
  {
    const int count = nlocal + (newton_pair ? nghost : 0);

    if (count > DualViewHelper<Space>::view(k_dpdThetaLocal).extent(0)) {
      memoryKK->destroy_kokkos (k_dpdThetaLocal, dpdThetaLocal);
      memoryKK->create_kokkos (k_dpdThetaLocal, dpdThetaLocal, count, "FixRxKokkos::dpdThetaLocal");
      this->d_dpdThetaLocal = DualViewHelper<Space>::view(k_dpdThetaLocal);
      this->h_dpdThetaLocal = k_dpdThetaLocal.h_view;
    }

    const int neighflag = lmp->kokkos->neighflag;

#define _template_switch(_wtflag, _localTempFlag) { \
       if (neighflag == HALF) \
          if (newton_pair) \
             computeLocalTemperature<_wtflag, _localTempFlag, true , HALF> (); \
          else \
             computeLocalTemperature<_wtflag, _localTempFlag, false, HALF> (); \
       else if (neighflag == HALFTHREAD) \
          if (newton_pair) \
             computeLocalTemperature<_wtflag, _localTempFlag, true , HALFTHREAD> (); \
          else \
             computeLocalTemperature<_wtflag, _localTempFlag, false, HALFTHREAD> (); \
       else if (neighflag == FULL) \
          if (newton_pair) \
             computeLocalTemperature<_wtflag, _localTempFlag, true , FULL> (); \
          else \
             computeLocalTemperature<_wtflag, _localTempFlag, false, FULL> (); \
    }

    // Are there is no other options than wtFlag = (0)LUCY and localTempFlag = NONE : HARMONIC?
    if (localTempFlag == HARMONIC) {
       _template_switch(LUCY, HARMONIC)
    }
    else {
       _template_switch(LUCY, NONE)
    }
#undef _template_switch
  }

  TimerType timer_localTemperature = getTimeStamp();

  // Total counters from the ODE solvers.
  CounterType TotalCounters;

  // Set data needed in the operators.
  // ...

  // Local references to the atomKK objects.
  //typename AT::t_float_1d d_dpdTheta = DualViewHelper<Space>::view(atomKK->k_dpdTheta);
  //typename AT::t_float_2d  d_dvector  = DualViewHelper<Space>::view(atomKK->k_dvector);
  //typename AT::t_int_1d    d_mask     = DualViewHelper<Space>::view(atomKK->k_mask);
  this->d_dpdTheta = DualViewHelper<Space>::view(atomKK->k_dpdTheta);
  this->d_dvector  = DualViewHelper<Space>::view(atomKK->k_dvector);
  this->d_mask     = DualViewHelper<Space>::view(atomKK->k_mask);

  // Get up-to-date data.
  atomKK->sync(execution_space, MASK_MASK | DVECTOR_MASK | DPDTHETA_MASK);

  // Set some constants outside of the parallel_for
  //const KK_FLOAT boltz = force->boltz;
  //const KK_FLOAT t_stop = update->dt; // DPD time-step and integration length.
  this->boltz = force->boltz;
  this->t_stop = update->dt; // DPD time-step and integration length.

  // Average DPD volume. Used in the RHS function.
  this->VDPD = domain->xprd * domain->yprd * domain->zprd / atom->natoms;

  if (odeIntegrationFlag == ODE_LAMMPS_RKF45 && diagnosticFrequency == 1)
  {
    memoryKK->create_kokkos (k_diagnosticCounterPerODEnSteps, diagnosticCounterPerODEnSteps, nlocal, "FixRxKokkos::diagnosticCounterPerODEnSteps");
    memoryKK->create_kokkos (k_diagnosticCounterPerODEnFuncs, diagnosticCounterPerODEnFuncs, nlocal, "FixRxKokkos::diagnosticCounterPerODEnFuncs");

    d_diagnosticCounterPerODEnSteps = DualViewHelper<Space>::view(k_diagnosticCounterPerODEnSteps);
    d_diagnosticCounterPerODEnFuncs = DualViewHelper<Space>::view(k_diagnosticCounterPerODEnFuncs);

    Kokkos::parallel_for (Kokkos::RangePolicy<DeviceType, Tag_FixRxKokkos_zeroCounterViews>(0,nlocal), *this);
    //Kokkos::parallel_for (nlocal,
    //      LAMMPS_LAMBDA(const int i)
    //      {
    //         d_diagnosticCounterPerODEnSteps(i) = 0;
    //         d_diagnosticCounterPerODEnFuncs(i) = 0;
    //      }
    //   );
  }

  // Error flag for any failures.
  //DAT::tdual_int_scalar k_error_flag("pair:error_flag");

  // Initialize and sync the device flag.
  k_error_flag.h_view() = 0;
  k_error_flag.modify_host();
  DualViewHelper<Space>::sync(k_error_flag);

  // Create scratch array space.
  //const size_t scratchSpaceSize = (8*nspecies + 2*nreactions);
  this->scratchSpaceSize = (8*nspecies + 2*nreactions);
  //KK_FLOAT *scratchSpace = new KK_FLOAT[ scratchSpaceSize * nlocal ];

  //typename AT::t_float_1d d_scratchSpace("d_scratchSpace", scratchSpaceSize * nlocal);
  if (nlocal*scratchSpaceSize > d_scratchSpace.extent(0)) {
    memoryKK->destroy_kokkos (d_scratchSpace);
    memoryKK->create_kokkos (d_scratchSpace, nlocal*scratchSpaceSize, "FixRxKokkos::d_scratchSpace");
  }

#if 0
  Kokkos::parallel_reduce(nlocal, LAMMPS_LAMBDA(int i, CounterType &counter)
    {
      if (d_mask(i) & groupbit)
      {
        //KK_FLOAT *y = new KK_FLOAT[8*nspecies];
        //KK_FLOAT *rwork = y + nspecies;

        //StridedArrayType<KK_FLOAT,1> _y(y);
        //StridedArrayType<KK_FLOAT,1> _rwork(rwork);

        StridedArrayType<KK_FLOAT,1> y(d_scratchSpace.data() + scratchSpaceSize * i);
        StridedArrayType<KK_FLOAT,1> rwork(&y[nspecies]);

        //UserRHSData userData;
        //userData.kFor = new KK_FLOAT[nreactions];
        //userData.rxnRateLaw = new KK_FLOAT[nreactions];

        //UserRHSDataKokkos<1> userDataKokkos;
        //userDataKokkos.kFor.m_data = userData.kFor;
        //userDataKokkos.rxnRateLaw.m_data = userData.rxnRateLaw;

        UserRHSDataKokkos<1> userData;
        userData.kFor.m_data = &(rwork[7*nspecies]);
        userData.rxnRateLaw.m_data = &(userData.kFor[ nreactions ]);

        CounterType counter_i;

        const KK_FLOAT theta = (localTempFlag) ? d_dpdThetaLocal(i) : d_dpdTheta(i);

        //Compute the reaction rate constants
        for (int irxn = 0; irxn < nreactions; irxn++)
        {
          if (setRatesToZero)
            userData.kFor[irxn] = 0.0;
          else
          {
            userData.kFor[irxn] = d_kineticsData.Arr(irxn) *
                                   pow(theta, d_kineticsData.nArr(irxn)) *
                                   exp(-d_kineticsData.Ea(irxn) / boltz / theta);
          }
        }

        // Update ConcOld and initialize the ODE solution vector y[].
        for (int ispecies = 0; ispecies < nspecies; ispecies++)
        {
          const KK_FLOAT tmp = d_dvector(ispecies, i);
          d_dvector(ispecies+nspecies, i) = tmp;
          y[ispecies] = tmp;
        }

        // Solver the ODE system.
        if (odeIntegrationFlag == ODE_LAMMPS_RK4)
        {
          k_rk4(t_stop, y, rwork, userData);
        }
        else if (odeIntegrationFlag == ODE_LAMMPS_RKF45)
        {
          k_rkf45(nspecies, t_stop, y, rwork, userData, counter_i);

          if (diagnosticFrequency == 1)
          {
            d_diagnosticCounterPerODEnSteps(i) = counter_i.nSteps;
            d_diagnosticCounterPerODEnFuncs(i) = counter_i.nFuncs;
          }
        }

        // Store the solution back in dvector.
        for (int ispecies = 0; ispecies < nspecies; ispecies++)
        {
          if (y[ispecies] < -1.0e-10)
          {
            //error->one(FLERR,"Computed concentration in RK solver is < -1.0e-10");
            DualViewHelper<Space>::view(k_error_flag)() = 2;
            // This should be an atomic update.
          }
          else if (y[ispecies] < MY_EPSILON)
            y[ispecies] = 0.0;

          d_dvector(ispecies,i) = y[ispecies];
        }

        //delete [] y;
        //delete [] userData.kFor;
        //delete [] userData.rxnRateLaw;

        // Update the iteration statistics counter. Is this unique for each iteration?
        counter += counter_i;

      } // if
    } // parallel_for lambda-body

    , TotalCounters // reduction value for all iterations.
  );
#else
  if (setRatesToZero)
    Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, Tag_FixRxKokkos_solveSystems<true > >(0,nlocal), *this, TotalCounters);
  else
    Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, Tag_FixRxKokkos_solveSystems<false> >(0,nlocal), *this, TotalCounters);
#endif

  TimerType timer_ODE = getTimeStamp();

  // Check the error flag for any failures.
  DualViewHelper<Space>::modify(k_error_flag);
  k_error_flag.sync_host();
  if (k_error_flag.h_view() == 2)
    error->one(FLERR,"Computed concentration in RK solver is < -1.0e-10");

  // Signal that dvector has been modified on this execution space.
  atomKK->modified(execution_space, DVECTOR_MASK);

  // Communicate the updated species data to all nodes
  atomKK->sync (Host, DVECTOR_MASK);

  comm->forward_comm_fix(this);

  atomKK->modified (Host, DVECTOR_MASK);

  TimerType timer_stop = getTimeStamp();

  KK_FLOAT time_ODE = getElapsedTime(timer_localTemperature, timer_ODE);

  //printf("me= %d kokkos total= %g temp= %g ode= %g comm= %g nlocal= %d nfc= %d %d\n", comm->me,
  //                       getElapsedTime(timer_start, timer_stop),
  //                       getElapsedTime(timer_start, timer_localTemperature),
  //                       getElapsedTime(timer_localTemperature, timer_ODE),
  //                       getElapsedTime(timer_ODE, timer_stop), nlocal, TotalCounters.nFuncs, TotalCounters.nSteps);

  // Warn the user if a failure was detected in the ODE solver.
  if (TotalCounters.nFails > 0){
    char sbuf[128];
    sprintf(sbuf,"in FixRX::pre_force, ODE solver failed for %d atoms.", TotalCounters.nFails);
    error->warning(FLERR, sbuf);
  }

  // Compute and report ODE diagnostics, if requested.
  if (odeIntegrationFlag == ODE_LAMMPS_RKF45 && diagnosticFrequency != 0)
  {
    // Update the counters.
    diagnosticCounter[StepSum] += TotalCounters.nSteps;
    diagnosticCounter[FuncSum] += TotalCounters.nFuncs;
    diagnosticCounter[TimeSum] += time_ODE;
    diagnosticCounter[AtomSum] += nlocal;
    diagnosticCounter[numDiagnosticCounters-1] ++;

    if ((diagnosticFrequency > 0 &&
               ((update->ntimestep - update->firststep) % diagnosticFrequency) == 0) ||
         (diagnosticFrequency < 0 && update->ntimestep == update->laststep))
      this->odeDiagnostics();
  }

  copymode = 0;
}

/* ---------------------------------------------------------------------- */

template <ExecutionSpace Space>
void FixRxKokkos<Space>::odeDiagnostics(void)
{
  TimerType timer_start = getTimeStamp();

  // Compute:
  // 1) Average # of ODE integrator steps and RHS evaluations per atom globally.
  // 2) RMS     # of  ...
  // 3) Average # of ODE steps and RHS evaluations per MPI task.
  // 4) RMS     # of ODE steps and RHS evaluations per MPI task.
  // 5) MAX     # of ODE steps and RHS evaluations per MPI task.
  //
  // ... 1,2 are for ODE control diagnostics.
  // ... 3-5 are for load balancing diagnostics.
  //
  // To do this, we'll need to
  // a) Allreduce (sum) the sum of nSteps / nFuncs. Dividing by atom->natoms
  //    gives the avg # of steps/funcs per atom globally.
  // b) Reduce (sum) to root the sum of squares of the differences.
  //    i) Sum_i (steps_i - avg_steps_global)^2
  //   ii) Sum_i (funcs_i - avg_funcs_global)^2
  //  iii) (avg_steps_local - avg_steps_global)^2
  //   iv) (avg_funcs_local - avg_funcs_global)^2

  const int numCounters = numDiagnosticCounters-1;

  // # of time-steps for averaging.
  const int nTimes = this->diagnosticCounter[numDiagnosticCounters-1];

  // # of ODE's per time-step (on average).
  //const int nODEs  = this->diagnosticCounter[AtomSum] / nTimes;

  // Sum up the sums from each task.
  KK_FLOAT sums[numCounters];
  KK_FLOAT my_vals[numCounters];
  KK_FLOAT max_per_proc[numCounters];
  KK_FLOAT min_per_proc[numCounters];

  // Compute counters per dpd time-step.
  for (int i = 0; i < numCounters; ++i){
    my_vals[i] = this->diagnosticCounter[i] / nTimes;
    //printf("my sum[%d] = %f %d\n", i, my_vals[i], comm->me);
  }

  MPI_Allreduce (my_vals, sums, numCounters, MPI_DOUBLE, MPI_SUM, world);

  MPI_Reduce (my_vals, max_per_proc, numCounters, MPI_DOUBLE, MPI_MAX, 0, world);
  MPI_Reduce (my_vals, min_per_proc, numCounters, MPI_DOUBLE, MPI_MIN, 0, world);

  const KK_FLOAT nODEs = sums[numCounters-1];

  KK_FLOAT avg_per_atom[numCounters], avg_per_proc[numCounters];

  // Averages per-ODE and per-proc per time-step.
  for (int i = 0; i < numCounters; ++i){
    avg_per_atom[i] = sums[i] / nODEs;
    avg_per_proc[i] = sums[i] / comm->nprocs;
  }

  // Sum up the differences from each task.
  KK_FLOAT sum_sq[2*numCounters];
  KK_FLOAT my_sum_sq[2*numCounters];
  for (int i = 0; i < numCounters; ++i){
    KK_FLOAT diff_i = my_vals[i] - avg_per_proc[i];
    my_sum_sq[i] = diff_i * diff_i;
  }

  KK_FLOAT max_per_ODE[numCounters], min_per_ODE[numCounters];

  // Process the per-ODE RMS of the # of steps/funcs
  if (diagnosticFrequency == 1)
  {
    h_diagnosticCounterPerODEnSteps = k_diagnosticCounterPerODEnSteps.h_view;
    h_diagnosticCounterPerODEnFuncs = k_diagnosticCounterPerODEnFuncs.h_view;

    Kokkos::deep_copy(h_diagnosticCounterPerODEnSteps, d_diagnosticCounterPerODEnSteps);
    Kokkos::deep_copy(h_diagnosticCounterPerODEnFuncs, d_diagnosticCounterPerODEnFuncs);

    KK_FLOAT my_max[numCounters], my_min[numCounters];

    //const int nlocal = atom->nlocal;
    nlocal = atom->nlocal;
    HAT::t_int_1d h_mask = atomKK->k_mask.h_view;

    for (int i = 0; i < numCounters; ++i)
    {
      my_sum_sq[i+numCounters] = 0;
      my_max[i] = 0;
      my_min[i] = DBL_MAX;
    }

    for (int j = 0; j < nlocal; ++j)
      if (h_mask(j) & groupbit)
      {
        int nSteps = h_diagnosticCounterPerODEnSteps(j);
        KK_FLOAT diff_nSteps = KK_FLOAT(nSteps) - avg_per_atom[StepSum];
        my_sum_sq[StepSum+numCounters] += diff_nSteps*diff_nSteps;
        my_max[StepSum] = std::max(my_max[StepSum], (KK_FLOAT)nSteps);
        my_min[StepSum] = std::min(my_min[StepSum], (KK_FLOAT)nSteps);

        int nFuncs = h_diagnosticCounterPerODEnFuncs(j);
        KK_FLOAT diff_nFuncs = KK_FLOAT(nFuncs) - avg_per_atom[FuncSum];
        my_sum_sq[FuncSum+numCounters] += diff_nFuncs*diff_nFuncs;

        my_max[FuncSum] = std::max(my_max[FuncSum], (KK_FLOAT)nFuncs);
        my_min[FuncSum] = std::min(my_min[FuncSum], (KK_FLOAT)nFuncs);
      }

    memoryKK->destroy_kokkos(k_diagnosticCounterPerODEnSteps, diagnosticCounterPerODEnSteps);
    memoryKK->destroy_kokkos(k_diagnosticCounterPerODEnFuncs, diagnosticCounterPerODEnFuncs);

    MPI_Reduce (my_sum_sq, sum_sq, 2*numCounters, MPI_DOUBLE, MPI_SUM, 0, world);

    MPI_Reduce (my_max, max_per_ODE, numCounters, MPI_DOUBLE, MPI_MAX, 0, world);
    MPI_Reduce (my_min, min_per_ODE, numCounters, MPI_DOUBLE, MPI_MIN, 0, world);
  }
  else
    MPI_Reduce (my_sum_sq, sum_sq, numCounters, MPI_DOUBLE, MPI_SUM, 0, world);

  TimerType timer_stop = getTimeStamp();
  KK_FLOAT time_local = getElapsedTime(timer_start, timer_stop);

  if (comm->me == 0){
    char smesg[128];

#define print_mesg(smesg) {\
    if (screen)  fprintf(screen,"%s\n", smesg); \
    if (logfile) fprintf(logfile,"%s\n", smesg); }

    sprintf(smesg, "FixRX::ODE Diagnostics:  # of iters  |# of rhs evals| run-time (sec) | # atoms");
    print_mesg(smesg);

    sprintf(smesg, "         AVG per ODE  : %-12.5g | %-12.5g | %-12.5g", avg_per_atom[0], avg_per_atom[1], avg_per_atom[2]);
    print_mesg(smesg);

    // only valid for single time-step!
    if (diagnosticFrequency == 1){
      KK_FLOAT rms_per_ODE[numCounters];
      for (int i = 0; i < numCounters; ++i)
        rms_per_ODE[i] = sqrt(sum_sq[i+numCounters] / nODEs);

      sprintf(smesg, "         RMS per ODE  : %-12.5g | %-12.5g ", rms_per_ODE[0], rms_per_ODE[1]);
      print_mesg(smesg);

      sprintf(smesg, "         MAX per ODE  : %-12.5g | %-12.5g ", max_per_ODE[0], max_per_ODE[1]);
      print_mesg(smesg);

      sprintf(smesg, "         MIN per ODE  : %-12.5g | %-12.5g ", min_per_ODE[0], min_per_ODE[1]);
      print_mesg(smesg);
    }

    sprintf(smesg, "         AVG per Proc : %-12.5g | %-12.5g | %-12.5g | %-12.5g", avg_per_proc[StepSum], avg_per_proc[FuncSum], avg_per_proc[TimeSum], avg_per_proc[AtomSum]);
    print_mesg(smesg);

    if (comm->nprocs > 1){
      KK_FLOAT rms_per_proc[numCounters];
      for (int i = 0; i < numCounters; ++i)
        rms_per_proc[i] = sqrt(sum_sq[i] / comm->nprocs);

      sprintf(smesg, "         RMS per Proc : %-12.5g | %-12.5g | %-12.5g | %-12.5g", rms_per_proc[0], rms_per_proc[1], rms_per_proc[2], rms_per_proc[AtomSum]);
      print_mesg(smesg);

      sprintf(smesg, "         MAX per Proc : %-12.5g | %-12.5g | %-12.5g | %-12.5g", max_per_proc[0], max_per_proc[1], max_per_proc[2], max_per_proc[AtomSum]);
      print_mesg(smesg);

      sprintf(smesg, "         MIN per Proc : %-12.5g | %-12.5g | %-12.5g | %-12.5g", min_per_proc[0], min_per_proc[1], min_per_proc[2], min_per_proc[AtomSum]);
      print_mesg(smesg);
    }

    sprintf(smesg, "  AVG'd over %d time-steps", nTimes);
    print_mesg(smesg);
    sprintf(smesg, "  AVG'ing took %g sec", time_local);
    print_mesg(smesg);

#undef print_mesg

  }

  // Reset the counters.
  for (int i = 0; i < numDiagnosticCounters; ++i)
    diagnosticCounter[i] = 0;

  return;
}

/* ---------------------------------------------------------------------- */

template <ExecutionSpace Space>
  KOKKOS_INLINE_FUNCTION
void FixRxKokkos<Space>::operator()(Tag_FixRxKokkos_zeroTemperatureViews, const int& i) const
{
  d_sumWeights(i) = 0.0;
  d_dpdThetaLocal(i) = 0.0;
}

/* ---------------------------------------------------------------------- */

template <ExecutionSpace Space>
  template <int WT_FLAG, bool NEWTON_PAIR, int NEIGHFLAG>
  KOKKOS_INLINE_FUNCTION
void FixRxKokkos<Space>::operator()(Tag_FixRxKokkos_firstPairOperator<WT_FLAG,NEWTON_PAIR,NEIGHFLAG>, const int& ii) const
{
  // Create an atomic view of sumWeights and dpdThetaLocal. Only needed
  // for Half/thread scenarios.
  typedef Kokkos::View<typename AT::t_float_1d::data_type, typename AT::t_float_1d::array_layout, typename KKDevice<DeviceType>::value, Kokkos::MemoryTraits< AtomicF< NEIGHFLAG >::value> > AtomicViewType;

  AtomicViewType a_dpdThetaLocal = d_dpdThetaLocal;
  AtomicViewType a_sumWeights    = d_sumWeights;

  // Local scalar accumulators.
  KK_FLOAT i_dpdThetaLocal = 0.0;
  KK_FLOAT i_sumWeights    = 0.0;

  const int i = d_ilist(ii);

  const KK_FLOAT xtmp = d_x(i,0);
  const KK_FLOAT ytmp = d_x(i,1);
  const KK_FLOAT ztmp = d_x(i,2);
  const int itype = d_type(i);

  const int jnum = d_numneigh(i);

  for (int jj = 0; jj < jnum; jj++)
  {
    const int j = (d_neighbors(i,jj) & NEIGHMASK);
    const int jtype = d_type(j);

    const KK_FLOAT delx = xtmp - d_x(j,0);
    const KK_FLOAT dely = ytmp - d_x(j,1);
    const KK_FLOAT delz = ztmp - d_x(j,2);
    const KK_FLOAT rsq = delx*delx + dely*dely + delz*delz;

    const KK_FLOAT cutsq_ij = d_cutsq(itype,jtype);

    if (rsq < cutsq_ij)
    {
      const KK_FLOAT rcut = sqrt(cutsq_ij);
      KK_FLOAT rij = sqrt(rsq);
      KK_FLOAT ratio = rij/rcut;

      KK_FLOAT wij = 0.0;

      // Lucy's Weight Function
      if (WT_FLAG == LUCY)
      {
        wij = (1.0+3.0*ratio) * (1.0-ratio)*(1.0-ratio)*(1.0-ratio);
        i_dpdThetaLocal += wij / d_dpdTheta(j);
        if ((NEIGHFLAG==HALF || NEIGHFLAG==HALFTHREAD) && (NEWTON_PAIR || j < nlocal))
          a_dpdThetaLocal(j) += wij / d_dpdTheta(i);
      }

      i_sumWeights += wij;
      if ((NEIGHFLAG==HALF || NEIGHFLAG==HALFTHREAD) && (NEWTON_PAIR || j < nlocal))
        a_sumWeights(j) += wij;
    }
  }

  // Update, don't assign, the array value (because another iteration may have hit it).
  a_dpdThetaLocal(i) += i_dpdThetaLocal;
  a_sumWeights(i) += i_sumWeights;
}

/* ---------------------------------------------------------------------- */

template <ExecutionSpace Space>
  template <int WT_FLAG, int LOCAL_TEMP_FLAG>
  KOKKOS_INLINE_FUNCTION
void FixRxKokkos<Space>::operator()(Tag_FixRxKokkos_2ndPairOperator<WT_FLAG,LOCAL_TEMP_FLAG>, const int& i) const
{
  KK_FLOAT wij = 0.0;

  // Lucy Weight Function
  if (WT_FLAG == LUCY)
  {
    wij = 1.0;
    d_dpdThetaLocal(i) += wij / d_dpdTheta(i);
  }
  d_sumWeights(i) += wij;

  // Normalized local temperature
  d_dpdThetaLocal(i) = d_dpdThetaLocal(i) / d_sumWeights(i);

  if (LOCAL_TEMP_FLAG == HARMONIC)
    d_dpdThetaLocal(i) = 1.0 / d_dpdThetaLocal(i);
}

/* ---------------------------------------------------------------------- */

template <ExecutionSpace Space>
  template <int WT_FLAG, int LOCAL_TEMP_FLAG, bool NEWTON_PAIR, int NEIGHFLAG>
void FixRxKokkos<Space>::computeLocalTemperature()
{
  //typename AT::t_float_1d_3_randomread d_x        = DualViewHelper<Space>::view(atomKK->k_x);
  //typename AT::t_int_1d_randomread  d_type     = DualViewHelper<Space>::view(atomKK->k_type);
  //typename AT::t_float_1d          d_dpdTheta = DualViewHelper<Space>::view(atomKK->k_dpdTheta);
  d_x        = DualViewHelper<Space>::view(atomKK->k_x);
  d_type     = DualViewHelper<Space>::view(atomKK->k_type);
  d_dpdTheta = DualViewHelper<Space>::view(atomKK->k_dpdTheta);

  atomKK->sync(execution_space, X_MASK | TYPE_MASK | DPDTHETA_MASK);

  //const int nlocal = atom->nlocal;
  nlocal = atom->nlocal;
  const int nghost = atom->nghost;

  //printf("Inside FixRxKokkos::computeLocalTemperature: %d %d %d %d %d %d %d\n", WT_FLAG, LOCAL_TEMP_FLAG, NEWTON_PAIR, (int)lmp->kokkos->neighflag, NEIGHFLAG, nlocal, nghost);

  // Pull from pairDPDE. The pairDPDEKK objects are protected so recreate here for now.
  //DualViewHelper<Space>::sync(pairDPDEKK->k_cutsq);
  //typename AT::t_float_2d d_cutsq = DualViewHelper<Space>::view(pairDPDEKK->k_cutsq);

  //!< Copies pulled from pairDPDE for local use since pairDPDEKK's objects are protected.
  //DAT::tdual_float_2d k_cutsq;
  //typename AT::t_float_2d     d_cutsq;
  //KK_FLOAT **h_cutsq;

  {
    const int ntypes = atom->ntypes;

    //memoryKK->create_kokkos (k_cutsq, h_cutsq, ntypes+1, ntypes+1, "pair:cutsq");
    if (ntypes+1 > k_cutsq.extent(0)) {
      memoryKK->destroy_kokkos (k_cutsq);
      memoryKK->create_kokkos (k_cutsq, ntypes+1, ntypes+1, "FixRxKokkos::k_cutsq");
      d_cutsq = DualViewHelper<Space>::view(k_cutsq);
    }

    for (int i = 1; i <= ntypes; ++i)
      for (int j = i; j <= ntypes; ++j)
      {
        k_cutsq.h_view(i,j) = pairDPDE->cutsq[i][j];
        k_cutsq.h_view(j,i) = k_cutsq.h_view(i,j);
      }

    k_cutsq.modify_host();
    DualViewHelper<Space>::sync(k_cutsq);
  }

  // Initialize the local temperature weight array
  int sumWeightsCt = nlocal + (NEWTON_PAIR ? nghost : 0);

  //memoryKK->create_kokkos (k_sumWeights, sumWeights, sumWeightsCt, "FixRxKokkos::sumWeights");
  if (sumWeightsCt > DualViewHelper<Space>::view(k_sumWeights).extent(0)) {
    memoryKK->destroy_kokkos(k_sumWeights, sumWeights);
    memoryKK->create_kokkos (k_sumWeights, sumWeightsCt, "FixRxKokkos::sumWeights");
    d_sumWeights = DualViewHelper<Space>::view(k_sumWeights);
    h_sumWeights = k_sumWeights.h_view;
  }

  // Initialize the accumulator to zero ...
  //Kokkos::parallel_for (sumWeightsCt,
  //      LAMMPS_LAMBDA(const int i)
  //      {
  //         d_sumWeights(i) = 0.0;
  //      }
  //   );

  Kokkos::parallel_for (Kokkos::RangePolicy<DeviceType, Tag_FixRxKokkos_zeroTemperatureViews>(0, sumWeightsCt), *this);

  // Local list views. (This isn't working!)
  NeighListKokkos<Space>* k_list = static_cast<NeighListKokkos<Space>*>(list);
  if (not(list->kokkos))
     error->one(FLERR,"list is not a Kokkos list\n");

  //typename AT::t_neighbors_2d d_neighbors = k_list->d_neighbors;
  //typename AT::t_int_1d       d_ilist     = k_list->d_ilist;
  //typename AT::t_int_1d       d_numneigh  = k_list->d_numneigh;
  d_neighbors = k_list->d_neighbors;
  d_ilist     = k_list->d_ilist;
  d_numneigh  = k_list->d_numneigh;

  const int inum = list->inum;

  // loop over neighbors of my atoms
#if 0
  Kokkos::parallel_for (inum,
        LAMMPS_LAMBDA(const int ii)
        {
          // Create an atomic view of sumWeights and dpdThetaLocal. Only needed
          // for Half/thread scenarios.
          //typedef Kokkos::View< KK_FLOAT*, typename AT::t_float_1d::array_layout, typename KKDevice<DeviceType>::value, Kokkos::MemoryTraits< AtomicF< NEIGHFLAG >::value> > AtomicViewType;
          typedef Kokkos::View< KK_FLOAT*, typename AT::t_float_1d::array_layout, typename KKDevice<DeviceType>::value, Kokkos::MemoryTraits< AtomicF< NEIGHFLAG >::value> > AtomicViewType;

          AtomicViewType a_dpdThetaLocal = d_dpdThetaLocal;
          AtomicViewType a_sumWeights    = d_sumWeights;

          // Local scalar accumulators.
          KK_FLOAT i_dpdThetaLocal = 0.0;
          KK_FLOAT i_sumWeights    = 0.0;

          const int i = d_ilist(ii);

          const KK_FLOAT xtmp = d_x(i,0);
          const KK_FLOAT ytmp = d_x(i,1);
          const KK_FLOAT ztmp = d_x(i,2);
          const int itype = d_type(i);

          const int jnum = d_numneigh(i);

          for (int jj = 0; jj < jnum; jj++)
          {
            const int j = (d_neighbors(i,jj) & NEIGHMASK);
            const int jtype = d_type(j);

            const KK_FLOAT delx = xtmp - d_x(j,0);
            const KK_FLOAT dely = ytmp - d_x(j,1);
            const KK_FLOAT delz = ztmp - d_x(j,2);
            const KK_FLOAT rsq = delx*delx + dely*dely + delz*delz;

            const KK_FLOAT cutsq_ij = d_cutsq(itype,jtype);

            if (rsq < cutsq_ij)
            {
              const KK_FLOAT rcut = sqrt(cutsq_ij);
              KK_FLOAT rij = sqrt(rsq);
              KK_FLOAT ratio = rij/rcut;

              KK_FLOAT wij = 0.0;

              // Lucy's Weight Function
              if (WT_FLAG == LUCY)
              {
                wij = (1.0+3.0*ratio) * (1.0-ratio)*(1.0-ratio)*(1.0-ratio);
                i_dpdThetaLocal += wij / d_dpdTheta(j);
                if (NEWTON_PAIR || j < nlocal)
                  a_dpdThetaLocal(j) += wij / d_dpdTheta(i);
              }

              i_sumWeights += wij;
              if (NEWTON_PAIR || j < nlocal)
                a_sumWeights(j) += wij;
            }
          }

          // Update, don't assign, the array value (because another iteration may have hit it).
          a_dpdThetaLocal(i) += i_dpdThetaLocal;
          a_sumWeights(i) += i_sumWeights;
        }
     );
#else
  Kokkos::parallel_for (Kokkos::RangePolicy<DeviceType, Tag_FixRxKokkos_firstPairOperator<WT_FLAG, NEWTON_PAIR, NEIGHFLAG> >(0, inum), *this);
#endif

  // Signal that dpdThetaLocal and sumWeights have been modified.
  DualViewHelper<Space>::modify(k_dpdThetaLocal);
  DualViewHelper<Space>::modify(k_sumWeights);

  // Communicate the sum dpdTheta and the weights on the host.
  if (NEWTON_PAIR) comm->reverse_comm_fix(this);

  // Update the device view in case they got changed.
  DualViewHelper<Space>::sync(k_dpdThetaLocal);
  DualViewHelper<Space>::sync(k_sumWeights);

  // self-interaction for local temperature
#if 0
  Kokkos::parallel_for (nlocal,
        LAMMPS_LAMBDA(const int i)
        {
          KK_FLOAT wij = 0.0;

          // Lucy Weight Function
          if (WT_FLAG == LUCY)
          {
            wij = 1.0;
            d_dpdThetaLocal(i) += wij / d_dpdTheta(i);
          }
          d_sumWeights(i) += wij;

          // Normalized local temperature
          d_dpdThetaLocal(i) = d_dpdThetaLocal(i) / d_sumWeights(i);

          if (LOCAL_TEMP_FLAG == HARMONIC)
            d_dpdThetaLocal(i) = 1.0 / d_dpdThetaLocal(i);
        }
     );
#else
  Kokkos::parallel_for (Kokkos::RangePolicy<DeviceType, Tag_FixRxKokkos_2ndPairOperator<WT_FLAG, LOCAL_TEMP_FLAG> >(0, nlocal), *this);
#endif

}

/* ---------------------------------------------------------------------- */

template <ExecutionSpace Space>
int FixRxKokkos<Space>::pack_forward_comm(int n, int *list, double *buf, int pbc_flag, int *pbc)
{
  //printf("inside FixRxKokkos::pack_forward_comm %d\n", comm->me);

  HAT::t_float_2d h_dvector = atomKK->k_dvector.h_view;

  int m = 0;
  for (int ii = 0; ii < n; ii++) {
    const int jj = list[ii];
    for(int ispecies = 0; ispecies < nspecies; ispecies++){
      buf[m++] = h_dvector(ispecies,jj);
      buf[m++] = h_dvector(ispecies+nspecies,jj);
    }
  }

  //printf("done with FixRxKokkos::pack_forward_comm %d\n", comm->me);

  return m;
}

/* ---------------------------------------------------------------------- */

template <ExecutionSpace Space>
void FixRxKokkos<Space>::unpack_forward_comm(int n, int first, double *buf)
{
  //printf("inside FixRxKokkos::unpack_forward_comm %d\n", comm->me);

  HAT::t_float_2d h_dvector = atomKK->k_dvector.h_view;

  const int last = first + n ;
  int m = 0;
  for (int ii = first; ii < last; ii++){
    for (int ispecies = 0; ispecies < nspecies; ispecies++){
      h_dvector(ispecies,ii) = buf[m++];
      h_dvector(ispecies+nspecies,ii) = buf[m++];
    }
  }

  //printf("done with FixRxKokkos::unpack_forward_comm %d\n", comm->me);
}

/* ---------------------------------------------------------------------- */

template <ExecutionSpace Space>
int FixRxKokkos<Space>::pack_reverse_comm(int n, int first, double *buf)
{
  //printf("inside FixRxKokkos::pack_reverse_comm %d %d %d\n", comm->me, first, n);
  // Sync the host view.
  k_dpdThetaLocal.sync_host();
  k_sumWeights.   sync_host();

  const int last = first + n;
  int m = 0;
  for (int i = first; i < last; ++i)
  {
    buf[m++] = h_dpdThetaLocal(i);
    buf[m++] = h_sumWeights(i);
  }
  //printf("done with FixRxKokkos::pack_reverse_comm %d\n", comm->me);

  return m;
}

/* ---------------------------------------------------------------------- */

template <ExecutionSpace Space>
void FixRxKokkos<Space>::unpack_reverse_comm(int n, int *list, double *buf)
{
  // printf("inside FixRxKokkos::unpack_reverse_comm %d\n", comm->me);
  int m = 0;
  for (int i = 0; i < n; i++) {
    const int j = list[i];

    h_dpdThetaLocal(j) += buf[m++];
    h_sumWeights(j) += buf[m++];
  }

  // Signal that the host view has been modified.
  k_dpdThetaLocal.modify_host();
  k_sumWeights.   modify_host();

  // printf("done with FixRxKokkos::unpack_reverse_comm %d\n", comm->me);
}

namespace LAMMPS_NS {
template class FixRxKokkos<Device>;
template class FixRxKokkos<Host>;
}

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

#include <stdio.h>
#include <string.h>
#include "fix_rx_kokkos.h"
#include "atom_masks.h"
#include "atom_kokkos.h"
#include "force.h"
#include "memory.h"
#include "update.h"
#include "respa.h"
#include "modify.h"
#include "neighbor.h"
#include "neigh_list_kokkos.h"
#include "neigh_request.h"
#include "error.h"
#include "math_special.h"

#include <float.h> // DBL_EPSILON

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathSpecial;

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

typedef double TimerType;
TimerType getTimeStamp(void) { return MPI_Wtime(); }
double getElapsedTime( const TimerType &t0, const TimerType &t1) { return t1-t0; }

} // end namespace

/* ---------------------------------------------------------------------- */

template <typename DeviceType>
FixRxKokkos<DeviceType>::FixRxKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixRX(lmp, narg, arg),
  pairDPDEKK(NULL),
  update_kinetics_data(true)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = EMPTY_MASK;
  datamask_modify = EMPTY_MASK;

  printf("Inside FixRxKokkos::FixRxKokkos\n");
}

template <typename DeviceType>
FixRxKokkos<DeviceType>::~FixRxKokkos()
{
  printf("Inside FixRxKokkos::~FixRxKokkos\n");
}

/* ---------------------------------------------------------------------- */

template <typename DeviceType>
void FixRxKokkos<DeviceType>::post_constructor()
{
  // Run the parents and then reset one value.
  FixRX::post_constructor();

  // Need a copy of this
  this->my_restartFlag = modify->fix[modify->nfix-1]->restart_reset;
}

/* ---------------------------------------------------------------------- */

template <typename DeviceType>
void FixRxKokkos<DeviceType>::init()
{
  printf("Inside FixRxKokkos::init\n");

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
    if (strcmp(modify->fix[i]->style,"eos/table/rx") == 0) eos_flag = true;
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
    kokkos_host = Kokkos::Impl::is_same<DeviceType,LMPHostType>::value &&
    !Kokkos::Impl::is_same<DeviceType,LMPDeviceType>::value;
  neighbor->requests[irequest]->
    kokkos_device = Kokkos::Impl::is_same<DeviceType,LMPDeviceType>::value;
 
  if (neighflag == FULL) {
    neighbor->requests[irequest]->full = 1;
    neighbor->requests[irequest]->half = 0;
  } else { //if (neighflag == HALF || neighflag == HALFTHREAD)
    neighbor->requests[irequest]->full = 0;
    neighbor->requests[irequest]->half = 1;
  }
}

/* ---------------------------------------------------------------------- */

template <typename DeviceType>
void FixRxKokkos<DeviceType>::init_list(int, class NeighList* ptr)
{
  printf("Inside FixRxKokkos::init_list\n");
  this->list = ptr;
}

/* ---------------------------------------------------------------------- */

template <typename DeviceType>
void FixRxKokkos<DeviceType>::rk4(const double t_stop, double *y, double *rwork, void* v_params) const
{
  double *k1 = rwork;
  double *k2 = k1 + nspecies;
  double *k3 = k2 + nspecies;
  double *k4 = k3 + nspecies;
  double *yp = k4 + nspecies;

  const int numSteps = minSteps;

  const double h = t_stop / double(numSteps);

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

template <typename DeviceType>
void FixRxKokkos<DeviceType>::rkf45_step (const int neq, const double h, double y[], double y_out[], double rwk[], void* v_param) const
{
   const double c21=0.25;
   const double c31=0.09375;
   const double c32=0.28125;
   const double c41=0.87938097405553;
   const double c42=-3.2771961766045;
   const double c43=3.3208921256258;
   const double c51=2.0324074074074;
   const double c52=-8.0;
   const double c53=7.1734892787524;
   const double c54=-0.20589668615984;
   const double c61=-0.2962962962963;
   const double c62=2.0;
   const double c63=-1.3816764132554;
   const double c64=0.45297270955166;
   const double c65=-0.275;
   const double a1=0.11574074074074;
   const double a3=0.54892787524366;
   const double a4=0.5353313840156;
   const double a5=-0.2;
   const double b1=0.11851851851852;
   const double b3=0.51898635477583;
   const double b4=0.50613149034201;
   const double b5=-0.18;
   const double b6=0.036363636363636;

   // local dependent variables (5 total)
   double* f1 = &rwk[    0];
   double* f2 = &rwk[  neq];
   double* f3 = &rwk[2*neq];
   double* f4 = &rwk[3*neq];
   double* f5 = &rwk[4*neq];
   double* f6 = &rwk[5*neq];

   // scratch for the intermediate solution.
   //double* ytmp = &rwk[6*neq];
   double* ytmp = y_out;

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
      //const double f6 = h * ydot[k];
      f6[k] *= h;

      // 5th-order solution.
      const double r5 = b1*f1[k] + b3*f3[k] + b4*f4[k] + b5*f5[k] + b6*f6[k];

      // 4th-order solution.
      const double r4 = a1*f1[k] + a3*f3[k] + a4*f4[k] + a5*f5[k];

      // Truncation error: difference between 4th and 5th-order solutions.
      rwk[k] = fabs(r5 - r4);

      // Update solution.
    //y_out[k] = y[k] + r5; // Local extrapolation
      y_out[k] = y[k] + r4;
   }

   return;
}

template <typename DeviceType>
int FixRxKokkos<DeviceType>::rkf45_h0
                    (const int neq, const double t, const double t_stop,
                     const double hmin, const double hmax,
                     double& h0, double y[], double rwk[], void* v_params) const
{
   // Set lower and upper bounds on h0, and take geometric mean as first trial value.
   // Exit with this value if the bounds cross each other.

   // Adjust upper bound based on ydot ...
   double hg = sqrt(hmin*hmax);

   //if (hmax < hmin)
   //{
   //   h0 = hg;
   //   return;
   //}

   // Start iteration to find solution to ... {WRMS norm of (h0^2 y'' / 2)} = 1

   double *ydot  = rwk;
   double *y1    = ydot + neq;
   double *ydot1 = y1 + neq;

   const int max_iters = 10;
   bool hnew_is_ok = false;
   double hnew = hg;
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
      double yddnrm = 0.0;
      for (int k = 0; k < neq; k++){
         double ydd = (ydot1[k] - ydot[k]) / hg;
         double wterr = ydd / (relTol * fabs( y[k] ) + absTol);
         yddnrm += wterr * wterr;
      }

      yddnrm = sqrt( yddnrm / double(neq) );

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
      double hrat = hnew / hg;

      // Accept this value ... the bias factor should bring it within range.
      if ( (hrat > 0.5) && (hrat < 2.0) )
         hnew_is_ok = true;

      // If y'' is still bad after a few iterations, just accept h and give up.
      if ( (iter > 1) && hrat > 2.0 ) {
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

template <typename DeviceType>
void FixRxKokkos<DeviceType>::rkf45(const int neq, const double t_stop, double *y, double *rwork, void *v_param, CounterType& counter) const
{
  // Rounding coefficient.
  const double uround = DBL_EPSILON;

  // Adaption limit (shrink or grow)
  const double adaption_limit = 4.0;

  // Safety factor on the adaption. very specific but not necessary .. 0.9 is common.
  const double hsafe = 0.840896415;

  // Time rounding factor.
  const double tround = t_stop * uround;

  // Counters for diagnostics.
  int nst = 0; // # of steps (accepted)
  int nit = 0; // # of iterations total
  int nfe = 0; // # of RHS evaluations

  // Min/Max step-size limits.
  const double h_min = 100.0 * tround;
  const double h_max = (minSteps > 0) ? t_stop / double(minSteps) : t_stop;

  // Set the initial step-size. 0 forces an internal estimate ... stable Euler step size.
  double h = (minSteps > 0) ? t_stop / double(minSteps) : 0.0;

  double t = 0.0;

  if (h < h_min){
    //fprintf(stderr,"hin not implemented yet\n");
    //exit(-1);
    nfe = rkf45_h0 (neq, t, t_stop, h_min, h_max, h, y, rwork, v_param);
  }

  //printf("t= %e t_stop= %e h= %e\n", t, t_stop, h);

  // Integrate until we reach the end time.
  while (fabs(t - t_stop) > tround){
    double *yout = rwork;
    double *eout = yout + neq;

    // Take a trial step.
    rkf45_step (neq, h, y, yout, eout, v_param);

    // Estimate the solution error.
      // ... weighted 2-norm of the error.
      double err2 = 0.0;
      for (int k = 0; k < neq; k++){
        const double wterr = eout[k] / (relTol * fabs( y[k] ) + absTol);
        err2 += wterr * wterr;
      }

    double err = fmax( uround, sqrt( err2 / double(nspecies) ));

    // Accept the solution?
    if (err <= 1.0 || h <= h_min){
      t += h;
      nst++;

      for (int k = 0; k < neq; k++)
        y[k] = yout[k];
    }

    // Adjust h for the next step.
    double hfac = hsafe * sqrt( sqrt( 1.0 / err ) );

    // Limit the adaption.
    hfac = fmax( hfac, 1.0 / adaption_limit );
    hfac = fmin( hfac,       adaption_limit );

    // Apply the adaption factor...
    h *= hfac;

    // Limit h.
    h = fmin( h, h_max );
    h = fmax( h, h_min );

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

template <typename DeviceType>
int FixRxKokkos<DeviceType>::rhs(double t, const double *y, double *dydt, void *params) const
{
  // Use the sparse format instead.
  if (useSparseKinetics)
    return this->rhs_sparse( t, y, dydt, params);
  else
    return this->rhs_dense ( t, y, dydt, params);
}

/* ---------------------------------------------------------------------- */

template <typename DeviceType>
int FixRxKokkos<DeviceType>::rhs_dense(double t, const double *y, double *dydt, void *params) const
{
  UserRHSData *userData = (UserRHSData *) params;

  double *rxnRateLaw = userData->rxnRateLaw;
  double *kFor       = userData->kFor;

  //const double VDPD = domain->xprd * domain->yprd * domain->zprd / atom->natoms;
  //const int nspecies = atom->nspecies_dpd;

  for(int ispecies=0; ispecies<nspecies; ispecies++)
    dydt[ispecies] = 0.0;

  // Construct the reaction rate laws
  for(int jrxn=0; jrxn<nreactions; jrxn++){
    double rxnRateLawForward = kFor[jrxn];

    for(int ispecies=0; ispecies<nspecies; ispecies++){
      const double concentration = y[ispecies]/VDPD;
      rxnRateLawForward *= pow( concentration, d_kineticsData.stoichReactants(jrxn,ispecies) );
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

template <typename DeviceType>
int FixRxKokkos<DeviceType>::rhs_sparse(double t, const double *y, double *dydt, void *v_params) const
{
   UserRHSData *userData = (UserRHSData *) v_params;

   //const double VDPD = domain->xprd * domain->yprd * domain->zprd / atom->natoms;

   #define kFor         (userData->kFor)
   #define kRev         (NULL)
   #define rxnRateLaw   (userData->rxnRateLaw)
   #define conc         (dydt)
   #define maxReactants (this->sparseKinetics_maxReactants)
   #define maxSpecies   (this->sparseKinetics_maxSpecies)
   #define nuk          (this->d_kineticsData.nuk)
   #define nu           (this->d_kineticsData.nu)
   #define inu          (this->d_kineticsData.inu)
   #define isIntegral(idx) ( SparseKinetics_enableIntegralReactions \
                             && this->d_kineticsData.isIntegral(idx) )

   for (int k = 0; k < nspecies; ++k)
      conc[k] = y[k] / VDPD;

   // Construct the reaction rate laws
   for (int i = 0; i < nreactions; ++i)
   {
      double rxnRateLawForward;
      if (isIntegral(i)){
         rxnRateLawForward = kFor[i] * powint( conc[ nuk(i,0) ], inu(i,0) );
         for (int kk = 1; kk < maxReactants; ++kk){
            const int k = nuk(i,kk);
            if (k == SparseKinetics_invalidIndex) break;
            //if (k != SparseKinetics_invalidIndex)
               rxnRateLawForward *= powint( conc[k], inu(i,kk) );
         }
      } else {
         rxnRateLawForward = kFor[i] * pow( conc[ nuk(i,0) ], nu(i,0) );
         for (int kk = 1; kk < maxReactants; ++kk){
            const int k = nuk(i,kk);
            if (k == SparseKinetics_invalidIndex) break;
            //if (k != SparseKinetics_invalidIndex)
               rxnRateLawForward *= pow( conc[k], nu(i,kk) );
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

/*template <typename DeviceType>
  template <typename SolverType>
    KOKKOS_INLINE_FUNCTION
void FixRxKokkos<DeviceType>::operator()(SolverType, const int &i) const
{
  if (atom->mask[i] & groupbit)
  {
    double *rwork = new double[8*nspecies];

    UserRHSData userData;
    userData.kFor = new double[nreactions];
    userData.rxnRateLaw = new double[nreactions];

    int ode_counter[4] = { 0 };

    const double theta = (localTempFlag) ? dpdThetaLocal[i] : atom->dpdTheta[i];

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

template <typename DeviceType>
void FixRxKokkos<DeviceType>::create_kinetics_data(void)
{
  printf("Inside FixRxKokkos::create_kinetics_data\n");

  memory->create_kokkos( d_kineticsData.Arr, h_kineticsData.Arr, nreactions, "KineticsType::Arr");
  memory->create_kokkos( d_kineticsData.nArr, h_kineticsData.nArr, nreactions, "KineticsType::nArr");
  memory->create_kokkos( d_kineticsData.Ea, h_kineticsData.Ea, nreactions, "KineticsType::Ea");

  for (int i = 0; i < nreactions; ++i)
  {
    h_kineticsData.Arr[i]  = Arr[i];
    h_kineticsData.nArr[i] = nArr[i];
    h_kineticsData.Ea[i]   = Ea[i];
  }

  Kokkos::deep_copy( d_kineticsData.Arr, h_kineticsData.Arr );
  Kokkos::deep_copy( d_kineticsData.nArr, h_kineticsData.nArr );
  Kokkos::deep_copy( d_kineticsData.Ea, h_kineticsData.Ea );

  if (useSparseKinetics)
  {

    memory->create_kokkos( d_kineticsData.nu , h_kineticsData.nu , nreactions, sparseKinetics_maxSpecies, "KineticsType::nu");
    memory->create_kokkos( d_kineticsData.nuk, h_kineticsData.nuk, nreactions, sparseKinetics_maxSpecies, "KineticsType::nuk");

    for (int i = 0; i < nreactions; ++i)
      for (int k = 0; k < sparseKinetics_maxSpecies; ++k)
      {
        h_kineticsData.nu (i,k) = sparseKinetics_nu [i][k];
        h_kineticsData.nuk(i,k) = sparseKinetics_nuk[i][k];
      }

    Kokkos::deep_copy( d_kineticsData.nu, h_kineticsData.nu );
    Kokkos::deep_copy( d_kineticsData.nuk, h_kineticsData.nuk );

    if (SparseKinetics_enableIntegralReactions)
    {
      memory->create_kokkos( d_kineticsData.inu, h_kineticsData.inu, nreactions, sparseKinetics_maxSpecies, "KineticsType::inu");
      memory->create_kokkos( d_kineticsData.isIntegral, h_kineticsData.isIntegral, nreactions, "KineticsType::isIntegral");

      for (int i = 0; i < nreactions; ++i)
      {
        h_kineticsData.isIntegral(i) = sparseKinetics_isIntegralReaction[i];

        for (int k = 0; k < sparseKinetics_maxSpecies; ++k)
          h_kineticsData.inu(i,k) = sparseKinetics_inu[i][k];
      }

      Kokkos::deep_copy( d_kineticsData.inu, h_kineticsData.inu );
      Kokkos::deep_copy( d_kineticsData.isIntegral, h_kineticsData.isIntegral );
    }
  }

  //else
  //{

    // Dense option
    memory->create_kokkos( d_kineticsData.stoich, h_kineticsData.stoich, nreactions, nspecies, "KineticsType::stoich");
    memory->create_kokkos( d_kineticsData.stoichReactants, h_kineticsData.stoichReactants, nreactions, nspecies, "KineticsType::stoichReactants");
    memory->create_kokkos( d_kineticsData.stoichProducts, h_kineticsData.stoichProducts, nreactions, nspecies, "KineticsType::stoichProducts");

    for (int i = 0; i < nreactions; ++i)
      for (int k = 0; k < nspecies; ++k)
      {
        h_kineticsData.stoich(i,k) = stoich[i][k];
        h_kineticsData.stoichReactants(i,k) = stoichReactants[i][k];
        h_kineticsData.stoichProducts(i,k) = stoichProducts[i][k];
      }

    Kokkos::deep_copy( d_kineticsData.stoich, h_kineticsData.stoich );
    Kokkos::deep_copy( d_kineticsData.stoichReactants, h_kineticsData.stoichReactants );
    Kokkos::deep_copy( d_kineticsData.stoichProducts, h_kineticsData.stoichProducts );

  //}

  update_kinetics_data = false;
}

/* ---------------------------------------------------------------------- */

template <typename DeviceType>
void FixRxKokkos<DeviceType>::setup_pre_force(int vflag)
{
  printf("Inside FixRxKokkos<DeviceType>::setup_pre_force restartFlag= %d\n", my_restartFlag);

  if (my_restartFlag)
    my_restartFlag = 0;
  else
    this->solve_reactions( vflag, false );
}

/* ---------------------------------------------------------------------- */

template <typename DeviceType>
void FixRxKokkos<DeviceType>::pre_force(int vflag)
{
  printf("Inside FixRxKokkos<DeviceType>::pre_force localTempFlag= %d\n", localTempFlag);

  this->solve_reactions( vflag, true );
}
/* ---------------------------------------------------------------------- */

template <typename DeviceType>
void FixRxKokkos<DeviceType>::solve_reactions(const int vflag, const bool isPreForce)
{
  printf("Inside FixRxKokkos<DeviceType>::solve_reactions localTempFlag= %d isPreForce= %s\n", localTempFlag, isPreForce ? "True" : "false");

  if (update_kinetics_data)
    create_kinetics_data();

  TimerType timer_start = getTimeStamp();

  const int nlocal = atom->nlocal;
  const int nghost = atom->nghost;
  const int newton_pair = force->newton_pair;

  // Set the forward rates to zero if acting as setup_pre_force.
  const bool setRatesToZero = (isPreForce == false);

  if (localTempFlag)
  {
    const int count = nlocal + (newton_pair ? nghost : 0);

    memory->create_kokkos (k_dpdThetaLocal, dpdThetaLocal, count, "FixRxKokkos::dpdThetaLocal");
    d_dpdThetaLocal = k_dpdThetaLocal.d_view;
    h_dpdThetaLocal = k_dpdThetaLocal.h_view;

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
  typename ArrayTypes<DeviceType>::t_efloat_1d d_dpdTheta = atomKK->k_dpdTheta.view<DeviceType>(); 
  typename ArrayTypes<DeviceType>::t_float_2d  d_dvector  = atomKK->k_dvector.view<DeviceType>();
  typename ArrayTypes<DeviceType>::t_int_1d    d_mask     = atomKK->k_mask.view<DeviceType>();

  // Get up-to-date data.
  atomKK->sync( execution_space, MASK_MASK | DVECTOR_MASK | DPDTHETA_MASK );

  // Set some constants outside of the parallel_for
  const double boltz = force->boltz;
  const double t_stop = update->dt; // DPD time-step and integration length.

  // Average DPD volume. Used in the RHS function.
  this->VDPD = domain->xprd * domain->yprd * domain->zprd / atom->natoms;

  /*if (odeIntegrationFlag == ODE_LAMMPS_RKF45 && diagnosticFrequency == 1)
  {
    memory->create( diagnosticCounterPerODE[StepSum], nlocal, "FixRX::diagnosticCounterPerODE");
    memory->create( diagnosticCounterPerODE[FuncSum], nlocal, "FixRX::diagnosticCounterPerODE");
  }*/

  Kokkos::parallel_reduce( nlocal, LAMMPS_LAMBDA(int i, CounterType &counter)
    {
      if (d_mask(i) & groupbit)
      {
        double *y = new double[8*nspecies];
        double *rwork = y + nspecies;

        UserRHSData userData;
        userData.kFor = new double[nreactions];
        userData.rxnRateLaw = new double[nreactions];

        CounterType counter_i;

        const double theta = (localTempFlag) ? d_dpdThetaLocal(i) : d_dpdTheta(i);

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
          const double tmp = d_dvector(ispecies, i);
          d_dvector(ispecies+nspecies, i) = tmp;
          y[ispecies] = tmp;
        }

        // Solver the ODE system.
        if (odeIntegrationFlag == ODE_LAMMPS_RK4)
        {
          rk4(t_stop, y, rwork, &userData);
        }
        else if (odeIntegrationFlag == ODE_LAMMPS_RKF45)
        {
          rkf45(nspecies, t_stop, y, rwork, &userData, counter_i);

          //if (diagnosticFrequency == 1 && diagnosticCounterPerODE[StepSum] != NULL)
          //if (diagnosticCounterPerODE[StepSum] != NULL)
          //{
          //  diagnosticCounterPerODE[StepSum][i] = counter_i.nSteps;
          //  diagnosticCounterPerODE[FuncSum][i] = counter_i.nFuncs;
          //}
        }

        // Store the solution back in dvector.
        for (int ispecies = 0; ispecies < nspecies; ispecies++)
        {
          if (y[ispecies] < -MY_EPSILON)
            error->one(FLERR,"Computed concentration in RK solver is < -10*DBL_EPSILON");
          else if (y[ispecies] < MY_EPSILON)
            y[ispecies] = 0.0;

          d_dvector(ispecies,i) = y[ispecies];
        }

        delete [] y;
        delete [] userData.kFor;
        delete [] userData.rxnRateLaw;

        // Update the iteration statistics counter. Is this unique for each iteration?
        counter += counter_i;

      } // if
    } // parallel_for lambda-body

    , TotalCounters // reduction value for all iterations.
  );

  TimerType timer_ODE = getTimeStamp();

  // Signal that dvector has been modified on this execution space.
  atomKK->modified( execution_space, DVECTOR_MASK );

  // Communicate the updated species data to all nodes
  atomKK->sync ( Host, DVECTOR_MASK );

  comm->forward_comm_fix(this);

  atomKK->modified ( Host, DVECTOR_MASK );

  if (localTempFlag)
  {
    //delete [] dpdThetaLocal;
    memory->destroy_kokkos(k_dpdThetaLocal, dpdThetaLocal);
  }

  TimerType timer_stop = getTimeStamp();

  double time_ODE = getElapsedTime(timer_localTemperature, timer_ODE);

  printf("me= %d kokkos total= %g temp= %g ode= %g comm= %g nlocal= %d nfc= %d %d\n", comm->me,
                         getElapsedTime(timer_start, timer_stop),
                         getElapsedTime(timer_start, timer_localTemperature),
                         getElapsedTime(timer_localTemperature, timer_ODE),
                         getElapsedTime(timer_ODE, timer_stop), nlocal, TotalCounters.nFuncs, TotalCounters.nSteps);

  // Warn the user if a failure was detected in the ODE solver.
  if (TotalCounters.nFails > 0){
    char sbuf[128];
    sprintf(sbuf,"in FixRX::pre_force, ODE solver failed for %d atoms.", TotalCounters.nFails);
    error->warning(FLERR, sbuf);
  }

/*
  // Compute and report ODE diagnostics, if requested.
  if (odeIntegrationFlag == ODE_LAMMPS_RKF45 && diagnosticFrequency != 0){
    // Update the counters.
    diagnosticCounter[StepSum] += nSteps;
    diagnosticCounter[FuncSum] += nFuncs;
    diagnosticCounter[TimeSum] += time_ODE;
    diagnosticCounter[AtomSum] += nlocal;
    diagnosticCounter[numDiagnosticCounters-1] ++;

    if ( (diagnosticFrequency > 0 &&
               ((update->ntimestep - update->firststep) % diagnosticFrequency) == 0) ||
         (diagnosticFrequency < 0 && update->ntimestep == update->laststep) )
      this->odeDiagnostics();

    for (int i = 0; i < numDiagnosticCounters; ++i)
      if (diagnosticCounterPerODE[i])
        memory->destroy( diagnosticCounterPerODE[i] );
  } */
}

/* ---------------------------------------------------------------------- */

template <typename DeviceType>
  template <int WT_FLAG, int LOCAL_TEMP_FLAG, bool NEWTON_PAIR, int NEIGHFLAG>
void FixRxKokkos<DeviceType>::computeLocalTemperature()
{
  typename ArrayTypes<DeviceType>::t_x_array_randomread d_x        = atomKK->k_x.view<DeviceType>();
  typename ArrayTypes<DeviceType>::t_int_1d_randomread  d_type     = atomKK->k_type.view<DeviceType>();
  typename ArrayTypes<DeviceType>::t_efloat_1d          d_dpdTheta = atomKK->k_dpdTheta.view<DeviceType>(); 

  atomKK->sync(execution_space, X_MASK | TYPE_MASK | DPDTHETA_MASK );

  const int nlocal = atom->nlocal;
  const int nghost = atom->nghost;

  printf("Inside FixRxKokkos::computeLocalTemperature: %d %d %d %d %d %d %d\n", WT_FLAG, LOCAL_TEMP_FLAG, NEWTON_PAIR, (int)lmp->kokkos->neighflag, NEIGHFLAG, nlocal, nghost);

  // Pull from pairDPDE. The pairDPDEKK objects are protected so recreate here for now.
  //pairDPDEKK->k_cutsq.template sync<DeviceType>();
  //typename ArrayTypes<DeviceType>::t_ffloat_2d d_cutsq = pairDPDEKK->k_cutsq.template view<DeviceType();

  //!< Copies pulled from pairDPDE for local use since pairDPDEKK's objects are protected.
  typename ArrayTypes<DeviceType>::tdual_ffloat_2d k_cutsq;
  typename ArrayTypes<DeviceType>::t_ffloat_2d     d_cutsq;
  double **h_cutsq;

  {
    const int ntypes = atom->ntypes;

    memory->create_kokkos (k_cutsq, h_cutsq, ntypes+1, ntypes+1, "pair:cutsq");
    d_cutsq = k_cutsq.template view<DeviceType>();

    for (int i = 1; i <= ntypes; ++i)
      for (int j = i; j <= ntypes; ++j)
      {
        k_cutsq.h_view(i,j) = pairDPDE->cutsq[i][j];
        k_cutsq.h_view(j,i) = k_cutsq.h_view(i,j);
      }

    k_cutsq.template modify<LMPHostType>();
    k_cutsq.template sync<DeviceType>();
  }

  // Initialize the local temperature weight array
  int sumWeightsCt = nlocal + (NEWTON_PAIR ? nghost : 0);

  memory->create_kokkos (k_sumWeights, sumWeights, sumWeightsCt, "FixRxKokkos::sumWeights");
  d_sumWeights = k_sumWeights.d_view;
  h_sumWeights = k_sumWeights.h_view;

  // Initialize the accumulator to zero ...
  Kokkos::parallel_for (sumWeightsCt,
        LAMMPS_LAMBDA(const int i)
        {
           d_sumWeights(i) = 0.0;
        }
     );

  // Local list views. (This isn't working!)
  NeighListKokkos<DeviceType>* k_list = static_cast<NeighListKokkos<DeviceType>*>(list);
  if (not(list->kokkos))
     error->one(FLERR,"list is not a Kokkos list\n");

  typename ArrayTypes<DeviceType>::t_neighbors_2d d_neighbors = k_list->d_neighbors;
  typename ArrayTypes<DeviceType>::t_int_1d       d_ilist     = k_list->d_ilist;
  typename ArrayTypes<DeviceType>::t_int_1d       d_numneigh  = k_list->d_numneigh;

  const int inum = list->inum;

  // loop over neighbors of my atoms
  Kokkos::parallel_for ( inum,
        LAMMPS_LAMBDA(const int ii)
        {
          // Create an atomic view of sumWeights and dpdThetaLocal. Only needed
          // for Half/thread scenarios.
          typedef Kokkos::View< E_FLOAT*, typename DAT::t_efloat_1d::array_layout, DeviceType, Kokkos::MemoryTraits< AtomicF< NEIGHFLAG >::value> > AtomicViewType;

          AtomicViewType a_dpdThetaLocal = d_dpdThetaLocal;
          AtomicViewType a_sumWeights    = d_sumWeights;

          // Local scalar accumulators.
          double i_dpdThetaLocal = 0.0;
          double i_sumWeights    = 0.0;

          const int i = d_ilist(ii);
 
          const double xtmp = d_x(i,0);
          const double ytmp = d_x(i,1);
          const double ztmp = d_x(i,2);
          const int itype = d_type(i);

          const int jnum = d_numneigh(i);

          for (int jj = 0; jj < jnum; jj++)
          {
            const int j = (d_neighbors(i,jj) & NEIGHMASK);
            const int jtype = d_type(j);

            const double delx = xtmp - d_x(j,0);
            const double dely = ytmp - d_x(j,1);
            const double delz = ztmp - d_x(j,2);
            const double rsq = delx*delx + dely*dely + delz*delz;

            const double cutsq_ij = d_cutsq(itype,jtype);

            if (rsq < cutsq_ij)
            {
              const double rcut = sqrt( cutsq_ij );
              double rij = sqrt(rsq);
              double ratio = rij/rcut;

              double wij = 0.0;

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

  // Signal that dpdThetaLocal and sumWeights have been modified.
  k_dpdThetaLocal.template modify<DeviceType>();
  k_sumWeights.   template modify<DeviceType>();

  // Communicate the sum dpdTheta and the weights on the host.
  if (NEWTON_PAIR) comm->reverse_comm_fix(this);

  // Update the device view in case they got changed.
  k_dpdThetaLocal.template sync<DeviceType>();
  k_sumWeights.   template sync<DeviceType>();

  // self-interaction for local temperature
  Kokkos::parallel_for ( nlocal,
        LAMMPS_LAMBDA(const int i)
        {
          double wij = 0.0;

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

  // Clean up the local kokkos data.
  memory->destroy_kokkos(k_cutsq, h_cutsq);
  memory->destroy_kokkos(k_sumWeights, sumWeights);
}

/* ---------------------------------------------------------------------- */

template <typename DeviceType>
int FixRxKokkos<DeviceType>::pack_forward_comm(int n, int *list, double *buf, int pbc_flag, int *pbc)
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

template <typename DeviceType>
void FixRxKokkos<DeviceType>::unpack_forward_comm(int n, int first, double *buf)
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

template <typename DeviceType>
int FixRxKokkos<DeviceType>::pack_reverse_comm(int n, int first, double *buf)
{
  //printf("inside FixRxKokkos::pack_reverse_comm %d %d %d\n", comm->me, first, n);
  // Sync the host view.
  k_dpdThetaLocal.template sync<LMPHostType>();
  k_sumWeights.   template sync<LMPHostType>();

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

template <typename DeviceType>
void FixRxKokkos<DeviceType>::unpack_reverse_comm(int n, int *list, double *buf)
{
  // printf("inside FixRxKokkos::unpack_reverse_comm %d\n", comm->me);
  int m = 0;
  for (int i = 0; i < n; i++) {
    const int j = list[i];

    h_dpdThetaLocal(j) += buf[m++];
    h_sumWeights(j) += buf[m++];
  }

  // Signal that the host view has been modified.
  k_dpdThetaLocal.template modify<LMPHostType>();
  k_sumWeights.   template modify<LMPHostType>();

  // printf("done with FixRxKokkos::unpack_reverse_comm %d\n", comm->me);
}

namespace LAMMPS_NS {
template class FixRxKokkos<LMPDeviceType>;
#ifdef KOKKOS_HAVE_CUDA
template class FixRxKokkos<LMPHostType>;
#endif
}

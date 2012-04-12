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
   Author:  Todd Plantenga (SNL)
   Sources: "Numerical Optimization", Nocedal and Wright, 2nd Ed, p170
            "Parallel Unconstrained Min", Plantenga, SAND98-8201
------------------------------------------------------------------------- */

#include "math.h"
#include "string.h"
#include "atom.h"
#include "fix_minimize.h"
#include "min_hftn.h"
#include "modify.h"
#include "output.h"
#include "pair.h"
#include "update.h"
#include "timer.h"

using namespace LAMMPS_NS;

/* ----------------------------------------------------------------------
 * This class performs Hessian-free truncated Newton minimization on an
 * unconstrained molecular potential.  The algorithm avoids computing the
 * Hessian matrix, but obtains a near-quadratic rate of convergence.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   File local data
------------------------------------------------------------------------- */

//---- CONSTANTS MAP TO stopstrings DECLARED IN Min.run (min.cpp).

static const int  STOP_MAX_ITER = 0;          //-- MAX ITERATIONS EXCEEDED
static const int  STOP_MAX_FORCE_EVALS = 1;   //-- MAX FORCE EVALUATIONS EXCEEDED
static const int  STOP_ENERGY_TOL = 2;        //-- STEP DID NOT CHANGE ENERGY
static const int  STOP_FORCE_TOL = 3;         //-- CONVERGED TO DESIRED FORCE TOL
static const int  STOP_TR_TOO_SMALL = 8;      //-- TRUST REGION TOO SMALL
static const int  STOP_ERROR = 9;             //-- INTERNAL ERROR

static const int  NO_CGSTEP_BECAUSE_F_TOL_SATISFIED = 0;
static const int  CGSTEP_NEWTON                     = 1;
static const int  CGSTEP_TO_TR                      = 2;
static const int  CGSTEP_TO_DMAX                    = 3;
static const int  CGSTEP_NEGATIVE_CURVATURE         = 4;
static const int  CGSTEP_MAX_INNER_ITERS            = 5;
static const int  CGSTEP_UNDETERMINED               = 6;

//---- WHEN TESTING ENERGY_TOL, THE ENERGY MAGNITUDE MUST BE AT LEAST THIS BIG.

static const double  MIN_ETOL_MAG = 1.0e-8;

//---- MACHINE PRECISION IS SOMETIMES DEFINED BY THE C RUNTIME.

#ifdef DBL_EPSILON
  #define MACHINE_EPS  DBL_EPSILON
#else
  #define MACHINE_EPS  2.220446049250313e-16
#endif

/* ----------------------------------------------------------------------
   Constructor
------------------------------------------------------------------------- */

MinHFTN::MinHFTN(LAMMPS *lmp) : Min(lmp)
{
  searchflag = 1;

  for (int  i = 1; i < NUM_HFTN_ATOM_BASED_VECTORS; i++)
    _daExtraGlobal[i] = NULL;
  for (int  i = 0; i < NUM_HFTN_ATOM_BASED_VECTORS; i++)
    _daExtraAtom[i] = NULL;
  
  _fpPrint = NULL;
  
  return;
}

/* ----------------------------------------------------------------------
   Destructor
------------------------------------------------------------------------- */

MinHFTN::~MinHFTN (void)
{
  for (int  i = 1; i < NUM_HFTN_ATOM_BASED_VECTORS; i++)
    if (_daExtraGlobal[i] != NULL)
      delete [] _daExtraGlobal[i];
  for (int  i = 0; i < NUM_HFTN_ATOM_BASED_VECTORS; i++)
    if (_daExtraAtom[i] != NULL)
      delete [] _daExtraAtom[i];
  
  return;
}

/* ----------------------------------------------------------------------
   Public method init
------------------------------------------------------------------------- */

void MinHFTN::init()
{
  Min::init();

  for (int  i = 1; i < NUM_HFTN_ATOM_BASED_VECTORS; i++) {
    if (_daExtraGlobal[i] != NULL)
      delete [] _daExtraGlobal[i];
    _daExtraGlobal[i] = NULL;
  }
  for (int  i = 0; i < NUM_HFTN_ATOM_BASED_VECTORS; i++) {
    if (_daExtraAtom[i] != NULL)
      delete [] _daExtraAtom[i];
    _daExtraAtom[i] = NULL;
  }
  
  return;
}

/* ----------------------------------------------------------------------
   Public method setup_style
------------------------------------------------------------------------- */

void MinHFTN::setup_style()
{
  //---- ALLOCATE MEMORY FOR ATOMIC DEGREES OF FREEDOM.
  for (int  i = 0; i < NUM_HFTN_ATOM_BASED_VECTORS; i++)
    fix_minimize->add_vector(3);
  
  //---- ALLOCATE MEMORY FOR EXTRA GLOBAL DEGREES OF FREEDOM.
  //---- THE FIX MODULE TAKES CARE OF THE FIRST VECTOR, X0 (XK).
  if (nextra_global) {
    for (int  i = 1; i < NUM_HFTN_ATOM_BASED_VECTORS; i++)
      _daExtraGlobal[i] = new double[nextra_global];
  }

  //---- ALLOCATE MEMORY FOR EXTRA PER-ATOM DEGREES OF FREEDOM.
  if (nextra_atom) {
    for (int  i = 0; i < NUM_HFTN_ATOM_BASED_VECTORS; i++)
      _daExtraAtom[i] = new double*[nextra_atom];
    
    for (int m = 0; m < nextra_atom; m++) {
      for (int  i = 0; i < NUM_HFTN_ATOM_BASED_VECTORS; i++)
	fix_minimize->add_vector (extra_peratom[m]);
    }
  }
  
  return;
}

/* ----------------------------------------------------------------------
   Public method reset_vectors
   After an energy/force calculation, atoms may migrate from one processor
   to another.  Any local vector correlated with atom positions or forces
   must also be migrated.  This is accomplished by a subclass of Fix.
   This method updates local pointers to the latest Fix copies.
------------------------------------------------------------------------- */

void MinHFTN::reset_vectors()
{
  nvec = 3 * atom->nlocal;
  
  //---- ATOMIC DEGREES OF FREEDOM.
  if (nvec > 0) {
    xvec = atom->x[0];
    fvec = atom->f[0];
  }
  for (int  i = 0; i < NUM_HFTN_ATOM_BASED_VECTORS; i++)
    _daAVectors[i] = fix_minimize->request_vector (i);
  
  //---- EXTRA PER-ATOM DEGREES OF FREEDOM.
  if (nextra_atom) {
    int  n = NUM_HFTN_ATOM_BASED_VECTORS;
    for (int m = 0; m < nextra_atom; m++) {
      extra_nlen[m] = extra_peratom[m] * atom->nlocal;
      requestor[m]->min_xf_pointers(m,&xextra_atom[m],&fextra_atom[m]);
      for (int  i = 0; i < NUM_HFTN_ATOM_BASED_VECTORS; i++)
	_daExtraAtom[i][m] = fix_minimize->request_vector (n++);
    }
  }
  
  return;
}

/* ----------------------------------------------------------------------
   Public method iterate
   Upon entry, Min::setup() and Min::run have executed, and energy has
   already been evaluated at the initial point.  Return an integer code
   that maps to a stop condition in min.cpp.
------------------------------------------------------------------------- */

int MinHFTN::iterate(int)
{
  //---- TURN THIS ON TO GENERATE AN OPTIMIZATION PROGRESS FILE.
  bool  bPrintProgress = false;
  
  if (bPrintProgress)
    open_hftn_print_file_();
  
  double  dFinalEnergy = 0.0;
  double  dFinalFnorm2 = 0.0;
  modify->min_clearstore();
  int  nStopCode = execute_hftn_ (bPrintProgress,
				  einitial,
				  fnorm2_init,
				  dFinalEnergy,
				  dFinalFnorm2);
  modify->min_clearstore();
  if (bPrintProgress)
    close_hftn_print_file_();
  
  return( nStopCode );
}

/* ----------------------------------------------------------------------
   Private method execute_hftn_
   @param[in] bPrintProgress - if true then print progress to a file
   @param[in] dInitialEnergy - energy at input x
   @param[in] dInitialForce2 - |F|_2 at input x
   @param[out] dFinalEnergy  - energy at output x
   @param[out] dFinalForce2  - |F|_2 at output x

   Return stop code described in the enumeration at the top of this file,
   and the following:
     atom->x  - positions at output x
     atom->f  - forces evaluated at output x
------------------------------------------------------------------------- */
int MinHFTN::execute_hftn_(const bool      bPrintProgress,
			   const double    dInitialEnergy,
			   const double    dInitialForce2,
			   double &  dFinalEnergy,
			   double &  dFinalForce2)
{
  //---- DEFINE OUTPUTS PRINTED BY "Finish".
  eprevious = dInitialEnergy;
  alpha_final = 0.0;
  dFinalEnergy = dInitialEnergy;
  dFinalForce2 = dInitialForce2;
  
  if (dInitialForce2 < update->ftol)
    return( STOP_FORCE_TOL );
  
  //---- SAVE ATOM POSITIONS BEFORE AN ITERATION.
  fix_minimize->store_box();
  for (int  i = 0; i < nvec; i++)
    _daAVectors[VEC_XK][i] = xvec[i];
  if (nextra_atom) {
    for (int  m = 0; m < nextra_atom; m++) {
      double *  xatom  = xextra_atom[m];
      double *  xkAtom = _daExtraAtom[VEC_XK][m];
      int  n = extra_nlen[m];
      for (int  i = 0; i < n; i++)
	xkAtom[i] = xatom[i];
    }
  }
  if (nextra_global)
    modify->min_store();
  
  double  dXInf = calc_xinf_using_mpi_();
  
  //---- FIND THE NUMBER OF UNKNOWNS.
  int  nLocalNumUnknowns = nvec + nextra_atom;
  MPI_Allreduce (&nLocalNumUnknowns, &_nNumUnknowns,
		 1, MPI_INT, MPI_SUM, world);
  
  //---- INITIALIZE THE TRUST RADIUS BASED ON THE GRADIENT.
  double  dTrustRadius = 1.5 * dInitialForce2;
  
  //---- TRUST RADIUS MUST KEEP STEPS FROM LETTING ATOMS MOVE SO FAR THEY
  //---- VIOLATE PHYSICS OR JUMP BEYOND A PARALLEL PROCESSING DOMAIN.
  //---- LINE SEARCH METHODS DO THIS BY RESTRICTING THE LARGEST CHANGE
  //---- OF ANY ATOM'S COMPONENT TO dmax.  AN EXACT CHECK IS MADE LATER,
  //---- BUT THIS GUIDES DETERMINATION OF A MAX TRUST RADIUS.
  double  dMaxTrustRadius = dmax * sqrt((double) _nNumUnknowns);
  
  dTrustRadius = MIN (dTrustRadius, dMaxTrustRadius);
  double  dLastNewtonStep2 = dMaxTrustRadius;
  
  if (bPrintProgress)
    hftn_print_line_ (false, -1, neval, dInitialEnergy, dInitialForce2,
		      -1, dTrustRadius, 0.0, 0.0, 0.0);
  
  bool    bHaveEvaluatedAtX = true;
  double  dCurrentEnergy    = dInitialEnergy;
  double  dCurrentForce2    = dInitialForce2;
  for (niter = 0; niter < update->nsteps; niter++) {
    (update->ntimestep)++;
    
    //---- CALL THE INNER LOOP TO GET THE NEXT TRUST REGION STEP.
    
    double  dCgForce2StopTol = MIN ((dCurrentForce2 / 2.0), 0.1 / (niter+1));
    dCgForce2StopTol = MAX (dCgForce2StopTol, update->ftol);
    
    double  dNewEnergy;
    double  dNewForce2;
    int     nStepType;
    double  dStepLength2;
    double  dStepLengthInf;
    if (compute_inner_cg_step_ (dTrustRadius,
				dCgForce2StopTol,
				update->max_eval,
				bHaveEvaluatedAtX,
				dCurrentEnergy, dCurrentForce2,
				dNewEnergy, dNewForce2,
				nStepType,
				dStepLength2, dStepLengthInf) == false) {
      //---- THERE WAS AN ERROR.  RESTORE TO LAST ACCEPTED STEP.
      if (nextra_global)
	modify->min_step (0.0, _daExtraGlobal[VEC_CG_P]);
      for (int i = 0; i < nvec; i++)
	xvec[i] = _daAVectors[VEC_XK][i];
      if (nextra_atom) {
	for (int  m = 0; m < nextra_atom; m++) {
	  double *  xatom  = xextra_atom[m];
	  double *  xkAtom = _daExtraAtom[VEC_XK][m];
	  int  n = extra_nlen[m];
	  for (int  i = 0; i < n; i++)
	    xatom[i] = xkAtom[i];
	  requestor[m]->min_x_set(m);
	}
      }
      dFinalEnergy = energy_force (0);
      neval++;
      dFinalForce2 = sqrt (fnorm_sqr());
      return( STOP_ERROR );
    }
    
    //---- STOP IF THE CURRENT POSITION WAS FOUND TO BE ALREADY GOOD ENOUGH.
    //---- IN THIS CASE THE ENERGY AND FORCES ARE ALREADY COMPUTED.
    if (nStepType == NO_CGSTEP_BECAUSE_F_TOL_SATISFIED) {
      if (bPrintProgress)
	hftn_print_line_ (true, niter+1, neval, dNewEnergy, dNewForce2,
			  nStepType, dTrustRadius, dStepLength2,
			  0.0, 0.0);
      dFinalEnergy = dNewEnergy;
      dFinalForce2 = dNewForce2;
      return( STOP_FORCE_TOL );
    }
    
    //---- COMPUTE THE DIRECTIONAL DERIVATIVE H(x_k) p.
    bool  bUseForwardDiffs = (dCurrentForce2 > 1000.0 * sqrt (MACHINE_EPS));
    evaluate_dir_der_ (bUseForwardDiffs,
		       VEC_CG_P,
		       VEC_CG_HD,
		       true,
		       dCurrentEnergy);
    
    //---- COMPUTE p^T grad(x_k) AND SAVE IT FOR PRED.
    double  dGradDotP = calc_grad_dot_v_using_mpi_ (VEC_CG_P);
    
    //---- MOVE TO THE NEW POINT AND EVALUATE ENERGY AND FORCES.
    //---- THIS IS THE PLACE WHERE energy_force IS ALLOWED TO RESET.
    for (int i = 0; i < nvec; i++)
      xvec[i] = _daAVectors[VEC_XK][i] + _daAVectors[VEC_CG_P][i];
    if (nextra_atom) {
      for (int  m = 0; m < nextra_atom; m++) {
	double *  xatom  = xextra_atom[m];
	double *  xkAtom = _daExtraAtom[VEC_XK][m];
	double *  pAtom  = _daExtraAtom[VEC_CG_P][m];
	int  n = extra_nlen[m];
	for (int  i = 0; i < n; i++)
	  xatom[i] = xkAtom[i] + pAtom[i];
	requestor[m]->min_x_set(m);
      }
    }
    if (nextra_global)
      modify->min_step (1.0, _daExtraGlobal[VEC_CG_P]);
    dNewEnergy = energy_force (1);
    neval++;
    
    dNewForce2 = sqrt (fnorm_sqr());
    
    double  dAred = dCurrentEnergy - dNewEnergy;
    
    //---- STOP IF THE FORCE TOLERANCE IS MET.
    if (dNewForce2 < update->ftol) {
      if (bPrintProgress)
	hftn_print_line_ (true, niter+1, neval, dNewEnergy, dNewForce2,
			  nStepType, dTrustRadius, dStepLength2,
			  dAred, -1.0);
      //---- (IMPLICITLY ACCEPT THE LAST STEP TO THE NEW POINT.)
      dFinalEnergy = dNewEnergy;
      dFinalForce2 = dNewForce2;
      return( STOP_FORCE_TOL );
    }
    
    //---- STOP IF THE ACTUAL ENERGY REDUCTION IS TINY.
    if (nStepType != CGSTEP_TO_DMAX) {
      double  dMag = 0.5 * (fabs (dCurrentEnergy) + fabs (dNewEnergy));
      dMag = MAX (dMag, MIN_ETOL_MAG);
      if (   (fabs (dAred) < (update->etol * dMag))
	     || (dStepLengthInf == 0.0) ) {
	if (bPrintProgress)
	  hftn_print_line_ (true, niter+1, neval,
			    dNewEnergy, dNewForce2,
			    nStepType, dTrustRadius, dStepLength2,
			    dAred, -1.0);
	//---- (IMPLICITLY ACCEPT THE LAST STEP TO THE NEW POINT.)
	dFinalEnergy = dNewEnergy;
	dFinalForce2 = dNewForce2;
	return( STOP_ENERGY_TOL );
      }
    }
    
    //---- COMPUTE THE PREDICTED REDUCTION  - p^T grad - 0.5 p^T Hp
    double  dPHP = calc_dot_prod_using_mpi_ (VEC_CG_P, VEC_CG_HD);
    double  dPred = - dGradDotP - (0.5 * dPHP);
    
    //---- ACCEPT OR REJECT THE STEP PROPOSED BY THE INNER CG LOOP.
    //---- WHEN NEAR A SOLUTION, THE FORCE NORM IS PROBABLY MORE ACCURATE,
    //---- SO DON'T ACCEPT A STEP THAT REDUCES ENERGY SOME TINY AMOUNT
    //---- WHILE INCREASING THE FORCE NORM.
    bool  bStepAccepted = (dAred > 0.0)
      && (   (dNewForce2 < dCurrentForce2)
	     || (dCurrentForce2 > 1.0e-6));
    if (bStepAccepted) {
      //---- THE STEP IS ACCEPTED.
      if (bPrintProgress)
	hftn_print_line_ (true, niter+1, neval, dNewEnergy, dNewForce2,
			  nStepType, dTrustRadius, dStepLength2,
			  dAred, dPred);
      
      fix_minimize->store_box();
      modify->min_clearstore();
      for (int  i = 0; i < nvec; i++)
	_daAVectors[VEC_XK][i] = xvec[i];
      if (nextra_atom) {
	for (int  m = 0; m < nextra_atom; m++) {
	  double *  xatom  = xextra_atom[m];
	  double *  xkAtom = _daExtraAtom[VEC_XK][m];
	  int  n = extra_nlen[m];
	  for (int  i = 0; i < n; i++)
	    xkAtom[i] = xatom[i];
	}
      }
      if (nextra_global)
	modify->min_store();
      
      if (niter > 0)
	eprevious = dCurrentEnergy;
      dCurrentEnergy = dNewEnergy;
      dCurrentForce2 = dNewForce2;
      bHaveEvaluatedAtX = true;
      
      if (nStepType == CGSTEP_NEWTON)
	dLastNewtonStep2 = dStepLength2;
      
      //---- UPDATE THE TRUST REGION BASED ON AGREEMENT BETWEEN
      //---- THE ACTUAL REDUCTION AND THE PREDICTED MODEL REDUCTION.
      if ((dAred > 0.75 * dPred) && (dStepLength2 >= 0.99 * dTrustRadius))
	dTrustRadius = 2.0 * dTrustRadius;
      dTrustRadius = MIN (dTrustRadius, dMaxTrustRadius);
      
      //---- DMAX VIOLATIONS TRUNCATE THE CG STEP WITHOUT COMPARISONS;
      //---- BETTER TO ADJUST THE TRUST REGION SO DMAX STOPS HAPPENING.
      if (nStepType == CGSTEP_TO_DMAX) {
	if (dStepLength2 <= MACHINE_EPS)
	  dTrustRadius = 0.1 * dTrustRadius;
	else
	  dTrustRadius = MIN (dTrustRadius, 2.0 * dStepLength2);
      }
    }
    else {
      //---- THE STEP IS REJECTED.
      if (bPrintProgress)
	hftn_print_line_ (false, niter+1, neval,
			  dCurrentEnergy, dCurrentForce2,
			  nStepType, dTrustRadius, dStepLength2,
			  dAred, dPred);
      
      //---- RESTORE THE LAST X_K POSITION.
      if (nextra_global)
	modify->min_step (0.0, _daExtraGlobal[VEC_CG_P]);
      for (int  i = 0; i < nvec; i++)
	xvec[i] = _daAVectors[VEC_XK][i];
      if (nextra_atom) {
	for (int  m = 0; m < nextra_atom; m++) {
	  double *  xatom  = xextra_atom[m];
	  double *  xkAtom = _daExtraAtom[VEC_XK][m];
	  int  n = extra_nlen[m];
	  for (int  i = 0; i < n; i++)
	    xatom[i] = xkAtom[i];
	  requestor[m]->min_x_set(m);
	}
      }
      modify->min_clearstore();
      bHaveEvaluatedAtX = false;
      
      //---- UPDATE THE TRUST REGION.
      //---- EXPERIMENTS INDICATE NEGATIVE CURVATURE CAN TAKE A BAD
      //---- STEP A LONG WAY, SO BE MORE AGGRESSIVE IN THIS CASE.
      //---- ALSO, IF NEAR A SOLUTION AND DONE WITH NEWTON STEPS,
      //---- THEN REDUCE TO SOMETHING NEAR THE LAST GOOD NEWTON STEP.
      if ((nStepType == CGSTEP_NEGATIVE_CURVATURE) && (-dAred > dPred))
	dTrustRadius = 0.10 * MIN (dTrustRadius, dStepLength2);
      else if (   (nStepType == CGSTEP_TO_DMAX)
		  && (dStepLength2 <= MACHINE_EPS))
	dTrustRadius = 0.10 * dTrustRadius;
      else if (-dAred > dPred)
	dTrustRadius = 0.20 * MIN (dTrustRadius, dStepLength2);
      else
	dTrustRadius = 0.25 * MIN (dTrustRadius, dStepLength2);
      
      if (   (nStepType != CGSTEP_NEWTON)
	     && (dCurrentForce2 < sqrt (MACHINE_EPS)))
	dTrustRadius = MIN (dTrustRadius, 2.0 * dLastNewtonStep2);
      
      dLastNewtonStep2 = dMaxTrustRadius;
      
      //---- STOP IF THE TRUST RADIUS IS TOO SMALL TO CONTINUE.
      if (   (dTrustRadius <= 0.0)
	     || (dTrustRadius <= MACHINE_EPS * MAX (1.0, dXInf))) {
	dFinalEnergy = dCurrentEnergy;
	dFinalForce2 = dCurrentForce2;
	return( STOP_TR_TOO_SMALL );
      }
    }
    
    //---- OUTPUT FOR thermo, dump, restart FILES.
    if (output->next == update->ntimestep) {
      //---- IF THE LAST STEP WAS REJECTED, THEN REEVALUATE ENERGY AND
      //---- FORCES AT THE OLD POINT SO THE OUTPUT DOES NOT DISPLAY
      //---- THE INCREASED ENERGY OF THE REJECTED STEP.
      if (bStepAccepted == false) {
	dCurrentEnergy = energy_force (1);
	neval++;
      }
      timer->stamp();
      output->write (update->ntimestep);
      timer->stamp (TIME_OUTPUT);
    }
    
    //---- RETURN IF NUMBER OF EVALUATIONS EXCEEDED.
    if (neval >= update->max_eval) {
      dFinalEnergy = dCurrentEnergy;
      dFinalForce2 = dCurrentForce2;
      return( STOP_MAX_FORCE_EVALS );
    }
    
  }     //-- END for LOOP OVER niter
  
  dFinalEnergy = dCurrentEnergy;
  dFinalForce2 = dCurrentForce2;
  return( STOP_MAX_ITER );
}

/* ----------------------------------------------------------------------
   Private method compute_inner_cg_step_
   Execute CG using Hessian-vector products approximated by finite difference
     directional derivatives.

   On input these must be defined:
     atom->x   - positions at x
     atom->f   - ignored
     VEC_XK    - positions at x
   On output these are defined:
     atom->x   - unchanged
     atom->f   - forces evaluated at x, but only if nStepType == NO_CGSTEP
     VEC_XK    - unchanged
     VEC_CG_P  - step from VEC_XK to new positions
   During processing these are modified:
     VEC_CG_D  - conjugate gradient inner loop step
     VEC_CG_HD - Hessian-vector product
     VEC_CG_R  - residual of inner loop step
     VEC_DIF1  - temp storage
     VEC_DIF2  - temp storage

   @param[in] dTrustRadius     - trust region radius for this subiteration
   @param[in] dForceTol        - stop tolerance on |F|_2 for this subiteration
   @param[in] nMaxEvals        - total energy/force evaluations allowed
   @param[in] bHaveEvalAtXin   - true if forces are valid at input x
   @param[in] dEnergyAtXin     - energy at input x, if bHaveEvalAtXin is true
   @param[in] dForce2AtXin     - |F|_2 at input x, if bHaveEvalAtXin is true
   @param[out] dEnergyAtXout   - energy at output x, if NO_CGSTEP (see below)
   @param[out] dForce2AtXout   - |F|_2 at output x, if NO_CGSTEP (see below)
   @param[out] nStepType       - step type for hftn_print_line_()
   @param[out] dStepLength2    - |step|_2
   @param[out] dStepLengthInf  - |step|_inf

   Return false if there was a fatal error.
   If nStepType equals NO_CGSTEP_BECAUSE_F_TOL_SATISFIED, then the energy
   and forces are evaluated and returned in dEnergyAtXout, dForce2AtXout;
   else energy and forces are not evaluated.
------------------------------------------------------------------------- */

bool MinHFTN::compute_inner_cg_step_(const double    dTrustRadius,
				     const double    dForceTol,
				     const int       nMaxEvals,
				     const bool      bHaveEvalAtXin,
				     const double    dEnergyAtXin,
				     const double    dForce2AtXin,
				     double &  dEnergyAtXout,
				     double &  dForce2AtXout,
				     int    &  nStepType,
				     double &  dStepLength2,
				     double &  dStepLengthInf)
{
  //---- SET  p_0 = 0.
  if (nextra_global) {
    for (int  i = 0; i < nextra_global; i++)
      _daExtraGlobal[VEC_CG_P][i] = 0.0;
  }
  for (int  i = 0; i < nvec; i++)
    _daAVectors[VEC_CG_P][i] = 0.0;
  if (nextra_atom) {
    for (int  m = 0; m < nextra_atom; m++) {
      double *  pAtom = _daExtraAtom[VEC_CG_P][m];
      int  n = extra_nlen[m];
      for (int  i = 0; i < n; i++)
	pAtom[i] = 0.0;
    }
  }
  double  dPP = 0.0;
  
  //---- OBTAIN THE ENERGY AND FORCES AT THE INPUT POSITION.
  double  dEnergyAtX = dEnergyAtXin;
  double  dForce2AtX = dForce2AtXin;
  if (bHaveEvalAtXin == false) {
    dEnergyAtX = energy_force (0);
    neval++;
    dForce2AtX = sqrt (fnorm_sqr());
  }
  
  //---- RETURN IMMEDIATELY IF THE FORCE TOLERANCE IS ALREADY MET.
  //---- THE STEP TYPE INFORMS THE CALLER THAT ENERGY AND FORCES HAVE
  //---- BEEN EVALUATED.
  if (dForce2AtX <= dForceTol) {
    dEnergyAtXout = dEnergyAtX;
    dForce2AtXout = dForce2AtX;
    nStepType = NO_CGSTEP_BECAUSE_F_TOL_SATISFIED;
    dStepLength2 = 0.0;
    dStepLengthInf = 0.0;
    return( true );
  }
  
  //---- r_0 = -grad  (FIRST SEARCH DIRECTION IS STEEPEST DESCENT)
  //---- d_0 = r_0
  //---- REMEMBER THAT FORCES = -GRADIENT.
  if (nextra_global) {
    for (int  i = 0; i < nextra_global; i++) {
      _daExtraGlobal[VEC_CG_R][i] = fextra[i];
      _daExtraGlobal[VEC_CG_D][i] = fextra[i];
    }
  }
  for (int  i = 0; i < nvec; i++) {
    _daAVectors[VEC_CG_R][i] = fvec[i];
    _daAVectors[VEC_CG_D][i] = fvec[i];
  }
  if (nextra_atom) {
    for (int  m = 0; m < nextra_atom; m++) {
      double *  fatom = fextra_atom[m];
      double *  rAtom = _daExtraAtom[VEC_CG_R][m];
      double *  dAtom = _daExtraAtom[VEC_CG_D][m];
      int  n = extra_nlen[m];
      for (int  i = 0; i < n; i++) {
	rAtom[i] = fatom[i];
	dAtom[i] = fatom[i];
      }
    }
  }
  double  dRR = dForce2AtX * dForce2AtX;
  double  dR0norm2 = sqrt (dRR);
  
  //---- LIMIT THE NUMBER OF INNER CG ITERATIONS.
  //---- BASE IT ON THE NUMBER OF UNKNOWNS, OR MAXIMUM EVALUATIONS ASSUMING
  //---- FORWARD DIFFERENCES ARE USED.
  //---- NOTE THAT SETTING MAX=1 GIVES STEEPEST DESCENT.
  int  nLimit1 = _nNumUnknowns / 5;
  if (nLimit1 < 100)
    nLimit1 = MIN (_nNumUnknowns, 100);
  int  nLimit2 = (nMaxEvals - neval) / 2;
  int  nMaxInnerIters = MIN (nLimit1, nLimit2);
  
  //---- FURTHER LIMIT ITERATIONS IF NEAR MACHINE ROUNDOFF.
  //---- THE METHOD CAN WASTE A LOT EVALUATIONS WITH LITTLE PAYOFF PROSPECT.
  if (dForce2AtX < (sqrt (MACHINE_EPS) * MAX (1.0, fabs (dEnergyAtX))) )
    nMaxInnerIters = MIN (nMaxInnerIters, _nNumUnknowns / 20);
  
  bool  bUseForwardDiffs = (dForce2AtX > 1000.0 * sqrt (MACHINE_EPS));
  
  //---- MAIN CG LOOP.
  for (int  nInnerIter = 0; nInnerIter < nMaxInnerIters; nInnerIter++) {
    //---- COMPUTE HESSIAN-VECTOR PRODUCT:  H(x_k) d_i.
    double  dDummyEnergy;
    evaluate_dir_der_ (bUseForwardDiffs,
		       VEC_CG_D,
		       VEC_CG_HD,
		       false,
		       dDummyEnergy);
    
    //---- CALCULATE  d_i^T H d_i AND d_i^T d_i.
    double  dDHD;
    double  dDD;
    calc_dhd_dd_using_mpi_ (dDHD, dDD);
    
    //---- HANDLE NEGATIVE CURVATURE.
    if (dDHD <= (MACHINE_EPS * dDD)) {
      //---- PROJECT BOTH DIRECTIONS TO THE TRUST RADIUS AND DECIDE
      //---- WHICH MAKES A BETTER PREDICTED REDUCTION.
      //---- p_i^T H(x_k) d_i AND grad_i^T d_i.
      
      double  dPdotD  = calc_dot_prod_using_mpi_ (VEC_CG_P, VEC_CG_D);
      double  dPdotHD = calc_dot_prod_using_mpi_ (VEC_CG_P, VEC_CG_HD);
      
      //---- MOVE TO X_K AND COMPUTE ENERGY AND FORCES.
      if (nextra_global)
	modify->min_step (0.0, _daExtraGlobal[VEC_CG_P]);
      for (int  i = 0; i < nvec; i++)
	xvec[i] = _daAVectors[VEC_XK][i];
      if (nextra_atom) {
	for (int  m = 0; m < nextra_atom; m++) {
	  double *  xatom  = xextra_atom[m];
	  double *  xkAtom = _daExtraAtom[VEC_XK][m];
	  int  n = extra_nlen[m];
	  for (int  i = 0; i < n; i++)
	    xatom[i] = xkAtom[i];
	  requestor[m]->min_x_set(m);
	}
      }
      dEnergyAtX = energy_force (0);
      neval++;
      
      double  dGradDotD = calc_grad_dot_v_using_mpi_ (VEC_CG_D);
      
      double  tau = compute_to_tr_ (dPP, dPdotD, dDD, dTrustRadius,
				    true, dDHD, dPdotHD, dGradDotD);
      
      //---- MOVE THE POINT.
      if (nextra_global) {
	double *  pGlobal = _daExtraGlobal[VEC_CG_P];
	double *  dGlobal = _daExtraGlobal[VEC_CG_D];
	for (int  i = 0; i < nextra_global; i++) {
	  pGlobal[i] += tau * dGlobal[i];
	}
      }
      for (int  i = 0; i < nvec; i++)
	_daAVectors[VEC_CG_P][i] += tau * _daAVectors[VEC_CG_D][i];
      if (nextra_atom) {
	for (int  m = 0; m < nextra_atom; m++) {
	  double *  pAtom = _daExtraAtom[VEC_CG_P][m];
	  double *  dAtom = _daExtraAtom[VEC_CG_D][m];
	  int  n = extra_nlen[m];
	  for (int  i = 0; i < n; i++)
	    pAtom[i] += tau * dAtom[i];
	}
      }
      
      nStepType = CGSTEP_NEGATIVE_CURVATURE;
      calc_plengths_using_mpi_ (dStepLength2, dStepLengthInf);
      return( true );
    }
    
    //---- COMPUTE THE OPTIMAL STEP LENGTH BASED ON THE QUADRATIC CG MODEL.
    double  dAlpha = dRR / dDHD;
    
    //---- MIGHT WANT TO ENABLE THIS TO DEBUG INTERNAL CG STEPS.
    //fprintf (_fpPrint, "     alpha = %11.8f  neval=%4d\n", dAlpha, neval);
    
    //---- p_i+1 = p_i + alpha_i d_i
    //---- (SAVE THE CURRENT p_i IN CASE THE STEP HAS TO BE SHORTENED.)
    if (nextra_global) {
      double *  pGlobal  = _daExtraGlobal[VEC_CG_P];
      double *  dGlobal  = _daExtraGlobal[VEC_CG_D];
      double *  d1Global = _daExtraGlobal[VEC_DIF1];
      for (int  i = 0; i < nextra_global; i++) {
	d1Global[i] = pGlobal[i];
	pGlobal[i] += dAlpha * dGlobal[i];
      }
    }
    for (int  i = 0; i < nvec; i++) {
      _daAVectors[VEC_DIF1][i] = _daAVectors[VEC_CG_P][i];
      _daAVectors[VEC_CG_P][i] += dAlpha * _daAVectors[VEC_CG_D][i];
    }
    if (nextra_atom) {
      for (int  m = 0; m < nextra_atom; m++) {
	double *  pAtom  = _daExtraAtom[VEC_CG_P][m];
	double *  dAtom  = _daExtraAtom[VEC_CG_D][m];
	double *  d1Atom = _daExtraAtom[VEC_DIF1][m];
	int  n = extra_nlen[m];
	for (int  i = 0; i < n; i++) {
	  d1Atom[i] = pAtom[i];
	  pAtom[i] += dAlpha * dAtom[i];
	}
      }
    }
    
    //---- COMPUTE VECTOR PRODUCTS  p_i+1^T p_i+1 AND p_i^T d_i.
    double  dPnewDotPnew;
    double  dPoldDotD;
    calc_ppnew_pdold_using_mpi_ (dPnewDotPnew, dPoldDotD);
    
    nStepType = CGSTEP_UNDETERMINED;
    
    //---- IF STEP LENGTH IS TOO LARGE, THEN REDUCE IT AND RETURN.
    double  tau;
    if (step_exceeds_TR_ (dTrustRadius, dPP, dPoldDotD, dDD, tau)) {
      adjust_step_to_tau_ (tau);
      nStepType = CGSTEP_TO_TR;
    }
    if (step_exceeds_DMAX_()) {
      adjust_step_to_tau_ (0.0);
      nStepType = CGSTEP_TO_DMAX;
    }
    if ((nStepType == CGSTEP_TO_TR) || (nStepType == CGSTEP_TO_DMAX)) {
      calc_plengths_using_mpi_ (dStepLength2, dStepLengthInf);
      return( true );
    }
    
    dStepLength2 = sqrt (dPnewDotPnew);
    
    //---- r_i+1 = r_i - alpha * H d_i
    if (nextra_global) {
      double *  rGlobal  = _daExtraGlobal[VEC_CG_R];
      double *  hdGlobal = _daExtraGlobal[VEC_CG_HD];
      for (int  i = 0; i < nextra_global; i++)
	rGlobal[i] -= dAlpha * hdGlobal[i];
    }
    for (int  i = 0; i < nvec; i++)
      _daAVectors[VEC_CG_R][i] -= dAlpha * _daAVectors[VEC_CG_HD][i];
    if (nextra_atom) {
      for (int  m = 0; m < nextra_atom; m++) {
	double *  rAtom  = _daExtraAtom[VEC_CG_R][m];
	double *  hdAtom = _daExtraAtom[VEC_CG_HD][m];
	int  n = extra_nlen[m];
	for (int  i = 0; i < n; i++)
	  rAtom[i] -= dAlpha * hdAtom[i];
      }
    }
    double  dRnewDotRnew = calc_dot_prod_using_mpi_ (VEC_CG_R, VEC_CG_R);
    
    //---- IF RESIDUAL IS SMALL ENOUGH, THEN RETURN THE CURRENT STEP.
    if (sqrt (dRnewDotRnew) < dForceTol * dR0norm2) {
      nStepType = CGSTEP_NEWTON;
      calc_plengths_using_mpi_ (dStepLength2, dStepLengthInf);
      return( true );
    }
    
    //---- beta = r_i+1^T r_i+1 / r_i^T r_i
    //---- d_i+1 = r_i+1 + beta d_i
    double  dBeta = dRnewDotRnew / dRR;
    if (nextra_global) {
      double *  rGlobal = _daExtraGlobal[VEC_CG_R];
      double *  dGlobal = _daExtraGlobal[VEC_CG_D];
      for (int  i = 0; i < nextra_global; i++)
	dGlobal[i] = rGlobal[i] + dBeta * dGlobal[i];
    }
    for (int  i = 0; i < nvec; i++)
      _daAVectors[VEC_CG_D][i] = _daAVectors[VEC_CG_R][i]
	+ dBeta * _daAVectors[VEC_CG_D][i];
    if (nextra_atom) {
      for (int  m = 0; m < nextra_atom; m++) {
	double *  rAtom = _daExtraAtom[VEC_CG_R][m];
	double *  dAtom = _daExtraAtom[VEC_CG_D][m];
	int  n = extra_nlen[m];
	for (int  i = 0; i < n; i++)
	  dAtom[i] = rAtom[i] + dBeta * dAtom[i];
      }
    }
    
    //---- CONTINUE THE LOOP.
    dRR = dRnewDotRnew;
    dPP = dPnewDotPnew;
  }
  
  nStepType = CGSTEP_MAX_INNER_ITERS;
  calc_plengths_using_mpi_ (dStepLength2, dStepLengthInf);
  return( true );    
}

/* ----------------------------------------------------------------------
   Private method calc_xinf_using_mpi_
------------------------------------------------------------------------- */

double MinHFTN::calc_xinf_using_mpi_(void) const
{
  double dXInfLocal = 0.0;
  for (int  i = 0; i < nvec; i++)
    dXInfLocal = MAX(dXInfLocal,fabs(xvec[i]));
  
  double  dXInf;
  MPI_Allreduce (&dXInfLocal, &dXInf, 1, MPI_DOUBLE, MPI_MAX, world);
  
  if (nextra_atom) {
    for (int  m = 0; m < nextra_atom; m++) {
      double *  xatom = xextra_atom[m];
      int  n = extra_nlen[m];
      double  dXInfLocalExtra = 0.0;
      for (int  i = 0; i < n; i++)
	dXInfLocalExtra = MAX (dXInfLocalExtra, fabs (xatom[i]));
      double  dXInfExtra;
      MPI_Allreduce (&dXInfLocalExtra, &dXInfExtra,
		     1, MPI_DOUBLE, MPI_MAX, world);
      dXInf = MAX (dXInf, dXInfExtra);
    }
  }
  
  return( dXInf );
}


/* ----------------------------------------------------------------------
   Private method calc_dot_prod_using_mpi_
------------------------------------------------------------------------- */

double MinHFTN::calc_dot_prod_using_mpi_(const int  nIx1,
					  const int  nIx2) const
{
  double dDotLocal = 0.0;
  for (int  i = 0; i < nvec; i++)
    dDotLocal += _daAVectors[nIx1][i] * _daAVectors[nIx2][i];
  if (nextra_atom) {
    for (int  m = 0; m < nextra_atom; m++) {
      double *  i1Atom = _daExtraAtom[nIx1][m];
      double *  i2Atom = _daExtraAtom[nIx2][m];
      int  n = extra_nlen[m];
      for (int  i = 0; i < n; i++)
	dDotLocal += i1Atom[i] * i2Atom[i];
    }
  }
  
  double  dDot;
  MPI_Allreduce (&dDotLocal, &dDot, 1, MPI_DOUBLE, MPI_SUM, world);
  
  if (nextra_global) {
    for (int  i = 0; i < nextra_global; i++) {
      double *  i1Global = _daExtraGlobal[nIx1];
      double *  i2Global = _daExtraGlobal[nIx2];
      dDot += i1Global[i] * i2Global[i];
    }
  }
  
  return( dDot );
}

/* ----------------------------------------------------------------------
   Private method calc_grad_dot_v_using_mpi_
------------------------------------------------------------------------- */

double MinHFTN::calc_grad_dot_v_using_mpi_(const int  nIx) const
{
  //---- ASSUME THAT FORCES HAVE BEEN EVALUATED AT DESIRED ATOM POSITIONS.
  //---- REMEMBER THAT FORCES = -GRADIENT.
  
  double dGradDotVLocal = 0.0;
  for (int i = 0; i < nvec; i++)
    dGradDotVLocal += - _daAVectors[nIx][i] * fvec[i];
  if (nextra_atom) {
    for (int m = 0; m < nextra_atom; m++) {
      double *  fatom = fextra_atom[m];
      double *  iAtom = _daExtraAtom[nIx][m];
      int  n = extra_nlen[m];
      for (int i = 0; i < n; i++)
	dGradDotVLocal += - iAtom[i] * fatom[i];
    }
  }
  
  double  dGradDotV;
  MPI_Allreduce (&dGradDotVLocal, &dGradDotV, 1, MPI_DOUBLE, MPI_SUM, world);
  
  if (nextra_global) {
    for (int  i = 0; i < nextra_global; i++) {
      double *  iGlobal = _daExtraGlobal[nIx];
      dGradDotV += - iGlobal[i] * fextra[i];
    }
  }
  
  return( dGradDotV );
}

/* ----------------------------------------------------------------------
   Private method calc_dhd_dd_using_mpi_
------------------------------------------------------------------------- */

void MinHFTN::calc_dhd_dd_using_mpi_(double &  dDHD,
				     double &  dDD) const
{
  double dDHDLocal = 0.0;
  double dDDLocal  = 0.0;
  for (int  i = 0; i < nvec; i++) {
    dDHDLocal += _daAVectors[VEC_CG_D][i] * _daAVectors[VEC_CG_HD][i];
    dDDLocal  += _daAVectors[VEC_CG_D][i] * _daAVectors[VEC_CG_D][i];
  }
  if (nextra_atom) {
    for (int  m = 0; m < nextra_atom; m++) {
      double *  dAtom  = _daExtraAtom[VEC_CG_D][m];
      double *  hdAtom = _daExtraAtom[VEC_CG_HD][m];
      int  n = extra_nlen[m];
      for (int  i = 0; i < n; i++) {
	dDHDLocal += dAtom[i] * hdAtom[i];
	dDDLocal += dAtom[i] * dAtom[i];
      }
    }
  }
  
  double  daDotsLocal[2];
  daDotsLocal[0] = dDHDLocal;
  daDotsLocal[1] = dDDLocal;
  double  daDots[2];
  MPI_Allreduce (daDotsLocal, daDots, 2, MPI_DOUBLE, MPI_SUM, world);
  
  if (nextra_global) {
    double *  dGlobal  = _daExtraGlobal[VEC_CG_D];
    double *  hdGlobal = _daExtraGlobal[VEC_CG_HD];
    for (int  i = 0; i < nextra_global; i++) {
      daDots[0] += dGlobal[i] * hdGlobal[i];
      daDots[1] += dGlobal[i] * dGlobal[i];
    }
  }
  dDHD = daDots[0];
  dDD  = daDots[1];
  
  return;
}

/* ----------------------------------------------------------------------
   Private method calc_ppnew_pdold_using_mpi_
------------------------------------------------------------------------- */

void MinHFTN::calc_ppnew_pdold_using_mpi_(double &  dPnewDotPnew,
					  double &  dPoldDotD) const
{
  double dPnewDotPnewLocal = 0.0;
  double dPoldDotDLocal    = 0.0;
  for (int  i = 0; i < nvec; i++) {
    dPnewDotPnewLocal
      += _daAVectors[VEC_CG_P][i] * _daAVectors[VEC_CG_P][i];
    dPoldDotDLocal
      += _daAVectors[VEC_DIF1][i] * _daAVectors[VEC_CG_D][i];
  }
  if (nextra_atom) {
    for (int  m = 0; m < nextra_atom; m++) {
      double *  dAtom  = _daExtraAtom[VEC_CG_D][m];
      double *  pAtom  = _daExtraAtom[VEC_CG_P][m];
      double *  d1Atom = _daExtraAtom[VEC_DIF1][m];
      int  n = extra_nlen[m];
      for (int  i = 0; i < n; i++) {
	dPnewDotPnewLocal += pAtom[i] * pAtom[i];
	dPoldDotDLocal += d1Atom[i] * dAtom[i];
      }
    }
  }
  
  double  daDotsLocal[2];
  daDotsLocal[0] = dPnewDotPnewLocal;
  daDotsLocal[1] = dPoldDotDLocal;
  double  daDots[2];
  MPI_Allreduce (daDotsLocal, daDots, 2, MPI_DOUBLE, MPI_SUM, world);
  
  if (nextra_global) {
    for (int  i = 0; i < nextra_global; i++) {
      double *  dGlobal  = _daExtraGlobal[VEC_CG_D];
      double *  pGlobal  = _daExtraGlobal[VEC_CG_P];
      double *  d1Global = _daExtraGlobal[VEC_DIF1];
      daDots[0] += pGlobal[i] * pGlobal[i];
      daDots[1] += d1Global[i] * dGlobal[i];
    }
  }
  
  dPnewDotPnew = daDots[0];
  dPoldDotD    = daDots[1];
  
  return;
}

/* ----------------------------------------------------------------------
   Private method calc_plengths_using_mpi_
------------------------------------------------------------------------- */

void MinHFTN::calc_plengths_using_mpi_(double &  dStepLength2,
				       double &  dStepLengthInf) const
{
  double dPPLocal   = 0.0;
  double dPInfLocal = 0.0;
  for (int  i = 0; i < nvec; i++) {
    dPPLocal += _daAVectors[VEC_CG_P][i] * _daAVectors[VEC_CG_P][i];
    dPInfLocal = MAX (dPInfLocal, fabs (_daAVectors[VEC_CG_P][i]));
  }
  if (nextra_atom) {
    for (int  m = 0; m < nextra_atom; m++) {
      double *  pAtom = _daExtraAtom[VEC_CG_P][m];
      int  n = extra_nlen[m];
      for (int  i = 0; i < n; i++) {
	dPPLocal += pAtom[i] * pAtom[i];
	dPInfLocal = MAX (dPInfLocal, fabs (pAtom[i]));
      }
    }
  }
  
  double  dPP;
  MPI_Allreduce (&dPPLocal, &dPP, 1, MPI_DOUBLE, MPI_SUM, world);
  
  double  dPInf;
  MPI_Allreduce (&dPInfLocal, &dPInf, 1, MPI_DOUBLE, MPI_MAX, world);
  
  if (nextra_global) {
    for (int  i = 0; i < nextra_global; i++) {
      double *  pGlobal = _daExtraGlobal[VEC_CG_P];
      dPP += pGlobal[i] * pGlobal[i];
      dPInf = MAX (dPInf, fabs (pGlobal[i]));
    }
  }
  
  dStepLength2 = sqrt (dPP);
  dStepLengthInf = dPInf;
  return;
}

/* ----------------------------------------------------------------------
   Private method step_exceeds_TR_
------------------------------------------------------------------------- */

bool MinHFTN::step_exceeds_TR_(const double    dTrustRadius,
			       const double    dPP,
			       const double    dPD,
			       const double    dDD,
			       double &  dTau) const
{
  double  dPnewNorm2;
  double  dPnewNormInf;
  calc_plengths_using_mpi_ (dPnewNorm2, dPnewNormInf);
  
  if (dPnewNorm2 > dTrustRadius) {
    dTau = compute_to_tr_ (dPP, dPD, dDD, dTrustRadius,
			   false, 0.0, 0.0, 0.0);
    return( true );
  }
  
  //---- STEP LENGTH IS NOT TOO LONG.
  dTau = 0.0;
  return( false );
}

/* ----------------------------------------------------------------------
   Private method step_exceeds_DMAX_

   Check that atoms do not move too far:
     for atom coordinates:    limit is dmax
     for extra per-atom DOF:  limit is extra_max[]
     for extra global DOF:    limit is set by modify->max_alpha,
                              which calls fix_box_relax->max_alpha
------------------------------------------------------------------------- */

bool MinHFTN::step_exceeds_DMAX_(void) const
{
  double  dAlpha = dmax * sqrt((double) _nNumUnknowns);
  
  double  dPInfLocal = 0.0;
  for (int  i = 0; i < nvec; i++)
    dPInfLocal = MAX (dPInfLocal, fabs (_daAVectors[VEC_CG_P][i]));
  double  dPInf;
  MPI_Allreduce (&dPInfLocal, &dPInf, 1, MPI_DOUBLE, MPI_MAX, world);
  if (dPInf > dmax)
    return( true );
  if (dPInf > MACHINE_EPS)
    dAlpha = MIN (dAlpha, dmax / dPInf);
  
  if (nextra_atom) {
    for (int  m = 0; m < nextra_atom; m++) {
      double *  pAtom = _daExtraAtom[VEC_CG_P][m];
      dPInfLocal = 0.0;
      int  n = extra_nlen[m];
      for (int  i = 0; i < n; i++)
	dPInfLocal = MAX (dPInfLocal, fabs (pAtom[i]));
      MPI_Allreduce (&dPInfLocal, &dPInf, 1, MPI_DOUBLE, MPI_MAX, world);
      if (dPInf > extra_max[m])
	return( true );
      if (dPInf > MACHINE_EPS)
	dAlpha = MIN (dAlpha, extra_max[m] / dPInf);
    }
  }
  
  if (nextra_global) {
    //---- IF THE MAXIMUM DISTANCE THAT THE GLOBAL BOX CONSTRAINT WILL
    //---- ALLOW IS SMALLER THAN THE PROPOSED DISTANCE, THEN THE STEP
    //---- IS TOO LONG.  PROPOSED DISTANCE IS ESTIMATED BY |P|_INF.
    double  dAlphaExtra = modify->max_alpha (_daExtraGlobal[VEC_CG_P]);
    if (dAlphaExtra < dAlpha)
      return( true );
  }
  
  //---- STEP LENGTH IS NOT TOO LONG.
  return( false );
}

/* ----------------------------------------------------------------------
   Private method adjust_step_to_tau_
   Adjust the step so that VEC_CG_P = VEC_DIF1 + tau * VEC_CG_D.
------------------------------------------------------------------------- */

void MinHFTN::adjust_step_to_tau_(const double tau)
{
  if (nextra_global) {
    double *  pGlobal  = _daExtraGlobal[VEC_CG_P];
    double *  dGlobal  = _daExtraGlobal[VEC_CG_D];
    double *  d1Global = _daExtraGlobal[VEC_DIF1];
    for (int  i = 0; i < nextra_global; i++)
      pGlobal[i] = d1Global[i] + (tau * dGlobal[i]);
  }
  for (int  i = 0; i < nvec; i++) {
    _daAVectors[VEC_CG_P][i] = _daAVectors[VEC_DIF1][i]
      + (tau * _daAVectors[VEC_CG_D][i]);
  }
  if (nextra_atom) {
    for (int  m = 0; m < nextra_atom; m++) {
      double *  pAtom  = _daExtraAtom[VEC_CG_P][m];
      double *  dAtom  = _daExtraAtom[VEC_CG_D][m];
      double *  d1Atom = _daExtraAtom[VEC_DIF1][m];
      int  n = extra_nlen[m];
      for (int  i = 0; i < n; i++)
	pAtom[i] = d1Atom[i] + (tau * dAtom[i]);
    }
  }
  return;
}

/* ----------------------------------------------------------------------
   Private method compute_to_tr_
   Return the value tau that solves
     || p + tau*d ||_2 = dTrustRadius

   If both roots are considered, the TR method chooses the one that minimizes
     grad^T (p + tau*d) + 0.5 (p + tau*d)^T H (p + tau*d)

   @param[in] dPP                - p^T p
   @param[in] dPD                - p^T d
   @param[in] dDD                - d^T d
   @param[in] dTrustRadius       - distance to match
   @param[in] bConsiderBothRoots - evaluate both roots, or return the positive
   @param[in] dDHD               - d^T H d
   @param[in] dPdotHD            - p^T H d
   @param[in] dGradDotD          - grad(x_k)^T d
------------------------------------------------------------------------- */

double MinHFTN::compute_to_tr_(const double  dPP,
			       const double  dPD,
			       const double  dDD,
			       const double  dTrustRadius,
			       const bool    bConsiderBothRoots,
			       const double  dDHD,
			       const double  dPdotHD,
			       const double  dGradDotD) const
{
  //---- SOLVE A QUADRATIC EQUATION FOR TAU.
  //---- THE COEFFICIENTS ARE SUCH THAT THERE ARE ALWAYS TWO REAL ROOTS,
  //---- ONE POSITIVE AND ONE NEGATIVE.
  
  //---- CHECK FOR ERRONEOUS DATA.
  if (   (dDD <= 0.0) || (dPP < 0.0) || (dTrustRadius < 0.0)
	 || (dTrustRadius * dTrustRadius < dPP) ) {
    printf ("HFTN internal error - bad data given to compute_to_tr_()\n");
    return( 0.0 );
  }
  
  double  dTRsqrd = dTrustRadius * dTrustRadius;
  double  dDiscr = (dPD * dPD) - (dDD * (dPP - dTRsqrd));
  dDiscr = MAX (0.0, dDiscr);    //-- SHOULD NEVER BE NEGATIVE
  dDiscr = sqrt (dDiscr);
  
  double  dRootPos = (-dPD + dDiscr) / dDD;
  double  dRootNeg = (-dPD - dDiscr) / dDD;
  
  if (bConsiderBothRoots == false)
    return( dRootPos );
  
  //---- EVALUATE THE CG OBJECTIVE FUNCTION FOR EACH ROOT.
  double  dTmpTerm = dGradDotD + dPdotHD;
  double  dCgRedPos = (dRootPos * dTmpTerm) + (0.5 * dRootPos*dRootPos * dDHD);
  double  dCgRedNeg = (dRootNeg * dTmpTerm) + (0.5 * dRootNeg*dRootNeg * dDHD);
  
  if ((-dCgRedPos) > (-dCgRedNeg))
    return( dRootPos );
  else
    return( dRootNeg );
}

/* ----------------------------------------------------------------------
   Private method evaluate_dir_der_
   Compute the directional derivative using a finite difference approximation.
   This is equivalent to the Hessian times direction vector p.
   As a side effect, the method computes the energy and forces at x.
 
   On input these must be defined:
     atom->x   - positions at x
     atom->f   - ignored
     nIxDir    - VEC_ index of the direction p
     nIxResult - ignored
   On output these are defined:
     atom->x   - unchanged
     atom->f   - forces evaluated at x, only if bEvaluateAtX is true
     nIxDir    - unchanged
     nIxResult - directional derivative Hp
   During processing these are modified:
     VEC_DIF1
     VEC_DIF2
 
   @param[in] bUseForwardDiffs - if true use forward difference approximation,
                                 else use central difference
   @param[in] nIxDir           - VEC_ index of the direction
   @param[in] nIxResult        - VEC_ index to place the result
                                 (it is acceptable for nIxDir = nIxResult)
   @param[in] bEvaluateAtX     - if true, then evaluate at x before returning
   @param[out] dNewEnergy      - energy at x, if bEvaluateAtX is true
   @param[out] dNewForce2      - |F|_2 at x, if bEvaluateAtX is true
------------------------------------------------------------------------- */

void MinHFTN::evaluate_dir_der_(const bool      bUseForwardDiffs,
				const int       nIxDir,
				const int       nIxResult,
				const bool      bEvaluateAtX,
				double &  dNewEnergy)
{
  //---- COMPUTE THE MAGNITUDE OF THE DIRECTION VECTOR:  |p|_2.
  double dDirNorm2SqrdLocal = 0.0;
  for (int  i = 0; i < nvec; i++)
    dDirNorm2SqrdLocal
      += _daAVectors[nIxDir][i] * _daAVectors[nIxDir][i];
  if (nextra_atom) {
    for (int  m = 0; m < nextra_atom; m++) {
      double *  iAtom = _daExtraAtom[nIxDir][m];
      int  n = extra_nlen[m];
      for (int  i = 0; i < n; i++)
	dDirNorm2SqrdLocal += iAtom[i] * iAtom[i];
    }
  }
  double  dDirNorm2Sqrd = 0.0;
  MPI_Allreduce (&dDirNorm2SqrdLocal, &dDirNorm2Sqrd,
		 1, MPI_DOUBLE, MPI_SUM, world);
  if (nextra_global) {
    for (int  i = 0; i < nextra_global; i++) {
      double *  iGlobal = _daExtraGlobal[nIxDir];
      dDirNorm2Sqrd += iGlobal[i] * iGlobal[i];
    }
  }
  double  dDirNorm2 = sqrt (dDirNorm2Sqrd);
  
  //---- IF THE STEP IS TOO SMALL, RETURN ZERO FOR THE DERIVATIVE.
  if (dDirNorm2 == 0.0) {
    for (int  i = 0; i < nvec; i++)
      _daAVectors[nIxResult][i] = 0.0;
    if (nextra_atom) {
      for (int  m = 0; m < nextra_atom; m++) {
	double *  iAtom = _daExtraAtom[nIxDir][m];
	int  n = extra_nlen[m];
	for (int  i = 0; i < n; i++)
	  iAtom[i] = 0;
      }
    }
    if (nextra_global) {
      for (int  i = 0; i < nextra_global; i++)
	_daExtraGlobal[nIxDir][i] = 0.0;
    }
    
    if (bEvaluateAtX) {
      dNewEnergy = energy_force (0);
      neval++;
    }
    
    return;
  }
  
  //---- FORWARD DIFFERENCES ARE LESS ACCURATE THAN CENTRAL DIFFERENCES,
  //---- BUT REQUIRE ONLY 2 ENERGY+FORCE EVALUATIONS VERSUS 3 EVALUATIONS.
  //---- STORAGE REQUIREMENTS ARE THE SAME.
  
  if (bUseForwardDiffs) {
    //---- EQUATION IS FROM THE OLD LAMMPS VERSION, SAND98-8201.
    double  dEps = 2.0 * sqrt (1000.0 * MACHINE_EPS) / dDirNorm2;
    
    //---- SAVE A COPY OF x.
    fix_minimize->store_box();
    for (int  i = 0; i < nvec; i++)
      _daAVectors[VEC_DIF1][i] = xvec[i];
    if (nextra_atom) {
      for (int  m = 0; m < nextra_atom; m++) {
	double *  xatom  = xextra_atom[m];
	double *  d1Atom = _daExtraAtom[VEC_DIF1][m];
	int  n = extra_nlen[m];
	for (int  i = 0; i < n; i++)
	  d1Atom[i] = xatom[i];
      }
    }
    if (nextra_global) {
      modify->min_pushstore();
      modify->min_store();
    }
    
    //---- EVALUATE FORCES AT x + eps*p.
    if (nextra_global)
      modify->min_step (dEps, _daExtraGlobal[nIxDir]);
    for (int  i = 0; i < nvec; i++)
      xvec[i] += dEps * _daAVectors[nIxDir][i];
    if (nextra_atom) {
      for (int  m = 0; m < nextra_atom; m++) {
	double *  xatom = xextra_atom[m];
	double *  iAtom = _daExtraAtom[nIxDir][m];
	int  n = extra_nlen[m];
	for (int  i = 0; i < n; i++)
	  xatom[i] += dEps * iAtom[i];
	requestor[m]->min_x_set(m);
      }
    }
    energy_force (0);
    neval++;
    
    //---- STORE THE FORCE IN DIF2.
    if (nextra_global) {
      for (int  i = 0; i < nextra_global; i++)
	_daExtraGlobal[VEC_DIF2][i] = fextra[i];
    }
    for (int  i = 0; i < nvec; i++)
      _daAVectors[VEC_DIF2][i] = fvec[i];
    if (nextra_atom) {
      for (int  m = 0; m < nextra_atom; m++) {
	double *  fatom  = fextra_atom[m];
	double *  d2Atom = _daExtraAtom[VEC_DIF2][m];
	int  n = extra_nlen[m];
	for (int  i = 0; i < n; i++)
	  d2Atom[i] = fatom[i];
      }
    }
    
    //---- MOVE BACK TO x AND EVALUATE FORCES.
    if (nextra_global) {
      modify->min_step (0.0, _daExtraGlobal[VEC_DIF1]);
      modify->min_popstore();
    }
    for (int  i = 0; i < nvec; i++)
      xvec[i] = _daAVectors[VEC_DIF1][i];
    if (nextra_atom) {
      for (int  m = 0; m < nextra_atom; m++) {
	double *  xatom  = xextra_atom[m];
	double *  d1Atom = _daExtraAtom[VEC_DIF1][m];
	int  n = extra_nlen[m];
	for (int  i = 0; i < n; i++)
	  xatom[i] += d1Atom[i];
	requestor[m]->min_x_set(m);
      }
    }
    dNewEnergy = energy_force (0);
    neval++;
    
    //---- COMPUTE THE DIFFERENCE VECTOR:  [grad(x + eps*p) - grad(x)] / eps.
    //---- REMEMBER THAT FORCES = -GRADIENT.
    for (int  i = 0; i < nvec; i++)
      _daAVectors[nIxResult][i] = (fvec[i] - _daAVectors[VEC_DIF2][i]) / dEps;
    if (nextra_atom) {
      for (int  m = 0; m < nextra_atom; m++) {
	double *  iAtom  = _daExtraAtom[nIxResult][m];
	double *  d2Atom = _daExtraAtom[VEC_DIF2][m];
	int  n = extra_nlen[m];
	for (int  i = 0; i < n; i++)
	  iAtom[i] = (fextra_atom[m][i] - d2Atom[i]) / dEps;
      }
    }
    if (nextra_global) {
      for (int  i = 0; i < nextra_global; i++)
	_daExtraGlobal[nIxResult][i]
	  = (fextra[i] - _daExtraGlobal[VEC_DIF2][i]) / dEps;
    }
  }
  
  else {    //-- bUseForwardDiffs == false
    //---- EQUATION IS FROM THE OLD LAMMPS VERSION, SAND98-8201.
    double  dEps = pow (3000.0 * MACHINE_EPS, 0.33333333) / dDirNorm2;
    
    //---- SAVE A COPY OF x.
    fix_minimize->store_box();
    for (int  i = 0; i < nvec; i++)
      _daAVectors[VEC_DIF1][i] = xvec[i];
    if (nextra_atom) {
      for (int  m = 0; m < nextra_atom; m++) {
	double *  xatom  = xextra_atom[m];
	double *  d1Atom = _daExtraAtom[VEC_DIF1][m];
	int  n = extra_nlen[m];
	for (int  i = 0; i < n; i++)
	  d1Atom[i] = xatom[i];
      }
    }
    if (nextra_global) {
      modify->min_pushstore();
      modify->min_store();
    }
    
    //---- EVALUATE FORCES AT x + eps*p.
    if (nextra_global)
      modify->min_step (dEps, _daExtraGlobal[nIxDir]);
    for (int  i = 0; i < nvec; i++)
      xvec[i] += dEps * _daAVectors[nIxDir][i];
    if (nextra_atom) {
      for (int  m = 0; m < nextra_atom; m++) {
	double *  xatom = xextra_atom[m];
	double *  iAtom = _daExtraAtom[nIxDir][m];
	int  n = extra_nlen[m];
	for (int  i = 0; i < n; i++)
	  xatom[i] += dEps * iAtom[i];
	requestor[m]->min_x_set(m);
      }
    }
    energy_force (0);
    neval++;
    
    //---- STORE THE FORCE IN DIF2.
    if (nextra_global) {
      for (int  i = 0; i < nextra_global; i++)
	_daExtraGlobal[VEC_DIF2][i] = fextra[i];
    }
    for (int  i = 0; i < nvec; i++)
      _daAVectors[VEC_DIF2][i] = fvec[i];
    if (nextra_atom) {
      for (int  m = 0; m < nextra_atom; m++) {
	double *  fatom  = fextra_atom[m];
	double *  d2Atom = _daExtraAtom[VEC_DIF2][m];
	int  n = extra_nlen[m];
	for (int  i = 0; i < n; i++)
	  d2Atom[i] = fatom[i];
      }
    }
    
    //---- EVALUATE FORCES AT x - eps*p.
    if (nextra_global)
      modify->min_step (-dEps, _daExtraGlobal[nIxDir]);
    for (int  i = 0; i < nvec; i++)
      xvec[i] = _daAVectors[VEC_DIF1][i]
	- dEps * _daAVectors[nIxDir][i];
    if (nextra_atom) {
      for (int  m = 0; m < nextra_atom; m++) {
	double *  xatom  = xextra_atom[m];
	double *  iAtom  = _daExtraAtom[nIxDir][m];
	double *  d1Atom = _daExtraAtom[VEC_DIF1][m];
	int  n = extra_nlen[m];
	for (int  i = 0; i < n; i++)
	  xatom[i] = d1Atom[i] - dEps * iAtom[i];
	requestor[m]->min_x_set(m);
      }
    }
    energy_force (0);
    neval++;
    
    //---- COMPUTE THE DIFFERENCE VECTOR:
    //----     [grad(x + eps*p) - grad(x - eps*p)] / 2*eps.
    //---- REMEMBER THAT FORCES = -GRADIENT.
    if (nextra_global) {
      double *  iGlobal  = _daExtraGlobal[nIxResult];
      double *  d2Global = _daExtraGlobal[VEC_DIF2];
      for (int  i = 0; i < nextra_global; i++)
	iGlobal[i] = (fextra[i] - d2Global[i]) / (2.0 + dEps);
    }
    for (int  i = 0; i < nvec; i++)
      _daAVectors[nIxResult][i] = 
	(fvec[i] - _daAVectors[VEC_DIF2][i]) / (2.0 * dEps);
    if (nextra_atom) {
      for (int  m = 0; m < nextra_atom; m++) {
	double *  fatom  = fextra_atom[m];
	double *  iAtom  = _daExtraAtom[nIxResult][m];
	double *  d2Atom = _daExtraAtom[VEC_DIF2][m];
	int  n = extra_nlen[m];
	for (int  i = 0; i < n; i++)
	  iAtom[i] = (fatom[i] - d2Atom[i]) / (2.0 + dEps);
      }
    }
    
    if (bEvaluateAtX) {
      //---- EVALUATE FORCES AT x.
      if (nextra_global) {
	modify->min_step (0.0, _daExtraGlobal[VEC_DIF1]);
	modify->min_popstore();
      }
      for (int  i = 0; i < nvec; i++)
	xvec[i] = _daAVectors[VEC_DIF1][i];
      if (nextra_atom) {
	for (int  m = 0; m < nextra_atom; m++) {
	  double *  xatom  = xextra_atom[m];
	  double *  d1Atom = _daExtraAtom[VEC_DIF1][m];
	  int  n = extra_nlen[m];
	  for (int  i = 0; i < n; i++)
	    xatom[i] = d1Atom[i];
	  requestor[m]->min_x_set(m);
	}
      }
      dNewEnergy = energy_force (0);
      neval++;
    }
  }
  
  return;
}

/* ----------------------------------------------------------------------
   Private method open_hftn_print_file_
------------------------------------------------------------------------- */

void MinHFTN::open_hftn_print_file_(void)
{
  int  nMyRank;
  MPI_Comm_rank (world, &nMyRank);
  
  char  szTmp[50];
  sprintf (szTmp, "progress_MinHFTN_%d.txt", nMyRank);
  _fpPrint = fopen (szTmp, "w");
  if (_fpPrint == NULL) {
    printf ("*** MinHFTN cannot open file '%s'\n", szTmp);
    printf ("*** continuing...\n");
    return;
  }
  
  fprintf (_fpPrint, "  Iter   Evals      Energy         |F|_2"
	   "    Step   TR used    |step|_2      ared      pred\n");
  return;
}

/* ----------------------------------------------------------------------
   Private method hftn_print_line_
   Step types:
   1 - Nw   (inner iteration converged like a Newton step)
   2 - TR   (inner iteration reached the trust region boundary)
   3 - Neg  (inner iteration ended with negative curvature)
------------------------------------------------------------------------- */

void MinHFTN::hftn_print_line_(const bool    bIsStepAccepted,
			       const int     nIteration,
			       const int     nTotalEvals,
			       const double  dEnergy,
			       const double  dForce2,
			       const int     nStepType,
			       const double  dTrustRadius,
			       const double  dStepLength2,
			       const double  dActualRed,
			       const double  dPredictedRed) const
{
  const char  sFormat1[]
    = "  %4d   %5d  %14.8f  %11.5e\n";
  const char  sFormatA[]
    = "  %4d   %5d  %14.8f  %11.5e  %3s  %9.3e   %8.2e  %10.3e %10.3e\n";
  const char  sFormatR[]
    = "r %4d   %5d  %14.8f  %11.5e  %3s  %9.3e   %8.2e  %10.3e %10.3e\n";
  
  if (_fpPrint == NULL)
    return;
  
  char  sStepType[4];
  if (nStepType == NO_CGSTEP_BECAUSE_F_TOL_SATISFIED)
    strcpy (sStepType, " - ");
  else if (nStepType == CGSTEP_NEWTON)
    strcpy (sStepType, "Nw ");
  else if (nStepType == CGSTEP_TO_TR)
    strcpy (sStepType, "TR ");
  else if (nStepType == CGSTEP_TO_DMAX)
    strcpy (sStepType, "dmx");
  else if (nStepType == CGSTEP_NEGATIVE_CURVATURE)
    strcpy (sStepType, "Neg");
  else if (nStepType == CGSTEP_MAX_INNER_ITERS)
    strcpy (sStepType, "its");
  else
    strcpy (sStepType, "???");
  
  if (nIteration == -1) {
    fprintf (_fpPrint, sFormat1,
	     0, nTotalEvals, dEnergy, dForce2);
  }
  else {
    if (bIsStepAccepted)
      fprintf (_fpPrint, sFormatA,
	       nIteration, nTotalEvals, dEnergy, dForce2,
	       sStepType, dTrustRadius, dStepLength2,
	       dActualRed, dPredictedRed);
    else
      fprintf (_fpPrint, sFormatR,
	       nIteration, nTotalEvals, dEnergy, dForce2,
	       sStepType, dTrustRadius, dStepLength2,
	       dActualRed, dPredictedRed);
  }
  
  fflush (_fpPrint);
  return;
}

/* ----------------------------------------------------------------------
   Private method close_hftn_print_file_
------------------------------------------------------------------------- */

void MinHFTN::close_hftn_print_file_(void)
{
  if (_fpPrint != NULL) fclose (_fpPrint);
  return;
}

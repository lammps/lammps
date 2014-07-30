/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef MINIMIZE_CLASS

MinimizeStyle(hftn,MinHFTN)

#else

#ifndef LMP_MIN_HFTN_H
#define LMP_MIN_HFTN_H

#include "min.h"

namespace LAMMPS_NS
{
  //---- THE ALGORITHM NEEDS TO STORE THIS MANY ATOM-BASED VECTORS,
  //---- IN ADDITION TO ATOM POSITIONS AND THE FORCE VECTOR.

  static const int  NUM_HFTN_ATOM_BASED_VECTORS = 7;
  static const int  VEC_XK    = 0;   //-- ATOM POSITIONS AT SUBITER START
  static const int  VEC_CG_P  = 1;   //-- STEP p IN CG SUBITER
  static const int  VEC_CG_D  = 2;   //-- DIRECTION d IN CG SUBITER
  static const int  VEC_CG_HD = 3;   //-- HESSIAN-VECTOR PRODUCT Hd
  static const int  VEC_CG_R  = 4;   //-- RESIDUAL r IN CG SUBITER
  static const int  VEC_DIF1  = 5;   //-- FOR FINITE DIFFERENCING
  static const int  VEC_DIF2  = 6;   //-- FOR FINITE DIFFERENCING

class MinHFTN : public Min
{
  public:

  MinHFTN (LAMMPS *);
  ~MinHFTN (void);
  void init();
  void setup_style();
  void reset_vectors();
  int iterate (int);

  private:

  //---- ATOM-BASED STORAGE VECTORS.
  double *   _daAVectors[NUM_HFTN_ATOM_BASED_VECTORS];
  double **  _daExtraAtom[NUM_HFTN_ATOM_BASED_VECTORS];

  //---- GLOBAL DOF STORAGE.  ELEMENT [0] IS X0 (XK), NOT USED IN THIS ARRAY.
  double *   _daExtraGlobal[NUM_HFTN_ATOM_BASED_VECTORS];

  int     _nNumUnknowns;
  FILE *  _fpPrint;

  int   execute_hftn_ (const bool      bPrintProgress,
                       const double    dInitialEnergy,
                       const double    dInitialForce2,
                       double &  dFinalEnergy,
                       double &  dFinalForce2);
  bool    compute_inner_cg_step_ (const double    dTrustRadius,
                                  const double    dForceTol,
                                  const int       nMaxEvals,
                                  const bool      bHaveEvalAtXin,
                                  const double    dEnergyAtXin,
                                  const double    dForce2AtXin,
                                  double &  dEnergyAtXout,
                                  double &  dForce2AtXout,
                                  int    &  nStepType,
                                  double &  dStepLength2,
                                  double &  dStepLengthInf);
  double  calc_xinf_using_mpi_ (void) const;
  double  calc_dot_prod_using_mpi_ (const int  nIx1,
                                    const int  nIx2) const;
  double  calc_grad_dot_v_using_mpi_ (const int  nIx) const;
  void    calc_dhd_dd_using_mpi_ (double &  dDHD,
                                  double &  dDD) const;
  void    calc_ppnew_pdold_using_mpi_ (double &  dPnewDotPnew,
                                       double &  dPoldDotD) const;
  void    calc_plengths_using_mpi_ (double &  dStepLength2,
                                    double &  dStepLengthInf) const;
  bool    step_exceeds_TR_ (const double    dTrustRadius,
                            const double    dPP,
                            const double    dPD,
                            const double    dDD,
                            double &  dTau) const;
  bool    step_exceeds_DMAX_ (void) const;
  void    adjust_step_to_tau_ (const double  tau);
  double  compute_to_tr_ (const double  dPP,
                          const double  dPD,
                          const double  dDD,
                          const double  dTrustRadius,
                          const bool    bConsiderBothRoots,
                          const double  dDHD,
                          const double  dPdotHD,
                          const double  dGradDotD) const;
  void    evaluate_dir_der_ (const bool      bUseForwardDiffs,
                             const int       nIxDir,
                             const int       nIxResult,
                             const bool      bEvaluateAtX,
                             double &  dNewEnergy);
  void  open_hftn_print_file_ (void);
  void  hftn_print_line_ (const bool    bIsStepAccepted,
                          const int     nIteration,
                          const int     nTotalEvals,
                          const double  dEnergy,
                          const double  dForce2,
                          const int     nStepType,
                          const double  dTrustRadius,
                          const double  dStepLength2,
                          const double  dActualRed,
                          const double  dPredictedRed) const;
  void  close_hftn_print_file_ (void);
};

}

#endif
#endif

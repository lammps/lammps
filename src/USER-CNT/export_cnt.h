/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
   
   Contributing author: Maxim Shugaev (UVA), mvs9t@virginia.edu
------------------------------------------------------------------------- */

#ifdef __cplusplus
extern "C" {
#endif
  // see ExportCNT.f90 in lib/cnt for function details
  void TPBInit();
  void TPMInit(const int& M, const int& N);
  void SetTablePath(const char* TPMSSTPFile, const int& N1, 
    const char* TPMAFile, const int& N2);
    
  void InitCNTPotModule(const int& STRModel, const int& STRParams,
    const int& YMType, const int& BNDModel, const double& Rref);
    
  double get_R();

  void TubeStretchingForceField(double& U1, double& U2, double* F1,
    double* F2, double* S1, double* S2, const double* X1, const double* X2,
    const double& R12, const double& L12);

  void TubeBendingForceField(double& U1, double& U2, double& U3, double* F1,
    double* F2, double* F3, double* S1, double* S2, double* S3, 
    const double* X1, const double* X2, const double* X3, const double& R123,
    const double& L123, int& BBF2 );
    
  void SegmentTubeForceField(double& U1, double& U2, double *U, double* F1,
    double* F2, double* F, double* Fe, double* S1, double* S2, double* S,
    double* Se, const double* X1, const double* X2,  const double& R12,
    const int& N, const double* X, const double* Xe, const int* BBF, 
    const double& R, const int& E1, const int& E2, const int& Ee, 
    const int& TPMType);
      
#ifdef __cplusplus
}
#endif

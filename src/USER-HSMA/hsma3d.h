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

/* ----------------------------------------------------------------------
   Contributing author: Jiuyang Liang (liangjiuyang@sjtu.edu.cn)
------------------------------------------------------------------------- */

#ifdef KSPACE_CLASS

KSpaceStyle(HSMA3D,HSMA3D)

#else

#ifndef LMP_HSMA3D_H
#define LMP_HSMA3D_H

#include "kspace.h"

namespace LAMMPS_NS {

class HSMA3D : public KSpace {
public:
  HSMA3D(class LAMMPS *);
  void settings(int, char **);
  void setup() {}
  void compute(int, int);
  double memory_usage();
  void init();

  void CalculateRDMultipoleExpansion(double* Q, int p, double x, double y, double z);//Q is a p*p vector.  Set multi-pole expansion coefficient
  void CalculateLocalRDMultipoleExpansion(double* Q, int p, double x, double y, double z, double Rs);//Q is a p*p vector.  Set multi-pole expansion coefficient
  void CalculateMultipoleExpansion(double* Q, int p, double x, double y, double z);//Q is a p*p vector.  Set multi-pole expansion coefficient
  void CalculateZDerivativeMultipoleExpansion(double* Q, int p, double x, double y, double z);//Multipole bases with Z-direction derivative
  void CalculateXDMultipoleExpansion(double* Q, int p, double x, double y, double z);
  void CalculateYDMultipoleExpansion(double* Q, int p, double x, double y, double z);

  void SetFibonacci(double Fibonacci[][4], double F, double Fp, int Np, double Rs, double PI);//Construct Fibonacci Points And There Local/Multipole Expansion.  Fibonacci[Np][4] storage the positions and the weights.
  void SetImageCharge(double ImageCharge[][4], int* ImageNumber, int TotalNumber, double Source[][3], double* Q, int NSource, double Rs, double Lx, double Ly, double Lz, int lx, int ly, int lz);//Find the image charge
  void AdjustParticle_Float(float Particle[][3], int N, float Lx, float Ly, float Lz);//Image Correction
  void AdjustParticle_Double(double Particle[][3], int N, double Lx, double Ly, double Lz);
  void CalculateNearFieldAndZD(double* Near, double ImageCharge[][4], int ImageNumber, int Nw, int IF_FMM_RightTerm, double Lx, double Ly, double Lz, double **PointSum, double **QuizSum,double tolerance);//Calculate near field potential
  void SolveLeastSquareProblem(double* C, double** A, double* Near, int p, int Nw);
  double FinalCalculateEnergyAndForce(double Force[][3], double* Pot, double Source[][3], double* Q, int NSource, double ImageCharge[][4], int ImageNumber, double **Fibonacci, double** QRD, double** QLocalRD, double* C, int p, double Fp, double F, double Rs, double PI, int IF_FMM_FinalPotential, double tolerance);
  void CalculateNearFieldAndZD_Single(double* Near, double ImageCharge[][4], int ImageNumber, int Nw, int IF_FMM_RightTerm, double Lx, double Ly, double Lz, double **PointSum, double **QuizSum, double Source[][3], double Force[][3], double* Pot, int NSource, double* Q, double tolerance);
  double FinalCalculateEnergyAndForce_Single(double Force[][3], double* Pot, double Source[][3], double* Q, int NSource, double ImageCharge[][4], int ImageNumber, double **Fibonacci, double** QRD, double** QLocalRD, double* C, int p, double Fp, double F, double Rs, double PI, int IF_FMM_FinalPotential, double tolerance);

 private:
  int me,RankID;
  double Lx, Ly, Lz;
  double Lambda;
  double Fp, F;
  int p;
  int Nw;
  int IF_FMM_RightTerm;
  int IF_FMM_FinalPotential;
  double **PointSum, **QuizSum;
  double **A, * PointSumMultipleExpansionMatrix, * QuizSumMultipleExpansionMatrix;
  double pi;
  double R, Rs;
  double **Fibonacci;
  double ** QRD, ** QLocalRD;
  int Np;
  int maxatom;

  double tolerance;
  int Step;
  float* Time;
};

}

#endif
#endif

/* ERROR/WARNING messages:

*/

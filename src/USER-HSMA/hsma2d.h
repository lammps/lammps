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

KSpaceStyle(HSMA2D, HSMA2D)

#else

#ifndef LMP_HSMA2D_H
#define LMP_HSMA2D_H

#include "kspace.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "math.h"
#include "memory.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>

namespace LAMMPS_NS {

	class HSMA2D : public KSpace {
	public:
		HSMA2D(class LAMMPS*);
		void settings(int, char**);
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

		void Gauss_Grule3D(double * Point, double* Weight, int N);
		void ConstructLeftTerm(double ** LeftTermReal, double** LeftTermImag, double*** IntegralTop, double*** IntegralDown, int Nw, int NJKBound, int p, double EU, double mu, double Lx, double Ly, double Lz, int S);
		void SetFibonacci(double **Fibonacci, double F, double Fp, int Np, double Rs, double PI);//Construct Fibonacci Points And There Local/Multipole Expansion.  Fibonacci[Np][4] storage the positions and the weights.
		void AdjustParticle_Double(double Particle[][3], int N, double Lx, double Ly, double Lz);
		double SetImageCharge(double ImageCharge[][5], int* ImageNumber, int TotalNumber, double Source[][3], double* Q, int NSource, double Rs, double Lx, double Ly, double Lz, int lx, int ly, int lz);
		void CalculateNearFieldAndZD(double** Top, double** TopZD, double** Down, double** DownZD, double ImageCharge[][5], int ImageNumber, double*** IntegralTop, double*** IntegralDown, int Nw, double Gamma, int IF_FMM_RightTerm, int S, double* AR, double Lx, double Ly, double Lz, double tolerance);
		void ConstructRightTerm(double* RightTermReal, double* RightTermImag, double** TopNear, double** TopZDNear, double** DownNear, double** DownZDNear, double*** IntegralTop, double*** IntegralDown, int Nw, int NJKBound, double EU, double mu, double Lx, double Ly, double Lz, int S, double tolerance);
		void SolveLeastSquareProblem(double* C, double** LeftTermReal, double** LeftTermImag, double* RightTermReal, double* RightTermImag, int p, int row);
		double FinalCalculateEnergyAndForce(double Force[][3], double* Pot, double Source[][3], double* Q, int NSource, double ImageCharge[][5], int ImageNumber, double** Fibonacci, double** QRD, double** QLocalRD, double Gamma, double* C, int p, double Fp, double F, double Rs, double PI, int IF_FMM_FinalPotential, double tolerance);
		double fac(double t);

	private:
		int me, RankID;
		double Lx, Ly, Lz;
		double Lambda;
		double Gamma;
		double Fp, F;
		int p;
		int Nw;
		double w;
		int IF_FMM_RightTerm;
		int IF_FMM_FinalPotential;

		double *Gauss, *Weight;
		double PI,mu,EU;
		int S;
		double *** IntegralTop, *** IntegralDown;
		double Rs;
		int IndexReal;
		int IndexImag, NJKBound, NJK;
		double** LeftTermReal, ** LeftTermImag;

		int Np;
		double** Fibonacci;
		double** QRD, ** QLocalRD;
		int maxatom;
		
		double ** TopNear, ** TopZDNear;
		double ** DownNear, ** DownZDNear;
		double *AR;

		double tolerance;
		int Step;
		float* Time;
	};

}

#endif
#endif

/* ERROR/WARNING messages:

*/

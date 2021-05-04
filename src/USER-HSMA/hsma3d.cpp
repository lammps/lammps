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
   Contributing author: Jiuyang Liang (liangjiuyang@sjtu.edu.cn)
------------------------------------------------------------------------- */

#include "hsma3d.h"

#include "atom.h"
#include "comm.h"
#include "complex.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "math.h"
#include "math_special.h"
#include "memory.h"
#include "update.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include<iomanip>
#include<iostream>
#include <immintrin.h>

#if defined(_OPENMP)
	#include<omp.h>
#endif

extern "C" {void lfmm3d_t_c_g_(double *eps, int *nsource,double *source, double *charge, int *nt, double *targ, double *pottarg, double *gradtarg, int *ier);}
extern int fab(int n);
extern int isfab(int m);

using namespace LAMMPS_NS;
using namespace std;

/* ---------------------------------------------------------------------- */
HSMA3D::HSMA3D(LAMMPS *lmp) : KSpace(lmp)
{
  maxatom = atom->natoms;
  me = comm->me;
  RankID = comm->nprocs;
  Lx = domain->xprd;
  Ly = domain->yprd;
  Lz = domain->zprd * slab_volfactor;
  pi = 3.141592653589793;
}

/* ---------------------------------------------------------------------- */
void HSMA3D::settings(int narg, char **arg)
{
  if (narg != 8) error->all(FLERR,"Illegal kspace_style FMM command");
  tolerance = fabs(utils::numeric(FLERR,arg[0],false,lmp));
  Lambda = utils::numeric(FLERR, arg[1], false, lmp);
  p = utils::numeric(FLERR, arg[2], false, lmp);
  Nw = utils::numeric(FLERR, arg[3], false, lmp);
  Fp = utils::numeric(FLERR, arg[4], false, lmp);
  F = utils::numeric(FLERR, arg[5], false, lmp);
  IF_FMM_RightTerm = utils::numeric(FLERR, arg[6], false, lmp);
  IF_FMM_FinalPotential = utils::numeric(FLERR, arg[7], false, lmp);
  Step = 0;
  Time = new float[2000];

  if(Lambda<0)
	  error->all(FLERR, "Lambda shoule be >0.");
  if (Lambda > 20 || Lambda < 0.2) {
	  error->warning(FLERR, fmt::format("The Lambda is too big or too small! Please use an approximate range of Lambda. Set Lambda to the default value.",
		  update->ntimestep));
	  Lambda = 1.3;
  }
  if (p<1)
	  error->all(FLERR, "p shoule be >=1.");
  if (p > 50){
	  error->warning(FLERR, fmt::format("The p is too big! Please use a smaller p. Set p to the default value.",
		  update->ntimestep));
	  p = 8;
  }
  if((IF_FMM_RightTerm)!=0&& (IF_FMM_RightTerm!=1))
	  error->all(FLERR, "Using wrong value of IF_FMM_RightTerm.");
  if ((IF_FMM_FinalPotential) != 0 && (IF_FMM_FinalPotential != 1))
	  error->all(FLERR, "Using wrong value of IF_FMM_FinalPotential.");
  if (Fp > F || Fp < 0 || F < 0 || isfab(int(Fp)) != 1 || isfab(int(F)) != 1) {
	  error->warning(FLERR, fmt::format("Using wrong value of Fp and F. Set Fp and F to the default value.",
		  update->ntimestep));
	  Fp = 89.0;
	  F = 144.0;
  }
  if (Nw <= p * p)
  {
	  error->warning(FLERR, fmt::format("Nw is too small! Set Nw to the default value.",
		  update->ntimestep));
	  Nw = int(2 * p * p);
  }
  if (domain->dimension == 2)
	  error->all(FLERR, "Cannot use HSMA3D in a 2d simulation. Please use HSMA2D instead.");

}

void HSMA3D::init()
{
  printf("Setting up HSMA implemented by Jiuyang Liang (Release 1.0.0)\n");
  PointSum = new double * [Nw];
  QuizSum = new double * [Nw];
  for (int i = 0; i < Nw; i++)
  {
      PointSum[i] = new double[3];
      QuizSum[i] = new double[3];
  }

  R= sqrt((Lx / 2) * (Lx / 2) + (Ly / 2) * (Ly / 2) + (Lz / 2) * (Lz / 2));
  Rs = Lambda* sqrt((Lx / 2) * (Lx / 2) + (Ly / 2) * (Ly / 2) + (Lz / 2) * (Lz / 2));//The radius of the cut-off sphere B
  for (int i = 0; i < Nw; i++)
  {
	  PointSum[i][2] = (2 * (i + 1) - 1) / (Nw + 0.00) - 1;
	  PointSum[i][0] = (sqrt(1 - PointSum[i][2] * PointSum[i][2]) * cos(2 * pi * (i + 1) * 0.618)) * R;
	  PointSum[i][1] = (sqrt(1 - PointSum[i][2] * PointSum[i][2]) * sin(2 * pi * (i + 1) * 0.618)) * R;
	  PointSum[i][2] = PointSum[i][2] * R;


	  if (abs(PointSum[i][0]) >= (Lx / 2))
	  {
		  QuizSum[i][0] = (abs(PointSum[i][0]) - Lx * int((abs(PointSum[i][0]) + Lx / 2) / Lx)) * PointSum[i][0] / abs(PointSum[i][0]);
	  }
	  else
	  {
		  QuizSum[i][0] = PointSum[i][0];
	  }
	  if (abs(PointSum[i][1]) >= (Ly / 2))
	  {
		  QuizSum[i][1] = (abs(PointSum[i][1]) - Ly * int((abs(PointSum[i][1]) + Ly / 2) / Ly)) * PointSum[i][1] / abs(PointSum[i][1]);
	  }
	  else
	  {
		  QuizSum[i][1] = PointSum[i][1];
	  }
	  if (abs(PointSum[i][2]) >= (Lz / 2))
	  {
		  QuizSum[i][2] = (abs(PointSum[i][2]) - Lz * int((abs(PointSum[i][2]) + Lz / 2) / Lz)) * PointSum[i][2] / abs(PointSum[i][2]);
	  }
	  else
	  {
		  QuizSum[i][2] = PointSum[i][2];
	  }
  }

  A = new double* [Nw];
  PointSumMultipleExpansionMatrix = new double[p * p];
  QuizSumMultipleExpansionMatrix = new double[p * p];
  for (int i = 0; i < Nw; i++)
  {
	  A[i] = new double[p*p-1];
  }
  for (int i = 0; i < Nw; i++)
  {
	  CalculateMultipoleExpansion(PointSumMultipleExpansionMatrix, p, PointSum[i][0], PointSum[i][1], PointSum[i][2]);
	  CalculateMultipoleExpansion(QuizSumMultipleExpansionMatrix, p, QuizSum[i][0], QuizSum[i][1], QuizSum[i][2]);
	  for (int j = 0; j < p * p - 1; j++)
	  {
		  A[i][j] = PointSumMultipleExpansionMatrix[j + 1] - QuizSumMultipleExpansionMatrix[j + 1];
	  }
  }

  //Construct Fibonacci Points And There Local/Multipole Expansion.  Fibonacci[Np][4] storage the positions and the weights.
  Np = int(2 * F + 2.0);
  Fibonacci = new double* [Np];
  for (int i = 0; i < Np; i++)
  {
	  Fibonacci[i] = new double[4];
  }
  double Fibonacci_New[Np][4];
  SetFibonacci(Fibonacci_New, F, Fp, Np, Rs, pi);
  for (int i = 0; i < Np; i++)
	  for (int j = 0; j < 4; j++)
		  Fibonacci[i][j] = Fibonacci_New[i][j];

  QRD = new double* [Np]; QLocalRD= new double* [Np];
  for (int i = 0; i < Np; i++)
  {
	  QRD[i] = new double[p*p];
	  QLocalRD[i]= new double[p * p];
  }
  #pragma omp parallel
  {
	  double MulQ[p * p], MulLocalQ[p * p];
      #pragma omp for schedule(static) private(MulQ,MulLocalQ)
	  for (int i = 0; i < Np; i++)
	  {
		  CalculateRDMultipoleExpansion(MulQ, p, Fibonacci[i][0], Fibonacci[i][1], Fibonacci[i][2]);
		  CalculateLocalRDMultipoleExpansion(MulLocalQ, p, Fibonacci[i][0], Fibonacci[i][1], Fibonacci[i][2], Rs);
		  for (int j = 0; j < p * p; j++)
		  {
			  QRD[i][j] = MulQ[j];
			  QLocalRD[i][j] = MulLocalQ[j];
		  }
	  }
  }
  cout << Lx << "=Lx   " << Ly << "=Ly    " << Lz << "=Lz   " << Lambda << "=Lambda    " << p << "=p   " << Nw << "=Nw    "  << Fp << "=Fp   " << F << "=F   " << IF_FMM_RightTerm << "=IF_FMM_RightTerm   " << IF_FMM_FinalPotential << "=IF_FMM_FinalPotential   " << endl;
}

void HSMA3D::compute(int eflag, int vflag)
{
	// set energy/virial flags
	ev_init(eflag, vflag);
	// if atom count has changed, update qsum and qsqsum
	if (atom->natoms != natoms_original) {
		qsum_qsq();
		natoms_original = atom->natoms;
	}
	// return if there are no charges
	if (qsqsum == 0.0) return;

  //Set interfaces
  double **x = atom->x;
  double *q = atom->q;
  int nlocal = atom->nlocal;
  qqrd2e = force->qqrd2e;
  double **f = atom->f;
  double boxlo[3] = { domain->boxlo[0] ,domain->boxlo[1],domain->boxlo[2]};
  double boxhi[3] = { domain->boxhi[0] ,domain->boxhi[1],domain->boxhi[2]};

  if (RankID == 1)
  {
	  double X[nlocal][3], Q[nlocal], Force[nlocal][3], Pot[nlocal];
	  for (int i = 0; i < nlocal; i++)
	  {
		  X[i][0] = x[i][0] - (boxhi[0] - Lx / 2);
		  X[i][1] = x[i][1] - (boxhi[1] - Ly / 2);
		  X[i][2] = x[i][2] - (boxhi[2] - Lz / 2);
		  Q[i] = q[i];
		  Force[i][0] = 0.00;
		  Force[i][1] = 0.00;
		  Force[i][2] = 0.00;
		  Pot[i] = 0.00;
	  }
	  AdjustParticle_Double(X, nlocal, Lx, Ly, Lz);

	  //Find the image charge
	  int lx = ceil((Rs - Lx / 2) / Lx), ly = ceil((Rs - Ly / 2) / Ly), lz = ceil((Rs - Lz / 2) / Lz);
	  int TotalNumber = ceil(nlocal * (2 * lx + 1) * (2 * ly + 1) * (2 * lz + 1) / 2);
	  double ImageCharge[TotalNumber][4];
	  int ImageNumber;
	  SetImageCharge(ImageCharge, &ImageNumber, (2 * lx + 1) * (2 * ly + 1) * (2 * lz + 1), X, Q, nlocal, Rs, Lx, Ly, Lz, lx, ly, lz);

	  //Calculate near field potential which contains all the contributions in the cut-off sphere B at the monitoring points
	  double Near[Nw];
	  CalculateNearFieldAndZD_Single(Near, ImageCharge, ImageNumber, Nw, IF_FMM_RightTerm, Lx, Ly, Lz, PointSum, QuizSum, X, Force, Pot, nlocal, Q,tolerance);

	  //Solve The Least Square Problem
	  double C[p * p];
	  SolveLeastSquareProblem(C, A, Near, p, Nw);

	  //Compute final force, potential and energy
	  double Energy_HSMA;
	  Energy_HSMA = FinalCalculateEnergyAndForce_Single(Force, Pot, X, Q, nlocal, ImageCharge, ImageNumber, Fibonacci, (double**)QRD, (double**)QLocalRD, C, p, Fp, F, Rs, pi, IF_FMM_FinalPotential,tolerance);

	  scale = 1.0;
	  const double qscale = qqrd2e * scale;

	  //Assemble results
	  for (int i = 0; i < nlocal; i++)
	  {
		  f[i][0] += Force[i][0] * qscale;
		  f[i][1] += Force[i][1] * qscale;
		  f[i][2] += Force[i][2] * qscale;
	  }

	  if (eflag_global) {
		  energy += Energy_HSMA * qscale / 2;
	  }

	  if (evflag_atom) {
		  if (eflag_atom) {
			  for (int i = 0; i < nlocal; i++) {
				  eatom[i] += Pot[i];
				  eatom[i] *= qscale;
			  }
		  }
		  if (vflag_atom)
			  for (int i = 0; i < nlocal; i++)
				  for (int j = 0; j < 6; j++) vatom[i][j] *= q[i] * qscale;
	  }
  }
  else if (RankID > 1) {

	  double AllSource[maxatom][3], AllQ[maxatom];
	  int nlocal_All[RankID], nlocal_All_Q[RankID];
	  MPI_Allgather(&nlocal, 1, MPI_INT, nlocal_All, 1, MPI_INT, world);
	  int Size_All[RankID], Size_All_Q[RankID];
	  for (int i = 0; i < RankID; i++)
	  {
		  nlocal_All_Q[i] = nlocal_All[i];
		  nlocal_All[i] = nlocal_All[i] * 3;
		  if (i == 0)
		  {
			  Size_All[i] = 0;
			  Size_All_Q[i] = 0;
		  }
		  else
		  {
			  Size_All[i] = Size_All[i - 1] + nlocal_All[i - 1];
			  Size_All_Q[i] = Size_All_Q[i - 1] + nlocal_All_Q[i - 1];
		  }
	  }

	  MPI_Request request, request_Q; MPI_Status status;

	  double X[nlocal][3], Q[nlocal];
	  for (int i = 0; i < nlocal; i++)
	  {
		  X[i][0] = x[i][0] - (boxhi[0] - Lx / 2);
		  X[i][1] = x[i][1] - (boxhi[1] - Ly / 2);
		  X[i][2] = x[i][2] - (boxhi[2] - Lz / 2);
		  Q[i] = q[i];
	  }

	  AdjustParticle_Double(X, nlocal, Lx, Ly, Lz);

	  MPI_Iallgatherv((double*)X, nlocal * 3, MPI_DOUBLE, (double*)AllSource, nlocal_All, Size_All, MPI_DOUBLE, world, &request);
	  MPI_Iallgatherv((double*)Q, nlocal, MPI_DOUBLE, (double*)AllQ, nlocal_All_Q, Size_All_Q, MPI_DOUBLE, world, &request_Q);

	  //Find the image charge
	  int lx = ceil((Rs - Lx / 2) / Lx), ly = ceil((Rs - Ly / 2) / Ly), lz = ceil((Rs - Lz / 2) / Lz);
	  int TotalNumber = ceil(nlocal * (2 * lx + 1) * (2 * ly + 1) * (2 * lz + 1) / 2);
	  double ImageCharge[TotalNumber][4];
	  int ImageNumber;
	  SetImageCharge(ImageCharge, &ImageNumber, (2 * lx + 1) * (2 * ly + 1) * (2 * lz + 1), X, Q, nlocal, Rs, Lx, Ly, Lz, lx, ly, lz);

	  //Calculate near field potential which contains all the contributions in the cut-off sphere B at the monitoring points
	  double Near[Nw], Near_All[Nw];

	  //CalculateNearFieldAndZD(Near, ImageCharge, ImageNumber, Nw, IF_FMM_RightTerm, Lx, Ly, Lz, PointSum_New, QuizSum_New);
	  CalculateNearFieldAndZD(Near, ImageCharge, ImageNumber, Nw, IF_FMM_RightTerm, Lx, Ly, Lz, PointSum, QuizSum,tolerance);
	  MPI_Reduce(Near, Near_All, Nw, MPI_DOUBLE, MPI_SUM, 0, world);

	  //Solve The Least Square Problem
	  double C[p * p];
	  if (me == 0) { SolveLeastSquareProblem(C, A, Near_All, p, Nw); }
	  MPI_Bcast(C, p * p, MPI_DOUBLE, 0, world);
	  MPI_Wait(&request, &status);
	  MPI_Wait(&request_Q, &status);

	  //Compute final force, potential and energy
	  double Energy_HSMA;
	  double Force[nlocal][3], Pot[nlocal];
	  TotalNumber = ceil(maxatom * (2 * lx + 1) * (2 * ly + 1) * (2 * lz + 1) / 2);//¼õÉÙ´æ´¢Á¿ÏûºÄ
	  double ImageCharge_All[TotalNumber][4];
	  SetImageCharge(ImageCharge_All, &ImageNumber, (2 * lx + 1) * (2 * ly + 1) * (2 * lz + 1), AllSource, AllQ, maxatom, Rs, Lx, Ly, Lz, lx, ly, lz);

	  Energy_HSMA = FinalCalculateEnergyAndForce(Force, Pot, X, Q, nlocal, ImageCharge_All, ImageNumber, Fibonacci, (double**)QRD, (double**)QLocalRD, C, p, Fp, F, Rs, pi, IF_FMM_FinalPotential,tolerance);

	  scale = 1.0;
	  const double qscale = qqrd2e * scale;

	  //Assemble results
	  for (int i = 0; i < nlocal; i++)
	  {
		  f[i][0] += Force[i][0] * qscale;
		  f[i][1] += Force[i][1] * qscale;
		  f[i][2] += Force[i][2] * qscale;
	  }

	  double Energy_HSMA_All;
	  MPI_Allreduce(&Energy_HSMA, &Energy_HSMA_All, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	  if (eflag_global) {
		  energy += Energy_HSMA_All * qscale / 2;
	  }

	  if (evflag_atom) {
		  if (eflag_atom) {
			  for (int i = 0; i < nlocal; i++) {
				  eatom[i] += Pot[i];
				  eatom[i] *= qscale;
			  }
		  }
		  if (vflag_atom)
			  for (int i = 0; i < nlocal; i++)
				  for (int j = 0; j < 6; j++) vatom[i][j] *= q[i] * qscale;
	  }
  }
}

double HSMA3D::memory_usage()
{
  double bytes = 0.0;
  bytes += maxatom * sizeof(double);
  bytes += 3*maxatom * sizeof(double);
  return bytes;
}

double fac(double t)//calculate factorial
{
	return factorial(int(t));
}

void HSMA3D::CalculateRDMultipoleExpansion(double* Q, int p, double x, double y, double z)
{
	Q[0] = 1.0;
	Q[1] = y / 2;
	Q[2] = (-z);
	Q[3] = (-x / 2);
	int t = 4;
	int m, n;
	for (int i = 2; i < p; i++)
	{
		while (t < (i + 1) * (i + 1))
		{
			m = t - i - i * i;
			n = i;
			if (m == -n)
			{
				Q[t] = (y * Q[n * n - 1] - x * Q[n * n - 2 * n + 1]) / (2 * n + 0.00);
			}
			else if (m == n)
			{
				Q[t] = (-(x * Q[n * n - 1] + y * Q[n * n - 2 * n + 1]) / (2 * n + 0.00));
			}
			else if (n - abs(m) == 1)
			{
				Q[t] = (-z) * Q[n * n - n + m];
			}
			else if ((n - abs(m)) > 1)
			{
				Q[t] = (-((2 * n - 1.0) * z * Q[n * n - n + m] + (x * x + y * y + z * z) * Q[n * n - 3 * n + m + 2]) / ((n - abs(m) + 0.0) * (n + abs(m) + 0.0)));
			}
			t++;
		}
	}

	double rr = sqrt(x * x + y * y + z * z);
	t = 0;
	for (int i = 0; i < p; i++)
	{
		while (t < (i + 1) * (i + 1))
		{
			m = t - i - i * i;
			n = i;
			Q[t] = Q[t] * (n + 0.00) / rr;
			t++;
		}
	}

	t = 0;//normlization     Please do not normlize after every step!! That's wrong!
	for (int i = 0; i < p; i++)
	{
		while (t < (i + 1) * (i + 1))
		{
			m = t - i - i * i;
			n = i;
			Q[t] = Q[t] * sqrt((fac(n - m)) * fac(n + m));
			t++;
		}
	}
}

void HSMA3D::CalculateLocalRDMultipoleExpansion(double* Q, int p, double x, double y, double z, double Rs)
{
	Q[0] = 1.0;
	Q[1] = y / 2;
	Q[2] = (-z);
	Q[3] = (-x / 2);
	int t = 4;
	int m, n;
	for (int i = 2; i < p; i++)
	{
		while (t < (i + 1) * (i + 1))
		{
			m = t - i - i * i;
			n = i;
			if (m == -n)
			{
				Q[t] = (y * Q[n * n - 1] - x * Q[n * n - 2 * n + 1]) / (2 * n + 0.00);
			}
			else if (m == n)
			{
				Q[t] = (-(x * Q[n * n - 1] + y * Q[n * n - 2 * n + 1]) / (2 * n + 0.00));
			}
			else if (n - abs(m) == 1)
			{
				Q[t] = (-z) * Q[n * n - n + m];
			}
			else if ((n - abs(m)) > 1)
			{
				Q[t] = (-((2 * n - 1.0) * z * Q[n * n - n + m] + (x * x + y * y + z * z) * Q[n * n - 3 * n + m + 2]) / ((n - abs(m) + 0.0) * (n + abs(m) + 0.0)));
			}
			t++;
		}
	}

	double rr = sqrt(x * x + y * y + z * z);
	t = 0;
	for (int i = 0; i < p; i++)
	{
		while (t < (i + 1) * (i + 1))
		{
			m = t - i - i * i;
			n = i;
			Q[t] = Q[t] * (-n - 1) * pow(Rs, 2 * n + 1.0) / pow(rr, 2 * n + 2.0);
			t++;
		}
	}

	t = 0;//normlization     Please do not normlize after every step!! That's wrong!
	for (int i = 0; i < p; i++)
	{
		while (t < (i + 1) * (i + 1))
		{
			m = t - i - i * i;
			n = i;
			Q[t] = Q[t] * sqrt((fac(n - m)) * fac(n + m));
			t++;
		}
	}
}

void HSMA3D::CalculateMultipoleExpansion(double* Q, int p, double x, double y, double z)
{
	Q[0] = 1.0;
	Q[1] = y / 2;
	Q[2] = (-z);
	Q[3] = (-x / 2);
	int t = 4;
	int m, n;
	for (int i = 2; i < p; i++)
	{
		while (t < (i + 1) * (i + 1))
		{
			m = t - i - i * i;
			n = i;
			if (m == -n)
			{
				Q[t] = (y * Q[n * n - 1] - x * Q[n * n - 2 * n + 1]) / (2 * n + 0.00);
			}
			else if (m == n)
			{
				Q[t] = (-(x * Q[n * n - 1] + y * Q[n * n - 2 * n + 1]) / (2 * n + 0.00));
			}
			else if (n - abs(m) == 1)
			{
				Q[t] = (-z) * Q[n * n - n + m];
			}
			else if ((n - abs(m)) > 1)
			{
				Q[t] = (-((2 * n - 1.0) * z * Q[n * n - n + m] + (x * x + y * y + z * z) * Q[n * n - 3 * n + m + 2]) / ((n - abs(m) + 0.0) * (n + abs(m) + 0.0)));
			}
			t++;
		}
	}


	t = 0;//normlization Please do not normlize after every step!! That's wrong!
	for (int i = 0; i < p; i++)
	{
		while (t < (i + 1) * (i + 1))
		{
			m = t - i - i * i;
			n = i;
			Q[t] = Q[t] * sqrt((fac(n - m)) * fac(n + m));
			t++;
		}
	}
}

void HSMA3D::CalculateZDerivativeMultipoleExpansion(double* Q, int p, double x, double y, double z)//Q is a p*p vector.  Set multi-pole expansion coefficient
{
	double Qold[p * p];
	Qold[0] = 1.0;
	Qold[1] = y / 2;
	Qold[2] = (-1.0) * z;
	Qold[3] = (-1.0) * x / 2.0;
	int t = 4;
	int m, n;
	for (int i = 2; i < p; i++)
	{
		while (t < (i + 1) * (i + 1))
		{
			m = t - i - i * i;
			n = i;
			if (m == -n)
			{
				Qold[t] = (y * Qold[n * n - 1] - x * Qold[n * n - 2 * n + 1]) / (2 * n + 0.00);
			}
			else if (m == n)
			{
				Qold[t] = (-(x * Qold[n * n - 1] + y * Qold[n * n - 2 * n + 1]) / (2 * n + 0.00));
			}
			else if (n - abs(m) == 1)
			{
				Qold[t] = (-z) * Qold[n * n - n + m];
			}
			else if ((n - abs(m)) > 1)
			{
				Qold[t] = (-((2 * n - 1.0) * z * Qold[n * n - n + m] + (x * x + y * y + z * z) * Qold[n * n - 3 * n + m + 2]) / ((n - abs(m) + 0.0) * (n + abs(m) + 0.0)));
			}
			t++;
		}
	}

	Q[0] = 0.00;
	Q[1] = 0.00;
	Q[2] = -1.0;
	Q[3] = 0.00;
	t = 4;
	for (int i = 2; i < p; i++)
	{
		while (t < (i + 1) * (i + 1))
		{
			m = t - i - i * i;
			n = i;
			if (m == -n)
			{
				Q[t] = (y * Q[n * n - 1] - x * Q[n * n - 2 * n + 1]) / (2 * n + 0.00);
			}
			else if (m == n)
			{
				Q[t] = ((-1.0) * (x * Q[n * n - 1] + y * Q[n * n - 2 * n + 1]) / (2 * n + 0.00));
			}
			else if (n - abs(m) == 1)
			{
				Q[t] = (-1.0) * z * Q[n * n - n + m] - Qold[n * n - n + m];
			}
			else if ((n - abs(m)) > 1)
			{
				Q[t] = (-((2 * n - 1.0) * z * Q[n * n - n + m] + (2 * n - 1.0) * Qold[n * n - n + m] + (x * x + y * y + z * z) * Q[n * n - 3 * n + m + 2] + (2.0 * z) * Qold[n * n - 3 * n + m + 2]) / ((n - abs(m) + 0.0) * (n + abs(m) + 0.0)));
			}
			t++;
		}
	}

	t = 0;//normlization     Please do not normlize after every step!! That's wrong!
	for (int i = 0; i < p; i++)
	{
		while (t < (i + 1) * (i + 1))
		{
			m = t - i - i * i;
			n = i;
			Q[t] = Q[t] * sqrt((fac(n - m)) * fac(n + m));
			t++;
		}
	}

}

void HSMA3D::CalculateXDMultipoleExpansion(double* Q, int p, double x, double y, double z)//Q is a p*p vector.  Set multi-pole expansion coefficient
{
	double Qold[p * p];
	Qold[0] = 1.0;
	Qold[1] = y / 2;
	Qold[2] = (-1.0) * z;
	Qold[3] = (-1.0) * x / 2.0;
	int t = 4;
	int m, n;
	for (int i = 2; i < p; i++)
	{
		while (t < (i + 1) * (i + 1))
		{
			m = t - i - i * i;
			n = i;
			if (m == -n)
			{
				Qold[t] = (y * Qold[n * n - 1] - x * Qold[n * n - 2 * n + 1]) / (2 * n + 0.00);
			}
			else if (m == n)
			{
				Qold[t] = (-(x * Qold[n * n - 1] + y * Qold[n * n - 2 * n + 1]) / (2 * n + 0.00));
			}
			else if (n - abs(m) == 1)
			{
				Qold[t] = (-z) * Qold[n * n - n + m];
			}
			else if ((n - abs(m)) > 1)
			{
				Qold[t] = (-((2 * n - 1.0) * z * Qold[n * n - n + m] + (x * x + y * y + z * z) * Qold[n * n - 3 * n + m + 2]) / ((n - abs(m) + 0.0) * (n + abs(m) + 0.0)));
			}
			t++;
		}
	}
	Q[0] = 0.00;
	Q[1] = 0.00;
	Q[2] = 0.0;
	Q[3] = -1.0 / 2.0;
	t = 4;
	for (int i = 2; i < p; i++)
	{
		while (t < (i + 1) * (i + 1))
		{
			m = t - i - i * i;
			n = i;
			if (m == -n)
			{
				Q[t] = (y * Q[n * n - 1] - x * Q[n * n - 2 * n + 1] - Qold[n * n - 2 * n + 1]) / (2 * n + 0.00);
			}
			else if (m == n)
			{
				Q[t] = ((-1.0) * (x * Q[n * n - 1] + y * Q[n * n - 2 * n + 1] + Qold[n * n - 1]) / (2 * n + 0.00));
			}
			else if (n - abs(m) == 1)
			{
				Q[t] = (-1.0) * z * Q[n * n - n + m];
			}
			else if ((n - abs(m)) > 1)
			{
				Q[t] = (-((2 * n - 1.0) * z * Q[n * n - n + m] + (x * x + y * y + z * z) * Q[n * n - 3 * n + m + 2] + (2.0 * x) * Qold[n * n - 3 * n + m + 2]) / ((n - abs(m) + 0.0) * (n + abs(m) + 0.0)));
				//mm++;
			}
			t++;
		}
	}

	t = 0;//normlization     Please do not normlize after every step!! That's wrong!
	for (int i = 0; i < p; i++)
	{
		while (t < (i + 1) * (i + 1))
		{
			m = t - i - i * i;
			n = i;
			Q[t] = Q[t] * sqrt((fac(n - m)) * fac(n + m));
			t++;
		}
	}
}

void HSMA3D::CalculateYDMultipoleExpansion(double* Q, int p, double x, double y, double z)//Q is a p*p vector.  Set multi-pole expansion coefficient
{
	double Qold[p * p];
	Qold[0] = 1.0;
	Qold[1] = y / 2;
	Qold[2] = (-1.0) * z;
	Qold[3] = (-1.0) * x / 2.0;
	int t = 4;
	int m, n;
	for (int i = 2; i < p; i++)
	{
		while (t < (i + 1) * (i + 1))
		{
			m = t - i - i * i;
			n = i;
			if (m == -n)
			{
				Qold[t] = (y * Qold[n * n - 1] - x * Qold[n * n - 2 * n + 1]) / (2 * n + 0.00);
			}
			else if (m == n)
			{
				Qold[t] = (-(x * Qold[n * n - 1] + y * Qold[n * n - 2 * n + 1]) / (2 * n + 0.00));
			}
			else if (n - abs(m) == 1)
			{
				Qold[t] = (-z) * Qold[n * n - n + m];
			}
			else if ((n - abs(m)) > 1)
			{
				Qold[t] = (-((2 * n - 1.0) * z * Qold[n * n - n + m] + (x * x + y * y + z * z) * Qold[n * n - 3 * n + m + 2]) / ((n - abs(m) + 0.0) * (n + abs(m) + 0.0)));
			}
			t++;
		}
	}

	Q[0] = 0.00;
	Q[1] = 1.0 / 2.0;
	Q[2] = 0.00;
	Q[3] = 0.00;
	t = 4;
	for (int i = 2; i < p; i++)
	{
		while (t < (i + 1) * (i + 1))
		{
			m = t - i - i * i;
			n = i;
			if (m == -n)
			{
				Q[t] = (y * Q[n * n - 1] - x * Q[n * n - 2 * n + 1] + Qold[n * n - 1]) / (2 * n + 0.00);
			}
			else if (m == n)
			{
				Q[t] = ((-1.0) * (x * Q[n * n - 1] + y * Q[n * n - 2 * n + 1] + Qold[n * n - 2 * n + 1]) / (2 * n + 0.00));
			}
			else if (n - abs(m) == 1)
			{
				Q[t] = (-1.0) * z * Q[n * n - n + m];
			}
			else if ((n - abs(m)) > 1)
			{
				Q[t] = (-((2 * n - 1.0) * z * Q[n * n - n + m] + (x * x + y * y + z * z) * Q[n * n - 3 * n + m + 2] + (2.0 * y) * Qold[n * n - 3 * n + m + 2]) / ((n - abs(m) + 0.0) * (n + abs(m) + 0.0)));
				//mm++;
			}
			t++;
		}
	}

	t = 0;//normlization     Please do not normlize after every step!! That's wrong!
	for (int i = 0; i < p; i++)
	{
		while (t < (i + 1) * (i + 1))
		{
			m = t - i - i * i;
			n = i;
			Q[t] = Q[t] * sqrt((fac(n - m)) * fac(n + m));
			t++;
		}
	}
}

void HSMA3D::SetFibonacci(double Fibonacci[][4], double F, double Fp, int Np, double Rs, double PI)
{
	double zz, zt, zs, phia, phib;
	double deltaz = 2.0 / F;
	for (int i = 0; i <= int(F); i++)
	{
		zz = -1.0 + i * deltaz;
		zt = zz + sin(PI * zz) / PI;
		zs = sqrt(1 - zt * zt);
		phia = PI * (i + 0.00) * Fp / F;
		phib = PI + phia;

		Fibonacci[2 * i][0] = cos(phia) * zs * (Rs + 0.0);
		Fibonacci[2 * i + 1][0] = cos(phib) * zs * (Rs + 0.00);
		Fibonacci[2 * i][1] = sin(phia) * zs * (Rs + 0.00);
		Fibonacci[2 * i + 1][1] = sin(phib) * zs * (Rs + 0.00);
		Fibonacci[2 * i][2] = zt * (Rs + 0.00);
		Fibonacci[2 * i + 1][2] = zt * (Rs + 0.00);
		Fibonacci[2 * i][3] = PI * deltaz * (1.0 + cos(PI * zz)) * Rs * Rs;
		Fibonacci[2 * i + 1][3] = PI * deltaz * (1.0 + cos(PI * zz)) * Rs * Rs;
	}
}

void HSMA3D::AdjustParticle_Float(float Particle[][3], int N, float Lx, float Ly, float Lz)//Image Correction 0~Lx
{
	for (int i = 0; i < N; i++)
	{
		if (Particle[i][0] > Lx)
		{
			Particle[i][0] = Particle[i][0] - ceil((Particle[i][0] - Lx) / (Lx + 0.00)) * Lx;
		}
		if (Particle[i][0] < 0.00)
		{
			Particle[i][0] = Particle[i][0] + ceil(fabs(Particle[i][0] / Lx)) * Lx;
		}
		if (Particle[i][1] > Ly)
		{
			Particle[i][1] = Particle[i][1] - ceil((Particle[i][1] - Ly) / (Ly + 0.00)) * Ly;
		}
		if (Particle[i][1] < 0.00)
		{
			Particle[i][1] = Particle[i][1] + ceil(fabs(Particle[i][1] / Ly)) * Ly;
		}
		if (Particle[i][2] > Lz)
		{
			Particle[i][2] = Particle[i][2] - ceil((Particle[i][2] - Lz) / (Lz + 0.00)) * Lz;
		}
		if (Particle[i][2] < 0.00)
		{
			Particle[i][2] = Particle[i][2] + ceil(fabs(Particle[i][2] / Lz)) * Lz;
		}
	}
}

void HSMA3D::AdjustParticle_Double(double Particle[][3], int N, double Lx, double Ly, double Lz)//Image Correction -Lx/2~Lx/2
{
	for (int i = 0; i < N; i++)
	{
		if (Particle[i][0] > Lx/2)
		{
			Particle[i][0] = Particle[i][0] - ceil((Particle[i][0] - Lx/2.0) / (Lx + 0.00)) * Lx;
		}
		if (Particle[i][0] < -Lx / 2)
		{
			Particle[i][0] = Particle[i][0] + ceil(fabs((Particle[i][0] + Lx / 2.0)/Lx)) * Lx;
		}
		if (Particle[i][1] > Ly/2)
		{
			Particle[i][1] = Particle[i][1] - ceil((Particle[i][1] - Ly/2.0) / (Ly + 0.00)) * Ly;
		}
		if (Particle[i][1] < -Ly / 2)
		{
			Particle[i][1] = Particle[i][1] + ceil(fabs( (Particle[i][1] + Ly / 2.0) / Ly)) * Ly;
		}
		if (Particle[i][2] > Lz/2)
		{
			Particle[i][2] = Particle[i][2] - ceil((Particle[i][2] - Lz/2.0) / (Lz + 0.00)) * Lz;
		}
		if (Particle[i][2] < -Lz / 2)
		{
			Particle[i][2] = Particle[i][2] + ceil(fabs((Particle[i][2] + Lz / 2.0)/Lz)) * Lz;
		}
	}
}

void HSMA3D::SetImageCharge(double ImageCharge[][4], int* ImageNumber, int TotalNumber, double Source[][3], double* Q, int NSource, double Rs, double Lx, double Ly, double Lz, int lx, int ly, int lz)
{
	int total = 0;
	int number = 0;

	for (int i = -lx; i <= lx; i++)
	{
		double CX, CY, CZ;
		for (int j = -ly; j <= ly; j++)
			for (int k = -lz; k <= lz; k++)
			{

				for (int m = 0; m < NSource; m++)
				{

					CX = Source[m][0] + i * Lx;
					CY = Source[m][1] + j * Ly;
					CZ = Source[m][2] + k * Lz;
					if (CX * CX + CY * CY + CZ * CZ <= Rs * Rs)
					{
						ImageCharge[number][0] = CX;
						ImageCharge[number][1] = CY;
						ImageCharge[number][2] = CZ;
						ImageCharge[number][3] = Q[m];
						number++;
					}
				}

			}
	}
	*ImageNumber = number;
}

void HSMA3D::CalculateNearFieldAndZD(double* Near, double ImageCharge[][4], int ImageNumber, int Nw, int IF_FMM_RightTerm, double Lx, double Ly, double Lz, double **PointSum, double **QuizSum, double tolerance)
{
	if (IF_FMM_RightTerm)//Using FMM to calculate pairwise sum
	{
		double eps = tolerance;

		/*            Set FMM parameters           */
		int ns = ImageNumber;
		int nt = Nw * 2;
		double* source = (double*)malloc(3 * ns * sizeof(double));
		double* target = (double*)malloc(3 * nt * sizeof(double));
		double* charge = (double*)malloc(ns * sizeof(double));

		for (int i = 0; i < ns; i++)
		{
			source[3 * i] = ImageCharge[i][0];
			source[3 * i + 1] = ImageCharge[i][1];
			source[3 * i + 2] = ImageCharge[i][2];
			charge[i] = ImageCharge[i][3];
		}

		for (int i = 0; i < Nw; i++)
		{
			target[3 * i] = PointSum[i][0];
			target[3 * i + 1] = PointSum[i][1];
			target[3 * i + 2] = PointSum[i][2];

			target[3 * Nw + 3 * i] = QuizSum[i][0];
			target[3 * Nw + 3 * i + 1] = QuizSum[i][1];
			target[3 * Nw + 3 * i + 2] = QuizSum[i][2];
		}

		double* pottarg = (double*)malloc(nt * sizeof(double));
		double* gradtarg = (double*)malloc(3 * nt * sizeof(double));

		int ier = 0;
		lfmm3d_t_c_g_(&eps, &ns, source, charge, &nt, target, pottarg, gradtarg, &ier);
		for (int i = 0; i < Nw; i++)
		{
			Near[i] = pottarg[Nw + i] - pottarg[i];
		}

		free(source); free(target); free(charge); free(pottarg);
		source = NULL; target = NULL; charge = NULL; pottarg = NULL;
	}
	else
	{
		if (tolerance > 0.000001)
		{
			float Paramet1[int(ceil(ImageNumber / 16.0)) * 16];
			float Image_X[int(ceil(ImageNumber / 16.0)) * 16], Image_Y[int(ceil(ImageNumber / 16.0)) * 16], Image_Z[int(ceil(ImageNumber / 16.0)) * 16];
			for (int i = 0; i < int(ceil(ImageNumber / 16.0)) * 16; i++)
			{
				if (i < ImageNumber) {
					Image_X[i] = ImageCharge[i][0];
					Image_Y[i] = ImageCharge[i][1];
					Image_Z[i] = ImageCharge[i][2];
					Paramet1[i] = ImageCharge[i][3];
				}
				else {
					Image_X[i] = 0.00;
					Image_Y[i] = 0.00;
					Image_Z[i] = 0.00;
					Paramet1[i] = 0.00;
				}
			}

            #pragma omp parallel
			{
				int id = omp_get_thread_num();
				int size = comm->nthreads;

				int min_atom = id * floor( Nw / (size+0.00)) + 1, max_atom = (id + 1) * floor(Nw / (size+0.00));
				if (id == size - 1)max_atom = Nw - 1;
				if (id == 0)min_atom = 0;

				for (int i = min_atom; i <= max_atom; i++)
				{
					__m512 XP, YP, ZP, XQ, YQ, ZQ, X1, Y1, Z1, Q1, dx, dy, dz, dx1, dy1, dz1, square, square1;
					XP = _mm512_set1_ps(PointSum[i][0]);
					YP = _mm512_set1_ps(PointSum[i][1]);
					ZP = _mm512_set1_ps(PointSum[i][2]);
					XQ = _mm512_set1_ps(QuizSum[i][0]);
					YQ = _mm512_set1_ps(QuizSum[i][1]);
					ZQ = _mm512_set1_ps(QuizSum[i][2]);
					float pottarg = 0.00, pottarg1 = 0.00;
					for (int j = 0; j < ImageNumber; j=j+16)
					{
						X1 = _mm512_load_ps(&Image_X[j]);
						Y1 = _mm512_load_ps(&Image_Y[j]);
						Z1 = _mm512_load_ps(&Image_Z[j]);
						Q1 = _mm512_load_ps(&Paramet1[j]);

						dx = XP - X1;
						dy = YP - Y1;
						dz = ZP - Z1;
						square = dx * dx + dy * dy + dz * dz;
						pottarg += _mm512_reduce_add_ps(Q1 * _mm512_invsqrt_ps(square));

						dx1 = XQ - X1;
						dy1 = YQ - Y1;
						dz1 = ZQ - Z1;
						square1 = dx1 * dx1 + dy1 * dy1 + dz1 * dz1;
						pottarg1 += _mm512_reduce_add_ps(Q1 * _mm512_invsqrt_ps(square1));
					}
					Near[i] = pottarg1 - pottarg;
				}
			}
		}
		else
		{
			double Paramet1[int(ceil(ImageNumber / 8.0)) * 8];
			double Image_X[int(ceil(ImageNumber / 8.0)) * 8], Image_Y[int(ceil(ImageNumber / 8.0)) * 8], Image_Z[int(ceil(ImageNumber / 8.0)) * 8];
			for (int i = 0; i < int(ceil(ImageNumber / 8.0)) * 8; i++)
			{
				if (i < ImageNumber) {
					Image_X[i] = ImageCharge[i][0];
					Image_Y[i] = ImageCharge[i][1];
					Image_Z[i] = ImageCharge[i][2];
					Paramet1[i] = ImageCharge[i][3];
				}
				else {
					Image_X[i] = 0.00;
					Image_Y[i] = 0.00;
					Image_Z[i] = 0.00;
					Paramet1[i] = 0.00;
				}
			}

			#pragma omp parallel
			{
				int id = omp_get_thread_num();
				int size = comm->nthreads;

				int min_atom = id * floor(Nw / size) + 1, max_atom = (id + 1) * floor(Nw / size);
				if (id == size - 1)max_atom = Nw - 1;
				if (id == 0)min_atom = 0;

				for (int i = min_atom; i <= max_atom; i++)
				{
					__m512d XP, YP, ZP, XQ, YQ, ZQ, X1, Y1, Z1, Q1, dx, dy, dz, dx1, dy1, dz1, square, square1;
					XP = _mm512_set1_pd(PointSum[i][0]);
					YP = _mm512_set1_pd(PointSum[i][1]);
					ZP = _mm512_set1_pd(PointSum[i][2]);
					XQ = _mm512_set1_pd(QuizSum[i][0]);
					YQ = _mm512_set1_pd(QuizSum[i][1]);
					ZQ = _mm512_set1_pd(QuizSum[i][2]);
					double pottarg = 0.00, pottarg1 = 0.00;
					for (int j = 0; j < ImageNumber; j = j + 8)
					{
						X1 = _mm512_load_pd(&Image_X[j]);
						Y1 = _mm512_load_pd(&Image_Y[j]);
						Z1 = _mm512_load_pd(&Image_Z[j]);
						Q1 = _mm512_load_pd(&Paramet1[j]);

						dx = XP - X1;
						dy = YP - Y1;
						dz = ZP - Z1;
						square = dx * dx + dy * dy + dz * dz;
						pottarg += _mm512_reduce_add_pd(Q1 * _mm512_invsqrt_pd(square));

						dx1 = XQ - X1;
						dy1 = YQ - Y1;
						dz1 = ZQ - Z1;
						square1 = dx1 * dx1 + dy1 * dy1 + dz1 * dz1;
						pottarg1 += _mm512_reduce_add_pd(Q1 * _mm512_invsqrt_pd(square1));
					}
					Near[i] = pottarg1 - pottarg;
				}
			}
		}
	}
}

void HSMA3D::SolveLeastSquareProblem(double* C, double** A, double* Near, int p, int Nw)
{
	int rowAT, columnATA, columnAT;
	double alpha, beta;
	rowAT = p * p - 1; columnAT = Nw; columnATA = p * p - 1;
	double MatrixAT[rowAT * columnAT], MatrixATA[rowAT * columnATA], BB[Nw], ATB[p * p - 1], INV_ATA_ATB[p * p - 1];
	alpha = 1.0; beta = 0.00;
	for (int i = 0; i < rowAT; i++)
		for (int j = 0; j < columnAT; j++)
		{
			MatrixAT[i * columnAT + j] = A[j][i];
		}

	for (int i = 0; i < (rowAT * columnATA); i++) {
		MatrixATA[i] = 0.0;
	}
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, rowAT, columnATA, columnAT, alpha, MatrixAT, columnAT, MatrixAT, columnAT, beta, MatrixATA, columnATA);
	int InfoHelp;
	int VectorHelp[columnATA];
	for (int i = 0; i < columnATA; i++)
		VectorHelp[i] = 0;
	InfoHelp = LAPACKE_dgetrf(CblasRowMajor, columnATA, columnATA, MatrixATA, columnATA, VectorHelp);
	InfoHelp = LAPACKE_dgetri(CblasRowMajor, columnATA, MatrixATA, columnATA, VectorHelp);

	for (int i = 0; i < columnAT; i++)
	{
		BB[i] = Near[i];
	}
	for (int i = 0; i < columnATA; i++)
	{
		ATB[i] = 0.00;
		INV_ATA_ATB[i] = 0.00;
	}

	cblas_dgemv(CblasRowMajor, CblasNoTrans, rowAT, columnAT, alpha, MatrixAT, columnAT, BB, 1, beta, ATB, 1);
	cblas_dgemv(CblasRowMajor, CblasNoTrans, columnATA, columnATA, alpha, MatrixATA, columnATA, ATB, 1, beta, INV_ATA_ATB, 1);


	C[0] = 0.00;
	for (int i = 0; i < rowAT; i++)
	{
		C[i + 1] = INV_ATA_ATB[i];
	}
}

double HSMA3D::FinalCalculateEnergyAndForce(double Force[][3], double* Pot, double Source[][3], double* Q, int NSource, double ImageCharge[][4], int ImageNumber, double **Fibonacci, double** QRD, double** QLocalRD, double* C, int p, double Fp, double F, double Rs, double PI, int IF_FMM_FinalPotential, double tolerance)
{
	if (!IF_FMM_FinalPotential)
	{
		if (tolerance > 0.000001)
		{
			float EF[NSource], EFX[NSource], EFY[NSource], EFZ[NSource];
			float EN[NSource], ENX[NSource], ENY[NSource], ENZ[NSource];
			float Q_Image[int(ceil(ImageNumber / 16.0)) * 16];
			for (int j = 0; j < ImageNumber; j++)
			{
				Q_Image[j] = ImageCharge[j][3];
			}
			for (int j = ImageNumber; j<int(ceil(ImageNumber / 16.0)) * 16; j++)
			{
				Q_Image[j] = 0.00;
			}

			double C_New[int(ceil(p * p / 8.0)) * 8];
			for (int i = 0; i<int(ceil(p * p / 8.0)) * 8; i++)
			{
				if (i < p * p)
				{
					C_New[i] = C[i];
				}
				else
				{
					C_New[i] = 0.00;
				}
			}

			float Image_X[int(ceil(ImageNumber / 16.0)) * 16], Image_Y[int(ceil(ImageNumber / 16.0)) * 16], Image_Z[int(ceil(ImageNumber / 16.0)) * 16];
			for (int i = 0; i < int(ceil(ImageNumber / 16.0)) * 16; i++)
			{
				if (i < ImageNumber) {
					Image_X[i] = ImageCharge[i][0];
					Image_Y[i] = ImageCharge[i][1];
					Image_Z[i] = ImageCharge[i][2];
				}
				else {
					Image_X[i] = 0.00;
					Image_Y[i] = 0.00;
					Image_Z[i] = 0.00;
				}
			}

            #pragma omp parallel
			{
				int id = omp_get_thread_num();
				int size = comm->nthreads;

				int min_atom = id * floor(NSource / size) + 1, max_atom = (id + 1) * floor(NSource / size);
				if (id == size - 1)max_atom = NSource - 1;
				if (id == 0)min_atom = 0;

				for (int i = min_atom; i <= max_atom; i++)
				{
					double QF[int(ceil(p * p / 8.0)) * 8], QFX[int(ceil(p * p / 8.0)) * 8], QFY[int(ceil(p * p / 8.0)) * 8], QFZ[int(ceil(p * p / 8.0)) * 8];
					CalculateMultipoleExpansion(QF, p, Source[i][0], Source[i][1], Source[i][2]);
					CalculateZDerivativeMultipoleExpansion(QFZ, p, Source[i][0], Source[i][1], Source[i][2]);
					CalculateXDMultipoleExpansion(QFX, p, Source[i][0], Source[i][1], Source[i][2]);
					CalculateYDMultipoleExpansion(QFY, p, Source[i][0], Source[i][1], Source[i][2]);
					EF[i] = 0.00; EFX[i] = 0.00; EFY[i] = 0.00; EFZ[i] = 0.00;
					EN[i] = 0.00; ENX[i] = 0.00; ENY[i] = 0.00; ENZ[i] = 0.00;
					for (int ii = p * p; ii<int(ceil(p * p / 8.0)) * 8; ii++)
					{
						QF[ii] = 0.00;
						QFX[ii] = 0.00;
						QFY[ii] = 0.00;
						QFZ[ii] = 0.00;
					}
					__m512d qf, qfz, qfx, qfy, c;
					for (int j = 0; j < p * p; j = j + 8)
					{
						qf = _mm512_load_pd(&QF[j]);
						qfz = _mm512_load_pd(&QFZ[j]);
						qfx = _mm512_load_pd(&QFX[j]);
						qfy = _mm512_load_pd(&QFY[j]);
						c = _mm512_load_pd(&C_New[j]);
						EF[i] = EF[i] + _mm512_reduce_add_pd(qf * c);
						EFX[i] = EFX[i] + _mm512_reduce_add_pd(qfx * c);
						EFY[i] = EFY[i] + _mm512_reduce_add_pd(qfy * c);
						EFZ[i] = EFZ[i] + _mm512_reduce_add_pd(qfz * c);
					}
					__m512 X0, Y0, Z0, X1, Y1, Z1, Q1, deltax, deltay, deltaz, delta, square, judge, Zero, midterm, delta_square;
					__mmask16 k0;
					X0 = _mm512_set1_ps(Source[i][0]);
					Y0 = _mm512_set1_ps(Source[i][1]);
					Z0 = _mm512_set1_ps(Source[i][2]);
					judge = _mm512_set1_ps(0.000000000001);
					Zero = _mm512_set1_ps(0.00);
					for (int j = 0; j < ImageNumber; j = j + 16)
					{
						X1 = _mm512_load_ps(&Image_X[j]);
						Y1 = _mm512_load_ps(&Image_Y[j]);
						Z1 = _mm512_load_ps(&Image_Z[j]);
						Q1 = _mm512_load_ps(&Q_Image[j]);
						deltax = X1 - X0;
						deltay = Y1 - Y0;
						deltaz = Z1 - Z0;
						square = deltax * deltax + deltay * deltay + deltaz * deltaz;
						k0 = _mm512_cmp_ps_mask(square, judge, _MM_CMPINT_GT);
						delta = _mm512_mask_invsqrt_ps(Zero, k0, square);
						midterm = Q1 * delta;
						delta_square = delta * delta;
						EN[i] += _mm512_reduce_add_ps(midterm);
						ENX[i] += _mm512_reduce_add_ps(midterm * delta_square * deltax);
						ENY[i] += _mm512_reduce_add_ps(midterm * delta_square * deltay);
						ENZ[i] += _mm512_reduce_add_ps(midterm * delta_square * deltaz);

					}
				}

			}

			double Energy = 0.00;
			for (int i = 0; i < NSource; i++)
			{
				Pot[i] = EN[i] + EF[i];
				Energy = Energy + Q[i] * Pot[i];
				Force[i][0] = -(EFX[i] + ENX[i]) * Q[i];
				Force[i][1] = -(EFY[i] + ENY[i]) * Q[i];
				Force[i][2] = -(EFZ[i] + ENZ[i]) * Q[i];
			}

			return Energy;
		}
		else
		{
			double EF[NSource], EFX[NSource], EFY[NSource], EFZ[NSource];
			double EN[NSource], ENX[NSource], ENY[NSource], ENZ[NSource];
			double Q_Image[int(ceil(ImageNumber / 8.0)) * 8];

			for (int j = 0; j < ImageNumber; j++)
			{
				Q_Image[j] = ImageCharge[j][3] ;
			}
			for (int j = ImageNumber; j<int(ceil(ImageNumber / 8.0)) * 8; j++)
			{
				Q_Image[j] = 0.00;
			}

			double C_New[int(ceil(p * p / 8.0)) * 8];
			for (int i = 0; i<int(ceil(p * p / 8.0)) * 8; i++)
			{
				if (i < p * p)
				{
					C_New[i] = C[i];
				}
				else
				{
					C_New[i] = 0.00;
				}
			}

			double Image_X[int(ceil(ImageNumber / 8.0)) * 8], Image_Y[int(ceil(ImageNumber / 8.0)) * 8], Image_Z[int(ceil(ImageNumber / 8.0)) * 8];
			for (int i = 0; i < int(ceil(ImageNumber / 8.0)) * 8; i++)
			{
				if (i < ImageNumber) {
					Image_X[i] = ImageCharge[i][0];
					Image_Y[i] = ImageCharge[i][1];
					Image_Z[i] = ImageCharge[i][2];
				}
				else {
					Image_X[i] = 0.00;
					Image_Y[i] = 0.00;
					Image_Z[i] = 0.00;
				}
			}

             #pragma omp parallel
			{
				int id = omp_get_thread_num();
				int size = comm->nthreads;

				int min_atom = id * floor(NSource / size) + 1, max_atom = (id + 1) * floor(NSource / size);
				if (id == size - 1)max_atom = NSource - 1;
				if (id == 0)min_atom = 0;

				for (int i = min_atom; i <= max_atom; i++)
				{
					double QF[int(ceil(p * p / 8.0)) * 8], QFX[int(ceil(p * p / 8.0)) * 8], QFY[int(ceil(p * p / 8.0)) * 8], QFZ[int(ceil(p * p / 8.0)) * 8];
					CalculateMultipoleExpansion(QF, p, Source[i][0], Source[i][1], Source[i][2]);
					CalculateZDerivativeMultipoleExpansion(QFZ, p, Source[i][0], Source[i][1], Source[i][2]);
					CalculateXDMultipoleExpansion(QFX, p, Source[i][0], Source[i][1], Source[i][2]);
					CalculateYDMultipoleExpansion(QFY, p, Source[i][0], Source[i][1], Source[i][2]);
					EF[i] = 0.00; EFX[i] = 0.00; EFY[i] = 0.00; EFZ[i] = 0.00;
					EN[i] = 0.00; ENX[i] = 0.00; ENY[i] = 0.00; ENZ[i] = 0.00;
					for (int ii = p * p; ii<int(ceil(p * p / 8.0)) * 8; ii++)
					{
						QF[ii] = 0.00;
						QFX[ii] = 0.00;
						QFY[ii] = 0.00;
						QFZ[ii] = 0.00;
					}
					__m512d qf, qfz, qfx, qfy, c;
					for (int j = 0; j < p * p; j = j + 8)
					{
						qf = _mm512_load_pd(&QF[j]);
						qfz = _mm512_load_pd(&QFZ[j]);
						qfx = _mm512_load_pd(&QFX[j]);
						qfy = _mm512_load_pd(&QFY[j]);
						c = _mm512_load_pd(&C_New[j]);
						EF[i] = EF[i] + _mm512_reduce_add_pd(qf * c);
						EFX[i] = EFX[i] + _mm512_reduce_add_pd(qfx * c);
						EFY[i] = EFY[i] + _mm512_reduce_add_pd(qfy * c);
						EFZ[i] = EFZ[i] + _mm512_reduce_add_pd(qfz * c);
					}
					__m512d X0, Y0, Z0, X1, Y1, Z1, Q1, deltax, deltay, deltaz, delta, square, judge, Zero, midterm, delta_square;
					__mmask16 k0;
					X0 = _mm512_set1_pd(Source[i][0]);
					Y0 = _mm512_set1_pd(Source[i][1]);
					Z0 = _mm512_set1_pd(Source[i][2]);
					judge = _mm512_set1_pd(0.000000000001);
					Zero = _mm512_set1_pd(0.00);
					for (int j = 0; j < ImageNumber; j = j + 8)
					{
						X1 = _mm512_load_pd(&Image_X[j]);
						Y1 = _mm512_load_pd(&Image_Y[j]);
						Z1 = _mm512_load_pd(&Image_Z[j]);
						Q1 = _mm512_load_pd(&Q_Image[j]);
						deltax = X1 - X0;
						deltay = Y1 - Y0;
						deltaz = Z1 - Z0;
						square = deltax * deltax + deltay * deltay + deltaz * deltaz;
						k0 = _mm512_cmp_pd_mask(square, judge, _MM_CMPINT_GT);
						delta = _mm512_mask_invsqrt_pd(Zero, k0, square);
						midterm = Q1 * delta;
						delta_square = delta * delta;
						EN[i] += _mm512_reduce_add_pd(midterm);
						ENX[i] += _mm512_reduce_add_pd(midterm * delta_square * deltax);
						ENY[i] += _mm512_reduce_add_pd(midterm * delta_square * deltay);
						ENZ[i] += _mm512_reduce_add_pd(midterm * delta_square * deltaz);
					}
				}
			}

			double Energy = 0.00;
			for (int i = 0; i < NSource; i++)
			{
				Pot[i] = EN[i] + EF[i];
				Energy = Energy + Q[i] * Pot[i];
				Force[i][0] = -(EFX[i] + ENX[i]) * Q[i];
				Force[i][1] = -(EFY[i] + ENY[i]) * Q[i];
				Force[i][2] = -(EFZ[i] + ENZ[i]) * Q[i];
			}
			return Energy;
		}
	}
	else
	{
		double eps = tolerance;
		int ns = ImageNumber;
		int nt = NSource;
		double *source = (double *)malloc(3*ns*sizeof(double));
		double* target = (double*)malloc(3 * nt * sizeof(double));
		double *charge = (double *)malloc(ns*sizeof(double));

		for(int i=0;i<ns;i++)
		{
			source[3*i]=ImageCharge[i][0];
			source[3*i+1]=ImageCharge[i][1];
			source[3*i+2]=ImageCharge[i][2];
			charge[i]=ImageCharge[i][3];	
		}
		
		for (int i = 0; i < nt; i++)
		{
			target[3 * i] = Source[i][0];
			target[3 * i + 1] = Source[i][1];
			target[3 * i + 2] = Source[i][2];
		}

		double *pottarg = (double *)malloc(nt*sizeof(double));
		double *gradtarg = (double *)malloc(3*nt*sizeof(double));

		int ier1;
		lfmm3d_t_c_g_(&eps, &ns, source, charge, &nt, target, pottarg, gradtarg,&ier1);

		int KL = int(2 * F + 2);

		/*             BEGIN HSMA ALGORITHM                          */
		double EF[KL];
		double CenterPara;

		//#pragma omp parallel for
		for (int i = 0; i < KL; i++)
		{
			CenterPara = 0.00;
			for (int j = 0; j < p * p; j++)
			{
				CenterPara = CenterPara + (1 / (4 * PI)) * (QRD[i][j] - QLocalRD[i][j]) * C[j];
			}
			EF[i] = CenterPara * Fibonacci[i][3];
		}

		/*            Set FMM parameters           */
		ns = 2 * F + 2;
		double* sourceF = (double*)malloc(3 * ns * sizeof(double));
		double* chargeF = (double*)malloc(ns * sizeof(double));

		for (int i = 0; i < ns; i++)
		{
			sourceF[3 * i] = Fibonacci[i][0];
			sourceF[3 * i + 1] = Fibonacci[i][1];
			sourceF[3 * i + 2] = Fibonacci[i][2];
			chargeF[i] = EF[i];
		}

		double* pottargF = (double*)malloc(nt * sizeof(double));
		double* gradtargF = (double*)malloc(3 * nt * sizeof(double));

		int ier;
		lfmm3d_t_c_g_(&eps, &ns, sourceF, chargeF, &nt, target, pottargF, gradtargF,&ier);

		/*	Final Summation	*/
		double Energy = 0.00;
		for (int i = 0; i < NSource; i++)
		{
			Pot[i] = pottargF[i]+ pottarg[i];
			Energy = Energy + Pot[i] * Q[i];
			Force[i][0] = - (gradtargF[3 * i]+ gradtarg[3 * i]) * Q[i];
			Force[i][1] = - (gradtargF[3 * i + 1]+ gradtarg[3 * i+1]) * Q[i];
			Force[i][2] = - (gradtargF[3 * i + 2]+ gradtarg[3 * i+2]) * Q[i];
		}

		free(target);
		target = NULL;
		free(sourceF); free(chargeF); free(pottargF); free(gradtargF);
		sourceF = NULL; chargeF = NULL; pottargF = NULL; gradtargF = NULL;
		return Energy;
	}
}

void HSMA3D::CalculateNearFieldAndZD_Single(double* Near, double ImageCharge[][4], int ImageNumber, int Nw, int IF_FMM_RightTerm, double Lx, double Ly, double Lz, double **PointSum, double **QuizSum, double Source[][3], double Force[][3], double* Pot, int NSource, double* Q, double tolerance)
{
	if (IF_FMM_RightTerm)//Using FMM to calculate pairwise sum
	{
		double eps = tolerance;

		/*            Set FMM parameters           */
		int ns = ImageNumber;
		int nt = Nw * 2 + NSource;
		double* source = (double*)malloc(3 * ns * sizeof(double));
		double* target = (double*)malloc(3 * nt * sizeof(double));
		double* charge = (double*)malloc(ns * sizeof(double));

		for (int i = 0; i < ns; i++)
		{
			source[3 * i] = ImageCharge[i][0];
			source[3 * i + 1] = ImageCharge[i][1];
			source[3 * i + 2] = ImageCharge[i][2];
			charge[i] = ImageCharge[i][3];
		}

		for (int i = 0; i < Nw; i++)
		{
			target[3 * i] = PointSum[i][0];
			target[3 * i + 1] = PointSum[i][1];
			target[3 * i + 2] = PointSum[i][2];

			target[3 * Nw + 3 * i] = QuizSum[i][0];
			target[3 * Nw + 3 * i + 1] = QuizSum[i][1];
			target[3 * Nw + 3 * i + 2] = QuizSum[i][2];
		}

		for (int i = 0; i < NSource; i++)
		{
			target[6 * Nw + 3 * i] = Source[i][0];
			target[6 * Nw + 3 * i + 1] = Source[i][1];
			target[6 * Nw + 3 * i + 2] = Source[i][2];
		}

		double* pottarg = (double*)malloc(nt * sizeof(double));
		double* gradtarg = (double*)malloc(3 * nt * sizeof(double));

		int ier;
		lfmm3d_t_c_g_(&eps, &ns, source, charge, &nt, target, pottarg, gradtarg,&ier);
		for (int i = 0; i < Nw; i++)
		{
			Near[i] = pottarg[Nw + i] - pottarg[i];
		}
		for (int i = 0; i < NSource; i++)
		{
			Pot[i] = pottarg[2 * Nw + i];
			Force[i][0] = -gradtarg[6 * Nw + 3 * i] * Q[i];
			Force[i][1] = -gradtarg[6 * Nw + 3 * i + 1] * Q[i];
			Force[i][2] = -gradtarg[6 * Nw + 3 * i + 2] * Q[i];
		}


		free(source); free(target); free(charge); free(pottarg);
		source = NULL; target = NULL; charge = NULL; pottarg = NULL;
	}
	else
	{
		if (tolerance > 0.000001)
		{
			float Paramet1[int(ceil(ImageNumber / 16.0)) * 16];
			float Image_X[int(ceil(ImageNumber / 16.0)) * 16], Image_Y[int(ceil(ImageNumber / 16.0)) * 16], Image_Z[int(ceil(ImageNumber / 16.0)) * 16];
			for (int i = 0; i < int(ceil(ImageNumber / 16.0)) * 16; i++)
			{
				if (i < ImageNumber) {
					Image_X[i] = ImageCharge[i][0];
					Image_Y[i] = ImageCharge[i][1];
					Image_Z[i] = ImageCharge[i][2];
					Paramet1[i] = ImageCharge[i][3];
				}
				else {
					Image_X[i] = 0.00;
					Image_Y[i] = 0.00;
					Image_Z[i] = 0.00;
					Paramet1[i] = 0.00;
				}
			}

            #pragma omp parallel
			{
				int id = omp_get_thread_num();
				int size = comm->nthreads;

				int min_atom = id * floor(Nw / (size + 0.00)) + 1, max_atom = (id + 1) * floor(Nw / (size + 0.00));
				if (id == size - 1)max_atom = Nw - 1;
				if (id == 0)min_atom = 0;

				for (int i = min_atom; i <= max_atom; i++)
				{
					__m512 XP, YP, ZP, XQ, YQ, ZQ, X1, Y1, Z1, Q1, dx, dy, dz, dx1, dy1, dz1, square, square1;
					XP = _mm512_set1_ps(PointSum[i][0]);
					YP = _mm512_set1_ps(PointSum[i][1]);
					ZP = _mm512_set1_ps(PointSum[i][2]);
					XQ = _mm512_set1_ps(QuizSum[i][0]);
					YQ = _mm512_set1_ps(QuizSum[i][1]);
					ZQ = _mm512_set1_ps(QuizSum[i][2]);
					float pottarg = 0.00, pottarg1 = 0.00;
					for (int j = 0; j < ImageNumber; j = j + 16)
					{
						X1 = _mm512_load_ps(&Image_X[j]);
						Y1 = _mm512_load_ps(&Image_Y[j]);
						Z1 = _mm512_load_ps(&Image_Z[j]);
						Q1 = _mm512_load_ps(&Paramet1[j]);

						dx = XP - X1;
						dy = YP - Y1;
						dz = ZP - Z1;
						square = dx * dx + dy * dy + dz * dz;
						pottarg += _mm512_reduce_add_ps(Q1 * _mm512_invsqrt_ps(square));

						dx1 = XQ - X1;
						dy1 = YQ - Y1;
						dz1 = ZQ - Z1;
						square1 = dx1 * dx1 + dy1 * dy1 + dz1 * dz1;
						pottarg1 += _mm512_reduce_add_ps(Q1 * _mm512_invsqrt_ps(square1));
					}
					Near[i] = pottarg1 - pottarg;
				}
			}

            #pragma omp parallel
			{
				int id = omp_get_thread_num();
				int size = comm->nthreads;

				int min_atom = id * floor(NSource / size) + 1, max_atom = (id + 1) * floor(NSource / size);
				if (id == size - 1)max_atom = NSource - 1;
				if (id == 0)min_atom = 0;

				for (int i = min_atom; i <= max_atom; i++)
				{
					__m512 X0, Y0, Z0, X1, Y1, Z1, Q1, deltax, deltay, deltaz, delta, square, judge, Zero, midterm, delta_square;
					__mmask16 k0;
					X0 = _mm512_set1_ps(Source[i][0]);
					Y0 = _mm512_set1_ps(Source[i][1]);
					Z0 = _mm512_set1_ps(Source[i][2]);
					judge = _mm512_set1_ps(0.000000000001);
					Zero = _mm512_set1_ps(0.00);
					float a = 0.00, b = 0.00, c = 0.00, d = 0.00;
					for (int j = 0; j < ImageNumber; j = j + 16)
					{
						X1 = _mm512_load_ps(&Image_X[j]);
						Y1 = _mm512_load_ps(&Image_Y[j]);
						Z1 = _mm512_load_ps(&Image_Z[j]);
						Q1 = _mm512_load_ps(&Paramet1[j]);
						deltax = X1 - X0;
						deltay = Y1 - Y0;
						deltaz = Z1 - Z0;
						square = deltax * deltax + deltay * deltay + deltaz * deltaz;
						k0 = _mm512_cmp_ps_mask(square, judge, _MM_CMPINT_GT);
						delta = _mm512_mask_invsqrt_ps(Zero, k0, square);
						midterm = Q1 * delta;
						delta_square = delta * delta;

						a += _mm512_reduce_add_ps(midterm);
						b -= _mm512_reduce_add_ps(midterm * delta_square * deltax);
						c -= _mm512_reduce_add_ps(midterm * delta_square * deltay);
						d -= _mm512_reduce_add_ps(midterm * delta_square * deltaz);
					}
					Pot[i] = a;
					Force[i][0] = b * Q[i];
					Force[i][1] = c * Q[i];
					Force[i][2] = d * Q[i];
				}
			}
		}
		else
		{
			double Paramet1[int(ceil(ImageNumber / 8.0)) * 8];
			double Image_X[int(ceil(ImageNumber / 8.0)) * 8], Image_Y[int(ceil(ImageNumber / 8.0)) * 8], Image_Z[int(ceil(ImageNumber / 8.0)) * 8];
			for (int i = 0; i < int(ceil(ImageNumber / 8.0)) * 8; i++)
			{
				if (i < ImageNumber) {
					Image_X[i] = ImageCharge[i][0];
					Image_Y[i] = ImageCharge[i][1];
					Image_Z[i] = ImageCharge[i][2];
					Paramet1[i] = ImageCharge[i][3];
				}
				else {
					Image_X[i] = 0.00;
					Image_Y[i] = 0.00;
					Image_Z[i] = 0.00;
					Paramet1[i] = 0.00;
				}
			}

             #pragma omp parallel
			{
				int id = omp_get_thread_num();
				int size = comm->nthreads;

				int min_atom = id * floor(Nw / (size + 0.00)) + 1, max_atom = (id + 1) * floor(Nw / (size + 0.00));
				if (id == size - 1)max_atom = Nw - 1;
				if (id == 0)min_atom = 0;

				for (int i = min_atom; i <= max_atom; i++)
				{
					__m512d XP, YP, ZP, XQ, YQ, ZQ, X1, Y1, Z1, Q1, dx, dy, dz, dx1, dy1, dz1, square, square1;
					XP = _mm512_set1_pd(PointSum[i][0]);
					YP = _mm512_set1_pd(PointSum[i][1]);
					ZP = _mm512_set1_pd(PointSum[i][2]);
					XQ = _mm512_set1_pd(QuizSum[i][0]);
					YQ = _mm512_set1_pd(QuizSum[i][1]);
					ZQ = _mm512_set1_pd(QuizSum[i][2]);
					double pottarg = 0.00, pottarg1 = 0.00;

					for (int j = 0; j < ImageNumber; j = j + 8)
					{
						X1 = _mm512_load_pd(&Image_X[j]);
						Y1 = _mm512_load_pd(&Image_Y[j]);
						Z1 = _mm512_load_pd(&Image_Z[j]);
						Q1 = _mm512_load_pd(&Paramet1[j]);

						dx = XP - X1;
						dy = YP - Y1;
						dz = ZP - Z1;
						square = dx * dx + dy * dy + dz * dz;
						pottarg += _mm512_reduce_add_pd(Q1 * _mm512_invsqrt_pd(square));

						dx1 = XQ - X1;
						dy1 = YQ - Y1;
						dz1 = ZQ - Z1;
						square1 = dx1 * dx1 + dy1 * dy1 + dz1 * dz1;
						pottarg1 += _mm512_reduce_add_pd(Q1 * _mm512_invsqrt_pd(square1));
					}
					Near[i] = pottarg1 - pottarg;
				}
			}

            #pragma omp parallel
			{
				int id = omp_get_thread_num();
				int size = comm->nthreads();

				int min_atom = id * floor(NSource / size) + 1, max_atom = (id + 1) * floor(NSource / size);
				if (id == size - 1)max_atom = NSource - 1;
				if (id == 0)min_atom = 0;

				for (int i = min_atom; i <= max_atom; i++)
				{
					__m512d X0, Y0, Z0, X1, Y1, Z1, Q1, deltax, deltay, deltaz, delta, square, judge, Zero, midterm, delta_square;
					__mmask16 k0;
					X0 = _mm512_set1_pd(Source[i][0]);
					Y0 = _mm512_set1_pd(Source[i][1]);
					Z0 = _mm512_set1_pd(Source[i][2]);
					judge = _mm512_set1_pd(0.000000000001);
					Zero = _mm512_set1_pd(0.00);
					double a = 0.00, b = 0.00, c = 0.00, d = 0.00;
					for (int j = 0; j < ImageNumber; j = j + 8)
					{
						X1 = _mm512_load_pd(&Image_X[j]);
						Y1 = _mm512_load_pd(&Image_Y[j]);
						Z1 = _mm512_load_pd(&Image_Z[j]);
						Q1 = _mm512_load_pd(&Paramet1[j]);
						deltax = X1 - X0;
						deltay = Y1 - Y0;
						deltaz = Z1 - Z0;
						square = deltax * deltax + deltay * deltay + deltaz * deltaz;
						k0 = _mm512_cmp_pd_mask(square, judge, _MM_CMPINT_GT);
						delta = _mm512_mask_invsqrt_pd(Zero, k0, square);
						midterm = Q1 * delta;
						delta_square = delta * delta;

						a += _mm512_reduce_add_pd(midterm);
						b -= _mm512_reduce_add_pd(midterm * delta_square * deltax);
						c -= _mm512_reduce_add_pd(midterm * delta_square * deltay);
						d -= _mm512_reduce_add_pd(midterm * delta_square * deltaz);
					}
					Pot[i] = a;
					Force[i][0] = b * Q[i];
					Force[i][1] = c * Q[i];
					Force[i][2] = d * Q[i];
				}
			}
		}
	}
}

double HSMA3D::FinalCalculateEnergyAndForce_Single(double Force[][3], double* Pot, double Source[][3], double* Q, int NSource, double ImageCharge[][4], int ImageNumber, double **Fibonacci, double** QRD, double** QLocalRD, double* C, int p, double Fp, double F, double Rs, double PI, int IF_FMM_FinalPotential, double tolerance)
{
	if (!IF_FMM_FinalPotential)
	{
			double EF[NSource], EFX[NSource], EFY[NSource], EFZ[NSource];
			double C_New[int(ceil(p * p / 8.0)) * 8];
			for (int i = 0; i<int(ceil(p * p / 8.0)) * 8; i++)
			{
				if (i < p * p)
				{
					C_New[i] = C[i];
				}
				else
				{
					C_New[i] = 0.00;
				}
			}

            #pragma omp parallel
			{
				int id = omp_get_thread_num();
				int size = comm->nthreads;

				int min_atom = id * floor(NSource / size) + 1, max_atom = (id + 1) * floor(NSource / size);
				if (id == size - 1)max_atom = NSource - 1;
				if (id == 0)min_atom = 0;

				for (int i = min_atom; i <= max_atom; i++)
				{
					double QF[int(ceil(p * p / 8.0)) * 8], QFX[int(ceil(p * p / 8.0)) * 8], QFY[int(ceil(p * p / 8.0)) * 8], QFZ[int(ceil(p * p / 8.0)) * 8];
					CalculateMultipoleExpansion(QF, p, Source[i][0], Source[i][1], Source[i][2]);
					CalculateZDerivativeMultipoleExpansion(QFZ, p, Source[i][0], Source[i][1], Source[i][2]);
					CalculateXDMultipoleExpansion(QFX, p, Source[i][0], Source[i][1], Source[i][2]);
					CalculateYDMultipoleExpansion(QFY, p, Source[i][0], Source[i][1], Source[i][2]);
					float a = 0.00, b = 0.00, d = 0.00, e = 0.00;
					for (int ii = p * p; ii<int(ceil(p * p / 8.0)) * 8; ii++)
					{
						QF[ii] = 0.00;
						QFX[ii] = 0.00;
						QFY[ii] = 0.00;
						QFZ[ii] = 0.00;
					}
					__m512d qf, qfz, qfx, qfy, c;
					for (int j = 0; j < p * p; j = j + 8)
					{
						qf = _mm512_load_pd(&QF[j]);
						qfz = _mm512_load_pd(&QFZ[j]);
						qfx = _mm512_load_pd(&QFX[j]);
						qfy = _mm512_load_pd(&QFY[j]);
						c = _mm512_load_pd(&C_New[j]);
						a = a + _mm512_reduce_add_pd(qf * c);
						b = b + _mm512_reduce_add_pd(qfx * c);
						d = d + _mm512_reduce_add_pd(qfy * c);
						e = e + _mm512_reduce_add_pd(qfz * c);
					}
					EF[i] = a;
					EFX[i] = b;
					EFY[i] = d;
					EFZ[i] = e;
				}

			}

			double Energy = 0.00;
			for (int i = 0; i < NSource; i++)
			{
				Pot[i] += EF[i];
				Energy += Q[i] * Pot[i];
				Force[i][0] -= (EFX[i]) * Q[i];
				Force[i][1] -= (EFY[i]) * Q[i];
				Force[i][2] -= (EFZ[i]) * Q[i];
			}

			return Energy;
	}
	else
	{
		double eps = tolerance;
		int ns = ImageNumber;
		int nt = NSource;
		double* target = (double*)malloc(3 * nt * sizeof(double));

		for (int i = 0; i < nt; i++)
		{
			target[3 * i] = Source[i][0];
			target[3 * i + 1] = Source[i][1];
			target[3 * i + 2] = Source[i][2];
		}

		int KL = int(2 * F + 2);

		/*             BEGIN HSMA ALGORITHM                          */
		double EF[KL];
		double CenterPara;

		//#pragma omp parallel for
		for (int i = 0; i < KL; i++)
		{
			CenterPara = 0.00;
			for (int j = 0; j < p * p; j++)
			{
				CenterPara = CenterPara + (1 / (4 * PI)) * (QRD[i][j] - QLocalRD[i][j]) * C[j];//*((double*)QRD + i * p * p + j) - *((double*)QLocalRD + i * p * p + j)
			}
			EF[i] = CenterPara * Fibonacci[i][3];
		}

		/*            Set FMM parameters           */
		ns = 2 * F + 2;
		double* sourceF = (double*)malloc(3 * ns * sizeof(double));
		double* chargeF = (double*)malloc(ns * sizeof(double));

		for (int i = 0; i < ns; i++)
		{
			sourceF[3 * i] = Fibonacci[i][0];
			sourceF[3 * i + 1] = Fibonacci[i][1];
			sourceF[3 * i + 2] = Fibonacci[i][2];
			chargeF[i] = EF[i];
		}

		double* pottargF = (double*)malloc(nt * sizeof(double));
		double* gradtargF = (double*)malloc(3 * nt * sizeof(double));

		int ier;
		lfmm3d_t_c_g_(&eps, &ns, sourceF, chargeF, &nt, target, pottargF, gradtargF,&ier);

		/*	Final Summation	*/
		double Energy = 0.00;
		for (int i = 0; i < NSource; i++)
		{
			Pot[i] = Pot[i] + pottargF[i];
			Energy = Energy + Pot[i] * Q[i];
			Force[i][0] = Force[i][0] - (gradtargF[3 * i]) * Q[i];
			Force[i][1] = Force[i][1] - (gradtargF[3 * i + 1]) * Q[i];
			Force[i][2] = Force[i][2] - (gradtargF[3 * i + 2]) * Q[i];
		}

		free(target);
		target = NULL;
		free(sourceF); free(chargeF); free(pottargF); free(gradtargF);
		sourceF = NULL; chargeF = NULL; pottargF = NULL; gradtargF = NULL;

		return Energy;
	}
}

int fab(int n)
{
	if (n == 1)
		return 1;
	if (n == 2)
		return 1;
	if (n > 2)
		return fab(n - 1) + fab(n - 2);
}

int isfab(int m)
{
	int result = 0;
	for (int i = 0; fab(i) < m; i++)
	{
		if (fab(i + 1) == m)
		{
			result = 1;
			break;
		}
		else
			result = 0;

	}
	return result;
}

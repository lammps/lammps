// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors: Shawn Coleman & Douglas Spearot (Arkansas)
   Updated: 06/17/2015-2
------------------------------------------------------------------------- */

#include "compute_xrd.h"
#include "compute_xrd_consts.h"

#include "atom.h"
#include "citeme.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "group.h"
#include "math_const.h"
#include "memory.h"
#include "update.h"

#include <cmath>
#include <cstring>

#include "omp_compat.h"
using namespace LAMMPS_NS;
using namespace MathConst;

static const char cite_compute_xrd_c[] =
  "compute xrd command: doi:10.1088/0965-0393/21/5/055020\n\n"
  "@Article{Coleman13,\n"
  " author = {S. P. Coleman and D. E. Spearot and L. Capolungo},\n"
  " title = {Virtual Diffraction Analysis of {Ni} [010] Symmetric Tilt Grain Boundaries},\n"
  " journal = {Modelling and Simulation in Materials Science and Engineering},\n"
  " year =    2013,\n"
  " volume =  21,\n"
  " pages =   {055020}\n"
  "}\n\n";

/* ---------------------------------------------------------------------- */

ComputeXRD::ComputeXRD(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg), ztype(nullptr), store_tmp(nullptr)
{
  if (lmp->citeme) lmp->citeme->add(cite_compute_xrd_c);

  int ntypes = atom->ntypes;
  int natoms = group->count(igroup);
  int dimension = domain->dimension;
  int *periodicity = domain->periodicity;
  int triclinic = domain->triclinic;
  me = comm->me;

  // Checking errors
  if (dimension == 2)
     error->all(FLERR,"Compute XRD does not work with 2d structures");
  if (narg < 4+ntypes)
     error->all(FLERR,"Illegal Compute XRD Command");
  if (triclinic == 1)
     error->all(FLERR,"Compute XRD does not work with triclinic structures");

  array_flag = 1;
  extarray = 0;

  // Store radiation wavelength
  lambda = utils::numeric(FLERR,arg[3],false,lmp);
  if (lambda < 0)
    error->all(FLERR,"Compute SAED: Wavelength must be greater than zero");

  // Define atom types for atomic scattering factor coefficients
  int iarg = 4;
  ztype = new int[ntypes];
  for (int i = 0; i < ntypes; i++) {
    ztype[i] = XRDmaxType + 1;
  }
  for (int i = 0; i < ntypes; i++) {
    for (int j = 0; j < XRDmaxType; j++) {
      if (utils::lowercase(arg[iarg]) == utils::lowercase(XRDtypeList[j])) {
        ztype[i] = j;
       }
     }
    if (ztype[i] == XRDmaxType + 1)
        error->all(FLERR,"Compute XRD: Invalid ASF atom type");
    iarg++;
  }

  // Set defaults for optional args
  Min2Theta = 1;
  Max2Theta = 179;
  radflag = 1;
  c[0] = 1; c[1] = 1; c[2] = 1;
  LP = 1;
  manual = false;
  echo = false;

  // Process optional args
  while (iarg < narg) {
    if (strcmp(arg[iarg],"2Theta") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal Compute XRD Command");
      Min2Theta = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      Max2Theta = utils::numeric(FLERR,arg[iarg+2],false,lmp);
      iarg += 3;

    } else if (strcmp(arg[iarg],"c") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Illegal Compute XRD Command");
      c[0] = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      c[1] = utils::numeric(FLERR,arg[iarg+2],false,lmp);
      c[2] = utils::numeric(FLERR,arg[iarg+3],false,lmp);
      if (c[0] < 0 || c[1] < 0 || c[2] < 0)
        error->all(FLERR,"Compute XRD: c's must be greater than 0");
      iarg += 4;

    } else if (strcmp(arg[iarg],"LP") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal Compute XRD Command");
      LP = utils::numeric(FLERR,arg[iarg+1],false,lmp);

      if (!(LP == 1 || LP == 0))
         error->all(FLERR,"Compute XRD: LP must have value of 0 or 1");
      iarg += 2;

    } else if (strcmp(arg[iarg],"echo") == 0) {
      echo = true;
      iarg += 1;

    } else if (strcmp(arg[iarg],"manual") == 0) {
      manual = true;
      iarg += 1;

    } else error->all(FLERR,"Illegal Compute XRD Command");
  }

  // error check and process min/max 2Theta values
  Min2Theta /= 2.0;
  Max2Theta /= 2.0;
  if (Max2Theta > MY_PI) {
    Min2Theta = Min2Theta * MY_PI / 180;  // converting to radians if necessary
    Max2Theta = Max2Theta * MY_PI / 180;
    radflag = 0;
  }
  if (Min2Theta <= 0)
    error->all(FLERR,"Minimum 2Theta value must be greater than zero");
  if (Max2Theta >= MY_PI )
    error->all(FLERR,"Maximum 2Theta value must be less than 180 degrees");
  if (Max2Theta-Min2Theta <= 0)
    error->all(FLERR,"2Theta range must be greater than zero");
  Kmax = 2 * sin(Max2Theta) / lambda;

  // Calculating spacing between reciprocal lattice points
  // Using distance based on periodic repeating distance
  if (!manual) {
    if (!periodicity[0] && !periodicity[1] && !periodicity[2])
      error->all(FLERR,"Compute SAED must have at least one periodic boundary unless manual spacing specified");

    double *prd;
    double ave_inv = 0.0;
    prd = domain->prd;

    if (periodicity[0]) {
      prd_inv[0] = 1 / prd[0];
      ave_inv += prd_inv[0];
    }
    if (periodicity[1]) {
      prd_inv[1] = 1 / prd[1];
      ave_inv += prd_inv[1];
    }
    if (periodicity[2]) {
      prd_inv[2] = 1 / prd[2];
      ave_inv += prd_inv[2];
    }

    // Using the average inverse dimensions for non-periodic direction
    ave_inv = ave_inv / (periodicity[0] + periodicity[1] + periodicity[2]);
    if (!periodicity[0]) {
      prd_inv[0] = ave_inv;
    }
    if (!periodicity[1]) {
      prd_inv[1] = ave_inv;
    }
    if (!periodicity[2]) {
      prd_inv[2] = ave_inv;
    }
  }

  // Use manual mapping of reciprocal lattice
  if (manual) {
    for (int i=0; i<3; i++) {
      prd_inv[i] = 1.0;
    }
  }

  // Find reciprocal spacing and integer dimensions
  for (int i=0; i<3; i++) {
    dK[i] = prd_inv[i]*c[i];
    Knmax[i] = ceil(Kmax / dK[i]);
  }

  // Finding the intersection of the reciprocal space and Ewald sphere
  int nRows = 0;
  double dinv2= 0.0;
  double ang = 0.0;
  double K[3];

  // Procedure to determine how many rows are needed given the constraints on 2theta
  for (int i = -Knmax[0]; i <= Knmax[0]; i++) {
    for (int j = -Knmax[1]; j <= Knmax[1]; j++) {
      for (int k = -Knmax[2]; k <= Knmax[2]; k++) {
        K[0] = i * dK[0];
        K[1] = j * dK[1];
        K[2] = k * dK[2];
        dinv2 = (K[0] * K[0] + K[1] * K[1] + K[2] * K[2]);
        if  (4 >= dinv2 * lambda * lambda) {
          ang = asin(lambda * sqrt(dinv2) * 0.5);
          if ((ang <= Max2Theta) && (ang >= Min2Theta)) {
          nRows++;
                }
        }
      }
    }
  }


  size_array_rows = nRows;
  size_array_cols = 2;

  if (me == 0 && echo)
    utils::logmesg(lmp,"-----\nCompute XRD id:{}, # of atoms:{}, # of relp:{}\n"
                   "Reciprocal point spacing in k1,k2,k3 = {:.8} {:.8} {:.8}\n-----\n",
                   id,natoms,nRows,dK[0],dK[1],dK[2]);

  memory->create(array,size_array_rows,size_array_cols,"xrd:array");
  memory->create(store_tmp,3*size_array_rows,"xrd:store_tmp");
}

/* ---------------------------------------------------------------------- */

ComputeXRD::~ComputeXRD()
{
  memory->destroy(array);
  memory->destroy(store_tmp);
  delete[] ztype;
}

/* ---------------------------------------------------------------------- */

void ComputeXRD::init()
{

  int mmax = (2*Knmax[0]+1)*(2*Knmax[1]+1)*(2*Knmax[2]+1);
  double K[3];
  double dinv2 = 0.0;
  double ang = 0.0;

  double convf = 360 / MY_PI;
  if (radflag ==1) convf = 1;

  int n = 0;
  for (int m = 0; m < mmax; m++) {
    int k = m%(2*Knmax[2]+1);
    int j = (m%((2*Knmax[2]+1)*(2*Knmax[1]+1))-k)/(2*Knmax[2]+1);
    int i = (m-j*(2*Knmax[2]+1)-k)/((2*Knmax[2]+1)*(2*Knmax[1]+1))-Knmax[0];
    j = j-Knmax[1];
    k = k-Knmax[2];
    K[0] = i * dK[0];
    K[1] = j * dK[1];
    K[2] = k * dK[2];
    dinv2 = (K[0] * K[0] + K[1] * K[1] + K[2] * K[2]);
    if  (4 >= dinv2 * lambda * lambda) {
       ang = asin(lambda * sqrt(dinv2) * 0.5);
       if ((ang <= Max2Theta) && (ang >= Min2Theta)) {
          store_tmp[3*n] = k;
          store_tmp[3*n+1] = j;
          store_tmp[3*n+2] = i;
          array[n][0] = ang * convf;
          n++;
       }
    }
  }
 if (n != size_array_rows)
     error->all(FLERR,"Compute XRD compute_array() rows mismatch");

}

/* ---------------------------------------------------------------------- */

void ComputeXRD::compute_array()
{
  invoked_array = update->ntimestep;

  if (me == 0 && echo) utils::logmesg(lmp, "-----\nComputing XRD intensities");

  double t0 = platform::walltime();

  auto Fvec = new double[2*size_array_rows]; // Strct factor (real & imaginary)
  // -- Note: array rows correspond to different RELP

  ntypes = atom->ntypes;
  int nlocal = atom->nlocal;
  int *type  = atom->type;
  int natoms = group->count(igroup);
  int *mask = atom->mask;

  nlocalgroup = 0;
  for (int ii = 0; ii < nlocal; ii++) {
    if (mask[ii] & groupbit) {
     nlocalgroup++;
    }
  }

  auto xlocal = new double [3*nlocalgroup];
  int *typelocal = new int [nlocalgroup];

  nlocalgroup = 0;
  for (int ii = 0; ii < nlocal; ii++) {
    if (mask[ii] & groupbit) {
     xlocal[3*nlocalgroup+0] = atom->x[ii][0];
     xlocal[3*nlocalgroup+1] = atom->x[ii][1];
     xlocal[3*nlocalgroup+2] = atom->x[ii][2];
     typelocal[nlocalgroup]=type[ii];
     nlocalgroup++;
    }
  }

// Setting up OMP
#if defined(_OPENMP)
  if ((me == 0) && echo) utils::logmesg(lmp," using {} OMP threads\n",comm->nthreads);
#endif

  if ((me == 0) && echo) {
    utils::logmesg(lmp,"\n");

    if (LP == 1)
        utils::logmesg(lmp,"Applying Lorentz-Polarization Factor During XRD Calculation 2\n");
  }

  int m = 0;
  double frac = 0.1;

#if defined(_OPENMP)
#pragma omp parallel LMP_DEFAULT_NONE LMP_SHARED(typelocal,xlocal,Fvec,m,frac,ASFXRD)
#endif
  {
    auto f = new double[ntypes];    // atomic structure factor by type
    int n,typei = 0;

    double Fatom1 = 0.0;               // structure factor per atom (real)
    double Fatom2 = 0.0;               // structure factor per atom (imaginary)

    double K[3];
    double dinv2 = 0.0;
    double dinv  = 0.0;
    double SinTheta_lambda  = 0.0;
    double SinTheta = 0.0;
    double ang = 0.0;
    double Cos2Theta = 0.0;
    double CosTheta = 0.0;

    double inners = 0.0;
    double sqrt_lp = 0.0;

    if (LP == 1) {

#if defined(_OPENMP)
#pragma omp for
#endif
      for (n = 0; n < size_array_rows; n++) {
        int k = store_tmp[3*n];
        int j = store_tmp[3*n+1];
        int i = store_tmp[3*n+2];
        K[0] = i * dK[0];
        K[1] = j * dK[1];
        K[2] = k * dK[2];

        dinv2 = (K[0] * K[0] + K[1] * K[1] + K[2] * K[2]);
        dinv = sqrt(dinv2);
        SinTheta_lambda = 0.5*dinv;
        SinTheta = SinTheta_lambda * lambda;
        ang = asin( SinTheta );
        Cos2Theta = cos( 2 * ang);
        CosTheta = cos( ang );

        Fatom1 = 0.0;
        Fatom2 = 0.0;

        // Calculate the atomic structure factor by type
        for (int ii = 0; ii < ntypes; ii++) {
          f[ii] = 0;
          for (int C = 0; C < 8 ; C+=2) {
            f[ii] += ASFXRD[ztype[ii]][C] * exp(-1 * ASFXRD[ztype[ii]][C+1] * SinTheta_lambda * SinTheta_lambda );
          }
          f[ii] += ASFXRD[ztype[ii]][8];
        }

        // Evaluate the structure factor equation -- looping over all atoms
        for (int ii = 0; ii < nlocalgroup; ii++) {
          typei=typelocal[ii]-1;
          inners = 2 * MY_PI * (K[0] * xlocal[3*ii] + K[1] * xlocal[3*ii+1] +
                    K[2] * xlocal[3*ii+2]);
          Fatom1 += f[typei] * cos(inners);
          Fatom2 += f[typei] * sin(inners);
        }
        sqrt_lp = sqrt( (1 + Cos2Theta * Cos2Theta) /
             ( CosTheta * SinTheta * SinTheta) );
        Fvec[2*n] = Fatom1 * sqrt_lp;
        Fvec[2*n+1] = Fatom2 * sqrt_lp;

        // reporting progress of calculation
        if (echo) {
#if defined(_OPENMP)
          #pragma omp critical
#endif
          {
            if (m == round(frac * size_array_rows)) {
              if (me == 0) utils::logmesg(lmp," {:2.0f}% -",frac*100);
              frac += 0.1;
            }
            m++;
          }
        }
      } // End of pragma omp for region

    } else {
#if defined(_OPENMP)
#pragma omp for
#endif
      for (n = 0; n < size_array_rows; n++) {
        int k = store_tmp[3*n];
        int j = store_tmp[3*n+1];
        int i = store_tmp[3*n+2];
        K[0] = i * dK[0];
        K[1] = j * dK[1];
        K[2] = k * dK[2];

        dinv2 = (K[0] * K[0] + K[1] * K[1] + K[2] * K[2]);
        dinv = sqrt(dinv2);
        SinTheta_lambda = 0.5*dinv;

        Fatom1 = 0.0;
        Fatom2 = 0.0;

        // Calculate the atomic structure factor by type
        for (int ii = 0; ii < ntypes; ii++) {
          f[ii] = 0;
          for (int C = 0; C < 8 ; C+=2) {
            f[ii] += ASFXRD[ztype[ii]][C] * exp(-1 * ASFXRD[ztype[ii]][C+1] * SinTheta_lambda * SinTheta_lambda );
          }
          f[ii] += ASFXRD[ztype[ii]][8];
        }

        // Evaluate the structure factor equation -- looping over all atoms
        for (int ii = 0; ii < nlocalgroup; ii++) {
          typei=typelocal[ii]-1;
          inners = 2 * MY_PI * (K[0] * xlocal[3*ii] + K[1] * xlocal[3*ii+1] +
                    K[2] * xlocal[3*ii+2]);
          Fatom1 += f[typei] * cos(inners);
          Fatom2 += f[typei] * sin(inners);
        }
        Fvec[2*n] = Fatom1;
        Fvec[2*n+1] = Fatom2;

        // reporting progress of calculation
        if (echo) {
#if defined(_OPENMP)
          #pragma omp critical
#endif
          {
            if (m == round(frac * size_array_rows)) {
              if (me == 0) utils::logmesg(lmp," {:2.0f}% -",frac*100);
              frac += 0.1;
            }
            m++;
          }
        }
      } // End of pragma omp for region
    } // End of if LP=1 check
    delete [] f;
  } // End of pragma omp parallel region

  auto scratch = new double[2*size_array_rows];

  // Sum intensity for each ang-hkl combination across processors
  MPI_Allreduce(Fvec,scratch,2*size_array_rows,MPI_DOUBLE,MPI_SUM,world);

  for (int i = 0; i < size_array_rows; i++) {
    array[i][1] = (scratch[2*i] * scratch[2*i] + scratch[2*i+1] * scratch[2*i+1]) / natoms;
  }

  double t2 = platform::walltime();

  // compute memory usage per processor
  double bytes = memory_usage();

  if (me == 0 && echo)
    utils::logmesg(lmp," 100% \nTime elapsed during compute_xrd = {:.2f} sec "
                   "using {:.2f} Mbytes/processor\n-----\n", t2-t0, bytes/1024.0/1024.0);

  delete [] scratch;
  delete [] Fvec;
  delete [] xlocal;
  delete [] typelocal;
}

/* ----------------------------------------------------------------------
 memory usage of arrays
 ------------------------------------------------------------------------- */

double ComputeXRD::memory_usage()
{
  double bytes = (double)size_array_rows * size_array_cols * sizeof(double); //array
  bytes += (double) 4.0 * size_array_rows * sizeof(double); //Fvec1 & 2, scratch1 & 2
  bytes += (double)3.0 * nlocalgroup * sizeof(double); // xlocal
  bytes += (double)nlocalgroup * sizeof(int); // typelocal
  bytes += (double)ntypes * sizeof(double); // f
  bytes += (double)3.0 * size_array_rows * sizeof(int); // store_temp

  return bytes;
}

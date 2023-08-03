// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors: Shawn Coleman & Douglas Spearot (Arkansas)
------------------------------------------------------------------------- */

#include "compute_saed.h"

#include "atom.h"
#include "citeme.h"
#include "comm.h"
#include "compute_saed_consts.h"
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

static const char cite_compute_saed_c[] =
  "compute saed command: doi:10.1088/0965-0393/21/5/055020\n\n"
  "@Article{Coleman13,\n"
  " author = {S. P. Coleman and D. E. Spearot and L. Capolungo},\n"
  " title = {Virtual Diffraction Analysis of {Ni} [010] Symmetric Tilt Grain Boundaries},\n"
  " journal = {Modelling and Simulation in Materials Science and Engineering},\n"
  " year =    2013,\n"
  " volume =  21,\n"
  " pages =   {055020}\n"
  "}\n\n";

/* ---------------------------------------------------------------------- */

ComputeSAED::ComputeSAED(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg), ztype(nullptr), store_tmp(nullptr)
{
  if (lmp->citeme) lmp->citeme->add(cite_compute_saed_c);

  int ntypes = atom->ntypes;
  int natoms = group->count(igroup);
  int dimension = domain->dimension;
  int *periodicity = domain->periodicity;
  int triclinic = domain->triclinic;
  me = comm->me;

  // Checking errors specific to the compute
  if (dimension == 2)
    error->all(FLERR,"Compute SAED does not work with 2d structures");
  if (narg < 4+ntypes)
    error->all(FLERR,"Illegal Compute SAED Command");
  if (triclinic == 1)
    error->all(FLERR,"Compute SAED does not work with triclinic structures");

  vector_flag = 1;
  extvector = 0;

  // Store radiation wavelength
  lambda = utils::numeric(FLERR,arg[3],false,lmp);
  if (lambda < 0)
    error->all(FLERR,"Compute SAED: Wavelength must be greater than zero");

  // Define atom types for atomic scattering factor coefficients
  int iarg = 4;
  ztype = new int[ntypes];
  for (int i = 0; i < ntypes; i++) {
    ztype[i] = SAEDmaxType + 1;
  }
  for (int i=0; i<ntypes; i++) {
     for (int j = 0; j < SAEDmaxType; j++) {
       if (utils::lowercase(arg[iarg]) == utils::lowercase(SAEDtypeList[j])) {
         ztype[i] = j;
       }
     }
     if (ztype[i] == SAEDmaxType + 1)
       error->all(FLERR,"Compute SAED: Invalid ASF atom type");
    iarg++;
  }

  // Set defaults for optional args
  Kmax = 1.70;
  Zone[0] = 1; Zone[1] = 0; Zone[2] = 0;
  c[0] = 1; c[1] = 1; c[2] = 1;
  dR_Ewald = 0.01 / 2;
  manual = false;
  double manual_double=0;
  echo = false;

  // Process optional args
  while (iarg < narg) {

    if (strcmp(arg[iarg],"Kmax") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal Compute SAED Command");
      Kmax = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      if (Kmax / 2 < 0 || Kmax / 2 > 6)
        error->all(FLERR,"Compute SAED: |K|max/2 must be between 0 and 6 ");
      iarg += 2;

    } else if (strcmp(arg[iarg],"Zone") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Illegal Compute SAED Command");
      Zone[0] = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      Zone[1] = utils::numeric(FLERR,arg[iarg+2],false,lmp);
      Zone[2] = utils::numeric(FLERR,arg[iarg+3],false,lmp);
      iarg += 4;

    } else if (strcmp(arg[iarg],"c") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Illegal Compute SAED Command");
      c[0] = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      c[1] = utils::numeric(FLERR,arg[iarg+2],false,lmp);
      c[2] = utils::numeric(FLERR,arg[iarg+3],false,lmp);
      if (c[0] < 0 || c[1] < 0 || c[2] < 0)
        error->all(FLERR,"Compute SAED: dKs must be greater than 0");
      iarg += 4;

    } else if (strcmp(arg[iarg],"dR_Ewald") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal Compute SAED Command");
      dR_Ewald = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      if (dR_Ewald < 0)
        error->all(FLERR,"Compute SAED: dR_Ewald slice must be greater than 0");
      iarg += 2;

    } else if (strcmp(arg[iarg],"echo") == 0) {
      echo = true;
      iarg += 1;

    } else if (strcmp(arg[iarg],"manual") == 0) {
      manual = true;
      manual_double = 1;
      iarg += 1;

    } else error->all(FLERR,"Illegal Compute SAED Command");
  }

  // Zone flag to capture entire recrocal space volume
  if ( (Zone[0] == 0) && (Zone[1] == 0) && (Zone[2] == 0)) {
  } else {
      R_Ewald = (1 / lambda);
      double Rnorm = R_Ewald/ sqrt(Zone[0] * Zone[0] +
                     Zone[1] * Zone[1] +  Zone[2]* Zone[2]);
      Zone[0] = Zone[0] * Rnorm;
      Zone[1] = Zone[1] * Rnorm;
      Zone[2] = Zone[2] * Rnorm;
  }

  // Procedure to determine how many rows are needed given the constraints on 2theta
  // Calculating spacing between reciprical lattice points
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
  int n = 0;
  double dinv2, r2, EmdR2, EpdR2;
  double K[3];

  // Zone flag to capture entire reciprocal space volume
  if ((Zone[0] == 0) && (Zone[1] == 0) && (Zone[2] == 0)) {
    for (int k = -Knmax[2]; k <= Knmax[2]; k++) {
      for (int j = -Knmax[1]; j <= Knmax[1]; j++) {
        for (int i = -Knmax[0]; i <= Knmax[0]; i++) {
          K[0] = i * dK[0];
          K[1] = j * dK[1];
          K[2] = k * dK[2];
          dinv2 = (K[0] * K[0] + K[1] * K[1] + K[2] * K[2]);
          if (dinv2 < Kmax * Kmax) n++;
        }
      }
    }
  } else {
    for (int k = -Knmax[2]; k <= Knmax[2]; k++) {
      for (int j = -Knmax[1]; j <= Knmax[1]; j++) {
        for (int i = -Knmax[0]; i <= Knmax[0]; i++) {
          K[0] = i * dK[0];
          K[1] = j * dK[1];
          K[2] = k * dK[2];
          dinv2 = (K[0] * K[0] + K[1] * K[1] + K[2] * K[2]);
          if (dinv2 < Kmax * Kmax) {
            r2 = 0.0;
            for (int m=0; m<3; m++)
              r2 += pow(K[m] - Zone[m],2.0);
            EmdR2 = pow(R_Ewald - dR_Ewald,2);
            EpdR2 = pow(R_Ewald + dR_Ewald,2);
            if ((r2 > EmdR2) && (r2 < EpdR2)) {
              n++;
            }
          }
        }
      }
    }
  }

  if (me == 0 && echo)
    utils::logmesg(lmp,"-----\nCompute SAED id:{}, # of atoms:{}, # of relp:{}\n"
                   "Reciprocal point spacing in k1,k2,k3 = {:.8} {:.8} {:.8}\n-----\n",
                   id,natoms,n,dK[0],dK[1],dK[2]);

  nRows = n;
  size_vector = n;
  memory->create(vector,size_vector,"saed:vector");
  memory->create(store_tmp,3*size_vector,"saed:store_tmp");


  // Create vector of variables to be passed to fix ave/time/saed
  saed_var[0] = lambda;
  saed_var[1] = Kmax;
  saed_var[2] = Zone[0];
  saed_var[3] = Zone[1];
  saed_var[4] = Zone[2];
  saed_var[5] = c[0];
  saed_var[6] = c[1];
  saed_var[7] = c[2];
  saed_var[8] = dR_Ewald;
  saed_var[9] = manual_double;
}

/* ---------------------------------------------------------------------- */

ComputeSAED::~ComputeSAED()
{
  memory->destroy(vector);
  memory->destroy(store_tmp);
  delete[] ztype;
}

/* ---------------------------------------------------------------------- */

void ComputeSAED::init()
{
  double dinv2, r2, EmdR2, EpdR2;
  double K[3];
  int n = 0;

  // Zone 0 0 0 flag to capture entire recrocal space volume
  if ((Zone[0] == 0) && (Zone[1] == 0) && (Zone[2] == 0)) {
    for (int k = -Knmax[2]; k <= Knmax[2]; k++) {
      for (int j = -Knmax[1]; j <= Knmax[1]; j++) {
        for (int i = -Knmax[0]; i <= Knmax[0]; i++) {
          K[0] = i * dK[0];
          K[1] = j * dK[1];
          K[2] = k * dK[2];
          dinv2 = (K[0] * K[0] + K[1] * K[1] + K[2] * K[2]);
          if (dinv2 < Kmax * Kmax) {
            store_tmp[3*n] = i;
            store_tmp[3*n+1] = j;
            store_tmp[3*n+2] = k;
            n++;
          }
        }
      }
    }
  } else {
    for (int k = -Knmax[2]; k <= Knmax[2]; k++) {
      for (int j = -Knmax[1]; j <= Knmax[1]; j++) {
        for (int i = -Knmax[0]; i <= Knmax[0]; i++) {
          K[0] = i * dK[0];
          K[1] = j * dK[1];
          K[2] = k * dK[2];
          dinv2 = (K[0] * K[0] + K[1] * K[1] + K[2] * K[2]);
          if (dinv2 < Kmax * Kmax) {
            r2=0.0;
            for (int m=0; m<3; m++)
              r2 += pow(K[m] - Zone[m],2.0);
            EmdR2 = pow(R_Ewald - dR_Ewald,2);
            EpdR2 = pow(R_Ewald + dR_Ewald,2);
            if ((r2 > EmdR2) && (r2 < EpdR2)) {
              store_tmp[3*n] = i;
              store_tmp[3*n+1] = j;
              store_tmp[3*n+2] = k;
              n++;
            }
          }
        }
      }
    }
  }
  if (n != nRows)  error->all(FLERR,"Compute SAED Nrows inconsistent");

}

/* ---------------------------------------------------------------------- */

void ComputeSAED::compute_vector()
{
  invoked_vector = update->ntimestep;

  if (me == 0 && echo)
    utils::logmesg(lmp,"-----\nComputing SAED intensities");

  double t0 = platform::walltime();
  auto Fvec = new double[2*nRows]; // Strct factor (real & imaginary)
  // -- Note, vector entries correspond to different RELP

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

/*
  double *x = new double [3*nlocal];
  int nlocalgroup = 0;
  for (int ii = 0; ii < nlocal; ii++) {
    if (mask[ii] & groupbit) {
     x[3*ii+0] = atom->x[ii][0];
     x[3*ii+1] = atom->x[ii][1];
     x[3*ii+2] = atom->x[ii][2];
     nlocalgroup++;
    }
  }
*/


 // determining parameter set to use based on maximum S = sin(theta)/lambda
  double Smax = Kmax / 2;

  int offset = 0;                 // offset the ASFSAED matrix for appropriate value
  if (Smax <= 2) offset = 0;
  if (Smax > 2)  offset = 10;

  // Setting up OMP
#if defined(_OPENMP)
  if (me == 0 && echo) utils::logmesg(lmp," using {}OMP threads\n",comm->nthreads);
#endif

  if (me == 0 && echo) utils::logmesg(lmp,"\n");

  int m = 0;
  double frac = 0.1;

#if defined(_OPENMP)
#pragma omp parallel LMP_DEFAULT_NONE LMP_SHARED(offset,ASFSAED,typelocal,xlocal,Fvec,m,frac)
#endif
  {
    auto f = new double[ntypes];    // atomic structure factor by type
    int typei = 0;
    double Fatom1 = 0.0;               // structure factor per atom
    double Fatom2 = 0.0;               // structure factor per atom (imaginary)
    double K[3];
    double dinv2 = 0.0;
    double dinv  = 0.0;
    double SinTheta_lambda = 0.0;
    double inners = 0.0;

#if defined(_OPENMP)
#pragma omp for
#endif
    for (int n = 0; n < nRows; n++) {
      int i = store_tmp[3*n+0];
      int j = store_tmp[3*n+1];
      int k = store_tmp[3*n+2];
      K[0] = i * dK[0];
      K[1] = j * dK[1];
      K[2] = k * dK[2];

      dinv2 = (K[0] * K[0] + K[1] * K[1] + K[2] * K[2]);
      dinv = sqrt(dinv2);
      SinTheta_lambda = 0.5*dinv;

      Fatom1 = 0.0;
      Fatom2 = 0.0;

      // Calculate the atomic structure factor by type
      // determining parameter set to use based on S = sin(theta)/lambda <> 2
      for (int ii = 0; ii < ntypes; ii++) {
        f[ii] = 0;
        for (int C = 0; C < 5; C++) {
          int D = C + offset;
          f[ii] += ASFSAED[ztype[ii]][D] * exp(-1*ASFSAED[ztype[ii]][5+D] * SinTheta_lambda * SinTheta_lambda);
        }
      }

      // Evaluate the structure factor equation -- looping over all atoms
      for (int ii = 0; ii < nlocalgroup; ii++) {
        typei=typelocal[ii]-1;
        inners = 2 * MY_PI * (K[0] * xlocal[3*ii+0] + K[1] * xlocal[3*ii+1] +
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
      // TODO use VMD timer style incrementer
#endif
        {
          if (m == round(frac * nRows)) {
            if (me == 0) utils::logmesg(lmp," {:2.0f}% -",frac*100);
            frac += 0.1;
          }
          m++;
        }
      }
    } // End of pragma omp for region
    delete [] f;
  }

  auto scratch = new double[2*nRows];

  // Sum intensity for each ang-hkl combination across processors
  MPI_Allreduce(Fvec,scratch,2*nRows,MPI_DOUBLE,MPI_SUM,world);

  for (int i = 0; i < nRows; i++) {
    vector[i] = (scratch[2*i] * scratch[2*i] + scratch[2*i+1] * scratch[2*i+1]) / natoms;
  }

  double t2 = platform::walltime();

  // compute memory usage per processor
  double bytes = memory_usage();

  if (me == 0 && echo)
    utils::logmesg(lmp," 100% \nTime elapsed during compute_saed = {:.2f} sec "
                   "using {:.2f} Mbytes/processor\n-----\n", t2-t0, bytes/1024.0/1024.0);

  delete [] xlocal;
  delete [] typelocal;
  delete [] scratch;
  delete [] Fvec;
}

/* ----------------------------------------------------------------------
 memory usage of arrays
 ------------------------------------------------------------------------- */

double ComputeSAED::memory_usage()
{
  double bytes = nRows * sizeof(double); //vector
  bytes += (double) 4.0 * nRows * sizeof(double); //Fvec1 & 2, scratch1 & 2
  bytes += (double)ntypes * sizeof(double); // f
  bytes += (double)3.0 * nlocalgroup * sizeof(double); // xlocal
  bytes += (double)nlocalgroup * sizeof(int); // typelocal
  bytes += (double)3.0 * nRows * sizeof(int); // store_temp

  return bytes;
}


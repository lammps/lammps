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
using MathConst::MY_PI;

static const char cite_compute_xrd_c[] =
  "compute xrd command: https://doi.org/10.1088/0965-0393/21/5/055020\n\n"
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
  nlocalgroup = 0;
  if (lmp->citeme) lmp->citeme->add(cite_compute_xrd_c);

  ntypes = atom->ntypes;
  bigint natoms = group->count(igroup);
  int dimension = domain->dimension;
  triclinic = domain->triclinic;
  me = comm->me;

  // Checking errors
  if (dimension == 2)
     error->all(FLERR,"Compute XRD does not work with 2d structures");
  if (narg < 4+ntypes)
     error->all(FLERR,"Compute XRD: expected a wavelength and the chemical symbol of each of "
                "the {} atom types",ntypes);

  array_flag = 1;
  extarray = 0;

  // Store radiation wavelength
  lambda = utils::numeric(FLERR,arg[3],false,lmp);
  if (lambda <= 0)
    error->all(FLERR,"Compute XRD: Wavelength must be greater than zero");

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
        break;
       }
     }
    if (ztype[i] == XRDmaxType + 1) error->all(FLERR,"Compute XRD: Invalid ASF atom type {}", arg[iarg]);
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

  nunclaimed = 0;
  unclaimed = new int[narg];

  // Process optional args
  while (iarg < narg) {
    if (strcmp(arg[iarg],"2Theta") == 0) {
      if (iarg+3 > narg) utils::missing_cmd_args(FLERR,"compute xrd 2Theta",error);
      Min2Theta = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      Max2Theta = utils::numeric(FLERR,arg[iarg+2],false,lmp);
      iarg += 3;

    } else if (strcmp(arg[iarg],"c") == 0) {
      if (iarg+4 > narg) utils::missing_cmd_args(FLERR,"compute xrd c",error);
      c[0] = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      c[1] = utils::numeric(FLERR,arg[iarg+2],false,lmp);
      c[2] = utils::numeric(FLERR,arg[iarg+3],false,lmp);
      if ((c[0] <= 0.0) || (c[1] <= 0.0) || (c[2] <= 0.0))
        error->all(FLERR,"Compute XRD: c's must be greater than 0");
      iarg += 4;

    } else if (strcmp(arg[iarg],"LP") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR,"compute xrd LP",error);
      LP = utils::inumeric(FLERR,arg[iarg+1],false,lmp);

      if (LP != 1 && LP != 0)
         error->all(FLERR,"Compute XRD: LP must have value of 0 or 1");
      iarg += 2;

    } else if (strcmp(arg[iarg],"echo") == 0) {
      echo = true;
      iarg += 1;

    } else if (strcmp(arg[iarg],"manual") == 0) {
      manual = true;
      iarg += 1;

    // a word this style does not know may be a keyword of a style derived from
    // it.  record it and step over it: the values of such a keyword are
    // recorded in the same way, and no keyword of this style can be mistaken
    // for one of them, since they are all words rather than numbers.

    } else unclaimed[nunclaimed++] = iarg++;
  }

  // compute xrd/fft claims its own keywords in its constructor and rejects what
  // is left there.  nothing else derives from this class.

  if (strcmp(style,"xrd") == 0) reject_unclaimed(arg);

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
  // the boundary conditions are validated in set_spacing(), which is called
  // again whenever the box changes

  set_spacing();

  // which reciprocal lattice nodes are explored is fixed here and never
  // changes, because the number of rows of the output array cannot change.
  // rlv itself follows the box from here on.

  for (int i=0; i<6; i++) h_last[i] = domain->h[i];
  warned_range = 0;

  // index i of a node satisfies i = K.a_i/c_i, with a_i the box edge vector, so
  // |i| <= Kmax*|a_i|/c_i bounds the search.  for an orthogonal cell that is
  // the familiar Kmax/dK.

  if (!triclinic || manual) {
    for (int i=0; i<3; i++) Knmax[i] = (int) ceil(Kmax / dK[i]);
  } else {
    double *h = domain->h;
    double alen[3];
    alen[0] = h[0];
    alen[1] = sqrt(h[5]*h[5] + h[1]*h[1]);
    alen[2] = sqrt(h[4]*h[4] + h[3]*h[3] + h[2]*h[2]);
    for (int i=0; i<3; i++) Knmax[i] = (int) ceil(Kmax * alen[i] / c[i]);
  }

  // Finding the intersection of the reciprocal space and Ewald sphere
  bigint nRows = 0;
  double dinv2= 0.0;
  double ang = 0.0;
  double K[3];

  // Procedure to determine how many rows are needed given the constraints on 2theta
  for (int i = -Knmax[0]; i <= Knmax[0]; i++) {
    for (int j = -Knmax[1]; j <= Knmax[1]; j++) {
      for (int k = -Knmax[2]; k <= Knmax[2]; k++) {
        kvector(rlv,i,j,k,K);
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


  // store_tmp holds three integers per row, so that product must stay in range too

  if (3*nRows > MAXSMALLINT)
    error->all(FLERR,"Compute XRD: too many reciprocal lattice nodes ({}); reduce the "
               "2Theta range or increase the c values",nRows);

  size_array_rows = nRows;
  size_array_cols = 2;

  if (me == 0 && echo)
    utils::logmesg(lmp,"-----\nCompute XRD id:{}, # of atoms:{}, # of relp:{}\n"
                   "Reciprocal point spacing in k1,k2,k3 = {:.8} {:.8} {:.8}\n-----\n",
                   id,natoms,nRows,dK[0],dK[1],dK[2]);

  memory->create(array,size_array_rows,size_array_cols,"xrd:array");
  memory->create(store_tmp,3*size_array_rows,"xrd:store_tmp");

  // record which nodes those are.  the set never changes, because the number
  // of rows of the output array cannot, so this scan of the search box is done
  // once here rather than again at the start of every run: for a large cell it
  // costs more than the diffraction calculation itself.

  bigint nfill = 0;
  double convf = 360 / MY_PI;
  if (radflag == 1) convf = 2;

  for (int i = -Knmax[0]; i <= Knmax[0]; i++) {
    for (int j = -Knmax[1]; j <= Knmax[1]; j++) {
      for (int k = -Knmax[2]; k <= Knmax[2]; k++) {
        kvector(rlv,i,j,k,K);
        dinv2 = (K[0] * K[0] + K[1] * K[1] + K[2] * K[2]);
        if  (4 >= dinv2 * lambda * lambda) {
          ang = asin(lambda * sqrt(dinv2) * 0.5);
          if ((ang <= Max2Theta) && (ang >= Min2Theta)) {
            store_tmp[3*nfill] = k;
            store_tmp[3*nfill+1] = j;
            store_tmp[3*nfill+2] = i;
            array[nfill][0] = ang * convf;
            nfill++;
          }
        }
      }
    }
  }

  if (nfill != nRows) error->all(FLERR,"Compute XRD rows mismatch");
}

/* ---------------------------------------------------------------------- */

void ComputeXRD::reject_unclaimed(char **arg)
{
  if (nunclaimed)
    error->all(FLERR,"Unknown compute {} keyword: {}",style,arg[unclaimed[0]]);

  delete[] unclaimed;
  unclaimed = nullptr;
}

/* ---------------------------------------------------------------------- */

ComputeXRD::~ComputeXRD()
{
  // a KOKKOS style derived from this one is copied by value to be handed to a
  // device kernel, and the copy must not free what the original owns

  if (copymode) return;

  memory->destroy(array);
  memory->destroy(store_tmp);
  delete[] ztype;
  delete[] unclaimed;
}

/* ---------------------------------------------------------------------- */

void ComputeXRD::init()
{
  // the set of reciprocal lattice nodes was fixed when the compute was defined
  // and cannot change, so all that is needed here is to pick up a box change
  // that happened before this run started and refresh the angles from it

  set_spacing();
  for (int i = 0; i < 6; i++) h_last[i] = domain->h[i];
  refresh_angles();
}

/* ----------------------------------------------------------------------
   spacing of the reciprocal lattice nodes for the current box
------------------------------------------------------------------------- */

void ComputeXRD::set_spacing()
{
  // the cell may have acquired tilt, or a boundary may have changed, since the
  // compute was defined, so both are picked up here rather than remembered

  triclinic = domain->triclinic;

  if (manual) {
    for (int i = 0; i < 3; i++) prd_inv[i] = 1.0;

  } else {
    int *periodicity = domain->periodicity;
    double *prd = domain->prd;
    double ave_inv = 0.0;

    if (!periodicity[0] && !periodicity[1] && !periodicity[2])
      error->all(FLERR,"Compute XRD must have at least one periodic boundary unless "
                 "manual spacing is specified");

    // a non-periodic direction has no reciprocal lattice vector of its own and
    // is given the average of the periodic ones, which is only meaningful when
    // the reciprocal lattice is aligned with the coordinate axes

    if (triclinic && (!periodicity[0] || !periodicity[1] || !periodicity[2]))
      error->all(FLERR,"Compute XRD with a triclinic cell requires all boundaries to be "
                 "periodic unless manual spacing is specified");

    for (int i = 0; i < 3; i++)
      if (periodicity[i]) {
        prd_inv[i] = 1.0 / prd[i];
        ave_inv += prd_inv[i];
      }

    // Using the average inverse dimensions for non-periodic direction
    ave_inv = ave_inv / (periodicity[0] + periodicity[1] + periodicity[2]);
    for (int i = 0; i < 3; i++)
      if (!periodicity[i]) prd_inv[i] = ave_inv;
  }

  for (int i = 0; i < 3; i++) dK[i] = prd_inv[i]*c[i];

  // reciprocal lattice vectors of the cell, scaled by the c parameters.  the
  // rows of the inverse box matrix are the reciprocal lattice vectors, and they
  // are exactly zero off the diagonal for an orthogonal cell, so this reduces
  // to the diagonal spacing without changing a single bit.  manual spacing asks
  // for a fixed grid in absolute units, so it stays aligned with the axes.

  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) rlv[i][j] = 0.0;

  rlv[0][0] = dK[0];
  rlv[1][1] = dK[1];
  rlv[2][2] = dK[2];

  if (triclinic && !manual) {
    double *h_inv = domain->h_inv;
    rlv[0][1] = c[0]*h_inv[5];
    rlv[0][2] = c[0]*h_inv[4];
    rlv[1][2] = c[1]*h_inv[3];
  }
}

/* ----------------------------------------------------------------------
   rescale the reciprocal lattice if the box has changed since it was last
   set up, following the same approach as the kspace styles: the set of nodes
   is fixed, but the reciprocal lattice vectors are scaled with the cell.
   returns 1 if anything changed.
------------------------------------------------------------------------- */

int ComputeXRD::update_reciprocal()
{
  if (manual) return 0;

  // compare the whole box matrix, not only its diagonal: a pure shear changes
  // the tilt factors and leaves the three edge lengths untouched, and it
  // changes the reciprocal lattice just as a resize does

  int changed = 0;
  for (int i = 0; i < 6; i++)
    if (domain->h[i] != h_last[i]) changed = 1;
  if (!changed) return 0;

  set_spacing();
  for (int i = 0; i < 6; i++) h_last[i] = domain->h[i];
  refresh_angles();
  return 1;
}

/* ----------------------------------------------------------------------
   recompute the diffraction angle of every node from the current spacing.
   a node that has moved beyond the limit of the Ewald sphere no longer
   diffracts; it keeps its row but is flagged with an angle of zero, which is
   outside any valid 2Theta range.
------------------------------------------------------------------------- */

void ComputeXRD::refresh_angles()
{
  double convf = 360 / MY_PI;
  if (radflag == 1) convf = 2;

  bigint noutside = 0;

  for (int n = 0; n < size_array_rows; n++) {
    double K[3];
    kvector(rlv,store_tmp[3*n+2],store_tmp[3*n+1],store_tmp[3*n],K);
    double dinv2 = K[0]*K[0] + K[1]*K[1] + K[2]*K[2];

    if (4 >= dinv2 * lambda * lambda) {
      double ang = asin(lambda * sqrt(dinv2) * 0.5);
      array[n][0] = ang * convf;
      if ((ang > Max2Theta) || (ang < Min2Theta)) noutside++;
    } else {
      array[n][0] = 0.0;
      noutside++;
    }
  }

  // the node set was fixed when the compute was defined, so a large change of
  // box size moves part of it out of the requested range

  if (!warned_range && (100*noutside > size_array_rows) && (comm->me == 0)) {
    warned_range = 1;
    error->warning(FLERR,"{:.3}% of the reciprocal lattice nodes of compute {} have moved "
                   "outside its 2Theta range as the box changed.  The node set is fixed when "
                   "the compute is defined; define it at a representative box size, or use "
                   "manual spacing",100.0*noutside/size_array_rows,id);
  }
}

/* ---------------------------------------------------------------------- */

void ComputeXRD::compute_array()
{
  invoked_array = update->ntimestep;

  update_reciprocal();

  if (me == 0 && echo) utils::logmesg(lmp, "-----\nComputing XRD intensities");

  double t0 = platform::walltime();

  auto *Fvec = new double[2*size_array_rows]; // Strct factor (real & imaginary)
  // -- Note: array rows correspond to different RELP

  ntypes = atom->ntypes;
  int nlocal = atom->nlocal;
  int *type  = atom->type;
  bigint natoms = group->count(igroup);
  if (natoms == 0) natoms = 1;    // an empty group scatters nothing, rather than NaN
  int *mask = atom->mask;

  nlocalgroup = 0;
  for (int ii = 0; ii < nlocal; ii++) {
    if (mask[ii] & groupbit) {
     nlocalgroup++;
    }
  }

  auto *xlocal = new double[3*nlocalgroup];
  auto *typelocal = new int[nlocalgroup];

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
  if ((me == 0) && echo) utils::logmesg(lmp," using {} OMP thread(s)\n",comm->nthreads);
#endif

  if ((me == 0) && echo) {
    utils::logmesg(lmp,"\n");

    if (LP == 1)
        utils::logmesg(lmp,"Applying Lorentz-Polarization Factor During XRD Calculation\n");
  }

  int m = 0;
  double frac = 0.1;

#if defined(_OPENMP)
#pragma omp parallel LMP_DEFAULT_NONE LMP_SHARED(typelocal,xlocal,Fvec,m,frac,ASFXRD)
#endif
  {
    auto *f = new double[ntypes];    // atomic structure factor by type
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
        kvector(rlv,i,j,k,K);

        dinv2 = (K[0] * K[0] + K[1] * K[1] + K[2] * K[2]);
        dinv = sqrt(dinv2);
        SinTheta_lambda = 0.5*dinv;
        SinTheta = SinTheta_lambda * lambda;

        // a node driven past the limit of the Ewald sphere by a change of box
        // size no longer diffracts and has no diffraction angle.  it keeps its
        // row, with zero intensity, just as refresh_angles() gives it an angle
        // of zero, which lies outside any valid 2Theta range

        if (SinTheta > 1.0) sqrt_lp = 0.0;
        else {
          ang = asin( SinTheta );
          Cos2Theta = cos( 2 * ang);
          CosTheta = cos( ang );
          sqrt_lp = sqrt( (1 + Cos2Theta * Cos2Theta) /
               ( CosTheta * SinTheta * SinTheta) );
        }

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
        kvector(rlv,i,j,k,K);

        dinv2 = (K[0] * K[0] + K[1] * K[1] + K[2] * K[2]);
        dinv = sqrt(dinv2);
        SinTheta_lambda = 0.5*dinv;

        // see the branch above: a node past the limit of the Ewald sphere no
        // longer diffracts

        double scale = 1.0;
        if (SinTheta_lambda * lambda > 1.0) scale = 0.0;

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
        Fvec[2*n] = Fatom1 * scale;
        Fvec[2*n+1] = Fatom2 * scale;

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
    delete[] f;
  } // End of pragma omp parallel region

  auto *scratch = new double[2*size_array_rows];

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

  delete[] scratch;
  delete[] Fvec;
  delete[] xlocal;
  delete[] typelocal;
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

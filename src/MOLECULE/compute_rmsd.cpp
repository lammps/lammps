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

#include "compute_rmsd.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "group.h"
#include "math_eigen_impl.h"
#include "math_extra.h"
#include "math_special.h"
#include "memory.h"
#include "text_file_reader.h"
#include "universe.h"
#include "update.h"
#include "utils.h"

#include <cmath>

#include <iostream>


using namespace LAMMPS_NS;
using namespace MathExtra;
using namespace MathSpecial;

/* ---------------------------------------------------------------------- */

ComputeRmsd::ComputeRmsd(LAMMPS *lmp, int narg, char **arg) :
    Compute(lmp, narg, arg)
{

  scalar_flag = 1;
  extscalar = 0;

  if (narg < 3) utils::missing_cmd_args(FLERR, "compute rmsd", error);

  // xyz_file
  group_count = group->count(igroup);
  read_xyz(arg[3]);

  MathEigen::Alloc2D(group_count, 3, &x_group);
  MathEigen::Alloc2D(group_count, 3, &x_group_shifted);
  MathEigen::Alloc2D(group_count, 3, &ref_positions_shifted);
  memory->create(group_taglist,group_count,"compute_rmsd:group_taglist");
}

/* ---------------------------------------------------------------------- */

ComputeRmsd::~ComputeRmsd()
{
  if (copymode) return;
  MathEigen::Dealloc2D(&x_group);
  MathEigen::Dealloc2D(&x_group_shifted);
  MathEigen::Dealloc2D(&ref_positions_shifted);
  memory->destroy(ref_positions);
  memory->destroy(group_taglist);
}

/* ---------------------------------------------------------------------- */

void ComputeRmsd::init()
{

  int *mask = atom->mask;
  for( int i=0, j=0 ; i < atom->nlocal ; i++ )
    if (mask[i] & groupbit) {
      group_taglist[j] = atom->tag[i];
      //std::cerr << fmt::format("group_taglist[{}] {}\n", j,group_taglist[j]);
      j++;
    }

  utils::merge_sort(group_taglist, group_count, nullptr, idcompare);

}

int ComputeRmsd::idcompare(const int i, const int j, void *ptr)
{
  if (i < j) return -1;
  else if (i > j) return 1;
  else return 0;
}

/* ---------------------------------------------------------------------- */

double ComputeRmsd::compute_scalar()
{

  // Skip if already calculated on this timestep
  if (invoked_scalar == update->ntimestep) return scalar;

  invoked_scalar = update->ntimestep;
  double inverse_quat[4];
  double scalar_local = rmsd(inverse_quat);
  MPI_Allreduce(&scalar_local, &scalar, 1, MPI_DOUBLE, MPI_SUM, world);
  return scalar;

}


/* ---------------------------------------------------------------------- */

using std::fdim;
using std::sqrt;

double ComputeRmsd::rmsd( double inverse_quat[4] )
{

  // Find the center-of-geometry of each object:
  double **x = atom->x;
  imageint *image = atom->image;
  double aCenter_f[3] = {0.0};
  double aCenter_m[3] = {0.0};
  for (int n = 0; n < group_count; n++) {
    const int i = atom->map(group_taglist[n]);
    domain->unmap(x[i],image[i],x_group[n]);
    //std::cerr << fmt::format(" *** BEFORE_CG x_group[{}] {:.6} {:.6} {:.6}\n", n,x_group[n][0],x_group[n][1],x_group[n][2]);
    for (int d = 0; d < 3; d++) {
      aCenter_f[d] += x_group[n][d];
      aCenter_m[d] += ref_positions[n][d];
    }
  }

  for (int d = 0; d < 3; d++) {
    aCenter_f[d] /= group_count;
    aCenter_m[d] /= group_count;
  }

  //Subtract the centers-of-geometry from the original coordinates for each object
  for (int n = 0; n < group_count; n++) {
    for (int d = 0; d < 3; d++) {
      // shift the coordinates so that the new center of geometry is at the origin
      x_group_shifted[n][d] = x_group[n][d] - aCenter_f[d];
      ref_positions_shifted[n][d] = ref_positions[n][d] - aCenter_m[d];
    }
        //std::cerr << fmt::format(" *** AFTER_CG x_group[{}] {:.6} {:.6} {:.6}\n", n,x_group_shifted[n][0],x_group_shifted[n][1],x_group_shifted[n][2]);
  }

  // Calculate the "M" array from the Diamond paper (equation 16)
  double M[3][3];
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) M[i][j] = 0.0;

  for (size_t n = 0; n < group_count; n++) {
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) { M[i][j] += ref_positions_shifted[n][i] * x_group_shifted[n][j]; }
    }
  }

  //write3(M);

  // Calculate Q (equation 17)
  double traceM = 0.0;
  for (int i = 0; i < 3; i++) traceM += M[i][i];
  double Q[3][3];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      Q[i][j] = M[i][j] + M[j][i];
      if (i == j) Q[i][j] -= 2.0 * traceM;
    }
  }
  
  // Calculate V (equation 18)
  double V[3];
  V[0] = M[1][2] - M[2][1];
  V[1] = M[2][0] - M[0][2];
  V[2] = M[0][1] - M[1][0];

  // Calculate "P" (equation 22)
  // First we must allocate space for the P matrix.  It's not safe to declare:
  // double P[4][4];
  // ...because most matrix solvers expect arrays in pointer-to-pointer format.
  // (a different format).  Below I create a fixed size matrix P in this format.
  double _PF[4 * 4];             // Contiguous 1D array for storing contents of the 2D P array
  double *P[4];                  // This version of P has has ** (pointer-to-pointer) format.
  for (int i = 0; i < 4; i++)    // We must make sure that
    P[i] = &(_PF[4 * i]);        // P[i] points to the appropriate location in memory

  // Now fill the P array
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) P[i][j] = Q[i][j];
  P[0][3] = V[0];
  P[3][0] = V[0];
  P[1][3] = V[1];
  P[3][1] = V[1];
  P[2][3] = V[2];
  P[3][2] = V[2];
  P[3][3] = 0.0;

// The vector "p" contains the optimal rotation (backwards quaternion format)
  double p[4] = {0.0, 0.0, 0.0, 1.0};

  double Evl[4];                 // Store the eigenvalues of P here.
  double *Evc[4];                // Store the eigevectors here. This version has ** format.
  double _Evc[4 * 4];            // Contiguous 1D array for storing contents of "Evc" array
  for (int i = 0; i < 4; i++)    // We must make sure that
    Evc[i] = &(_Evc[4 * i]);     // Evc[i] points to the correct location in memory

  int ierror = MathEigen::Jacobi<double, double *, double **>(4).Diagonalize(P, Evl, Evc);

  if(ierror)
    error->all(FLERR, "compute rmsd: Too many iterations in jacobi diagonalization.\n"
      "This is usually the result of an ill-defined set of atoms for "
      "rotational alignment (RMSD, rotateReference, etc).\n");

  for (int i = 0; i < 4; i++)
    p[i] = Evc[0][i];    //copy eigenvector corresponding to this eigenvalue to p

  // Now normalize p
  double pnorm = 0.0;
  for (int i = 0; i < 4; i++) pnorm += p[i] * p[i];
  pnorm = sqrt(pnorm);
  for (int i = 0; i < 4; i++) p[i] /= pnorm;


  // Note: The "p" variable is not a quaternion in the
  //       conventional sense because its elements
  //       are in the wrong order.  I correct for that here.
  //       "q" is the quaternion correspond to rotation R.
  //       [Andrew Jewett (Scripps Research)]
  //q[0] = p[3];
  //q[1] = p[0];
  //q[2] = p[1];
  //q[3] = p[2];

  // BUT... we actually need the INVERSE ROTATION
  // [alphataubio (2024/08)]

  inverse_quat[0] = p[3];
  inverse_quat[1] = -p[0];
  inverse_quat[2] = -p[1];
  inverse_quat[3] = -p[2];

  // Finally compute the RMSD between the two coordinate sets:
  // First compute E0 from equation 24 of the paper
  double E0 = 0.0;

  for (size_t n = 0; n < group_count; n++) {
    double tmp[3];
    // quaternion rotation of vector: c = a*b*conj(a)
    quatrotvec(inverse_quat, x_group_shifted[n], tmp);
    //std::cerr << fmt::format(" *** AFTER_ROT inverse_quat {:.6} {:.6} {:.6} {:.6} x_group[{}] {:.6} {:.6} {:.6}\n", inverse_quat[0],inverse_quat[1],inverse_quat[2],inverse_quat[3], n,tmp[0],tmp[1],tmp[2]);
    for (int d = 0; d < 3; d++) {
      E0 += (square(x_group_shifted[n][d] - ref_positions_shifted[n][d]));
      x_group[n][d] = tmp[d]+aCenter_m[d];
    }
        //std::cerr << fmt::format(" *** x_group_shifted[{}] {:.6} {:.6} {:.6} ref_positions_shifted[{}] {:.6} {:.6} {:.6}\n", n, x_group_shifted[n][0],x_group_shifted[n][1],x_group_shifted[n][2], n, ref_positions_shifted[n][0],ref_positions_shifted[n][1],ref_positions_shifted[n][2] );

    //std::cerr << fmt::format(" *** AFTER_TRANSLATION x_group[{}] {:.6} {:.6} {:.6}\n", n,x_group[n][0],x_group[n][1],x_group[n][2]);
  }

  return sqrt(fdim(E0, 2.0 * Evl[0])/group_count); // Evl[0] = the maximum eigenvalue of P

}

/* ---------------------------------------------------------------------- */

void ComputeRmsd::read_xyz(char *filename)
{

  memory->create(ref_positions,group_count,3,"ComputeRmsd:ref_positions");
  if (comm->me != 0) return;
  FILE *fp = fopen(filename, "r");

  if (fp == nullptr)
    error->one(FLERR, "ComputeRmsd: Cannot open file {}: {}", filename, utils::getsyserror());

  TextFileReader reader(fp, "xyz");
  try {
    char *line = reader.next_line(1);
    int xyz_count = ValueTokenizer(line).next_int();

    if( xyz_count != group_count )
      error->all(FLERR, "ComputeRmsd: number {} of reference positions in {} does not match the number {} of atoms of group {}", xyz_count, filename, group_count, group->names[igroup]);

    reader.skip_line();
    for( int i=0 ; i<xyz_count ; i++ ) {
      double buffer[4];
      reader.next_dvector(buffer,4);
      ref_positions[i][0] = buffer[1];
      ref_positions[i][1] = buffer[2];
      ref_positions[i][2] = buffer[3];
    }
  } catch (std::exception &e) {
    error->all(FLERR, "ComputeRmsd: error reading xyz file {}: {}", filename, e.what());
  }
  fclose(fp);
}


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
 *    Contributing author:  Evangelos Voyiatzis (Royal DSM)
 * ------------------------------------------------------------------------- */

#include "compute_gyration_shape_chunk.h"

#include "error.h"
#include "math_eigen.h"
#include "math_special.h"
#include "memory.h"
#include "modify.h"
#include "update.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeGyrationShapeChunk::ComputeGyrationShapeChunk(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg), id_gyration_chunk(nullptr), shape_parameters(nullptr)
{
  if (narg != 4) error->all(FLERR,"Illegal compute gyration/shape/chunk command");

  // ID of compute gyration
  id_gyration_chunk = utils::strdup(arg[3]);

  ComputeGyrationShapeChunk::init();

  array_flag = 1;
  size_array_cols = 6;
  size_array_rows = 0;
  size_array_rows_variable = 1;
  extarray = 0;

  firstflag = 1;
  former_nchunks = 0;
  current_nchunks = 1;
  allocate();
}

/* ---------------------------------------------------------------------- */

ComputeGyrationShapeChunk::~ComputeGyrationShapeChunk()
{
  delete [] id_gyration_chunk;
  memory->destroy(shape_parameters);
}

/* ---------------------------------------------------------------------- */

void ComputeGyrationShapeChunk::init()
{
  // check that the compute gyration command exist
  int icompute = modify->find_compute(id_gyration_chunk);
  if (icompute < 0)
    error->all(FLERR,"Compute gyration/chunk ID does not exist for "
               "compute gyration/shape/chunk");

  // check the id_gyration_chunk corresponds really to a compute gyration/chunk command
  c_gyration_chunk = (Compute *) modify->compute[icompute];
  if (strcmp(c_gyration_chunk->style,"gyration/chunk") != 0)
    error->all(FLERR,"Compute gyration/shape/chunk does not point to "
               "gyration compute/chunk");

  // check the compute gyration/chunk command computes the whole gyration tensor
  if (c_gyration_chunk->array_flag == 0)
    error->all(FLERR,"Compute gyration/chunk where gyration/shape/chunk points to "
               "does not calculate the gyration tensor");

}

/* ---------------------------------------------------------------------- */

void ComputeGyrationShapeChunk::setup()
{
  // one-time calculation of per-chunk mass
  // done in setup, so that ComputeChunkAtom::setup() is already called

  if (firstflag) {
    compute_array();
    firstflag = 0;
  }
}

/* ----------------------------------------------------------------------
   compute shape parameters based on the eigenvalues of the
   gyration tensor of group of atoms
------------------------------------------------------------------------- */

void ComputeGyrationShapeChunk::compute_array()
{
  invoked_array = update->ntimestep;
  c_gyration_chunk->compute_array();

  current_nchunks = c_gyration_chunk->size_array_rows; // how to check for the number of chunks in the gyration/chunk?
  if (former_nchunks != current_nchunks) allocate();

  double **gyration_tensor = c_gyration_chunk->array;

  // call the function for the calculation of the eigenvalues
  double ione[3][3], evalues[3], evectors[3][3];

  for (int ichunk = 0; ichunk < current_nchunks; ichunk++) {

    ione[0][0] = gyration_tensor[ichunk][0];
    ione[1][1] = gyration_tensor[ichunk][1];
    ione[2][2] = gyration_tensor[ichunk][2];
    ione[0][1] = ione[1][0] = gyration_tensor[ichunk][3];
    ione[0][2] = ione[2][0] = gyration_tensor[ichunk][4];
    ione[1][2] = ione[2][1] = gyration_tensor[ichunk][5];

    int ierror = MathEigen::jacobi3(ione,evalues,evectors);
    if (ierror) error->all(FLERR, "Insufficient Jacobi rotations "
                         "for gyration/shape");

    // sort the eigenvalues according to their size with bubble sort
    double t;
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 2-i; j++) {
        if (fabs(evalues[j]) < fabs(evalues[j+1])) {
          t = evalues[j];
          evalues[j] = evalues[j+1];
          evalues[j+1] = t;
        }
      }
    }

    // compute the shape parameters of the gyration tensor
    double sq_eigen_x = MathSpecial::square(evalues[0]);
    double sq_eigen_y = MathSpecial::square(evalues[1]);
    double sq_eigen_z = MathSpecial::square(evalues[2]);

    double nominator = sq_eigen_x + sq_eigen_y + sq_eigen_z;
    double denominator = MathSpecial::square(evalues[0]+evalues[1]+evalues[2]);

    shape_parameters[ichunk][0] = evalues[0];
    shape_parameters[ichunk][1] = evalues[1];
    shape_parameters[ichunk][2] = evalues[2];
    shape_parameters[ichunk][3] = evalues[0] - 0.5*(evalues[1] + evalues[2]);
    shape_parameters[ichunk][4] = evalues[1] - evalues[2];
    shape_parameters[ichunk][5] = 1.5*nominator/denominator - 0.5;
  }
}

/* ----------------------------------------------------------------------
 *    calculate and return # of chunks = length of vector/array
 *    ------------------------------------------------------------------------- */

int ComputeGyrationShapeChunk::lock_length()
{
  int number_of_chunks = c_gyration_chunk->size_array_rows;
  return number_of_chunks;
}

/* ----------------------------------------------------------------------
 *    free and reallocate per-chunk arrays
 * ---------------------------------------------------------------------- */

void ComputeGyrationShapeChunk::allocate()
{
  memory->destroy(shape_parameters);
  former_nchunks = current_nchunks;
  memory->create(shape_parameters,current_nchunks,6,"gyration/shape/chunk:shape_parameters");
  array = shape_parameters;
  size_array_rows = current_nchunks;
}

/* ----------------------------------------------------------------------
 *    memory usage of local data
 * ---------------------------------------------------------------------- */

double ComputeGyrationShapeChunk::memory_usage()
{
  double bytes = (bigint) current_nchunks * 6 * sizeof(double);
  return bytes;
}

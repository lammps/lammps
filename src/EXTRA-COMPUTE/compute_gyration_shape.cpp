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
 *    Contributing author:  Evangelos Voyiatzis (Royal DSM)
 * ------------------------------------------------------------------------- */

#include "compute_gyration_shape.h"

#include "error.h"
#include "math_eigen.h"
#include "math_special.h"
#include "modify.h"
#include "update.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeGyrationShape::ComputeGyrationShape(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg), id_gyration(nullptr)
{
  if (narg != 4) error->all(FLERR,"Illegal compute gyration/shape command");

  vector_flag = 1;
  size_vector = 6;
  extscalar = 0;
  extvector = 0;

  // ID of compute gyration
  id_gyration = utils::strdup(arg[3]);

  ComputeGyrationShape::init();

  vector = new double[6];
}

/* ---------------------------------------------------------------------- */

ComputeGyrationShape::~ComputeGyrationShape()
{
  delete [] id_gyration;
  delete [] vector;
}

/* ---------------------------------------------------------------------- */

void ComputeGyrationShape::init()
{
  // check that the compute gyration command exist
  int icompute = modify->find_compute(id_gyration);
  if (icompute < 0)
    error->all(FLERR,"Compute gyration ID does not exist for "
               "compute gyration/shape");

  // check the id_gyration corresponds really to a compute gyration command
  c_gyration = (Compute *) modify->compute[icompute];
  if (strcmp(c_gyration->style,"gyration") != 0)
    error->all(FLERR,"Compute gyration compute ID does not point to "
               "gyration compute for compute gyration/shape");
}

/* ----------------------------------------------------------------------
   compute shape parameters based on the eigenvalues of the
   gyration tensor of group of atoms
------------------------------------------------------------------------- */

void ComputeGyrationShape::compute_vector()
{
  invoked_vector = update->ntimestep;
  c_gyration->compute_vector();
  double *gyration_tensor = c_gyration->vector;

  // call the function for the calculation of the eigenvalues
  double ione[3][3], evalues[3], evectors[3][3];

  ione[0][0] = gyration_tensor[0];
  ione[1][1] = gyration_tensor[1];
  ione[2][2] = gyration_tensor[2];
  ione[0][1] = ione[1][0] = gyration_tensor[3];
  ione[0][2] = ione[2][0] = gyration_tensor[4];
  ione[1][2] = ione[2][1] = gyration_tensor[5];

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
  double nominator = MathSpecial::square(evalues[0])
    + MathSpecial::square(evalues[1])
    + MathSpecial::square(evalues[2]);
  double denominator = MathSpecial::square(evalues[0]+evalues[1]+evalues[2]);

  vector[0] = evalues[0];
  vector[1] = evalues[1];
  vector[2] = evalues[2];
  vector[3] = evalues[0] - 0.5*(evalues[1] + evalues[2]);
  vector[4] = evalues[1] - evalues[2];
  vector[5] = 1.5*nominator/denominator - 0.5;
}

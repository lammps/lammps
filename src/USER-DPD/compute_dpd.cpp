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
   Contributing author: James Larentzos (U.S. Army Research Laboratory)
------------------------------------------------------------------------- */

#include "compute_dpd.h"
#include <mpi.h>
#include "atom.h"
#include "update.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeDpd::ComputeDpd(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg != 3) error->all(FLERR,"Illegal compute dpd command");

  vector_flag = 1;
  size_vector = 5;
  extvector = 0;

  vector = new double[size_vector];

  if (atom->dpd_flag != 1) error->all(FLERR,"compute dpd requires atom_style with internal temperature and energies (e.g. dpd)");
}

/* ---------------------------------------------------------------------- */

ComputeDpd::~ComputeDpd()
{

  delete [] vector;

}

/* ---------------------------------------------------------------------- */

void ComputeDpd::compute_vector()
{
  invoked_vector = update->ntimestep;

  double *uCond = atom->uCond;
  double *uMech = atom->uMech;
  double *uChem = atom->uChem;
  double *dpdTheta = atom->dpdTheta;
  int nlocal = atom->nlocal;
  int *mask = atom->mask;
  int natoms;

  dpdU = new double[size_vector];

  for (int i = 0; i < size_vector; i++) dpdU[i] = 0.0;

  for (int i = 0; i < nlocal; i++){
    if (mask[i] & groupbit){
      dpdU[0] += uCond[i];
      dpdU[1] += uMech[i];
      dpdU[2] += uChem[i];
      dpdU[3] += 1.0 / dpdTheta[i];
      dpdU[4]++;
    }
  }

  MPI_Allreduce(dpdU,vector,size_vector,MPI_DOUBLE,MPI_SUM,world);

  natoms = vector[4];
  vector[3] = natoms / vector[3];

  delete [] dpdU;

}

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

#include "string.h"
#include "compute_tally_stress.h"
#include "atom.h"
#include "pair.h"
#include "update.h"
#include "memory.h"
#include "error.h"
#include "force.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeTallyStress::ComputeTallyStress(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg < 3) error->all(FLERR,"Illegal compute tally/stress command");
  scalar_flag = 1;
  vector_flag = 1;
  peratom_flag = 1;
  timeflag = 1;

  comm_reverse = size_peratom_cols = 6;
  extscalar = 0;
  extvector = 0;
  size_vector = 6;
  peflag = 1;                   // we need Pair::ev_tally() to be run
  nmax = -1;
  stress = NULL;
  vector = new double[size_vector];
}

/* ---------------------------------------------------------------------- */

ComputeTallyStress::~ComputeTallyStress()
{
  if (force->pair) force->pair->del_tally_callback(this);
  delete[] vector;
}

/* ---------------------------------------------------------------------- */

void ComputeTallyStress::init()
{
  if (force->pair == NULL)
    error->all(FLERR,"Trying to use compute tally/stress with no pair style");
  else
    force->pair->add_tally_callback(this);

  did_compute = -1;
}


/* ---------------------------------------------------------------------- */
void ComputeTallyStress::pair_tally_callback(int i, int j, int nlocal, int newton,
                                             double, double, double fpair,
                                             double dx, double dy, double dz)
{
  const int ntotal = atom->nlocal + atom->nghost;

  // do setup work that needs to be done only once per timestep

  if (did_compute != update->ntimestep) {
    did_compute = update->ntimestep;
    
    // grow local stress array if necessary
    // needs to be atom->nmax in length

    if (atom->nmax > nmax) {
      memory->destroy(stress);
      nmax = atom->nmax;
      memory->create(stress,nmax,size_peratom_cols,"tally/stress:stress");
      array_atom = stress;
    }

    // clear storage as needed

    if (newton) {
      for (int i=0; i < ntotal; ++i)
        for (int j=0; j < 6; ++j)
          stress[i][j] = 0.0;
    } else {
      for (int i=0; i < atom->nlocal; ++i)
        for (int j=0; j < 6; ++j)
          stress[i][j] = 0.0;
    }

    for (int i=0; i < 6; ++i)
      vector[i] = virial[i] = 0.0;
  }

  fpair *= 0.5;
  const double v0 = dx*dx*fpair;
  const double v1 = dy*dy*fpair;
  const double v2 = dz*dz*fpair;
  const double v3 = dx*dy*fpair;
  const double v4 = dx*dz*fpair;
  const double v5 = dy*dz*fpair;
  

  if (newton || i < nlocal) {
    virial[0] += v0; stress[i][0] += v0;
    virial[1] += v1; stress[i][1] += v1;
    virial[2] += v2; stress[i][2] += v2;
    virial[3] += v3; stress[i][3] += v3;
    virial[4] += v4; stress[i][4] += v4;
    virial[5] += v5; stress[i][5] += v5;
  }
  if (newton || j < nlocal) {
    virial[0] += v0; stress[j][0] += v0;
    virial[1] += v1; stress[j][1] += v1;
    virial[2] += v2; stress[j][2] += v2;
    virial[3] += v3; stress[j][3] += v3;
    virial[4] += v4; stress[j][4] += v4;
    virial[5] += v5; stress[j][5] += v5;
  }
}

/* ---------------------------------------------------------------------- */

int ComputeTallyStress::pack_reverse_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    buf[m++] = stress[i][0];
    buf[m++] = stress[i][1];
    buf[m++] = stress[i][2];
    buf[m++] = stress[i][3];
    buf[m++] = stress[i][4];
    buf[m++] = stress[i][5];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void ComputeTallyStress::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    stress[j][0] += buf[m++];
    stress[j][1] += buf[m++];
    stress[j][2] += buf[m++];
    stress[j][3] += buf[m++];
    stress[j][4] += buf[m++];
    stress[j][5] += buf[m++];
  }
}

/* ---------------------------------------------------------------------- */

double ComputeTallyStress::compute_scalar()
{
  invoked_scalar = update->ntimestep;
  if (update->eflag_global != invoked_scalar)
    error->all(FLERR,"Energy was not tallied on needed timestep");

  compute_vector();
  scalar = (vector[0]+vector[1]+vector[2])/3.0;
  return scalar;
}

/* ---------------------------------------------------------------------- */

void ComputeTallyStress::compute_vector()
{
  invoked_vector = update->ntimestep;
  if (update->eflag_global != invoked_vector)
    error->all(FLERR,"Energy was not tallied on needed timestep");

  // sum accumulated virial across procs

  MPI_Allreduce(virial,vector,size_vector,MPI_DOUBLE,MPI_SUM,world);

  const double nktv2p = -force->nktv2p;
  for (int i=0; i < 6; ++i)
    vector[i] *= nktv2p;
}

/* ---------------------------------------------------------------------- */

void ComputeTallyStress::compute_peratom()
{
  invoked_peratom = update->ntimestep;
  if (update->eflag_global != invoked_peratom)
    error->all(FLERR,"Energy was not tallied on needed timestep");

  if (force->newton_pair)
    comm->reverse_comm_compute(this);

  // convert to stress*volume units = -pressure*volume

  const double nktv2p = -force->nktv2p;
  for (int i = 0; i < atom->nlocal; i++) {
    stress[i][0] *= nktv2p;
    stress[i][1] *= nktv2p;
    stress[i][2] *= nktv2p;
    stress[i][3] *= nktv2p;
    stress[i][4] *= nktv2p;
    stress[i][5] *= nktv2p;
  }
}


/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputeTallyStress::memory_usage()
{
  double bytes = nmax*6 * sizeof(double);
  return bytes;
}


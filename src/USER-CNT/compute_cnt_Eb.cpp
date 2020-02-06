/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
   
   Contributing author: Maxim Shugaev (UVA), mvs9t@virginia.edu
------------------------------------------------------------------------- */

#include "compute_cnt_Eb.h"
#include "pair_cnt_tpm.h"
#include <cstring>
#include "atom.h"
#include "update.h"
#include "comm.h"
#include "force.h"
#include "bond.h"
#include "modify.h"
#include "memory.h"
#include "error.h"
#include <typeinfo>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeCNT_Eb::ComputeCNT_Eb(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg),
  energy(NULL)
{
  if (narg < 3) error->all(FLERR,"Illegal compute cnt/Eb command");

  peratom_flag = 1;
  size_peratom_cols = 0;
  peatomflag = 1;
  timeflag = 1;
  comm_reverse = 1;
  nmax = 0;
}

/* ---------------------------------------------------------------------- */

ComputeCNT_Eb::~ComputeCNT_Eb()
{
  memory->destroy(energy);
}

/* ---------------------------------------------------------------------- */

void ComputeCNT_Eb::compute_peratom()
{
  int i;

  invoked_peratom = update->ntimestep;
  if (update->eflag_atom != invoked_peratom)
    error->all(FLERR,"Per-atom energy was not tallied on needed timestep");

  // grow local energy array if necessary
  // needs to be atom->nmax in length
  if (atom->nmax > nmax) {
    memory->destroy(energy);
    nmax = atom->nmax;
    memory->create(energy,nmax,"cnt_Eb:energy");
    vector_atom = energy;
  }

  // npair includes ghosts if newton_bond is set
  // ntotal includes ghosts if either newton flag is set
  int nlocal = atom->nlocal;
  int npair = nlocal;
  if (force->newton) npair += atom->nghost;
  int ntotal = nlocal;
  if (force->newton) ntotal += atom->nghost;

  // clear local energy array
  for (i = 0; i < ntotal; i++) energy[i] = 0.0;
  
  if (force->pair){
    try {
      PairCNTTPM* pair = dynamic_cast<PairCNTTPM*>(force->pair);
      for (i = 0; i < npair; i++) energy[i] += pair->eatom_b[i];
    }
    catch (std::bad_cast& bc){
      error->all(FLERR,"cnt/Eb is allowed only with cnt pair style");
    }
  }
  // communicate ghost energy between neighbor procs
  if (force->newton) comm->reverse_comm_compute(this);

  // zero energy of atoms not in group
  // only do this after comm since ghost contributions must be included
  int *mask = atom->mask;
  for (i = 0; i < nlocal; i++)
    if (!(mask[i] & groupbit)) energy[i] = 0.0;
}

/* ---------------------------------------------------------------------- */

int ComputeCNT_Eb::pack_reverse_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) buf[m++] = energy[i];
  return m;
}

/* ---------------------------------------------------------------------- */

void ComputeCNT_Eb::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    energy[j] += buf[m++];
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputeCNT_Eb::memory_usage()
{
  double bytes = nmax * sizeof(double);
  return bytes;
}

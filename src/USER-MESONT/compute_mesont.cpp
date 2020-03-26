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

#include "compute_mesont.h"
#include <cstring>
#include "atom.h"
#include "update.h"
#include "comm.h"
#include "force.h"
#include "modify.h"
#include "memory.h"
#include "error.h"
#include "pair.h"
#include <string>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeMesoNT::ComputeMesoNT(LAMMPS *lmp, int narg, char **arg) :
 Compute(lmp, narg, arg), energy(NULL) {
  if (narg != 4) error->all(FLERR,"Illegal compute mesont command");
  std::string ctype = arg[3];
  if (ctype == "estretch") compute_type = ES;
  else if (ctype == "estretch_tot") compute_type = ESTOT;
  else if (ctype == "ebend") compute_type = EB;
  else if (ctype == "ebend_tot") compute_type = EBTOT;
  else if (ctype == "etube") compute_type = ET;
  else if (ctype == "etube_tot") compute_type = ETTOT;
  else error->all(FLERR,"Illegal compute mesont command");

  if ((compute_type == ES) || (compute_type == EB) || (compute_type == ET)) {
    peratom_flag = 1;
    size_peratom_cols = 0;
    peatomflag = 1;
    timeflag = 1;
    comm_reverse = 1;
    nmax = 0;
  } else {
    timeflag = 1;
    extscalar = 1;
    scalar_flag = 1;
  }
}

/* ---------------------------------------------------------------------- */

ComputeMesoNT::~ComputeMesoNT() {
  memory->destroy(energy);
}

/* ---------------------------------------------------------------------- */

double ComputeMesoNT::compute_scalar() {
  invoked_scalar = update->ntimestep;
  if (update->eflag_global != invoked_scalar)
    error->all(FLERR,"Energy was not tallied on needed timestep");

  int i;
  double* ptr = NULL;
  if (compute_type == ESTOT)
   ptr = static_cast<double*>(force->pair->extract("mesonttpm_Es_tot",i));
  else if (compute_type == EBTOT)
   ptr = static_cast<double*>(force->pair->extract("mesonttpm_Eb_tot",i));
  else if (compute_type == ETTOT)
   ptr = static_cast<double*>(force->pair->extract("mesonttpm_Et_tot",i));
  else error->all(FLERR,"Illegal compute mesont command");

  if (!ptr) error->all(FLERR,
   "compute mesont is allowed only with mesont/tpm pair style");
  MPI_Allreduce(ptr,&scalar,1,MPI_DOUBLE,MPI_SUM,world);

  return scalar;
}

/* ---------------------------------------------------------------------- */

void ComputeMesoNT::compute_peratom() {
  invoked_peratom = update->ntimestep;
  if (update->eflag_atom != invoked_peratom)
    error->all(FLERR,"Per-atom energy was not tallied on needed timestep");

  // grow local energy array if necessary
  // needs to be atom->nmax in length
  if (atom->nmax > nmax) {
    memory->destroy(energy);
    nmax = atom->nmax;
    memory->create(energy,nmax,"mesont_Eb:energy");
    vector_atom = energy;
  }

  // npair includes ghosts if newton_bond is set
  // ntotal includes ghosts if either newton flag is set
  int nlocal = atom->nlocal;
  int npair = nlocal;
  if (force->newton) npair += atom->nghost;
  int ntotal = nlocal;
  if (force->newton) ntotal += atom->nghost;
  int i;
  // clear local energy array
  for (int i = 0; i < ntotal; i++) energy[i] = 0.0;
  double* ptr = NULL;
  if (compute_type == ES)
   ptr = static_cast<double*>(force->pair->extract("mesonttpm_Es",i));
  else if (compute_type == EB)
   ptr = static_cast<double*>(force->pair->extract("mesonttpm_Eb",i));
  else if (compute_type == ET)
   ptr = static_cast<double*>(force->pair->extract("mesonttpm_Et",i));
  else error->all(FLERR,"Illegal compute mesont command");

  if (ptr) for (i = 0; i < npair; i++) energy[i] += ptr[i];
  else error->all(FLERR,
   "compute mesont is allowed only with mesont/tpm pair style");

  // communicate ghost energy between neighbor procs
  if (force->newton) comm->reverse_comm_compute(this);

  // zero energy of atoms not in group
  // only do this after comm since ghost contributions must be included
  int *mask = atom->mask;
  for (int i = 0; i < nlocal; i++)
    if (!(mask[i] & groupbit)) energy[i] = 0.0;
}

/* ---------------------------------------------------------------------- */

int ComputeMesoNT::pack_reverse_comm(int n, int first, double *buf) {
  int m = 0;
  int last = first + n;
  for (int i = first; i < last; i++) buf[m++] = energy[i];
  return m;
}

/* ---------------------------------------------------------------------- */

void ComputeMesoNT::unpack_reverse_comm(int n, int *list, double *buf) {
  int m = 0;
  for (int i = 0; i < n; i++) {
    int j = list[i];
    energy[j] += buf[m++];
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputeMesoNT::memory_usage() {
  double bytes = nmax * sizeof(double);
  return bytes;
}

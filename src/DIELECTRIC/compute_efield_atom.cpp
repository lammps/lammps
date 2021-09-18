// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Trung Nguyen (Northwestern)
------------------------------------------------------------------------- */

#include "compute_efield_atom.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "kspace.h"
#include "memory.h"
#include "msm_dielectric.h"
#include "pair_coul_cut_dielectric.h"
#include "pair_coul_long_dielectric.h"
#include "pair_lj_cut_coul_cut_dielectric.h"
#include "pair_lj_cut_coul_long_dielectric.h"
#include "pair_lj_cut_coul_msm_dielectric.h"
#include "pppm_dielectric.h"
#include "update.h"

#include <cstring>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeEfieldAtom::ComputeEfieldAtom(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg), efield(nullptr)
{
  if (narg < 3) error->all(FLERR,"Illegal compute efield/atom command");

  peratom_flag = 1;
  size_peratom_cols = 3;
  timeflag = 1;
  comm_reverse = 3;

  pairflag = 0;
  kspaceflag = 0;

  if (narg == 3) {
    pairflag = 1;
    kspaceflag = 1;
  } else {
    int iarg = 3;
    while (iarg < narg) {
      if (strcmp(arg[iarg],"pair") == 0) pairflag = 1;
      else if (strcmp(arg[iarg],"kspace") == 0) kspaceflag = 1;
      else error->all(FLERR,"Illegal compute efield/atom command");
      iarg++;
    }
  }

  nmax = 0;

  comm_reverse = 1;
}

/* ---------------------------------------------------------------------- */

ComputeEfieldAtom::~ComputeEfieldAtom()
{
  memory->destroy(efield);
}

/* ---------------------------------------------------------------------- */

void ComputeEfieldAtom::init()
{
  if (!atom->q_flag) error->all(FLERR,"compute efield/atom requires atom attribute q");
  if (!force->kspace) kspaceflag = 0;
}

/* ---------------------------------------------------------------------- */

void ComputeEfieldAtom::setup()
{
  if (strcmp(force->pair_style,"lj/cut/coul/long/dielectric") == 0)
    efield_pair = ((PairLJCutCoulLongDielectric*)force->pair)->efield;
  else if (strcmp(force->pair_style,"lj/cut/coul/long/dielectric/omp") == 0)
    efield_pair = ((PairLJCutCoulMSMDielectric*)force->pair)->efield;
  else if (strcmp(force->pair_style,"lj/cut/coul/msm/dielectric") == 0)
    efield_pair = ((PairLJCutCoulMSMDielectric*)force->pair)->efield;
  else if (strcmp(force->pair_style,"lj/cut/coul/cut/dielectric") == 0)
    efield_pair = ((PairLJCutCoulCutDielectric*)force->pair)->efield;
  else if (strcmp(force->pair_style,"lj/cut/coul/cut/dielectric/omp") == 0)
    efield_pair = ((PairLJCutCoulCutDielectric*)force->pair)->efield;
  else if (strcmp(force->pair_style,"coul/long/dielectric") == 0)
    efield_pair = ((PairCoulLongDielectric*)force->pair)->efield;
  else if (strcmp(force->pair_style,"coul/cut/dielectric") == 0)
    efield_pair = ((PairCoulCutDielectric*)force->pair)->efield;
  else error->all(FLERR,"Compute efield/atom not supported by pair style");

  if (force->kspace) {
    if (strcmp(force->kspace_style,"pppm/dielectric") == 0)
      efield_kspace = ((PPPMDielectric*)force->kspace)->efield;
    else if (strcmp(force->kspace_style,"msm/dielectric") == 0)
      efield_kspace = ((MSMDielectric*)force->kspace)->efield;
    else error->all(FLERR,"Compute efield/atom not supported by kspace style");
    kspaceflag = 1;
  }

  if (!efield_pair && !efield_kspace)
    error->all(FLERR, "Compute efield/atom does not access to efield");
}

/* ---------------------------------------------------------------------- */

void ComputeEfieldAtom::compute_peratom()
{
  int i,j;

  invoked_peratom = update->ntimestep;
  if (update->vflag_atom != invoked_peratom)
    error->all(FLERR,"Per-atom virial was not tallied on needed timestep");

  // grow local stress array if necessary
  // needs to be atom->nmax in length

  if (atom->nmax > nmax) {
    memory->destroy(efield);
    nmax = atom->nmax;
    memory->create(efield,nmax,3,"stress/atom:efield");
    array_atom = efield;
  }

  // npair includes ghosts if either newton flag is set
  //   b/c some bonds/dihedrals call pair::ev_tally with pairwise info
  // nbond includes ghosts if newton_bond is set
  // ntotal includes ghosts if either newton flag is set
  // KSpace includes ghosts if tip4pflag is set

  double* q = atom->q;
  int nlocal = atom->nlocal;
  int npair = nlocal;
  int ntotal = nlocal;
  int nkspace = nlocal;
  if (force->newton) npair += atom->nghost;
  if (force->newton) ntotal += atom->nghost;
  if (force->kspace && force->kspace->tip4pflag) nkspace += atom->nghost;

  // clear local stress array

  for (i = 0; i < ntotal; i++)
    for (j = 0; j < 3; j++)
      efield[i][j] = 0.0;

  // add in per-atom contributions from each force

  if (pairflag && force->pair) {
    for (i = 0; i < npair; i++)
      for (j = 0; j < 3; j++) {
        if (q[i] != 0) efield[i][j] += efield_pair[i][j];
      }
  }

  if (kspaceflag && force->kspace) {
    for (i = 0; i < nkspace; i++)
      for (j = 0; j < 3; j++)
        efield[i][j] += efield_kspace[i][j];
  }

  // communicate ghost efield between neighbor procs

  if (force->newton || (force->kspace && force->kspace->tip4pflag))
    comm->reverse_comm_compute(this);

  // zero efield of atoms not in group
  // only do this after comm since ghost contributions must be included

  int *mask = atom->mask;

  for (i = 0; i < nlocal; i++)
    if (!(mask[i] & groupbit)) {
      efield[i][0] = 0.0;
      efield[i][1] = 0.0;
      efield[i][2] = 0.0;
    }
}


/* ---------------------------------------------------------------------- */

int ComputeEfieldAtom::pack_reverse_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    buf[m++] = efield[i][0];
    buf[m++] = efield[i][1];
    buf[m++] = efield[i][2];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void ComputeEfieldAtom::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    efield[j][0] += buf[m++];
    efield[j][1] += buf[m++];
    efield[j][2] += buf[m++];
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputeEfieldAtom::memory_usage()
{
  double bytes = nmax*3 * sizeof(double);
  return bytes;
}

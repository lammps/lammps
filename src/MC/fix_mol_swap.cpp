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

#include "fix_mol_swap.h"

#include "angle.h"
#include "atom.h"
#include "bond.h"
#include "comm.h"
#include "compute.h"
#include "dihedral.h"
#include "domain.h"
#include "error.h"
#include "fix.h"
#include "force.h"
#include "improper.h"
#include "kspace.h"
#include "modify.h"
#include "neighbor.h"
#include "pair.h"
#include "random_park.h"
#include "update.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixMolSwap::FixMolSwap(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg != 9) error->all(FLERR,"Illegal fix mol/swap command");

  vector_flag = 1;
  size_vector = 2;
  global_freq = 1;
  extvector = 0;
  restart_global = 1;
  time_depend = 1;

  // parse args

  nevery = utils::inumeric(FLERR,arg[3],false,lmp);
  ncycles = utils::inumeric(FLERR,arg[4],false,lmp);
  itype = utils::inumeric(FLERR,arg[5],false,lmp);
  jtype = utils::inumeric(FLERR,arg[6],false,lmp);
  seed = utils::inumeric(FLERR,arg[7],false,lmp);
  double temperature = utils::numeric(FLERR,arg[8],false,lmp);

  if (nevery <= 0) error->all(FLERR,"Illegal fix mol/swap command");
  if (ncycles < 0) error->all(FLERR,"Illegal fix mol/swap command");
  if (itype <= 0 || itype > atom->ntypes ||
      jtype <= 0 || jtype > atom->ntypes)
    error->all(FLERR,"Fix mol/swap atom types are invalid");
  if (seed <= 0) error->all(FLERR,"Illegal fix mol/swap command");
  if (temperature <= 0.0) error->all(FLERR,"Illegal fix mol/swap command");

  beta = 1.0/(force->boltz*temperature);

  // random number generator, same for all procs

  random = new RanPark(lmp,seed);

  // set up reneighboring

  force_reneighbor = 1;
  next_reneighbor = update->ntimestep + 1;

  // zero out counters

  nswap_attempts = 0.0;
  nswap_successes = 0.0;

  // set comm size needed by this Fix

  if (atom->q_flag) comm_forward = 2;
  else comm_forward = 1;
}

/* ---------------------------------------------------------------------- */

FixMolSwap::~FixMolSwap()
{
  delete random;
}

/* ---------------------------------------------------------------------- */

int FixMolSwap::setmask()
{
  int mask = 0;
  mask |= PRE_EXCHANGE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixMolSwap::init()
{
  // c_pe = compute used to calculate before/after potential energy

  char *id_pe = (char *) "thermo_pe";
  int ipe = modify->find_compute(id_pe);
  c_pe = modify->compute[ipe];

  // minmol = smallest molID with atoms of itype or jtype
  // maxmol = largest molID with atoms of itype or jtype
  // require: molID > 0, atoms in fix group

  int *mask = atom->mask;
  int *type = atom->type;
  tagint *molecule = atom->molecule;
  int nlocal = atom->nlocal;

  tagint minmol_me = MAXTAGINT;
  tagint maxmol_me = 0;

  for (int i = 0; i < nlocal; i++) {
    if (molecule[i] == 0) continue;
    if (!(mask[i] & groupbit)) continue;

    if (molecule[i] < minmol_me) {
      if (type[i] == itype || type[i] == jtype) minmol_me = molecule[i];
    }
    if (molecule[i] > maxmol_me) {
      if (type[i] == itype || type[i] == jtype) maxmol_me = molecule[i];
    }
  }

  MPI_Allreduce(&minmol_me,&minmol,1,MPI_LMP_TAGINT,MPI_MIN,world);
  MPI_Allreduce(&maxmol_me,&maxmol,1,MPI_LMP_TAGINT,MPI_MAX,world);

  // check if all cutoffs involving itype and jtype are the same
  // if not, reneighboring will be needed after swaps

  double **cutsq = force->pair->cutsq;

  unequal_cutoffs = false;
  for (int ktype = 1; ktype <= atom->ntypes; ktype++)
    if (cutsq[itype][ktype] != cutsq[jtype][ktype]) unequal_cutoffs = true;
}

/* ----------------------------------------------------------------------
   perform Ncycle Monte Carlo swaps
------------------------------------------------------------------------- */

void FixMolSwap::pre_exchange()
{
  // just return if should not be called on this timestep

  if (next_reneighbor != update->ntimestep) return;

  // insure current system is ready to compute energy

  if (domain->triclinic) domain->x2lamda(atom->nlocal);
  domain->pbc();
  comm->exchange();
  comm->borders();
  if (domain->triclinic) domain->lamda2x(atom->nlocal+atom->nghost);
  if (modify->n_pre_neighbor) modify->pre_neighbor();
  neighbor->build(1);

  // energy_stored = energy of current state
  // will be updated after accepted swaps

  energy_stored = energy_full();

  // attempt Ncycle molecule swaps

  int nsuccess = 0;
  for (int m = 0; m < ncycles; m++) nsuccess += attempt_swap();

  nswap_attempts += ncycles;
  nswap_successes += nsuccess;

  next_reneighbor = update->ntimestep + nevery;
}

/* ----------------------------------------------------------------------
   attempt a swap of atom types within a random molecule
   compare before/after energy and accept/reject the swap
------------------------------------------------------------------------- */

int FixMolSwap::attempt_swap()
{
  // pre-swap energy

  double energy_before = energy_stored;

  // pick a random molecule
  // swap all of its eligible itype & jtype atoms

  tagint molID = 
    minmol + static_cast<tagint> (random->uniform() * (maxmol-minmol+1));
  if (molID > maxmol) molID = maxmol;

  //printf("  Attempt %ld molID %d\n",update->ntimestep,molID);

  int *mask = atom->mask;
  int *type = atom->type;
  tagint *molecule = atom->molecule;
  int nlocal = atom->nlocal;

  // DEBUG
  int origitype = 0;
  int origjtype = 0;

  for (int i = 0; i < nlocal; i++) {
    if (molecule[i] != molID) continue;
    if (!(mask[i] & groupbit)) continue;
    if (type[i] == itype) {
      origitype++;
      type[i] = jtype;
    } else if (type[i] == jtype) {
      origjtype++;
      type[i] = itype;
    }
  }

  //printf("    ijtype %d %d orig %d %d count %d\n",
   //      itype,jtype,origitype,origjtype,origitype+origjtype);

  // if unequal_cutoffs, call comm->borders() and rebuild neighbor list
  // else communicate ghost atoms
  // call to comm->exchange() is a no-op but clears ghost atoms

  if (unequal_cutoffs) {
    if (domain->triclinic) domain->x2lamda(atom->nlocal);
    comm->exchange();
    comm->borders();
    if (domain->triclinic) domain->lamda2x(atom->nlocal+atom->nghost);
    if (modify->n_pre_neighbor) modify->pre_neighbor();
    neighbor->build(1);
  } else {
    comm->forward_comm_fix(this);
  }

  // post-swap energy

  double energy_after = energy_full();

  //printf("    eng before: %g, eng after %g\n",energy_before,energy_after);

  // swap accepted

  if (random->uniform() < exp(beta*(energy_before - energy_after))) {
    energy_stored = energy_after;
    //printf("    accept\n");
    return 1;

  // swap not accepted
  // restore the swapped itype & jtype atoms
  // do not need to re-call comm->borders() and rebuild neighbor list
  //   since will be done on next cycle or in Verlet when this fix is done

  } else {
    //printf("    reject\n");
    for (int i = 0; i < nlocal; i++) {
      if (molecule[i] != molID) continue;
      if (!(mask[i] & groupbit)) continue;
      if (type[i] == itype) type[i] = jtype;
      else if (type[i] == jtype) type[i] = itype;
    }
  }

  return 0;
}

/* ----------------------------------------------------------------------
   compute system potential energy before or after a swap
------------------------------------------------------------------------- */

double FixMolSwap::energy_full()
{
  int eflag = 1;
  int vflag = 0;

  if (modify->n_pre_force) modify->pre_force(vflag);

  if (force->pair) force->pair->compute(eflag,vflag);

  if (atom->molecular != Atom::ATOMIC) {
    if (force->bond) force->bond->compute(eflag,vflag);
    if (force->angle) force->angle->compute(eflag,vflag);
    if (force->dihedral) force->dihedral->compute(eflag,vflag);
    if (force->improper) force->improper->compute(eflag,vflag);
  }

  if (force->kspace) force->kspace->compute(eflag,vflag);

  if (modify->n_post_force) modify->post_force(vflag);

  update->eflag_global = update->ntimestep;
  double total_energy = c_pe->compute_scalar();

  return total_energy;
}

/* ---------------------------------------------------------------------- */

int FixMolSwap::pack_forward_comm(int n, int *list, double *buf, int /*pbc_flag*/, int * /*pbc*/)
{
  int i,j,m;

  int *type = atom->type;
  double *q = atom->q;

  m = 0;

  if (atom->q_flag) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = type[j];
      buf[m++] = q[j];
    }
  } else {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = type[j];
    }
  }

  return m;
}

/* ---------------------------------------------------------------------- */

void FixMolSwap::unpack_forward_comm(int n, int first, double *buf)
{
  int i,m,last;

  int *type = atom->type;
  double *q = atom->q;

  m = 0;
  last = first + n;

  if (atom->q_flag) {
    for (i = first; i < last; i++) {
      type[i] = static_cast<int> (buf[m++]);
      q[i] = buf[m++];
    }
  } else {
    for (i = first; i < last; i++)
      type[i] = static_cast<int> (buf[m++]);
  }
}

/* ----------------------------------------------------------------------
  return acceptance ratio
------------------------------------------------------------------------- */

double FixMolSwap::compute_vector(int n)
{
  if (n == 0) return nswap_attempts;
  if (n == 1) return nswap_successes;
  return 0.0;
}

/* ----------------------------------------------------------------------
   pack entire state of Fix into one write
------------------------------------------------------------------------- */

void FixMolSwap::write_restart(FILE *fp)
{
  int n = 0;
  double list[6];
  list[n++] = random->state();
  list[n++] = ubuf(next_reneighbor).d;
  list[n++] = nswap_attempts;
  list[n++] = nswap_successes;
  list[n++] = ubuf(update->ntimestep).d;

  if (comm->me == 0) {
    int size = n * sizeof(double);
    fwrite(&size,sizeof(int),1,fp);
    fwrite(list,sizeof(double),n,fp);
  }
}

/* ----------------------------------------------------------------------
   use state info from restart file to restart the Fix
------------------------------------------------------------------------- */

void FixMolSwap::restart(char *buf)
{
  int n = 0;
  double *list = (double *) buf;

  seed = static_cast<int> (list[n++]);
  random->reset(seed);

  next_reneighbor = (bigint) ubuf(list[n++]).i;

  nswap_attempts = static_cast<int>(list[n++]);
  nswap_successes = static_cast<int>(list[n++]);

  bigint ntimestep_restart = (bigint) ubuf(list[n++]).i;
  if (ntimestep_restart != update->ntimestep)
    error->all(FLERR,"Must not reset timestep when restarting fix mol/swap");
}

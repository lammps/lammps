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

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

#define BIG 1.0e20

/* ---------------------------------------------------------------------- */

FixMolSwap::FixMolSwap(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg), random(nullptr), c_pe(nullptr)
{
  if (narg < 9) error->all(FLERR,"Illegal fix mol/swap command");

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

  // optional args

  ke_flag = 1;

  int iarg = 9;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"ke") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix mol/swap command");
      ke_flag = utils::logical(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;
    } else error->all(FLERR,"Illegal fix mol/swap command");
  }

  // error check

  if (nevery <= 0) error->all(FLERR,"Illegal fix mol/swap command");
  if (ncycles < 0) error->all(FLERR,"Illegal fix mol/swap command");
  if (itype == jtype) error->all(FLERR,"Illegal fix mol/swap command");
  if (itype <= 0 || itype > atom->ntypes ||
      jtype <= 0 || jtype > atom->ntypes)
    error->all(FLERR,"Fix mol/swap atom types are invalid");
  if (seed <= 0) error->all(FLERR,"Illegal fix mol/swap command");
  if (temperature <= 0.0) error->all(FLERR,"Illegal fix mol/swap command");
  if (ke_flag && atom->rmass)
    error->all(FLERR,"Cannot conserve kinetic energy with fix mol/swap "
               "unless per-type masses");

  beta = 1.0/(force->boltz*temperature);

  // random number generator, same for all procs

  random = new RanPark(lmp,seed);

  // set up reneighboring

  force_reneighbor = 1;
  next_reneighbor = update->ntimestep + 1;

  // zero out counters

  nswap_attempt = 0.0;
  nswap_accept = 0.0;

  // qflag = 1 if charges are defined

  if (atom->q_flag) qflag = 1;
  else qflag = 0;

  // set comm size needed by this Fix

  if (qflag) comm_forward = 2;
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

  auto id_pe = (char *) "thermo_pe";
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

  // if ke_flag, check if itype/jtype masses are different
  // if yes, pre-calcuate velocity scale factors that keep KE constant
  // if no, turn ke_flag off

  if (ke_flag) {
    double *mass = atom->mass;
    if (mass[itype] != mass[jtype]) {
      i2j_vscale = sqrt(mass[itype]/mass[jtype]);
      j2i_vscale = sqrt(mass[jtype]/mass[itype]);
    } else ke_flag = 0;
  }

  // if qflag, check if all charges for itype are identical, ditto for jtype
  // if not, cannot swap charges, issue warning

  if (qflag) {
    double *q = atom->q;

    double iqone,jqone;
    iqone = jqone = -BIG;

    for (int i = 0; i < nlocal; i++) {
      if (molecule[i] == 0) continue;
      if (!(mask[i] & groupbit)) continue;
      if (type[i] == itype) iqone = q[i];
      if (type[i] == jtype) jqone = q[i];
    }

    MPI_Allreduce(&iqone,&iq,1,MPI_DOUBLE,MPI_MAX,world);
    MPI_Allreduce(&jqone,&jq,1,MPI_DOUBLE,MPI_MAX,world);

    int flag = 0;

    for (int i = 0; i < nlocal; i++) {
      if (molecule[i] == 0) continue;
      if (!(mask[i] & groupbit)) continue;
      if (type[i] == itype && q[i] != iq) flag = 1;
      if (type[i] == jtype && q[i] != jq) flag = 1;
    }

    int flagall;
    MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_MAX,world);
    if (flagall) qflag = 0;

    if (!qflag && comm->me == 0)
      error->warning(FLERR,"Cannot swap charges in fix mol/swap");
  }

  // check if all cutoffs involving itype and jtype are the same
  // if not, reneighboring will be needed after swaps

  double **cutsq = force->pair->cutsq;

  unequal_cutoffs = false;
  for (int ktype = 1; ktype <= atom->ntypes; ktype++)
    if (cutsq[itype][ktype] != cutsq[jtype][ktype]) unequal_cutoffs = true;
}

/* ----------------------------------------------------------------------
   perform ncycle Monte Carlo swaps
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

  int naccept = 0;
  for (int m = 0; m < ncycles; m++) naccept += attempt_swap();

  // udpate MC stats

  nswap_attempt += ncycles;
  nswap_accept += naccept;

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
  // swap properties of all its eligible itype & jtype atoms

  tagint molID =
    minmol + static_cast<tagint> (random->uniform() * (maxmol-minmol+1));
  if (molID > maxmol) molID = maxmol;

  int *mask = atom->mask;
  double **v = atom->v;
  double *q = atom->q;
  int *type = atom->type;
  tagint *molecule = atom->molecule;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (molecule[i] != molID) continue;
    if (!(mask[i] & groupbit)) continue;
    if (type[i] == itype) {
      type[i] = jtype;
      if (qflag) q[i] = jq;
      if (ke_flag) {
        v[i][0] *= i2j_vscale;
        v[i][1] *= i2j_vscale;
        v[i][2] *= i2j_vscale;
      }
    } else if (type[i] == jtype) {
      type[i] = itype;
      if (qflag) q[i] = iq;
      if (ke_flag) {
        v[i][0] *= j2i_vscale;
        v[i][1] *= j2i_vscale;
        v[i][2] *= j2i_vscale;
      }
    }
  }

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
    comm->forward_comm(this);
  }

  // post-swap energy

  double energy_after = energy_full();

  // swap accepted, return 1

  if (random->uniform() < exp(beta*(energy_before - energy_after))) {
    energy_stored = energy_after;
    return 1;
  }

  // swap not accepted, return 0
  // restore the swapped itype & jtype atoms
  // do not need to re-call comm->borders() and rebuild neighbor list
  //   since will be done on next cycle or in Verlet when this fix finishes

  for (int i = 0; i < nlocal; i++) {
    if (molecule[i] != molID) continue;
    if (!(mask[i] & groupbit)) continue;
    if (type[i] == itype) {
      type[i] = jtype;
      if (qflag) q[i] = jq;
      if (ke_flag) {
        v[i][0] *= i2j_vscale;
        v[i][1] *= i2j_vscale;
        v[i][2] *= i2j_vscale;
      }
    } else if (type[i] == jtype) {
      type[i] = itype;
      if (qflag) q[i] = iq;
      if (ke_flag) {
        v[i][0] *= j2i_vscale;
        v[i][1] *= j2i_vscale;
        v[i][2] *= j2i_vscale;
      }
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

  if (modify->n_post_force_any) modify->post_force(vflag);

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

  if (qflag) {
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

  if (qflag) {
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
  if (n == 0) return nswap_attempt;
  if (n == 1) return nswap_accept;
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
  list[n++] = nswap_attempt;
  list[n++] = nswap_accept;
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
  auto list = (double *) buf;

  seed = static_cast<int> (list[n++]);
  random->reset(seed);

  next_reneighbor = (bigint) ubuf(list[n++]).i;

  nswap_attempt = static_cast<int>(list[n++]);
  nswap_accept = static_cast<int>(list[n++]);

  bigint ntimestep_restart = (bigint) ubuf(list[n++]).i;
  if (ntimestep_restart != update->ntimestep)
    error->all(FLERR,"Must not reset timestep when restarting fix mol/swap");
}

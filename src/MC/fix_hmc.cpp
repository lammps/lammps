/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
   LAMMPS Development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors: Charlles R. A. Abreu (abreu@eq.ufrj.br)
                         Ana J. Silveira (asilveira@plapiqui.edu.ar)
                         Jack S. Draney (jacksdraney@gmail.com)
                         Louis E. S. Hoffenberg (lhoff@princeton.edu)
                         Baris E. Ugur (bu9134@princeton.edu)
------------------------------------------------------------------------- */

#include "fix_hmc.h"

#include "angle.h"
#include "atom.h"
#include "atom_vec.h"
#include "bond.h"
#include "comm.h"
#include "compute.h"
#include "dihedral.h"
#include "domain.h"
#include "error.h"
#include "fix.h"
#include "fix_rigid_small.h"
#include "force.h"
#include "group.h"
#include "improper.h"
#include "kspace.h"
#include "memory.h"
#include "modify.h"
#include "pair.h"
#include "random_park.h"
#include "update.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

static constexpr double BUFFACTOR = 1.5;
static constexpr auto SIX = sizeof(double) * 6;

/* ---------------------------------------------------------------------- */

FixHMC::FixHMC(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), buf_store(nullptr), id_rigid(nullptr), fix_rigid(nullptr), random(nullptr),
    random_equal(nullptr), eglobal(nullptr), eglobalptr(nullptr), vglobal(nullptr),
    vglobalptr(nullptr), pe(nullptr), ke(nullptr), peatom(nullptr), press(nullptr),
    pressatom(nullptr)
{
  // defaults

  if (narg < 5) utils::missing_cmd_args(FLERR, "fix hmc", error);

  // required arguments

  nevery = utils::inumeric(FLERR, arg[3], false, lmp);
  int seed = utils::inumeric(FLERR, arg[4], false, lmp);
  double temp = utils::numeric(FLERR, arg[5], false, lmp);

  if (seed <= 0) error->all(FLERR, 4, "Fix hmc seed must be > 0");
  if (temp <= 0) error->all(FLERR, 5, "Fix hmc temperature must be > 0.0");

  KT = force->boltz * temp / force->mvv2e;    // K*T in mvv units
  mbeta = -1.0 / (force->boltz * temp);       // -1/(K*T) in energy units

  // optional keywords

  mom_flag = 1;
  resample_on_accept_flag = 0;
  flag_rigid = 0;

  int iarg = 6;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "mom") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "hmc mom", error);
      mom_flag = utils::logical(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "resample") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "hmc resample", error);
      resample_on_accept_flag = utils::logical(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "rigid") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "hmc rigid", error);
      delete[] id_rigid;
      id_rigid = utils::strdup(arg[iarg + 1]);
      auto *ifix = modify->get_fix_by_id(id_rigid);
      if (!ifix) error->all(FLERR, iarg + 1, "Unknown rigid fix id {} for fix hmc", id_rigid);
      fix_rigid = dynamic_cast<FixRigidSmall *>(ifix);
      if (!fix_rigid ||
          (!utils::strmatch(ifix->style, "^rigid/small") &&
           !utils::strmatch(ifix->style, "^rigid/nve/small")))
        error->all(FLERR, Error::NOLASTLINE,
                   "Fix ID {} for fix hmc does not point to fix rigid/small or rigid/nve/small",
                   id_rigid);
      flag_rigid = 1;
      iarg += 2;
    } else {
      error->all(FLERR, iarg, "Unknown fix hmc keyword {}", arg[iarg]);
    }
  }

  // random = different RNG on each processor
  // random_equal = same RNG on each processor

  random = new RanPark(lmp, seed + comm->me);
  for (int i = 0; i < 100; i++) random->uniform();
  random_equal = new RanPark(lmp, seed);

  // register callback

  atom->add_callback(0);

  // add new computes for global and per-atom properties

  ke = modify->add_compute(fmt::format("hmc_ke_{} all ke", id));
  pe = modify->add_compute(fmt::format("hmc_pe_{} all pe", id));
  peatom = modify->add_compute(fmt::format("hmc_peatom_{} all pe/atom", id));
  press = modify->add_compute(fmt::format("hmc_press_{} all pressure NULL virial", id));
  pressatom = modify->add_compute(fmt::format("hmc_pressatom_{} all stress/atom NULL virial", id));

  // fix attributes

  global_freq = 1;
  scalar_flag = 1;
  extscalar = 0;
  vector_flag = 1;
  extvector = 0;
  size_vector = 5;
  force_reneighbor = 1;
  next_reneighbor = -1;

  first_init_complete = false;
  first_setup_complete = false;

  // initializations

  nattempts = 0;
  naccepts = 0;
  DeltaPE = 0.0;
  DeltaKE = 0.0;
}

/* ---------------------------------------------------------------------- */

FixHMC::~FixHMC()
{
  atom->delete_callback(id, 0);

  memory->destroy(eglobal);
  memory->destroy(vglobal);
  memory->destroy(buf_store);

  delete[] id_rigid;

  delete[] eglobalptr;
  if (vglobalptr)
    for (int m = 0; m < nv; m++) delete[] vglobalptr[m];
  delete[] vglobalptr;

  delete random;
  delete random_equal;

  modify->delete_compute(std::string("hmc_ke_") + id);
  modify->delete_compute(std::string("hmc_pe_") + id);
  modify->delete_compute(std::string("hmc_peatom_") + id);
  modify->delete_compute(std::string("hmc_press_") + id);
  modify->delete_compute(std::string("hmc_pressatom_") + id);
}

/* ---------------------------------------------------------------------- */

void FixHMC::setup_arrays_and_pointers()
{
  int i, m;
  int pair_flag;
  int bond_flag;
  int angle_flag;
  int dihedral_flag;
  int improper_flag;
  int kspace_flag;

  // determine which energy contributions must be computed

  ne = 0;
  if (force->pair) {
    pair_flag = 1;
    ne++;
  } else
    pair_flag = 0;
  if (force->bond) {
    bond_flag = 1;
    ne++;
  } else
    bond_flag = 0;
  if (force->angle) {
    angle_flag = 1;
    ne++;
  } else
    angle_flag = 0;
  if (force->dihedral) {
    dihedral_flag = 1;
    ne++;
  } else
    dihedral_flag = 0;
  if (force->improper) {
    improper_flag = 1;
    ne++;
  } else
    improper_flag = 0;
  if (force->kspace) {
    kspace_flag = 1;
    ne++;
  } else
    kspace_flag = 0;

  // initialize arrays for managing global energy terms

  neg = pair_flag ? ne + 1 : ne;
  memory->create(eglobal, neg, "fix_hmc:eglobal");
  delete[] eglobalptr;
  eglobalptr = new double *[neg];

  m = 0;
  if (pair_flag) {
    eglobalptr[m++] = &force->pair->eng_vdwl;
    eglobalptr[m++] = &force->pair->eng_coul;
  }

  if (bond_flag) eglobalptr[m++] = &force->bond->energy;
  if (angle_flag) eglobalptr[m++] = &force->angle->energy;
  if (dihedral_flag) eglobalptr[m++] = &force->dihedral->energy;
  if (improper_flag) eglobalptr[m++] = &force->improper->energy;
  if (kspace_flag) eglobalptr[m++] = &force->kspace->energy;

  // initialize arrays for managing global virial terms

  nv = ne;
  for (const auto &ifix : modify->get_fix_list())
    if (ifix->virial_global_flag) nv++;
  memory->create(vglobal, nv, 6, "fix_hmc:vglobal");
  if (vglobalptr)
    for (m = 0; m < nv; m++) delete[] vglobalptr[m];
  delete[] vglobalptr;
  vglobalptr = new double **[nv];

  for (m = 0; m < nv; m++) vglobalptr[m] = new double *[6];

  for (i = 0; i < 6; i++) {
    m = 0;
    if (pair_flag) vglobalptr[m++][i] = &force->pair->virial[i];
    if (bond_flag) vglobalptr[m++][i] = &force->bond->virial[i];
    if (angle_flag) vglobalptr[m++][i] = &force->angle->virial[i];
    if (dihedral_flag) vglobalptr[m++][i] = &force->dihedral->virial[i];
    if (improper_flag) vglobalptr[m++][i] = &force->improper->virial[i];
    if (kspace_flag) vglobalptr[m++][i] = &force->kspace->virial[i];
    for (const auto &ifix : modify->get_fix_list())
      if (ifix->virial_global_flag) vglobalptr[m++][i] = &ifix->virial[i];
  }
}

/* ---------------------------------------------------------------------- */

int FixHMC::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixHMC::init()
{
  int ntimestep = update->ntimestep;

  // check whether there is any fixes that change box size and/or shape

  for (const auto &ifix : modify->get_fix_list())
    if (ifix->box_change)
      error->all(FLERR, Error::NOLASTLINE,
                 "Fix hmc is incompatible with fixes that change box size or shape");

  // check whether fix rigid/small still exists

  if (flag_rigid) {
    auto *ifix = modify->get_fix_by_id(id_rigid);
    if (!ifix)
      error->all(FLERR, Error::NOLASTLINE, "Unknown rigid fix id {} for fix hmc", id_rigid);
    fix_rigid = dynamic_cast<FixRigidSmall *>(ifix);
    if (!fix_rigid ||
        (!utils::strmatch(ifix->style, "^rigid/small") &&
         !utils::strmatch(ifix->style, "^rigid/nve/small")))
      error->all(FLERR, Error::NOLASTLINE,
                 "Fix ID {} for fix hmc does not point to fix rigid/small or rigid/nve/small",
                 id_rigid);
  }

  // check if nve time integration fix exists

  for (const auto &ifix : modify->get_fix_list()) {
    if (ifix->time_integrate && !ifix->rigid_flag) {
      if (!utils::strmatch(ifix->style, "^nve"))
        if (comm->me == 0)
          error->warning(FLERR, "Non NVE time integration fix {} {} found", ifix->id, ifix->style);
    }
  }

  // check whether there are subsequent fixes with active virial_flag

  int past_this_fix = false;
  int past_rigid = !flag_rigid;
  for (const auto &ifix : modify->get_fix_list()) {
    if (past_this_fix && past_rigid && ifix->virial_global_flag) {
      if (comm->me == 0)
        utils::logmesg(lmp, "Fix {} with id {} defined after fix hmc.\n", ifix->style, ifix->id);
      error->all(FLERR, Error::NOLASTLINE,
                 "Fix hmc cannot precede fixes that contribute to the system pressure");
    }
    if (!strcmp(ifix->id, id)) past_this_fix = true;
    if (flag_rigid && !strcmp(ifix->id, fix_rigid->id)) past_rigid = true;
  }

  if (!first_init_complete) {

    // look for computes with active press_flag

    press_flag = 0;
    pressatom_flag = 0;
    peatom_flag = 0;
    for (const auto &icompute : modify->get_compute_list()) {
      if (utils::strmatch(icompute->id, "^hmc_")) {
        press_flag = press_flag | icompute->pressflag;
        pressatom_flag = pressatom_flag | icompute->pressatomflag;
        peatom_flag = peatom_flag | icompute->peatomflag;
      }
    }
    // initialize arrays and pointers for saving/restoring state

    setup_arrays_and_pointers();

    // count per-atom properties to be exchanged

    nvalues = 0;
    if (peatom_flag) nvalues += ne;
    if (pressatom_flag) nvalues += 6 * nv;

    // activate potential energy and other necessary calculations at setup

    pe->addstep(ntimestep);
    if (peatom_flag) peatom->addstep(ntimestep);
    if (press_flag) press->addstep(ntimestep);
    if (pressatom_flag) pressatom->addstep(ntimestep);

    first_init_complete = true;
  }
}

/* ----------------------------------------------------------------------
   initialize MC move, save current state, and activate computes
------------------------------------------------------------------------- */

void FixHMC::setup(int vflag)
{
  if (!first_setup_complete) {

    // initialize fix rigid first to avoid saving uninitialized state

    if (flag_rigid) fix_rigid->setup(vflag);

    // compute properties of initial state

    if (flag_rigid) {
      random_velocities();
      rigid_body_random_velocities();
    } else {
      random_velocities();
    }

    update->eflag_global = update->ntimestep;
    PE = pe->compute_scalar();
    KE = ke->compute_scalar();

    // activate potential energy and other necessary calculations

    int nextstep = update->ntimestep + nevery;
    pe->addstep(nextstep);
    if (peatom_flag) peatom->addstep(nextstep);
    if (press_flag) press->addstep(nextstep);
    if (pressatom_flag) pressatom->addstep(nextstep);

    // create buffer, store fixes for warnings

    int maxexchange_fix = 0;
    int maxexchange_atom = atom->avec->maxexchange;

    for (const auto &fix : modify->get_fix_list()) maxexchange_fix += fix->maxexchange;
    maxexchange = maxexchange_atom + maxexchange_fix;
    bufextra = maxexchange + Comm::BUFEXTRA;

    maxstore = Comm::BUFEXTRA;
    grow_store(maxstore, 2);
    save_current_state();
  }
}

/* ----------------------------------------------------------------------
   apply the Metropolis acceptance criterion
   restore saved system state if move is rejected
   activate computes for the next MC step
------------------------------------------------------------------------- */

void FixHMC::end_of_step()
{
  nattempts++;

  // compute potential and kinetic energy variations

  update->eflag_global = update->ntimestep;
  double newPE = pe->compute_scalar();
  double newKE = ke->compute_scalar();
  DeltaPE = newPE - PE;
  DeltaKE = newKE - KE;

  // apply the Metropolis criterion

  double DeltaE = DeltaPE + DeltaKE;
  int accept = DeltaE < 0.0;
  if (!accept) {
    accept = random_equal->uniform() <= exp(mbeta * DeltaE);
    MPI_Bcast(&accept, 1, MPI_INT, 0, world);
  }

  // for accept: update potential energy and save current state
  // for reject: restore saved state and trigger reneighboring on next step
  //   this will perform an exchange() and borders() to re-acquire ghost atoms
  // NOTE: assumption here for rejection
  //   after N steps an exchange() operation on next timestep will be able to
  //     migrate atoms with old coords back to their original owning procs
  //   if this is not the case, atoms will potentially be lost
  if (accept) {
    naccepts++;
    PE = newPE;
    KE = newKE;
    save_current_state();
  } else {
    restore_saved_state();
    next_reneighbor = update->ntimestep + 1;
  }

  // choose new velocities and compute kinetic energy

  if (!accept || resample_on_accept_flag) {
    if (flag_rigid)
      rigid_body_random_velocities();
    else
      random_velocities();
    KE = ke->compute_scalar();
  }

  // activate potential energy and other necessary calculations

  int nextstep = update->ntimestep + nevery;
  if (nextstep <= update->laststep) {
    pe->addstep(nextstep);
    if (peatom_flag) peatom->addstep(nextstep);
    if (press_flag) press->addstep(nextstep);
    if (pressatom_flag) pressatom->addstep(nextstep);
  }
}

/* ----------------------------------------------------------------------
   return acceptance fraction of proposed MC moves
------------------------------------------------------------------------- */

double FixHMC::compute_scalar()
{
  double acc_frac = naccepts;
  acc_frac /= MAX(1, nattempts);
  return acc_frac;
}

/* ----------------------------------------------------------------------
   return the acceptance fraction of proposed MC moves, or
   the total energy variation of the last proposed MC move, or
   the mean-square atom displacement in the last proposed MC move
------------------------------------------------------------------------- */

double FixHMC::compute_vector(int item)
{
  int n = item + 1;
  if (n == 1)
    return naccepts;
  else if (n == 2)
    return nattempts;
  else if (n == 3)
    return DeltaPE;
  else if (n == 4)
    return DeltaKE;
  else if (n == 5)
    return DeltaPE + DeltaKE;
  else
    return 0.0;
}

/* ----------------------------------------------------------------------
   realloc the size of the send buffer as needed with BUFFACTOR and bufextra
   flag = 0, don't need to realloc with copy, just free/malloc w/ BUFFACTOR
   flag = 1, realloc with BUFFACTOR
   flag = 2, free/malloc w/out BUFFACTOR
------------------------------------------------------------------------- */

void FixHMC::grow_store(int n, int flag)
{
  if (flag == 0) {
    maxstore = static_cast<int>(BUFFACTOR * n);
    memory->destroy(buf_store);
    memory->create(buf_store, maxstore + bufextra, "fix_hmc:buf_store");
  } else if (flag == 1) {
    maxstore = static_cast<int>(BUFFACTOR * n);
    memory->grow(buf_store, maxstore + bufextra, "fix_hmc:buf_store");
  } else {
    memory->grow(buf_store, maxstore + bufextra, "fix_hmc:buf_store");
  }
}

/* ----------------------------------------------------------------------
   save current system state for eventual restoration if move is rejected
------------------------------------------------------------------------- */

void FixHMC::save_current_state()
{
  int m;

  int nlocal = atom->nlocal;
  AtomVec *avec = atom->avec;
  nstore = 0;

  // store all needed info about owned atoms via pack_exchange()
  for (int i = 0; i < nlocal; i++) {
    if (nstore > maxstore) grow_store(nstore, 1);
    nstore += avec->pack_exchange(i, &buf_store[nstore]);
  }
  // save global energy terms

  for (m = 0; m < neg; m++) eglobal[m] = *eglobalptr[m];

  // save global virial terms

  if (press_flag)
    for (m = 0; m < nv; m++) memcpy(vglobal[m], *vglobalptr[m], SIX);
}

/* ----------------------------------------------------------------------
   restore system state saved at the beginning of the MC step
------------------------------------------------------------------------- */

void FixHMC::restore_saved_state()
{
  int i;
  AtomVec *avec = atom->avec;

  // clear atom map

  if (atom->map_style != Atom::MAP_NONE) atom->map_clear();

  // delete all owned + ghost atoms

  atom->nghost = 0;
  atom->avec->clear_bonus();
  atom->nlocal = 0;

  // unpack exchange buffer
  // this will restore atoms with coords and other properties from N steps ago
  // NOTE: atoms will not be re-assigned to correct procs until reneighboring on next step
  //       likewise ghost atoms will not be re-created until reneighboring on next step

  int m = 0;
  while (m < nstore) m += avec->unpack_exchange(&buf_store[m]);

  for (const auto &ifix : modify->get_fix_list()) ifix->pre_exchange();
  domain->pbc();
  domain->reset_box();
  comm->setup();
  comm->exchange();
  comm->borders();

  // reinit atom_map

  if (atom->map_style != Atom::MAP_NONE) {
    atom->map_clear();
    atom->map_init();
    atom->map_set();
  }

  // ensure fix_rigid images are OK

  if (flag_rigid) fix_rigid->pre_neighbor();

  // restore global energy terms

  for (i = 0; i < neg; i++) *eglobalptr[i] = eglobal[i];

  // restore global virial terms

  if (press_flag)
    for (i = 0; i < nv; i++) memcpy(*vglobalptr[i], vglobal[i], SIX);
}

/* ----------------------------------------------------------------------
   randomly choose velocities from a Maxwell-Boltzmann distribution
------------------------------------------------------------------------- */

void FixHMC::random_velocities()
{
  double **v = atom->v;
  int *type = atom->type;
  int *mask = atom->mask;

  double stdev;
  int nlocal, dimension;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  nlocal = atom->nlocal;
  dimension = domain->dimension;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      if (rmass)
        stdev = sqrt(KT / rmass[i]);
      else
        stdev = sqrt(KT / mass[type[i]]);
      for (int j = 0; j < dimension; j++) v[i][j] = stdev * random->gaussian();
    }

  if (mom_flag) {
    double vcm[3];
    group->vcm(igroup, group->mass(igroup), vcm);
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        for (int j = 0; j < dimension; j++) v[i][j] -= vcm[j];
      }
  }
}

/* ----------------------------------------------------------------------
   randomize VCMs of rigid bodies and their associated atoms
------------------------------------------------------------------------- */

void FixHMC::rigid_body_random_velocities()
{
  fix_rigid->resample_momenta(groupbit, mom_flag, random, KT);
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double FixHMC::memory_usage()
{
  double bytes = 0;
  bytes += (maxstore + bufextra) * sizeof(double);    // exchange
  return bytes;
}

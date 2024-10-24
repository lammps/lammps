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
#include "math_extra.h"
#include "memory.h"
#include "modify.h"
#include "neighbor.h"
#include "output.h"
#include "pair.h"
#include "random_park.h"
#include "update.h"

#include <cstdlib>
#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

static constexpr double EPSILON = 1.0e-7;

enum { ATOMS, VCM_OMEGA, XCM, ITENSOR, ROTATION, FORCE_TORQUE };

/* ---------------------------------------------------------------------- */

FixHMC::FixHMC(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), stored_tag(nullptr), stored_body(nullptr), stored_bodyown(nullptr),
    stored_bodytag(nullptr), stored_atom2body(nullptr), stored_xcmimage(nullptr),
    stored_displace(nullptr), stored_eflags(nullptr), stored_orient(nullptr),
    stored_dorient(nullptr), fix_rigid(nullptr), random(nullptr), random_equal(nullptr),
    rev_comm(nullptr), eatom(nullptr), eatomptr(nullptr), eglobal(nullptr), eglobalptr(nullptr),
    vglobal(nullptr), vglobalptr(nullptr), vatom(nullptr), vatomptr(nullptr), pe(nullptr),
    ke(nullptr), peatom(nullptr), press(nullptr), pressatom(nullptr)
{
  // set some defaults

  mom_flag = 1;
  resample_on_accept_flag = 0;
  rigid_flag = false;
  if (narg < 5) utils::missing_cmd_args(FLERR, "fix hmc", error);

  // Retrieve user-defined options:

  nevery = utils::numeric(FLERR, arg[3], false, lmp);         // Number of MD steps per MC step
  int seed = utils::numeric(FLERR, arg[4], false, lmp);       // Seed for random number generation
  double temp = utils::numeric(FLERR, arg[5], false, lmp);    // System temperature
  if (seed <= 0)
    error->all(FLERR, "Illegal fix hmc seed argument: {}. Seed must be greater than 0.0", seed);
  if (temp <= 0)
    error->all(FLERR, "Illegal fix hmc temp argument: {}. Temp must be greater than 0.0", temp);

  KT = force->boltz * temp / force->mvv2e;    // K*T in mvv units
  mbeta = -1.0 / (force->boltz * temp);       // -1/(K*T) in energy units

  // Check optional keywords:

  int iarg = 6;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "mom") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "hmc mom", error);
      mom_flag = utils::logical(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "ra") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "hmc ra", error);
      resample_on_accept_flag = utils::logical(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "rigid") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "hmc rigid", error);
      fix_rigid = (FixRigidSmall*) modify->get_fix_by_id(arg[iarg + 1]);
      if (fix_rigid == nullptr) error->all(FLERR, "Unknown rigid fix id {}", arg[iarg + 1]);
      rigid_flag = true;
      iarg += 2;
    } else {
      error->all(FLERR, "Unknown fix hmc keyword {}", arg[iarg]);
    }
  }

  // Initialize RNG with a different seed for each process:
  random = new RanPark(lmp, seed + comm->me);
  for (int i = 0; i < 100; i++) random->uniform();
  random_equal = new RanPark(lmp, seed);

  // Register callback:
  atom->add_callback(0);

  // Add new computes for global and per-atom properties:

  ke = modify->add_compute(fmt::format("hmc_ke_{} all ke", id));
  pe = modify->add_compute(fmt::format("hmc_pe_{} all pe", id));
  peatom = modify->add_compute(fmt::format("hmc_peatom_{} all pe/atom", id));
  press = modify->add_compute(fmt::format("hmc_press_{} all pressure NULL virial", id));
  pressatom = modify->add_compute(fmt::format("hmc_pressatom_{} all stress/atom NULL virial", id));

  // Define non-default fix attributes:

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

  // initialize quantities
  nattempts = 0;
  naccepts = 0;
  DeltaPE = 0;
  DeltaKE = 0;
}

/* ---------------------------------------------------------------------- */

FixHMC::~FixHMC()
{
  atom->delete_callback(id, 0);

  memory->destroy(eglobal);
  memory->destroy(vglobal);
  memory->destroy(eatom);
  memory->destroy(vatom);
  delete[] eglobalptr;
  if (vglobalptr)
    for (int m = 0; m < nv; m++) delete[] vglobalptr[m];
  delete[] vglobalptr;
  delete[] eatomptr;
  delete[] vatomptr;
  delete[] rev_comm;
  delete random;
  delete random_equal;
  modify->delete_compute(std::string("hmc_ke_") + id);
  modify->delete_compute(std::string("hmc_pe_") + id);
  modify->delete_compute(std::string("hmc_peatom_") + id);
  modify->delete_compute(std::string("hmc_press_") + id);
  modify->delete_compute(std::string("hmc_pressatom_") + id);

  memory->destroy(stored_tag);
  memory->destroy(stored_bodyown);
  memory->destroy(stored_bodytag);
  memory->destroy(stored_atom2body);
  memory->destroy(stored_xcmimage);
  memory->destroy(stored_displace);
  memory->destroy(stored_eflags);
  memory->destroy(stored_orient);
  memory->destroy(stored_dorient);

  for (Atom::PerAtom &stored_peratom_member : stored_peratom) {
    free(stored_peratom_member.address);
    free(stored_peratom_member.address_maxcols);
    free(stored_peratom_member.address_length);
  }
}

/* ---------------------------------------------------------------------- */

template <typename T>
void FixHMC::store_peratom_member(Atom::PerAtom &stored_peratom_member,
                                  Atom::PerAtom current_peratom_member, int ntotal, int nmax,
                                  int realloc)
{
  if (stored_peratom_member.name.compare(current_peratom_member.name)) {
    error->all(FLERR, "fix hmc tried to store incorrect peratom data");
  }

  // free old memory if reallocating and stored_peratom_member isn't a copy of current_peratom_member

  if (realloc && stored_peratom_member.address != current_peratom_member.address) {
    free(stored_peratom_member.address);
    stored_peratom_member.address = nullptr;
  }
  if (realloc && stored_peratom_member.address != current_peratom_member.address) {
    free(stored_peratom_member.address_maxcols);
    stored_peratom_member.address_maxcols = nullptr;
  }
  // peratom scalers
  if (current_peratom_member.cols == 0) {
    if (*(T **) current_peratom_member.address != nullptr) {
      if (realloc) stored_peratom_member.address = malloc(sizeof(T) * nmax);
      memcpy(stored_peratom_member.address, *(T **) current_peratom_member.address,
             ntotal * sizeof(T));
    } else {
      stored_peratom_member.address = nullptr;
    }
  } else {
    // peratom vectors
    int cols;
    if (current_peratom_member.cols < 0) {
      // variable column case
      cols = *(current_peratom_member.address_maxcols);
      if (realloc) stored_peratom_member.address_maxcols = (int *) malloc(sizeof(int));
      *(stored_peratom_member.address_maxcols) = *(current_peratom_member.address_maxcols);
    } else {
      // non-variable column case
      cols = current_peratom_member.cols;
    }
    if (*(T ***) current_peratom_member.address != nullptr) {
      if (realloc) stored_peratom_member.address = malloc(sizeof(T) * nmax * cols);
      for (int i = 0; i < ntotal; i++) {
        memcpy((T *) stored_peratom_member.address + i * cols,
               (**(T ***) current_peratom_member.address) + i * cols, sizeof(T) * cols);
      }
    } else {
      stored_peratom_member.address = nullptr;
    }
  }
  stored_peratom_member.cols = current_peratom_member.cols;
  stored_peratom_member.collength = current_peratom_member.collength;
  stored_peratom_member.address_length = nullptr;
}

/* ---------------------------------------------------------------------- */

template <typename T>
void FixHMC::restore_peratom_member(Atom::PerAtom stored_peratom_member,
                                    Atom::PerAtom &current_peratom_member, int nlocal)
{
  if (stored_peratom_member.name.compare(current_peratom_member.name)) {
    error->all(FLERR, "fix hmc tried to store incorrect peratom data");
  }
  if (stored_peratom_member.address == nullptr) return;
  int cols;
  // peratom scalers
  if (stored_peratom_member.cols == 0) {
    if (*(T **) current_peratom_member.address != nullptr) {
      memcpy(*(T **) current_peratom_member.address, stored_peratom_member.address,
             nlocal * sizeof(T));
    }
  } else {
    // peratom vectors
    if (stored_peratom_member.cols < 0) {
      // variable column case
      cols = *(stored_peratom_member.address_maxcols);
      *(current_peratom_member.address_maxcols) = *(stored_peratom_member.address_maxcols);
    } else {
      // non-variable column case
      cols = stored_peratom_member.cols;
    }
    if (*(T ***) current_peratom_member.address != nullptr) {
      for (int i = 0; i < nlocal; i++) {
        memcpy((**(T ***) current_peratom_member.address) + i * cols,
               (T *) stored_peratom_member.address + i * cols, sizeof(T) * cols);
      }
    }
  }
  current_peratom_member.cols = stored_peratom_member.cols;
  current_peratom_member.collength = stored_peratom_member.collength;
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

  // Per-atom scalar properties to be saved and restored:

  stored_nmax = 0;    // initialize so the memory gets allocated on first save
  stored_ntotal_body = 0;

  // Determine which energy contributions must be computed:
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

  // Initialize arrays for managing global energy terms:
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

  // Initialize arrays for managing global virial terms:
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

  // Determine which per-atom energy terms require reverse communication:
  rev_comm = new int[nv];
  m = 0;
  if (pair_flag) rev_comm[m++] = force->newton;
  if (bond_flag) rev_comm[m++] = force->newton_bond;
  if (angle_flag) rev_comm[m++] = force->newton_bond;
  if (dihedral_flag) rev_comm[m++] = force->newton_bond;
  if (improper_flag) rev_comm[m++] = force->newton_bond;
  if (kspace_flag) rev_comm[m++] = force->kspace->tip4pflag;
  for (i = ne; i < nv; i++) rev_comm[m++] = 0;

  // Initialize array of pointers to manage per-atom energies:
  delete[] eatomptr;
  eatomptr = new double **[ne];
  m = 0;
  if (pair_flag) eatomptr[m++] = &force->pair->eatom;
  if (bond_flag) eatomptr[m++] = &force->bond->eatom;
  if (angle_flag) eatomptr[m++] = &force->angle->eatom;
  if (dihedral_flag) eatomptr[m++] = &force->dihedral->eatom;
  if (improper_flag) eatomptr[m++] = &force->improper->eatom;
  if (kspace_flag) eatomptr[m++] = &force->kspace->eatom;

  // Initialize array of pointers to manage per-atom virials:
  delete[] vatomptr;
  vatomptr = new double ***[nv];
  m = 0;
  if (pair_flag) vatomptr[m++] = &force->pair->vatom;
  if (bond_flag) vatomptr[m++] = &force->bond->vatom;
  if (angle_flag) vatomptr[m++] = &force->angle->vatom;
  if (dihedral_flag) vatomptr[m++] = &force->dihedral->vatom;
  if (improper_flag) vatomptr[m++] = &force->improper->vatom;
  if (kspace_flag) vatomptr[m++] = &force->kspace->vatom;
  for (const auto &ifix : modify->get_fix_list())
    if (ifix->virial_peratom_flag) vatomptr[m++] = &ifix->vatom;

  // Determine the maximum and the actual numbers of per-atom variables in reverse
  // communications:
  // Note: fix-related virials do not communicate (thus 'ne' used instead of 'nv')
  comm_reverse = 0;
  ncommrev = 0;
  for (m = 0; m < ne; m++)
    if (rev_comm[m]) {
      comm_reverse += 7;    // 1 energy + 6 virials
      if (peatom_flag) ncommrev++;
      if (pressatom_flag) ncommrev += 6;
    }

  // Determine maximum number of per-atom variables in forward and reverse
  // communications when dealing with rigid bodies:
  if (rigid_flag) {
    comm_reverse = MAX(comm_reverse, 6);
    comm_forward = 12;
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

  if (atom->tag_enable == 0) error->all(FLERR, "Cannot use fix hmc unless atoms have IDs");

  // Check whether there is any fixes that change box size and/or shape:
  for (const auto &ifix : modify->get_fix_list())
    if (ifix->box_change)
      error->all(FLERR, "fix hmc is incompatible with fixes that change box size or shape");

  // Check whether there are subsequent fixes with active virial_flag:
  int past_this_fix = false;
  int past_rigid = !rigid_flag;
  for (const auto &ifix : modify->get_fix_list()) {
    if (past_this_fix && past_rigid && ifix->virial_global_flag) {
      if (comm->me == 0) utils::logmesg(lmp, "Fix {} defined after fix hmc.\n", ifix->style);
      error->all(FLERR, "fix hmc cannot precede fixes that modify the system pressure");
    }
    if (!strcmp(ifix->id, id)) past_this_fix = true;
     if (rigid_flag && !strcmp(ifix->id, fix_rigid->id)) past_rigid = true;
  }

  if (!first_init_complete)
  {
    // Look for computes with active peatomflag, press_flag, or pressatomflag:
    peatom_flag = 0;
    press_flag = 0;
    pressatom_flag = 0;
    for (const auto &icompute : modify->get_compute_list())
      if (strncmp(icompute->id, "hmc_", 4)) {
        peatom_flag = peatom_flag | icompute->peatomflag;
        press_flag = press_flag | icompute->pressflag;
        pressatom_flag = pressatom_flag | icompute->pressatomflag;
      }

    // Initialize arrays and pointers for saving/restoration of states:
    setup_arrays_and_pointers();

    // Count per-atom properties to be exchanged:
    nvalues = 0;
    if (peatom_flag) nvalues += ne;
    if (pressatom_flag) nvalues += 6 * nv;

    // Activate potential energy and other necessary calculations at setup:
    pe->addstep(ntimestep);
    if (peatom_flag) peatom->addstep(ntimestep);
    if (press_flag) press->addstep(ntimestep);
    if (pressatom_flag) pressatom->addstep(ntimestep);

    first_init_complete = true;
  }


  // (Re)allocate array of per-atom properties:
  grow_arrays(atom->nmax);
}

/* ----------------------------------------------------------------------
   Initialize MC move, save current state, and activate computes
------------------------------------------------------------------------- */

void FixHMC::setup(int vflag)
{
  if (!first_setup_complete) {
    // initialize rigid first to avoid saving uninitialized state
    if (rigid_flag) fix_rigid->setup(vflag);

    // Compute properties of the initial state:

    if (rigid_flag) {
      rigid_body_random_velocities();
    } else {
      random_velocities();
    }

    update->eflag_global = update->ntimestep;
    PE = pe->compute_scalar();
    KE = ke->compute_scalar();
    save_current_state();

    // Activate potential energy and other necessary calculations:
    int nextstep = update->ntimestep + nevery;
    pe->addstep(nextstep);
    if (peatom_flag) peatom->addstep(nextstep);
    if (press_flag) press->addstep(nextstep);
    if (pressatom_flag) pressatom->addstep(nextstep);
  }
}

/* ----------------------------------------------------------------------
   Apply the Metropolis acceptance criterion
   Restore saved system state if move is rejected
   Activate computes for the next MC step
------------------------------------------------------------------------- */

void FixHMC::end_of_step()
{
  nattempts++;

  // Compute potential and kinetic energy variations:
  update->eflag_global = update->ntimestep;
  double newPE = pe->compute_scalar();
  double newKE = ke->compute_scalar();
  DeltaPE = newPE - PE;
  DeltaKE = newKE - KE;

  // Apply the Metropolis criterion:
  double DeltaE = DeltaPE + DeltaKE;
  int accept = DeltaE < 0.0;
  if (!accept) {
    accept = random_equal->uniform() <= exp(mbeta * DeltaE);
    MPI_Bcast(&accept, 1, MPI_INT, 0, world);
  }
  if (accept) {
    // Update potential energy and save the current state:
    naccepts++;
    PE = newPE;
    save_current_state();
  } else {
    // Restore saved state and enforce check_distance/reneighboring in the next step:
    restore_saved_state();
    next_reneighbor = update->ntimestep + 1;
  }

  // Choose new velocities and compute kinetic energy:
  if (!accept || resample_on_accept_flag) {
    if (rigid_flag)
      rigid_body_random_velocities();
    else
      random_velocities();
    KE = ke->compute_scalar();
  }

  // Activate potential energy and other necessary calculations:
  int nextstep = update->ntimestep + nevery;
  if (nextstep <= update->laststep) {
    pe->addstep(nextstep);
    if (peatom_flag) peatom->addstep(nextstep);
    if (press_flag) press->addstep(nextstep);
    if (pressatom_flag) pressatom->addstep(nextstep);
  }
}

/* ----------------------------------------------------------------------
   Return the acceptance fraction of proposed MC moves
------------------------------------------------------------------------- */

double FixHMC::compute_scalar()
{
  double acc_frac = naccepts;
  acc_frac /= MAX(1, nattempts);
  return acc_frac;
}

/* ----------------------------------------------------------------------
   Return the acceptance fraction of proposed MC moves, or
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
   Save the system state for eventual restoration if move is rejected
------------------------------------------------------------------------- */
void FixHMC::save_current_state()
{
  int m;

  int nlocal = atom->nlocal;
  int ntotal = nlocal + atom->nghost;
  int nmax = atom->nmax;
  int reallocate_peratoms = false;
  current_peratom = atom->peratom;

  if (nmax > stored_nmax) {
    // reallocate tag array
    stored_nmax = nmax;
    memory->destroy(stored_tag);
    stored_tag = memory->create(stored_tag, stored_nmax, "hmc:stored_tag");
    reallocate_peratoms = true;
    // reallocate body peratom data
    if (rigid_flag) {
      memory->destroy(stored_bodyown);
      memory->destroy(stored_bodytag);
      memory->destroy(stored_atom2body);
      memory->destroy(stored_xcmimage);
      memory->destroy(stored_displace);
      memory->destroy(stored_eflags);
      memory->destroy(stored_orient);
      memory->destroy(stored_dorient);
      stored_bodyown = memory->create(stored_bodyown, stored_nmax, "hmc:stored_bodyown");
      stored_bodytag = memory->create(stored_bodytag, stored_nmax, "hmc:stored_bodytag");
      stored_atom2body = memory->create(stored_atom2body, stored_nmax, "hmc:stored_atom2body");
      stored_xcmimage = memory->create(stored_xcmimage, stored_nmax, "hmc:stored_xcmimage");
      stored_displace = memory->create(stored_displace, stored_nmax, 3, "hmc:stored_displace");
      if (fix_rigid->extended) {
        stored_eflags = memory->create(stored_eflags, nmax, "hmc:stored_eflags");
        if (fix_rigid->orientflag)
          stored_orient =
              memory->create(stored_orient, nmax, fix_rigid->orientflag, "hmc:stored_orient");
        if (fix_rigid->dorientflag)
          stored_dorient = memory->create(stored_dorient, nmax, 3, "hmc:stored_dorient");
      }
    }
  }

  // store tag array
  memcpy(stored_tag, atom->tag, ntotal * sizeof(tagint));

  // store body peratom data
  if (rigid_flag) {
    memcpy(stored_bodyown, fix_rigid->bodyown, ntotal * sizeof(int));
    memcpy(stored_bodytag, fix_rigid->bodytag, ntotal * sizeof(tagint));
    memcpy(stored_atom2body, fix_rigid->atom2body, ntotal * sizeof(int));
    memcpy(stored_xcmimage, fix_rigid->xcmimage, ntotal * sizeof(imageint));
    for (int i = 0; i < ntotal; i++)
      memcpy(stored_displace[i], fix_rigid->displace[i], 3 * sizeof(double));
    if (fix_rigid->extended) {
      memcpy(stored_eflags, fix_rigid->eflags, ntotal * sizeof(int));
      if (fix_rigid->orientflag)
        for (int i = 0; i < ntotal; i++)
          memcpy(stored_orient[i], fix_rigid->orient[i], fix_rigid->orientflag * sizeof(double));
      if (fix_rigid->dorientflag)
        for (int i = 0; i < ntotal; i++)
          memcpy(stored_dorient[i], fix_rigid->dorient[i], 3 * sizeof(double));
    }
  }

  // also reallocate if the number of peratoms has changed
  if (current_peratom.size() != stored_peratom.size()) reallocate_peratoms = true;

  // clear peratom data and store a new struct if reallocation, else just re-store
  if (reallocate_peratoms) {
    for (Atom::PerAtom &stored_peratom_member : stored_peratom) {
      free(stored_peratom_member.address);
      free(stored_peratom_member.address_maxcols);
    }
    stored_peratom.clear();
    Atom::PerAtom stored_peratom_member;
    for (Atom::PerAtom &current_peratom_member : current_peratom) {
      stored_peratom_member = current_peratom_member;
      switch (current_peratom_member.datatype) {
        case (Atom::INT):
          store_peratom_member<int>(stored_peratom_member, current_peratom_member, ntotal, nmax,
                                    reallocate_peratoms);
          break;
        case (Atom::DOUBLE):
          store_peratom_member<double>(stored_peratom_member, current_peratom_member, ntotal, nmax,
                                       reallocate_peratoms);
          break;
        case (Atom::BIGINT):
          store_peratom_member<bigint>(stored_peratom_member, current_peratom_member, ntotal, nmax,
                                       reallocate_peratoms);
          break;
      }
      stored_peratom.push_back(stored_peratom_member);
    }
  } else {
    for (Atom::PerAtom &stored_peratom_member : stored_peratom) {
      for (Atom::PerAtom &current_peratom_member : current_peratom) {
        if (stored_peratom_member.name.compare(current_peratom_member.name)) {
          continue;
        } else {
          switch (current_peratom_member.datatype) {
            case (Atom::INT):
              store_peratom_member<int>(stored_peratom_member, current_peratom_member, ntotal, nmax,
                                        reallocate_peratoms);
              break;
            case (Atom::DOUBLE):
              store_peratom_member<double>(stored_peratom_member, current_peratom_member, ntotal,
                                           nmax, reallocate_peratoms);
              break;
            case (Atom::BIGINT):
              store_peratom_member<bigint>(stored_peratom_member, current_peratom_member, ntotal,
                                           nmax, reallocate_peratoms);
              break;
          }
        }
      }
    }
  }

  // store totals
  stored_ntotal = ntotal;
  stored_nlocal = nlocal;
  stored_nmax = nmax;
  stored_nghost = atom->nghost;
  stored_nbonds = atom->nbonds;
  stored_nangles = atom->nangles;
  stored_ndihedrals = atom->ndihedrals;
  stored_nimpropers = atom->nimpropers;

  // store bodies
  if (rigid_flag) {
    stored_nlocal_body = fix_rigid->nlocal_body;
    stored_nghost_body = fix_rigid->nghost_body;
    stored_ntotal_body = stored_nlocal_body + stored_nghost_body;
    delete stored_body;
    stored_body = new RigidSmallBody[stored_ntotal_body];
    for (int i = 0; i < stored_ntotal_body; i++) stored_body[i] = fix_rigid->body[i];
  }

  // Save global energy terms:
  for (m = 0; m < neg; m++) eglobal[m] = *eglobalptr[m];

  // Save global virial terms:
  if (press_flag)
    for (m = 0; m < nv; m++) memcpy(vglobal[m], *vglobalptr[m], six);

  // Perform reverse communication to incorporate ghost atoms info:
  if (comm_reverse && (peatom_flag || pressatom_flag)) {
    comm_flag = ATOMS;
    comm->reverse_comm(this, ncommrev);
  }
}

/* ----------------------------------------------------------------------
   Restore the system state saved at the beginning of the MC step
------------------------------------------------------------------------- */

void FixHMC::restore_saved_state()
{
  int i;
  int nlocal = atom->nlocal;
  int ntotal = nlocal + atom->nghost;

  current_peratom = atom->peratom;

  // restore tag and peratom body data
  memcpy(atom->tag, stored_tag, stored_ntotal * sizeof(tagint));

  if (rigid_flag) {
    memcpy(fix_rigid->bodyown, stored_bodyown, ntotal * sizeof(int));
    memcpy(fix_rigid->bodytag, stored_bodytag, ntotal * sizeof(tagint));
    memcpy(fix_rigid->atom2body, stored_atom2body, ntotal * sizeof(int));
    memcpy(fix_rigid->xcmimage, stored_xcmimage, ntotal * sizeof(imageint));
    for (int i = 0; i < ntotal; i++)
      memcpy(fix_rigid->displace[i], stored_displace[i], 3 * sizeof(double));
    if (fix_rigid->extended) {
      memcpy(fix_rigid->eflags, stored_eflags, ntotal * sizeof(int));
      if (fix_rigid->orientflag)
        for (int i = 0; i < ntotal; i++)
          memcpy(fix_rigid->orient[i], stored_orient[i], fix_rigid->orientflag * sizeof(double));
      if (fix_rigid->dorientflag)
        for (int i = 0; i < ntotal; i++)
          memcpy(fix_rigid->dorient[i], stored_dorient[i], 3 * sizeof(double));
    }
  }

  if (stored_ntotal > atom->nlocal + atom->nghost) { atom->avec->grow(stored_ntotal); }

  // restore counts
  atom->nlocal = stored_nlocal;
  atom->nghost = stored_nghost;
  atom->nbonds = stored_nbonds;
  atom->nangles = stored_nangles;
  atom->ndihedrals = stored_ndihedrals;
  atom->nimpropers = stored_nimpropers;

  // restore peratom data
  for (Atom::PerAtom &stored_peratom_member : stored_peratom) {
    for (Atom::PerAtom &current_peratom_member : current_peratom) {
      if (stored_peratom_member.name.compare(current_peratom_member.name)) {
        continue;
      } else {
        switch (current_peratom_member.datatype) {
          case (Atom::INT):
            restore_peratom_member<int>(stored_peratom_member, current_peratom_member, ntotal);
            break;
          case (Atom::DOUBLE):
            restore_peratom_member<double>(stored_peratom_member, current_peratom_member, ntotal);
            break;
          case (Atom::BIGINT):
            restore_peratom_member<bigint>(stored_peratom_member, current_peratom_member, ntotal);
            break;
        }
        break;
      }
    }
  }

  // restore bodies
  if (rigid_flag) {
    fix_rigid->nlocal_body = stored_nlocal_body;
    fix_rigid->nghost_body = stored_nghost_body;
    for (int i = 0; i < stored_ntotal_body; i++) fix_rigid->body[i] = stored_body[i];
  }

  // reinit atom_map
  if (atom->map_style != Atom::MAP_NONE) {
    atom->map_clear();
    atom->map_init();
    atom->map_set();
  }

  // Restore global energy terms:
  for (i = 0; i < neg; i++) *eglobalptr[i] = eglobal[i];

  // Restore global virial terms:
  if (press_flag)
    for (i = 0; i < nv; i++) memcpy(*vglobalptr[i], vglobal[i], six);
}

/* ----------------------------------------------------------------------
   Randomly choose velocities from a Maxwell-Boltzmann distribution
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

------------------------------------------------------------------------- */

void FixHMC::rigid_body_random_velocities()
{
  RigidSmallBody *body = fix_rigid->body;
  int nlocal = fix_rigid->nlocal_body;
  int ntotal = nlocal + fix_rigid->nghost_body;
  int *mask = atom->mask;

  double stdev, bmass, wbody[3], mbody[3];
  double total_mass = 0;
  RigidSmallBody *b;
  double vcm[] = {0.0, 0.0, 0.0};

  for (int ibody = 0; ibody < nlocal; ibody++) {
    b = &body[ibody];
    if (mask[b->ilocal] & groupbit) {
      bmass = b->mass;
      stdev = sqrt(KT / bmass);
      total_mass += bmass;
      for (int j = 0; j < 3; j++) {
        b->vcm[j] = stdev * random->gaussian();
        vcm[j] += b->vcm[j] * bmass;
        if (b->inertia[j] > 0.0)
          wbody[j] = sqrt(KT / b->inertia[j]) * random->gaussian();
        else
          wbody[j] = 0.0;
      }
    }
    MathExtra::matvec(b->ex_space, b->ey_space, b->ez_space, wbody, b->omega);
  }

  if (mom_flag) {
    for (int j = 0; j < 3; j++) vcm[j] /= total_mass;
    for (int ibody = 0; ibody < nlocal; ibody++) {
      if (mask[b->ilocal] & groupbit) {
        b = &body[ibody];
        for (int j = 0; j < 3; j++) b->vcm[j] -= vcm[j];
      }
    }
  }

  // Forward communicate vcm and omega to ghost bodies:
  comm_flag = VCM_OMEGA;
  comm->forward_comm(this, 6);

  // Compute angular momenta of rigid bodies:
  for (int ibody = 0; ibody < ntotal; ibody++) {
    b = &body[ibody];
    MathExtra::omega_to_angmom(b->omega, b->ex_space, b->ey_space, b->ez_space, b->inertia,
                               b->angmom);
    MathExtra::transpose_matvec(b->ex_space, b->ey_space, b->ez_space, b->angmom, mbody);
    MathExtra::quatvec(b->quat, mbody, b->conjqm);
    b->conjqm[0] *= 2.0;
    b->conjqm[1] *= 2.0;
    b->conjqm[2] *= 2.0;
    b->conjqm[3] *= 2.0;
  }

  // Compute velocities of individual atoms:
  fix_rigid->set_v();
}

/* ----------------------------------------------------------------------
   Pack rigid body info for forward communication
------------------------------------------------------------------------- */

int FixHMC::pack_forward_comm(int n, int *list, double *buf, int /*pbc_flag*/, int * /*pbc*/)
{
  int *bodyown = fix_rigid->bodyown;
  RigidSmallBody *body = fix_rigid->body;

  int i, m, ibody;
  RigidSmallBody *b;

  m = 0;
  for (i = 0; i < n; i++) {
    ibody = bodyown[list[i]];
    if (ibody >= 0) {
      b = &body[ibody];
      if (comm_flag == VCM_OMEGA) {
        memcpy(&buf[m], b->vcm, three);
        memcpy(&buf[m + 3], b->omega, three);
        m += 6;
      } else if (comm_flag == XCM) {
        memcpy(&buf[m], b->xcm, three);
        m += 3;
      } else if (comm_flag == ROTATION) {
        memcpy(&buf[m], b->ex_space, three);
        memcpy(&buf[m + 3], b->ey_space, three);
        memcpy(&buf[m + 6], b->ez_space, three);
        memcpy(&buf[m + 9], b->quat, four);
        m += 12;
      }
    }
  }
  return m;
}

/* ----------------------------------------------------------------------
   Unpack rigid body info from forward communication
------------------------------------------------------------------------- */

void FixHMC::unpack_forward_comm(int n, int first, double *buf)
{
  int *bodyown = fix_rigid->bodyown;
  RigidSmallBody *body = fix_rigid->body;

  int i, m, last;
  RigidSmallBody *b;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++)
    if (bodyown[i] >= 0) {
      b = &body[bodyown[i]];
      if (comm_flag == VCM_OMEGA) {
        memcpy(b->vcm, &buf[m], three);
        memcpy(b->omega, &buf[m + 3], three);
        m += 6;
      }
    }
}

/* ----------------------------------------------------------------------
   Pack per-atom energies and/or virials or
   rigid body info for reverse communication
------------------------------------------------------------------------- */

int FixHMC::pack_reverse_comm(int n, int first, double *buf)
{
  int last = first + n;

  int i, k, m;

  m = 0;
  if (comm_flag == ATOMS) {
    for (i = first; i < last; i++)
      for (k = 0; k < ne; k++)
        if (rev_comm[k]) {
          if (peatom_flag) buf[m++] = eatom[k][i];
          if (pressatom_flag) {
            memcpy(&buf[m], vatom[k][i], six);
            m += 6;
          }
        }
  }
  return m;
}

/* ----------------------------------------------------------------------
   Unpack per-atom energies and/or virials or
   rigid body info for reverse communication
------------------------------------------------------------------------- */

void FixHMC::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i, j, k, m;

  m = 0;
  if (comm_flag == ATOMS) {
    for (j = 0; j < n; j++) {
      i = list[j];
      for (k = 0; k < ne; k++)
        if (rev_comm[k]) {
          if (peatom_flag) eatom[k][i] += buf[m++];
          if (pressatom_flag) vatom[k][i][0] += buf[m++];
          vatom[k][i][1] += buf[m++];
          vatom[k][i][2] += buf[m++];
          vatom[k][i][3] += buf[m++];
          vatom[k][i][4] += buf[m++];
          vatom[k][i][5] += buf[m++];
        }
    }
  }
}

/* ----------------------------------------------------------------------
   allocate atom-based arrays
------------------------------------------------------------------------- */

void FixHMC::grow_arrays(int nmax)
{
  memory->grow(eatom, ne, nmax, "fix_hmc:eatom");
  memory->grow(vatom, nv, nmax, 6, "fix_hmc:vatom");
}

/* ----------------------------------------------------------------------
   copy values within local atom-based arrays
------------------------------------------------------------------------- */

void FixHMC::copy_arrays(int i, int j, int delflag)
{
  int m;

  if (peatom_flag)
    for (m = 0; m < ne; m++) eatom[m][j] = eatom[m][i];

  if (pressatom_flag)
    for (m = 0; m < nv; m++) memcpy(vatom[m][j], vatom[m][i], six);
}

/* ----------------------------------------------------------------------
   pack values in local atom-based array for exchange with another proc
------------------------------------------------------------------------- */

int FixHMC::pack_exchange(int i, double *buf)
{
  int k, m = 0;

  if (peatom_flag)
    for (k = 0; k < ne; k++) buf[m++] = eatom[k][i];

  if (pressatom_flag)
    for (k = 0; k < nv; k++) {
      memcpy(&buf[m], vatom[k][i], six);
      m += 6;
    }

  return m;
}

/* ----------------------------------------------------------------------
   unpack values in local atom-based array from exchange with another proc
------------------------------------------------------------------------- */

int FixHMC::unpack_exchange(int i, double *buf)
{
  int k, m = 0;

  if (peatom_flag)
    for (k = 0; k < ne; k++) eatom[k][i] = buf[m++];

  if (pressatom_flag)
    for (k = 0; k < nv; k++) {
      memcpy(vatom[k][i], &buf[m], six);
      m += 6;
    }

  return m;
}

/* ----------------------------------------------------------------------
   memory usage of stored_peratom_member
------------------------------------------------------------------------- */

template <typename T>
double FixHMC::memory_usage_peratom_member(Atom::PerAtom &stored_peratom_member)
{
  int cols;
  if (stored_peratom_member.address == nullptr) return 0;
  if (stored_peratom_member.cols == 0) { return stored_ntotal * sizeof(T); };
  if (stored_peratom_member.cols < 0)
    cols = *(stored_peratom_member.address_maxcols);
  else
    cols = stored_peratom_member.cols;
  return sizeof(T) * stored_ntotal * cols;
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double FixHMC::memory_usage()
{
  double bytes = 0;
  bytes += stored_nmax * sizeof(int);    // tag
  if (rigid_flag) {
    bytes += stored_nmax * sizeof(int);           // stored_bodyown
    bytes += stored_nmax * sizeof(tagint);        // stored_bodytag
    bytes += stored_nmax * sizeof(int);           // stored_atom2body
    bytes += stored_nmax * sizeof(imageint);      // stored_xcmimage
    bytes += stored_nmax * 3 * sizeof(double);    // stored_displace
    if (fix_rigid->extended) {
      bytes += stored_nmax * sizeof(int);    // stored_eflags
      if (fix_rigid->orientflag)
        bytes += stored_nmax * fix_rigid->orientflag * sizeof(double);          // stored_orient
      if (fix_rigid->dorientflag) bytes += stored_nmax * 3 * sizeof(double);    // stored_dorient
    }
  }

  for (Atom::PerAtom &stored_peratom_member : stored_peratom) {
    switch (stored_peratom_member.datatype) {
      case (Atom::INT):
        bytes += memory_usage_peratom_member<int>(stored_peratom_member);
        break;
      case (Atom::DOUBLE):
        bytes += memory_usage_peratom_member<double>(stored_peratom_member);
        break;
      case (Atom::BIGINT):
        bytes += memory_usage_peratom_member<bigint>(stored_peratom_member);
        break;
    }
  }

  bytes += stored_ntotal_body * sizeof(RigidSmallBody);

  bytes += nvalues * atom->nmax * sizeof(double);
  return bytes;
}

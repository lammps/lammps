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

/* ----------------------------------------------------------------------
   Contributing author: Tine Curk (tcurk5@gmail.com) and Jiaxing Yuan (yuanjiaxing123@hotmail.com)
------------------------------------------------------------------------- */

#include "fix_charge_regulation.h"

#include "angle.h"
#include "atom.h"
#include "atom_vec.h"
#include "bond.h"
#include "citeme.h"
#include "comm.h"
#include "compute.h"
#include "dihedral.h"
#include "domain.h"
#include "error.h"
#include "fix.h"
#include "force.h"
#include "group.h"
#include "improper.h"
#include "input.h"
#include "kspace.h"
#include "math_const.h"
#include "math_special.h"
#include "memory.h"
#include "modify.h"
#include "neighbor.h"
#include "pair.h"
#include "random_park.h"
#include "update.h"
#include "variable.h"

#include <cmath>
#include <cstring>
#include <memory>

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;
using namespace MathSpecial;

static const char cite_fix_charge_regulation[] =
  "fix charge/regulation: doi:10.1063/5.0066432\n\n"
  "@Article{Curk22,\n"
  " author = {T. Curk and J. Yuan and E. Luijten},\n"
  " title = {Accelerated Simulation Method for Charge Regulation Effects},\n"
  " journal = {Journal of Chemical Physics},\n"
  " year = 2022,\n"
  " volume = 156\n"
  "}\n\n";

enum{CONSTANT,EQUAL}; // parsing input variables

// large energy value used to signal overlap
#define MAXENERGYSIGNAL 1.0e100
#define MAXENERGYTEST 1.0e50
#define SMALL 0.0000001
#define NA_RHO0 0.602214 // Avogadro's constant times reference concentration  (N_A * mol / liter)  [nm^-3]

/* ---------------------------------------------------------------------- */

FixChargeRegulation::FixChargeRegulation(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  ngroups(0), groupstrings(nullptr), ptype_ID(nullptr), pHstr(nullptr),
  c_pe(nullptr), random_equal(nullptr), random_unequal(nullptr),
  idftemp(nullptr)
{
  if (lmp->citeme) lmp->citeme->add(cite_fix_charge_regulation);

  // Region restrictions not yet implemented ..

  vector_flag = 1;
  size_vector = 8;
  global_freq = 1;
  extvector = 0;
  restart_global = 1;
  time_depend = 1;
  cr_nmax = 0;
  overlap_flag = 0;
  energy_stored = 0;

  // necessary to specify the free ion types
  cation_type = utils::inumeric(FLERR, arg[3], false, lmp);
  anion_type = utils::inumeric(FLERR, arg[4], false, lmp);

  // set defaults and read optional arguments
  options(narg - 5, &arg[5]);

  if ((nevery <= 0) || (nmc < 0) || (llength_unit_in_nm < 0.0)
      || (*target_temperature_tcp < 0.0) || (cation_type <= 0)
      || (anion_type <= 0) || (reaction_distance < 0.0)
      || (salt_charge[0] <= 0) || (salt_charge[1] >= 0))
    error->all(FLERR, "Illegal fix charge/regulation command");

  if (seed <= 0)
    error->all(FLERR, "Illegal fix charge/regulation command: "
               "Seed value (positive integer) must be provided ");
  if ((salt_charge[1] % salt_charge[0] != 0)
      && (salt_charge[0] % salt_charge[1] != 0))
    error->all(FLERR,"Illegal fix charge/regulation command, "
               "multivalent cation/anion charges are allowed, "
               "but must be divisible, e.g. (3,-1) is fine, "
               "but (3,-2) is not implemented");

  if (pmcmoves[0] < 0 || pmcmoves[1] < 0 || pmcmoves[2] < 0)
    error->all(FLERR, "Illegal fix charge/regulation command");
  if (acid_type < 0) pmcmoves[0] = 0;
  if (base_type < 0) pmcmoves[1] = 0;

  // normalize
  double psum = pmcmoves[0] + pmcmoves[1] + pmcmoves[2];
  if (psum <= 0) error->all(FLERR, "Illegal fix charge/regulation command");
  pmcmoves[0] /= psum;
  pmcmoves[1] /= psum;
  pmcmoves[2] /= psum;

  force_reneighbor = 1;
  next_reneighbor = update->ntimestep + 1;
  random_equal = new RanPark(lmp, seed);
  random_unequal = new RanPark(lmp, seed);
  nacid_attempts = 0;
  nacid_successes = 0;
  nbase_attempts = 0;
  nbase_successes = 0;
  nsalt_attempts = 0;
  nsalt_successes = 0;
}

FixChargeRegulation::~FixChargeRegulation() {

  memory->destroy(ptype_ID);

  delete random_equal;
  delete random_unequal;
  delete[] pHstr;
  delete[] idftemp;

  if (group) {
    int igroupall = group->find("all");
    neighbor->exclusion_group_group_delete(exclusion_group, igroupall);
  }
}

int FixChargeRegulation::setmask() {
  int mask = 0;
  mask |= PRE_EXCHANGE;
  return mask;
}

void FixChargeRegulation::init() {

  triclinic = domain->triclinic;
  int ipe = modify->find_compute("thermo_pe");
  c_pe = modify->compute[ipe];

  if (pHstr) {
    pHvar = input->variable->find(pHstr);
    if (pHvar < 0)
      error->all(FLERR,"Variable name for fix charge/regulation does not exist");
    if (input->variable->equalstyle(pHvar)) pHstyle = EQUAL;
    else error->all(FLERR,"Variable for fix charge/regulation is invalid style");

  }
  if (atom->molecule_flag) {

    int flag = 0;
    for (int i = 0; i < atom->nlocal; i++)
      if (atom->type[i] == cation_type || atom->type[i] == anion_type)
        if (atom->molecule[i]) flag = 1;
    int flagall = flag;

    MPI_Allreduce(&flag, &flagall, 1, MPI_INT, MPI_SUM, world);
    if (flagall && comm->me == 0)
      error->all(FLERR, "fix charge/regulation cannot exchange "
                 "individual atoms (ions) belonging to a molecule");
  }

  if (domain->dimension == 2)
    error->all(FLERR, "Cannot use fix charge/regulation in a 2d simulation");

  // create a new group for interaction exclusions
  // used for attempted atom deletions
  // skip if already exists from previous init()

  if (!exclusion_group_bit) {

    // create unique group name for atoms to be excluded

    auto group_id = fmt::format("FixChargeRegulation:exclusion_group:{}",id);
    group->assign(group_id + " subtract all all");
    exclusion_group = group->find(group_id);
    if (exclusion_group == -1)
      error->all(FLERR,"Could not find fix charge/regulation exclusion "
                 "group ID");
    exclusion_group_bit = group->bitmask[exclusion_group];

    // neighbor list exclusion setup
    // turn off interactions between group all and the exclusion group

    neighbor->modify_params(fmt::format("exclude group {} all",group_id));
  }

  // check that no deletable atoms are in atom->firstgroup
  // deleting such an atom would not leave firstgroup atoms first

  if (atom->firstgroup >= 0) {
    int *mask = atom->mask;
    int firstgroupbit = group->bitmask[atom->firstgroup];

    int flag = 0;
    for (int i = 0; i < atom->nlocal; i++)
      if ((mask[i] == groupbit) && (mask[i] && firstgroupbit)) flag = 1;

    int flagall;
    MPI_Allreduce(&flag, &flagall, 1, MPI_INT, MPI_SUM, world);

    if (flagall)
      error->all(FLERR, "Cannot use fix charge/regulation on atoms "
                 "in atom_modify first group");
  }

  // construct group bitmask for all new atoms
  // aggregated over all group keywords

  groupbitall = 1 | groupbit;

  for (int igroup = 0; igroup < ngroups; igroup++) {
    int jgroup = group->find(groupstrings[igroup]);
    if (jgroup == -1)
      error->all(FLERR, "Could not find fix charge/regulation group ID");
    groupbitall |= group->bitmask[jgroup];
  }
}

void FixChargeRegulation::pre_exchange() {

  if (next_reneighbor != update->ntimestep) return;
  xlo = domain->boxlo[0];
  xhi = domain->boxhi[0];
  ylo = domain->boxlo[1];
  yhi = domain->boxhi[1];
  zlo = domain->boxlo[2];
  zhi = domain->boxhi[2];

  if (triclinic) {
    sublo = domain->sublo_lamda;
    subhi = domain->subhi_lamda;
  } else {
    sublo = domain->sublo;
    subhi = domain->subhi;
  }
  volume = domain->xprd * domain->yprd * domain->zprd;
  if (triclinic) domain->x2lamda(atom->nlocal);
  domain->pbc();
  comm->exchange();
  atom->nghost = 0;
  comm->borders();

  if (triclinic) domain->lamda2x(atom->nlocal + atom->nghost);
  energy_stored = energy_full();
  if ((energy_stored > MAXENERGYTEST) && (comm->me == 0))
    error->warning(FLERR, "Energy of old configuration in fix "
                   "charge/regulation is > MAXENERGYTEST.");

  if ((reaction_distance > domain->prd_half[0]) ||
      (reaction_distance > domain->prd_half[1]) ||
      (reaction_distance > domain->prd_half[2])) {
    if (comm->me == 0)
      error->warning(FLERR,"reaction distance (rxd) is larger than "
                     "half the box dimension, resetting default: xrd = 0.");
    reaction_distance = 0;
  }
  // volume in units of (N_A * mol / liter)
  volume_rx = (xhi - xlo) * (yhi - ylo) * (zhi - zlo)
    * cube(llength_unit_in_nm) * NA_RHO0;
  if (reaction_distance < SMALL) {
    vlocal_xrd = volume_rx;
  } else {
    vlocal_xrd = 4.0 * MY_PI * cube(reaction_distance)
      / 3.0 * cube(llength_unit_in_nm) * NA_RHO0;
  }
  beta = 1.0 / (force->boltz * *target_temperature_tcp);

  if (pHstyle == EQUAL)
    pH = input->variable->compute_equal(pHvar);

  // pre-compute powers
  c10pH = pow(10.0,-pH); // dissociated ion (H+) activity
  c10pKa = pow(10.0,-pKa); // acid dissociation constant
  c10pKb = pow(10.0,-pKb); // base dissociation constant
  c10pOH = pow(10.0,-pKs + pH); // dissociated anion (OH-) activity
  c10pI_plus = pow(10.0,-pI_plus); // free cation activity
  c10pI_minus = pow(10.0,-pI_minus); // free anion activity

  // reinitialize counters
  nacid_neutral = particle_number(acid_type, 0);
  nacid_charged = particle_number(acid_type, -1);
  nbase_neutral = particle_number(base_type, 0);
  nbase_charged = particle_number(base_type, 1);
  ncation = particle_number(cation_type, salt_charge[0]);
  nanion = particle_number(anion_type, salt_charge[1]);


  // Attempt exchanges
  if (!only_salt_flag) {

    // Do charge regulation
    for (int i = 0; i < nmc; i++) {
      double rand_number = random_equal->uniform();
      if (rand_number < pmcmoves[0] / 2) {
        forward_acid();
        nacid_attempts++;
      } else if (rand_number < pmcmoves[0]) {
        backward_acid();
        nacid_attempts++;
      } else if (rand_number < pmcmoves[0] + pmcmoves[1] / 2) {
        forward_base();
        nbase_attempts++;
      } else if (rand_number < pmcmoves[0] + pmcmoves[1]) {
        backward_base();
        nbase_attempts++;
      } else if (rand_number < pmcmoves[0] + pmcmoves[1] + pmcmoves[2] / 2) {
        forward_ions();
        nsalt_attempts++;
      } else {
        backward_ions();
        nsalt_attempts++;
      }
    }
  } else {
    // do only ion insertion, multivalent cation/anions are implemented
    if (salt_charge[0] >= -salt_charge[1]) {
      salt_charge_ratio = -salt_charge[0] / salt_charge[1];
    } else {
      salt_charge_ratio = -salt_charge[1] / salt_charge[0];
    }
    for (int i = 0; i < nmc; i++) {
      double rand_number = random_equal->uniform();
      if (rand_number < 0.5) {
        forward_ions_multival();
        nsalt_attempts++;
      } else {
        backward_ions_multival();
        nsalt_attempts++;
      }
    }
  }

  // assign unique tags to newly inserted ions
  if (add_tags_flag && atom->tag_enable) assign_tags();

  if (triclinic) domain->x2lamda(atom->nlocal);
  domain->pbc();
  comm->exchange();
  atom->nghost = 0;
  comm->borders();
  if (triclinic) domain->lamda2x(atom->nlocal + atom->nghost);
  next_reneighbor = update->ntimestep + nevery;
}

void FixChargeRegulation::forward_acid() {

  double energy_before = energy_stored;
  double factor;
  double dummyp[3];
  double pos[3];
  pos[0] = 0;
  pos[1] = 0;
  pos[2] = 0; // acid/base particle position
  double pos_all[3];
  int m1 = -1, m2 = -1;

  m1 = get_random_particle(acid_type, 0, 0, dummyp);
  if (npart_xrd != nacid_neutral) error->all(FLERR, "fix charge/regulation acid count inconsistent");

  if (nacid_neutral > 0) {
    if (m1 >= 0) {
      atom->q[m1] = -1; // assign negative charge to acid
      pos[0] = atom->x[m1][0];
      pos[1] = atom->x[m1][1];
      pos[2] = atom->x[m1][2];
    }
    npart_xrd2 = ncation;
    if (reaction_distance >= SMALL) {
      pos_all[0] = pos[0];
      pos_all[1] = pos[1];
      pos_all[2] = pos[2];
      MPI_Allreduce(pos, pos_all, 3, MPI_DOUBLE, MPI_SUM, world);
      npart_xrd2 = particle_number_xrd(cation_type, 1, reaction_distance, pos_all);
    }
    m2 = insert_particle(cation_type, 1, reaction_distance, pos_all);
    factor = nacid_neutral * vlocal_xrd * c10pKa * c10pI_plus /
            (c10pH * (1 + nacid_charged) * (1 + npart_xrd2));

    if (force->kspace) force->kspace->qsum_qsq();
    if (force->pair->tail_flag) force->pair->reinit();
    double energy_after = energy_full();

    if (energy_after < MAXENERGYTEST &&
        random_equal->uniform() < factor * exp(beta * (energy_before - energy_after))) {
      energy_stored = energy_after;
      nacid_successes += 1;
      ncation++;
      nacid_charged++;
      nacid_neutral--;
    } else {
      energy_stored = energy_before;
      atom->natoms--;
      if (m2 >= 0) {
        atom->nlocal--;
      }
      if (m1 >= 0) {
        atom->q[m1] = 0;
      }
      if (force->kspace) force->kspace->qsum_qsq();
      if (force->pair->tail_flag) force->pair->reinit();
    }
  }
}

void FixChargeRegulation::backward_acid() {

  double energy_before = energy_stored;
  double factor;
  int mask_tmp;
  double dummyp[3];
  double pos[3];
  pos[0] = 0;
  pos[1] = 0;
  pos[2] = 0; // acid/base particle position
  double pos_all[3];
  int m1 = -1, m2 = -1;

  m1 = get_random_particle(acid_type, -1, 0, dummyp);
  if (npart_xrd != nacid_charged)
    error->all(FLERR, "fix charge/regulation acid count inconsistent");

  if (nacid_charged > 0) {
    if (m1 >= 0) {
      atom->q[m1] = 0;
      pos[0] = atom->x[m1][0];
      pos[1] = atom->x[m1][1];
      pos[2] = atom->x[m1][2];
    }
    if (reaction_distance >= SMALL) {
      pos_all[0] = pos[0];
      pos_all[1] = pos[1];
      pos_all[2] = pos[2];
      MPI_Allreduce(pos, pos_all, 3, MPI_DOUBLE, MPI_SUM, world);
    }
    m2 = get_random_particle(cation_type, 1, reaction_distance, pos_all);
    // note: npart_xrd changes everytime get_random_particle is called.

    if (npart_xrd > 0) {
      if (m2 >= 0) {
        atom->q[m2] = 0;
        mask_tmp = atom->mask[m2];  // remember group bits.
        atom->mask[m2] = exclusion_group_bit;
      }
      factor = (1 + nacid_neutral) * vlocal_xrd * c10pKa * c10pI_plus  /
              (c10pH * nacid_charged * npart_xrd);

      if (force->kspace) force->kspace->qsum_qsq();
      if (force->pair->tail_flag) force->pair->reinit();
      double energy_after = energy_full();

      if (energy_after < MAXENERGYTEST &&
          random_equal->uniform() < (1.0 / factor) * exp(beta * (energy_before - energy_after))) {
        nacid_successes += 1;
        atom->natoms--;
        energy_stored = energy_after;
        nacid_charged--;
        nacid_neutral++;
        ncation--;
        if (m2 >= 0) {
          atom->avec->copy(atom->nlocal - 1, m2, 1);
          atom->nlocal--;
        }
      } else {
        energy_stored = energy_before;
        if (m1 >= 0) {
          atom->q[m1] = -1;
        }
        if (m2 >= 0) {
          atom->q[m2] = 1;
          atom->mask[m2] = mask_tmp;
        }
        if (force->kspace) force->kspace->qsum_qsq();
        if (force->pair->tail_flag) force->pair->reinit();
      }
    } else {
      if (m1 >= 0) {
        atom->q[m1] = -1;
      }
    }
  }
}

void FixChargeRegulation::forward_base() {

  double energy_before = energy_stored;
  double factor;
  double dummyp[3];
  double pos[3];
  pos[0] = 0;
  pos[1] = 0;
  pos[2] = 0; // acid/base particle position
  double pos_all[3];
  int m1 = -1, m2 = -1;

  m1 = get_random_particle(base_type, 0, 0, dummyp);
  if (npart_xrd != nbase_neutral)
    error->all(FLERR, "fix charge/regulation acid count inconsistent");

  if (nbase_neutral > 0) {
    if (m1 >= 0) {
      atom->q[m1] = 1; // assign negative charge to acid
      pos[0] = atom->x[m1][0];
      pos[1] = atom->x[m1][1];
      pos[2] = atom->x[m1][2];
    }
    npart_xrd2 = nanion;
    if (reaction_distance >= SMALL) {
      pos_all[0] = pos[0];
      pos_all[1] = pos[1];
      pos_all[2] = pos[2];
      MPI_Allreduce(pos, pos_all, 3, MPI_DOUBLE, MPI_SUM, world);
      npart_xrd2 = particle_number_xrd(anion_type, -1, reaction_distance, pos_all);
    }
    factor = nbase_neutral * vlocal_xrd * c10pKb * c10pI_minus /
             (c10pOH * (1 + nbase_charged) * (1 + npart_xrd2));
    m2 = insert_particle(anion_type, -1, reaction_distance, pos_all);

    if (force->kspace) force->kspace->qsum_qsq();
    if (force->pair->tail_flag) force->pair->reinit();
    double energy_after = energy_full();
    if (energy_after < MAXENERGYTEST &&
        random_equal->uniform() < factor * exp(beta * (energy_before - energy_after))) {
      energy_stored = energy_after;
      nbase_successes += 1;
      nbase_charged++;
      nbase_neutral--;
      nanion++;
    } else {
      energy_stored = energy_before;
      atom->natoms--;
      if (m2 >= 0) {
        atom->nlocal--;
      }
      if (m1 >= 0) {
        atom->q[m1] = 0;
      }
      if (force->kspace) force->kspace->qsum_qsq();
      if (force->pair->tail_flag) force->pair->reinit();
    }
  }
}

void FixChargeRegulation::backward_base() {

  double energy_before = energy_stored;
  double factor;
  double dummyp[3];
  int mask_tmp;
  double pos[3];
  pos[0] = 0;
  pos[1] = 0;
  pos[2] = 0; // acid/base particle position
  double pos_all[3];
  int m1 = -1, m2 = -1;

  m1 = get_random_particle(base_type, 1, 0, dummyp);
  if (npart_xrd != nbase_charged)
    error->all(FLERR, "fix charge/regulation acid count inconsistent");

  if (nbase_charged > 0) {
    if (m1 >= 0) {
      atom->q[m1] = 0;
      pos[0] = atom->x[m1][0];
      pos[1] = atom->x[m1][1];
      pos[2] = atom->x[m1][2];
    }
    if (reaction_distance >= SMALL) {
      pos_all[0] = pos[0];
      pos_all[1] = pos[1];
      pos_all[2] = pos[2];
      MPI_Allreduce(pos, pos_all, 3, MPI_DOUBLE, MPI_SUM, world);
    }
    m2 = get_random_particle(anion_type, -1, reaction_distance, pos_all);

    if (npart_xrd > 0) {
      if (m2 >= 0) {
        atom->q[m2] = 0;
        mask_tmp = atom->mask[m2];  // remember group bits.
        atom->mask[m2] = exclusion_group_bit;
      }
      factor = (1 + nbase_neutral) * vlocal_xrd * c10pKb * c10pI_minus /
              (c10pOH * nbase_charged * npart_xrd);

      if (force->kspace) force->kspace->qsum_qsq();
      if (force->pair->tail_flag) force->pair->reinit();
      double energy_after = energy_full();

      if (energy_after < MAXENERGYTEST &&
          random_equal->uniform() < (1.0 / factor) * exp(beta * (energy_before - energy_after))) {
        nbase_successes += 1;
        atom->natoms--;
        energy_stored = energy_after;
        nbase_charged--;
        nbase_neutral++;
        nanion--;
        if (m2 >= 0) {
          atom->avec->copy(atom->nlocal - 1, m2, 1);
          atom->nlocal--;
        }
      } else {
        energy_stored = energy_before;
        if (m1 >= 0) {
          atom->q[m1] = 1;
        }
        if (m2 >= 0) {
          atom->q[m2] = -1;
          atom->mask[m2] = mask_tmp;
        }
        if (force->kspace) force->kspace->qsum_qsq();
        if (force->pair->tail_flag) force->pair->reinit();
      }
    } else {
      if (m1 >= 0) {
        atom->q[m1] = 1;
      }
    }
  }
}

void FixChargeRegulation::forward_ions() {

  double energy_before = energy_stored;
  double factor;
  double dummyp[3];
  int m1 = -1, m2 = -1;
  factor = volume_rx * volume_rx * c10pI_plus * c10pI_minus /
           ((1 + ncation) * (1 + nanion));

  m1 = insert_particle(cation_type, +1, 0, dummyp);
  m2 = insert_particle(anion_type, -1, 0, dummyp);
  if (force->kspace) force->kspace->qsum_qsq();
  if (force->pair->tail_flag) force->pair->reinit();
  double energy_after = energy_full();
  if (energy_after < MAXENERGYTEST &&
      random_equal->uniform() < factor * exp(beta * (energy_before - energy_after))) {
    energy_stored = energy_after;
    nsalt_successes += 1;
    ncation++;
    nanion++;
  } else {
    energy_stored = energy_before;
    atom->natoms--;
    if (m1 >= 0) {
      atom->nlocal--;
    }
    atom->natoms--;
    if (m2 >= 0) {
      atom->nlocal--;
    }
    if (force->kspace) force->kspace->qsum_qsq();
    if (force->pair->tail_flag) force->pair->reinit();
  }
}


void FixChargeRegulation::backward_ions() {

  double energy_before = energy_stored;
  double factor;
  int mask1_tmp = 0, mask2_tmp = 0;
  double dummyp[3];
  int m1 = -1, m2 = -1;

  m1 = get_random_particle(cation_type, +1, 0, dummyp);
  if (npart_xrd != ncation)
    error->all(FLERR, "fix charge/regulation salt count inconsistent");
  if (ncation > 0) {
    m2 = get_random_particle(anion_type, -1, 0, dummyp);
    if (npart_xrd != nanion)
      error->all(FLERR, "fix charge/regulation salt count inconsistent");
    if (nanion > 0) {

      // attempt deletion

      if (m1 >= 0) {
        atom->q[m1] = 0;
        mask1_tmp = atom->mask[m1];
        atom->mask[m1] = exclusion_group_bit;
      }
      if (m2 >= 0) {
        atom->q[m2] = 0;
        mask2_tmp = atom->mask[m2];
        atom->mask[m2] = exclusion_group_bit;
      }
      factor = volume_rx * volume_rx * c10pI_plus * c10pI_minus / (ncation * nanion);

      if (force->kspace) force->kspace->qsum_qsq();
      if (force->pair->tail_flag) force->pair->reinit();
      double energy_after = energy_full();
      if (energy_after < MAXENERGYTEST &&
          random_equal->uniform() < (1.0 / factor) * exp(beta * (energy_before - energy_after))) {
        energy_stored = energy_after;
        nsalt_successes += 1;
        ncation--;
        nanion--;

        // ions must be deleted in order, otherwise index m could change upon the first deletion
        atom->natoms -= 2;
        if (m1 > m2) {
          if (m1 >= 0) {
            atom->avec->copy(atom->nlocal - 1, m1, 1);
            atom->nlocal--;
          }
          if (m2 >= 0) {
            atom->avec->copy(atom->nlocal - 1, m2, 1);
            atom->nlocal--;
          }
        } else {
          if (m2 >= 0) {
            atom->avec->copy(atom->nlocal - 1, m2, 1);
            atom->nlocal--;
          }
          if (m1 >= 0) {
            atom->avec->copy(atom->nlocal - 1, m1, 1);
            atom->nlocal--;
          }
        }
      } else {
        energy_stored = energy_before;

        // reassign original charge and mask
        if (m1 >= 0) {
          atom->q[m1] = 1;
          atom->mask[m1] = mask1_tmp;
        }
        if (m2 >= 0) {
          atom->q[m2] = -1;
          atom->mask[m2] = mask2_tmp;
        }
        if (force->kspace) force->kspace->qsum_qsq();
        if (force->pair->tail_flag) force->pair->reinit();
      }
    }
  }
}

void FixChargeRegulation::forward_ions_multival() {

  double energy_before = energy_stored;
  double factor = 1;
  double dummyp[3];

  // particle ID array for all ions to be inserted
  auto mm = std::unique_ptr<int[]>(new int[salt_charge_ratio + 1]);

  if (salt_charge[0] <= -salt_charge[1]) {
    // insert one anion and (salt_charge_ratio) cations

    mm[0] = insert_particle(anion_type, salt_charge[1], 0, dummyp);
    factor *= volume_rx * c10pI_minus / (1 + nanion);
    for (int i = 0; i < salt_charge_ratio; i++) {
      mm[i + 1] = insert_particle(cation_type, salt_charge[0], 0, dummyp);
      factor *= volume_rx *c10pI_plus / (1 + ncation + i);
    }
  } else {
    // insert one cation and (salt_charge_ratio) anions

    mm[0] = insert_particle(cation_type, salt_charge[0], 0, dummyp);
    factor *= volume_rx * c10pI_plus / (1 + ncation);
    for (int i = 0; i < salt_charge_ratio; i++) {
      mm[i + 1] = insert_particle(anion_type, salt_charge[1], 0, dummyp);
      factor *= volume_rx * c10pI_minus / (1 + nanion + i);
    }
  }

  if (force->kspace) force->kspace->qsum_qsq();
  if (force->pair->tail_flag) force->pair->reinit();
  double energy_after = energy_full();
  if (energy_after < MAXENERGYTEST && random_equal->uniform() < factor * exp(beta * (energy_before - energy_after))) {
    energy_stored = energy_after;
    nsalt_successes += 1;

    if (salt_charge[0] <= -salt_charge[1]) {
      ncation += salt_charge_ratio;
      nanion++;
    } else {
      nanion += salt_charge_ratio;
      ncation++;
    }
  } else {
    energy_stored = energy_before;

    // delete inserted ions
    for (int i = 0; i < salt_charge_ratio + 1; i++) {
      atom->natoms--;
      if (mm[i] >= 0) {
        atom->nlocal--;
      }
    }
    if (force->kspace) force->kspace->qsum_qsq();
    if (force->pair->tail_flag) force->pair->reinit();
  }
}

void FixChargeRegulation::backward_ions_multival() {

  double energy_before = energy_stored;
  double factor = 1;
  double dummyp[3];  // dummy particle
  // particle ID array for all deleted ions
  auto mm = std::unique_ptr<int[]>(new int[salt_charge_ratio + 1]);
  // charge array for all deleted ions
  auto qq = std::unique_ptr<double[]>(new double[salt_charge_ratio + 1]);
  // temporary mask array
  auto mask_tmp = std::unique_ptr<int[]>(new int[salt_charge_ratio + 1]);

  if (salt_charge[0] <= -salt_charge[1]) {
    // delete one anion and (salt_charge_ratio) cations
    if (ncation < salt_charge_ratio || nanion < 1) return;

    mm[0] = get_random_particle(anion_type, salt_charge[1], 0, dummyp);
    if (npart_xrd != nanion)
      error->all(FLERR, "fix charge/regulation salt count inconsistent");
    factor *= volume_rx * c10pI_minus / (nanion);
    if (mm[0] >= 0) {
      qq[0] = atom->q[mm[0]];
      atom->q[mm[0]] = 0;
      mask_tmp[0] = atom->mask[mm[0]];
      atom->mask[mm[0]] = exclusion_group_bit;
    }
    for (int i = 0; i < salt_charge_ratio; i++) {
      mm[i + 1] = get_random_particle(cation_type, salt_charge[0], 0, dummyp);
      if (npart_xrd != ncation - i)
        error->all(FLERR, "fix charge/regulation salt count inconsistent");
      factor *= volume_rx * c10pI_plus / (ncation - i);
      if (mm[i + 1] >= 0) {
        qq[i + 1] = atom->q[mm[i + 1]];
        atom->q[mm[i + 1]] = 0;
        mask_tmp[i + 1] = atom->mask[mm[i + 1]];
        atom->mask[mm[i + 1]] = exclusion_group_bit;
      }
    }
  } else {
    // delete one cation and (salt_charge_ratio) anions

    if (nanion < salt_charge_ratio || ncation < 1) return;
    mm[0] = get_random_particle(cation_type, salt_charge[0], 0, dummyp);
    if (npart_xrd != ncation)
      error->all(FLERR, "fix charge/regulation salt count inconsistent");
    factor *= volume_rx * c10pI_plus / (ncation);
    if (mm[0] >= 0) {
      qq[0] = atom->q[mm[0]];
      atom->q[mm[0]] = 0;
      mask_tmp[0] = atom->mask[mm[0]];
      atom->mask[mm[0]] = exclusion_group_bit;
    }
    for (int i = 0; i < salt_charge_ratio; i++) {
      mm[i + 1] = get_random_particle(anion_type, salt_charge[1], 0, dummyp);
      if (npart_xrd != nanion - i)
        error->all(FLERR, "fix charge/regulation salt count inconsistent");
      if (mm[i + 1] >= 0) {
        qq[i + 1] = atom->q[mm[i + 1]];
        atom->q[mm[i + 1]] = 0;
        mask_tmp[i + 1] = atom->mask[mm[i + 1]];
        atom->mask[mm[i + 1]] = exclusion_group_bit;
      }
      factor *= volume_rx * c10pI_minus / (nanion - i);
    }
  }

  // attempt deletion

  if (force->kspace) force->kspace->qsum_qsq();
  if (force->pair->tail_flag) force->pair->reinit();
  double energy_after = energy_full();
  if (energy_after < MAXENERGYTEST &&
      random_equal->uniform() < (1.0 / factor) * exp(beta * (energy_before - energy_after))) {
    energy_stored = energy_after;

    atom->natoms -= 1 + salt_charge_ratio;
    // ions must be deleted in order, otherwise index m could change upon the first deletion
    for (int i = 0; i < salt_charge_ratio + 1; i++) {
      // get max mm value, poor N^2 scaling, but charge ratio is a SMALL number (2 or 3).
      int maxmm = -1, jmaxm = -1;
      for (int j = 0; j < salt_charge_ratio + 1; j++) {
        if (mm[j] > maxmm) {
          maxmm = mm[j];
          jmaxm = j;
        }
      }
      if (maxmm < 0) {
        break;  // already deleted all particles in this thread
      } else {
        // delete particle maxmm
        atom->avec->copy(atom->nlocal - 1, maxmm, 1);
        atom->nlocal--;
        mm[jmaxm] = -1;
      }
    }

    // update indices
    nsalt_successes += 1;
    if (salt_charge[0] <= -salt_charge[1]) {
      ncation -= salt_charge_ratio;
      nanion--;
    } else {
      nanion -= salt_charge_ratio;
      ncation--;
    }

  } else {
    energy_stored = energy_before;

    // reassign original charge and mask
    for (int i = 0; i < salt_charge_ratio + 1; i++) {
      if (mm[i] >= 0) {
        atom->q[mm[i]] = qq[i];
        atom->mask[mm[i]] = mask_tmp[i];
      }
    }
    if (force->kspace) force->kspace->qsum_qsq();
    if (force->pair->tail_flag) force->pair->reinit();
  }
}

int FixChargeRegulation::insert_particle(int ptype, double charge, double rd, double *target) {

  // insert a particle of type (ptype) with charge (charge) within distance (rd) of (target)

  double coord[3];
  int m = -1;
  if (rd < SMALL) {
    coord[0] = xlo + random_equal->uniform() * (xhi - xlo);
    coord[1] = ylo + random_equal->uniform() * (yhi - ylo);
    coord[2] = zlo + random_equal->uniform() * (zhi - zlo);
  } else {
    // get a random point inside a sphere with radius rd
    // simple rejection sampling, probably the fastest method
    double dxx=1,dyy=1,dzz=1;
    while (dxx * dxx + dyy * dyy + dzz * dzz > 1.0) {
      dxx = 2 * random_equal->uniform() - 1.0;
      dyy = 2 * random_equal->uniform() - 1.0;
      dzz = 2 * random_equal->uniform() - 1.0;
    }
    coord[0] = target[0] + rd * dxx;
    coord[1] = target[1] + rd * dyy;
    coord[2] = target[2] + rd * dzz;

    // Alternative way, but likely somewhat less efficient
    /*
    double radius = rd * pow(random_equal->uniform(), THIRD);
    double theta = acos(2 * random_equal->uniform() - 1);
    double phi = random_equal->uniform() * 2 * MY_PI;
    double sinphi = sin(phi);
    double cosphi = cos(phi);
    double sintheta = sin(theta);
    double costheta = cos(theta);
    coord[0] = target[0] + radius * sintheta * cosphi;
    coord[1] = target[1] + radius * sintheta * sinphi;
    coord[2] = target[2] + radius * costheta;
    */
    coord[0] = coord[0] - floor(1.0 * (coord[0] - xlo) / (xhi - xlo)) * (xhi - xlo);
    coord[1] = coord[1] - floor(1.0 * (coord[1] - ylo) / (yhi - ylo)) * (yhi - ylo);
    coord[2] = coord[2] - floor(1.0 * (coord[2] - zlo) / (zhi - zlo)) * (zhi - zlo);
  }

  if (coord[0] >= sublo[0] && coord[0] < subhi[0] &&
      coord[1] >= sublo[1] && coord[1] < subhi[1] &&
      coord[2] >= sublo[2] && coord[2] < subhi[2]) {
    atom->avec->create_atom(ptype, coord);
    m = atom->nlocal - 1;
    atom->mask[m] = groupbitall;

    sigma = sqrt(force->boltz * *target_temperature_tcp / atom->mass[ptype] / force->mvv2e);
    atom->v[m][0] = random_unequal->gaussian() * sigma;
    atom->v[m][1] = random_unequal->gaussian() * sigma;
    atom->v[m][2] = random_unequal->gaussian() * sigma;
    atom->q[m] = charge;
    modify->create_attribute(m);

  }
  atom->natoms++;
  atom->nghost = 0;
  if (atom->tag_enable) {
    if (atom->tag_enable) {
      atom->tag_extend();
      if (atom->map_style != Atom::MAP_NONE) atom->map_init();
    }
  }
  if (triclinic) domain->x2lamda(atom->nlocal);
  comm->borders();
  if (triclinic) domain->lamda2x(atom->nlocal+atom->nghost);
  return m;
}

int FixChargeRegulation::get_random_particle(int ptype, double charge, double rd, double *target) {

  // returns a randomly chosen particle of type (ptype) with charge (charge)
  // chosen among particles within distance (rd) of (target)

  int nlocal = atom->nlocal;

  // expand memory, if necessary
  if (atom->nmax > cr_nmax) {
    memory->sfree(ptype_ID);
    cr_nmax = atom->nmax;
    ptype_ID = (int *) memory->smalloc(cr_nmax * sizeof(int),
                                       "CR: local_atom_list");
  }

  int count_local, count_global, count_before;
  int m = -1;
  count_local = 0;
  count_global = 0;
  count_before = 0;

  if (rd < SMALL) {  //reaction_distance < SMALL: No constraint on random particle choice
    for (int i = 0; i < nlocal; i++) {
      if (atom->type[i] == ptype && fabs(atom->q[i] - charge) < SMALL &&
          atom->mask[i] != exclusion_group_bit) {
        ptype_ID[count_local] = i;
        count_local++;
      }
    }
  } else {
    double dx, dy, dz, distance_check;
    for (int i = 0; i < nlocal; i++) {
      dx = fabs(atom->x[i][0] - target[0]);
      dx -= static_cast<int>(1.0 * dx / (xhi - xlo) + 0.5) * (xhi - xlo);
      dy = fabs(atom->x[i][1] - target[1]);
      dy -= static_cast<int>(1.0 * dy / (yhi - ylo) + 0.5) * (yhi - ylo);
      dz = fabs(atom->x[i][2] - target[2]);
      dz -= static_cast<int>(1.0 * dz / (zhi - zlo) + 0.5) * (zhi - zlo);
      distance_check = dx * dx + dy * dy + dz * dz;
      if ((distance_check < rd * rd) && atom->type[i] == ptype &&
          fabs(atom->q[i] - charge) < SMALL && atom->mask[i] != exclusion_group_bit) {
        ptype_ID[count_local] = i;
        count_local++;
      }
    }
  }
  count_global = count_local;
  count_before = count_local;
  MPI_Allreduce(&count_local, &count_global, 1, MPI_INT, MPI_SUM, world);
  MPI_Scan(&count_local, &count_before, 1, MPI_INT, MPI_SUM, world);
  count_before -= count_local;

  npart_xrd = count_global; // save the number of particles, for use in MC acceptance ratio
  if (count_global > 0) {
    const int ID_global = floor(random_equal->uniform() * count_global);
    if ((ID_global >= count_before) && (ID_global < (count_before + count_local))) {
      const int ID_local = ID_global - count_before;
      m = ptype_ID[ID_local]; // local ID of the chosen particle
      return m;
    }
  }
  return -1;
}

double FixChargeRegulation::energy_full() {
  if (triclinic) domain->x2lamda(atom->nlocal);
  domain->pbc();
  comm->exchange();
  atom->nghost = 0;
  comm->borders();
  if (triclinic) domain->lamda2x(atom->nlocal + atom->nghost);
  if (modify->n_pre_neighbor) modify->pre_neighbor();
  neighbor->build(1);
  int eflag = 1;
  int vflag = 0;
  if (overlap_flag) {
    int overlaptestall;
    int overlaptest = 0;
    double delx, dely, delz, rsq;
    double **x = atom->x;
    int nall = atom->nlocal + atom->nghost;
    for (int i = 0; i < atom->nlocal; i++) {
      for (int j = i + 1; j < nall; j++) {
        delx = x[i][0] - x[j][0];
        dely = x[i][1] - x[j][1];
        delz = x[i][2] - x[j][2];
        rsq = delx * delx + dely * dely + delz * delz;
        if (rsq < overlap_cutoffsq) {
          overlaptest = 1;
          break;
        }
      }
      if (overlaptest) break;
    }
    overlaptestall = overlaptest;
    MPI_Allreduce(&overlaptest, &overlaptestall, 1,
                  MPI_INT, MPI_MAX, world);
    if (overlaptestall) return MAXENERGYSIGNAL;
  }
  size_t nbytes = sizeof(double) * (atom->nlocal + atom->nghost);
  if (nbytes) memset(&atom->f[0][0], 0, 3 * nbytes);

  if (modify->n_pre_force) modify->pre_force(vflag);

  if (force->pair) force->pair->compute(eflag, vflag);

  if (atom->molecular) {
    if (force->bond) force->bond->compute(eflag, vflag);
    if (force->angle) force->angle->compute(eflag, vflag);
    if (force->dihedral) force->dihedral->compute(eflag, vflag);
    if (force->improper) force->improper->compute(eflag, vflag);
  }

  if (force->kspace) force->kspace->compute(eflag, vflag);

  if (modify->n_pre_reverse) modify->pre_reverse(eflag,vflag);
  if (modify->n_post_force_any) modify->post_force(vflag);

  update->eflag_global = update->ntimestep;
  double total_energy = c_pe->compute_scalar();
  return total_energy;
}

int FixChargeRegulation::particle_number_xrd(int ptype, double charge, double rd, double *target) {

  int count = 0;
  if (rd < SMALL) {
    for (int i = 0; i < atom->nlocal; i++) {
      if (atom->type[i] == ptype && fabs(atom->q[i] - charge) < SMALL && atom->mask[i] != exclusion_group_bit)
        count++;
    }
  } else {
    double dx, dy, dz, distance_check;
    for (int i = 0; i < atom->nlocal; i++) {
      dx = fabs(atom->x[i][0] - target[0]);
      dx -= static_cast<int>(1.0 * dx / (xhi - xlo) + 0.5) * (xhi - xlo);
      dy = fabs(atom->x[i][1] - target[1]);
      dy -= static_cast<int>(1.0 * dy / (yhi - ylo) + 0.5) * (yhi - ylo);
      dz = fabs(atom->x[i][2] - target[2]);
      dz -= static_cast<int>(1.0 * dz / (zhi - zlo) + 0.5) * (zhi - zlo);
      distance_check = dx * dx + dy * dy + dz * dz;
      if ((distance_check < rd * rd) && atom->type[i] == ptype &&
          fabs(atom->q[i] - charge) < SMALL && atom->mask[i] != exclusion_group_bit) {
        count++;
      }
    }
  }
  int count_sum = count;
  MPI_Allreduce(&count, &count_sum, 1, MPI_INT, MPI_SUM, world);
  return count_sum;
}

int FixChargeRegulation::particle_number(int ptype, double charge) {

  int count = 0;
  for (int i = 0; i < atom->nlocal; i++) {
    if (atom->type[i] == ptype && fabs(atom->q[i] - charge) < SMALL && atom->mask[i] != exclusion_group_bit)
      count = count + 1;
  }
  int count_sum = count;
  MPI_Allreduce(&count, &count_sum, 1, MPI_INT, MPI_SUM, world);
  return count_sum;
}

double FixChargeRegulation::compute_vector(int n) {
  if (n == 0) {
    return nacid_attempts + nbase_attempts + nsalt_attempts;
  } else if (n == 1) {
    return nacid_successes + nbase_successes + nsalt_successes;
  } else if (n == 2) {
    return particle_number(acid_type, 0);
  } else if (n == 3) {
    return particle_number(acid_type, -1);
  } else if (n == 4) {
    return particle_number(base_type, 0);
  } else if (n == 5) {
    return particle_number(base_type, 1);
  } else if (n == 6) {
    return particle_number(cation_type, salt_charge[0]);
  } else if (n == 7) {
    return particle_number(anion_type, salt_charge[1]);
  }
  return 0.0;
}


/* ----------------------------------------------------------------------
   pack entire state of Fix into one write
------------------------------------------------------------------------- */

void FixChargeRegulation::write_restart(FILE *fp)
{
  int n = 0;
  double list[10];
  list[n++] = random_equal->state();
  list[n++] = random_unequal->state();
  list[n++] = nacid_attempts;
  list[n++] = nacid_successes;
  list[n++] = nbase_attempts;
  list[n++] = nbase_successes;
  list[n++] = nsalt_attempts;
  list[n++] = nsalt_successes;
  list[n++] = ubuf(next_reneighbor).d;
  list[n++] = ubuf(update->ntimestep).d;

  if (comm->me == 0) {
    int size = (int) sizeof(list);
    fwrite(&size,sizeof(int),1,fp);
    fwrite(list,sizeof(list),1,fp);
  }
}

/* ----------------------------------------------------------------------
   use state info from restart file to restart the Fix
------------------------------------------------------------------------- */

void FixChargeRegulation::restart(char *buf)
{
  int n = 0;
  auto list = (double *) buf;

  seed = static_cast<int> (list[n++]);
  random_equal->reset(seed);

  seed = static_cast<int> (list[n++]);
  random_unequal->reset(seed);

  nacid_attempts  = list[n++];
  nacid_successes = list[n++];
  nbase_attempts  = list[n++];
  nbase_successes = list[n++];
  nsalt_attempts  = list[n++];
  nsalt_successes = list[n++];

  next_reneighbor = (bigint) ubuf(list[n++]).i;
  bigint ntimestep_restart = (bigint) ubuf(list[n++]).i;
  if (ntimestep_restart != update->ntimestep)
    error->all(FLERR,"Must not reset timestep when restarting fix gcmc");
}

void FixChargeRegulation::setThermoTemperaturePointer() {
  int ifix = -1;
  ifix = modify->find_fix(idftemp);
  if (ifix == -1) {
    error->all(FLERR, "fix charge/regulation regulation could not find "
               "a temperature fix id provided by tempfixid\n");
  }
  Fix *temperature_fix = modify->fix[ifix];
  int dim;
  target_temperature_tcp = (double *) temperature_fix->extract("t_target", dim);

}

void FixChargeRegulation::assign_tags() {
  // Assign tags to ions with zero tags
  if (atom->tag_enable) {
    tagint *tag = atom->tag;
    tagint maxtag_all = 0;
    tagint maxtag = 0;
    for (int i = 0; i < atom->nlocal; i++) maxtag = MAX(maxtag, tag[i]);
    maxtag_all = maxtag;
    MPI_Allreduce(&maxtag, &maxtag_all, 1, MPI_LMP_TAGINT, MPI_MAX, world);
    if (maxtag_all >= MAXTAGINT) error->all(FLERR, "New atom IDs exceed maximum allowed ID");

    tagint notag = 0;
    tagint notag_all;
    for (int i = 0; i < atom->nlocal; i++)
      if (tag[i] == 0 && (atom->type[i] == cation_type || atom->type[i] == anion_type))notag++;
    notag_all = notag;
    MPI_Allreduce(&notag, &notag_all, 1, MPI_LMP_TAGINT, MPI_SUM, world);
    if (notag_all >= MAXTAGINT)
      error->all(FLERR, "New atom IDs exceed maximum allowed ID");

    tagint notag_sum = notag;
    MPI_Scan(&notag, &notag_sum, 1, MPI_LMP_TAGINT, MPI_SUM, world);
    // itag = 1st new tag that my untagged atoms should use

    tagint itag = maxtag_all + notag_sum - notag + 1;
    for (int i = 0; i < atom->nlocal; i++) {
      if (tag[i] == 0 && (atom->type[i] == cation_type || atom->type[i] == anion_type)) {
        tag[i] = itag++;
      }
    }
    if (atom->map_style) atom->map_init();
    atom->nghost = 0;
    comm->borders();
  }
}

/* ----------------------------------------------------------------------
   parse input options
------------------------------------------------------------------------- */

void FixChargeRegulation::options(int narg, char **arg) {
  if (narg < 0) error->all(FLERR, "Illegal fix charge/regulation command");

  // defaults

  pH = 7.0;
  pI_plus = 5;
  pI_minus = 5;
  acid_type = -1;
  base_type = -1;
  pKa = 100;
  pKb = 100;
  pKs = 14.0;
  nevery = 100;
  nmc = 100;
  pmcmoves[0] = pmcmoves[1] = pmcmoves[2] = THIRD;
  llength_unit_in_nm= 0.71; // Default set to Bjerrum length in water at 20 degrees C [nm]

  reservoir_temperature = 1.0;
  reaction_distance = 0;
  seed = 0;
  target_temperature_tcp = &reservoir_temperature;
  add_tags_flag = false;
  only_salt_flag = false;
  salt_charge[0] = 1; // cation charge
  salt_charge[1] = -1; // anion charge

  exclusion_group = 0;
  exclusion_group_bit = 0;
  ngroups = 0;
  int ngroupsmax = 0;
  groupstrings = nullptr;

  int iarg = 0;
  while (iarg < narg) {

    if (strcmp(arg[iarg], "lunit_nm") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal fix charge/regulation command");
      llength_unit_in_nm = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "acid_type") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal fix charge/regulation command");
      acid_type = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "base_type") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal fix charge/regulation command");
      base_type = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "pH") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal fix charge/regulation command");
      if (strstr(arg[iarg + 1],"v_") == arg[iarg + 1]) {
        pHstr = utils::strdup(&arg[iarg + 1][2]);
      } else {
        pH = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
        pHstyle = CONSTANT;
      }
      iarg += 2;
    } else if (strcmp(arg[iarg], "pIp") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal fix charge/regulation command");
      pI_plus = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "pIm") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal fix charge/regulation command");
      pI_minus = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "pKa") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal fix charge/regulation command");
      pKa = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "pKb") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal fix charge/regulation command");
      pKb = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;

    } else if (strcmp(arg[iarg], "temp") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal fix charge/regulation command");
      reservoir_temperature = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "pKs") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal fix charge/regulation command");
      pKs = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "tempfixid") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal fix charge/regulation command");
      idftemp = utils::strdup(arg[iarg+1]);
      setThermoTemperaturePointer();
      iarg += 2;
    } else if (strcmp(arg[iarg], "rxd") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal fix charge/regulation command");
      reaction_distance = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      if ((reaction_distance > domain->prd_half[0]) ||
          (reaction_distance > domain->prd_half[1]) ||
          (reaction_distance > domain->prd_half[2])) {
        if (comm->me == 0)
          error->warning(FLERR,"reaction distance (rxd) is larger than half "
                         "the box dimension, resetting default: xrd = 0.");
        reaction_distance = 0;
      }
      iarg += 2;
    } else if (strcmp(arg[iarg], "nevery") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal fix charge/regulation command");
      nevery = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "nmc") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal fix charge/regulation command");
      nmc = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "pmcmoves") == 0) {
      if (iarg + 4 > narg) error->all(FLERR, "Illegal fix charge/regulation command");
      pmcmoves[0] = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      pmcmoves[1] = utils::numeric(FLERR, arg[iarg + 2], false, lmp);
      pmcmoves[2] = utils::numeric(FLERR, arg[iarg + 3], false, lmp);
      iarg += 4;
    } else if (strcmp(arg[iarg], "seed") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal fix charge/regulation command");
      seed = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "tag") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal fix charge/regulation command");
      add_tags_flag = utils::logical(FLERR,arg[iarg+1],false,lmp) == 1;
      iarg += 2;
    } else if (strcmp(arg[iarg], "onlysalt") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal fix charge/regulation command");
      only_salt_flag = utils::logical(FLERR,arg[iarg+1],false,lmp) == 1;
      iarg += 2;
      if (only_salt_flag) {
        // need to specify salt charge
        if (iarg + 2 > narg) error->all(FLERR, "Illegal fix charge/regulation command");
        salt_charge[0] = utils::inumeric(FLERR, arg[iarg], false, lmp);
        salt_charge[1] = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
        iarg += 2;
      }
    } else if (strcmp(arg[iarg], "group") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal fix fix charge/regulation command");
      if (ngroups >= ngroupsmax) {
        ngroupsmax = ngroups + 1;
        groupstrings = (char **)
          memory->srealloc(groupstrings, ngroupsmax * sizeof(char *), "fix_charge_regulation:groupstrings");
      }
      groupstrings[ngroups] = utils::strdup(arg[iarg+1]);
      ngroups++;
      iarg += 2;
    } else error->all(FLERR, "Illegal fix charge/regulation command");
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double FixChargeRegulation::memory_usage() {
  return (double)cr_nmax * sizeof(int);
}

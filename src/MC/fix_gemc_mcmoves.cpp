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

/* ----------------------------------------------------------------------
   Contributing author: Andrew Hong, Aidan Thompson (SNL)
------------------------------------------------------------------------- */

#include "fix_gemc.h"

#include "atom.h"
#include "atom_vec.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "group.h"
#include "kspace.h"
#include "memory.h"
#include "modify.h"
#include "neighbor.h"
#include "pair.h"
#include "random_park.h"
#include "universe.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;

// this must be lower than MAXENERGYSIGNAL
// by a large amount, so that it is still
// less than total energy when negative
// energy contributions are added to MAXENERGYSIGNAL

static constexpr double MAXENERGYTEST = 1.0e50;

/* ----------------------------------------------------------------------
  Shrink/expand boxes (always requires full energy)
------------------------------------------------------------------------- */
void FixGEMC::attempt_volume_change_full()
{
  double dlogvolratio;

  nvolume_attempts++;

  // sample change in logvolratio
  // logvolratio = v_self/v_other
  // v1 = v_total/(1+exp(-logvolratio))
  // v2 = v_total/(1+exp(logvolratio))
  // - equal and opposite on both replicas
  // - never goes out of bounds, better sampling efficiency
  // - no communication required

  dlogvolratio = max_dlogvolratio*(2.0 * random_universe->uniform() - 1.0);

  // fvolume = vnew/vold

  double fvolume;
  if (myworld == 0)
    fvolume = (1.0+exp(-logvolratio))/(1.0+exp(-(logvolratio+dlogvolratio)));
  else
    fvolume = (1.0+exp(logvolratio))/(1.0+exp((logvolratio+dlogvolratio)));

  //  printf("fvolume %d %d %g %g\n",myworld, me, fvolume, logvolratio);

  double scale_length = pow(fvolume, 1.0/domain->dimension);

  // convert to lamda coords so they get scaled

  domain->x2lamda(atom->nlocal);
  for (auto &ifix : rfix) ifix->deform(0);

  // shrink box toward lower corner
  // lower box coordinates always same

  xhi_tmp = xlo + (xhi-xlo)*scale_length;
  yhi_tmp = ylo + (yhi-ylo)*scale_length;
  zhi_tmp = zlo + (zhi-zlo)*scale_length;

  // set temporarily

  domain->boxhi[0] = xhi_tmp;
  domain->boxhi[1] = yhi_tmp;
  domain->boxhi[2] = zhi_tmp;

  // reset box and subbox dimensions

  domain->set_global_box();
  domain->set_local_box();    // reassigns sub domains
  comm->setup();
  if (neighbor->style) neighbor->setup_bins();

  // positions are scaled now

  domain->lamda2x(atom->nlocal);
  for (auto &ifix : rfix) ifix->deform(1);

  // remap call

  domain->remap_all();

  // Frenkel & Smit, 3rd Ed. (2023), p. 221, Eq. (6.6.10)
  // prob = ((V, self, new) / (V, self, old))^(N+1) *
  //        ((V, other, new) / (V, other, old))^(N+1) * exp(-beta*dU)

  // change in potential due to volume change

  double dU_volume = (atom->natoms+1) * force->boltz * box_temp * log(fvolume);

  // current system energy

  double energy_before = energy_stored;

  // (possible) future system energy

  double energy_after = energy_full();

  // sum change in full energy across both boxes

  // both boxes must reach the same decision, so the energy change and
  // the overflow flag of both boxes are combined

  double dU = 0.0;
  if (me == 0) {
    double mine[2] = {energy_after - energy_before - dU_volume,
                      (energy_after < MAXENERGYTEST) ? 0.0 : 1.0};
    double other[2];
    MPI_Sendrecv(mine, 2, MPI_DOUBLE, 1 - myworld, 0,
                 other, 2, MPI_DOUBLE, 1 - myworld, 0,
                 comm_replica, MPI_STATUS_IGNORE);
    dU = mine[0] + other[0];
    if ((mine[1] != 0.0) || (other[1] != 0.0) || std::isnan(dU)) dU = INFINITY;
  }

  // bcast potential change to rest of my world

  MPI_Bcast(&dU, 1, MPI_DOUBLE, 0, world);

  // evaluate probability

  double prob = MIN(exp(-beta * dU), 1.0);
  double rf = random_universe->uniform();

  // volume change rejected -> revert atom positions

  if (prob < rf) {

    domain->x2lamda(atom->nlocal);
    for (auto &ifix : rfix) ifix->deform(0);

    domain->boxhi[0] = xhi;
    domain->boxhi[1] = yhi;
    domain->boxhi[2] = zhi;

    // reset box and subbox dimensions

    domain->set_global_box();
    domain->set_local_box();    // reassigns sub domains
    comm->setup();
    if (neighbor->style) neighbor->setup_bins();

    domain->lamda2x(atom->nlocal);
    for (auto &ifix : rfix) ifix->deform(1);

    // remap call (may lose atoms if no remap)

    domain->remap_all();

    // build neighbor list

    neighbor->build(1);

    // accept volume change

  } else {
    nvolume_successes += 1.0;
    logvolratio += dlogvolratio;

    // store new energy

    energy_stored = energy_after;

    // reacquire upper domain bounds

    xhi = domain->boxhi[0];
    yhi = domain->boxhi[1];
    zhi = domain->boxhi[2];

    // reacquire subdomain bounds

    if (triclinic_flag) {
      sublo = domain->sublo_lamda;
      subhi = domain->subhi_lamda;
    } else {
      sublo = domain->sublo;
      subhi = domain->subhi;
    }
  }
}

/* ----------------------------------------------------------------------
  Attempt atom exchange
------------------------------------------------------------------------- */

void FixGEMC::attempt_atomic_exchange_full()
{
  nexchange_attempts++;

  // Choose sender and receiver

  int sender;

  double drand = random_universe->uniform();
  if (drand > 0.5)
    if (myworld == 0)
      sender = 1;
    else
      sender = 0;
  else
    if (myworld == 0)
      sender = 0;
    else
      sender = 1;

  // number of group atoms in my box before the exchange

  int nold = natom_total;

  // atom to delete/insert

  int iatom = -1;
  int tmp_mask;
  double q_tmp;

  // pick atom to send

  if (sender) {

    // pick one atom randomly from all atoms in box
    // only one proc will actually delete atom

    iatom = pick_random_gas_atom();
    if (iatom >= 0) {

      // temporarily set mask to exclusion for full energy later

      tmp_mask = atom->mask[iatom];
      atom->mask[iatom] = exclusion_group_bit;

      // temporarily zero out charge for kspace later)

      if (q_flag) {
        q_tmp = atom->q[iatom];
        atom->q[iatom] = 0.0;
      }
    }
  }

  // pick random proc to place atom in

  int proc_flag = 0;
  tagint newtag = 0;
  if (!sender) {

    // sample random point in box

    double lamda[3], coord[3];
    if (me == 0) {
      if (triclinic_flag) {
        lamda[0] = random_proc->uniform();
        lamda[1] = random_proc->uniform();
        lamda[2] = random_proc->uniform();

        // wasteful, but necessary

        if (lamda[0] == 1.0) lamda[0] = 0.0;
        if (lamda[1] == 1.0) lamda[1] = 0.0;
        if (lamda[2] == 1.0) lamda[2] = 0.0;

        domain->lamda2x(lamda, coord);
      } else {
        coord[0] = xlo + random_proc->uniform() * (xhi - xlo);
        coord[1] = ylo + random_proc->uniform() * (yhi - ylo);
        coord[2] = zlo + random_proc->uniform() * (zhi - zlo);
      }
    }

    // find proc that contains coordinate

    MPI_Bcast(&coord, 3, MPI_DOUBLE, 0, world);
    if (triclinic_flag) {
      domain->x2lamda(coord, lamda);
      if (lamda[0] >= sublo[0] && lamda[0] < subhi[0] && lamda[1] >= sublo[1] &&
          lamda[1] < subhi[1] && lamda[2] >= sublo[2] && lamda[2] < subhi[2])
        proc_flag = 1;
    } else {
      domain->remap(coord);
      if (!domain->inside(coord)) error->universe_one(FLERR, "Fix gemc put atom outside box");
      if (coord[0] >= sublo[0] && coord[0] < subhi[0] && coord[1] >= sublo[1] &&
          coord[1] < subhi[1] && coord[2] >= sublo[2] && coord[2] < subhi[2])
        proc_flag = 1;
    }

    if (proc_flag) {

      // this treatment of exchange type needs to be generalized

      int ngemc_type = 1;
      atom->avec->create_atom(ngemc_type,coord);
      int m = atom->nlocal - 1;

      // add new atom to group all and this fix's group

      int groupbitall = 1 | groupbit;
      atom->mask[m] = groupbitall;

      atom->v[m][0] = random_proc->gaussian()*sigma;
      atom->v[m][1] = random_proc->gaussian()*sigma;
      atom->v[m][2] = random_proc->gaussian()*sigma;
      modify->create_attribute(m);
    }

    atom->natoms++;

    // the new atom has overwritten a ghost atom, ghost atoms are rebuilt later

    atom->nghost = 0;
    if (atom->tag_enable) {
      atom->tag_extend();
      tagint mytag = proc_flag ? atom->tag[atom->nlocal - 1] : 0;
      MPI_Allreduce(&mytag, &newtag, 1, MPI_LMP_TAGINT, MPI_MAX, world);
      if (atom->map_style != Atom::MAP_NONE) atom->map_init();
    }
  }

  // the tail correction must not count the atom to be deleted,
  // so its type is negated while the tail correction is updated

  if (force->pair->tail_flag) {
    if (iatom >= 0) atom->type[iatom] = -atom->type[iatom];
    force->pair->reinit();
    if (iatom >= 0) atom->type[iatom] = -atom->type[iatom];
  }

  update_gas_atoms_list();

  // evaluate probability for exchange

  double energy_before = energy_stored;
  double energy_after = energy_full();

  int success;
  double prob;

  if (me == 0) {

    // Frenkel & Smit, 3rd Ed. (2023), p. 221, Eq. (6.6.11)
    // prob = (V/N, receiver, new) / (V/N, sender, old) *  exp(-beta*dU)
    // natom_total = (N, sender, old)
    // natom_total = (N, receiver, new)

    double idU, jdU, all_dU;
    idU = energy_after - energy_before;

    double volume = (xhi - xlo) * (yhi - ylo) * (zhi - zlo);
    double logVN;

    // the donor box uses its number of atoms before the removal,
    // the receiving box its number of atoms after the insertion

    int overflow = (energy_after < MAXENERGYTEST) ? 0 : 1;
    if (sender) {
      if (nold > 0) {
        logVN = -log(volume / nold);
      } else {
        logVN = 0.0;
        overflow = 1;
      }
    } else {
      logVN = log(volume / natom_total);
    }
    idU += -box_temp * force->boltz * logVN;
    double mine[2] = {idU, (double) overflow};
    double other[2];
    MPI_Sendrecv(mine, 2, MPI_DOUBLE, 1 - myworld, 0,
                 other, 2, MPI_DOUBLE, 1 - myworld, 0,
                 comm_replica, MPI_STATUS_IGNORE);
    jdU = other[0];
    all_dU = idU + jdU;
    if ((overflow != 0) || (other[1] != 0.0) || std::isnan(all_dU))
      prob = 0.0;
    else
      prob = MIN(exp(-beta * all_dU), 1.0);

  }

  // both boxes draw the same random number and reach the same decision

  MPI_Bcast(&prob, 1, MPI_DOUBLE, 0, world);
  if (prob > random_universe->uniform())
    success = 1;
  else
    success = 0;

  // handle deletion/insertions or revert

  if (sender) {
    if (success) {
      nexchange_successes += 1.0;
      if (iatom >= 0) {
        atom->avec->copy(atom->nlocal - 1, iatom, 1);
        atom->nlocal--;
      }
      atom->natoms--;
      if (atom->map_style != Atom::MAP_NONE) atom->map_init();
      if (force->pair->tail_flag) force->pair->reinit();
      energy_stored = energy_after;

    } else {
      if (iatom >= 0) {
        atom->mask[iatom] = tmp_mask;
        if (q_flag) atom->q[iatom] = q_tmp;
      }
      if (force->kspace) force->kspace->qsum_qsq();
      if (force->pair->tail_flag) force->pair->reinit();
    }
  } else {
    if (success) {
      nexchange_successes += 1.0;
      energy_stored = energy_after;
    } else {
      atom->natoms--;
      if (newtag) {
        for (int k = 0; k < atom->nlocal; k++) {
          if (atom->tag[k] == newtag) {
            atom->avec->copy(atom->nlocal - 1, k, 1);
            atom->nlocal--;
            break;
          }
        }
        if (atom->map_style != Atom::MAP_NONE) atom->map_init();
      } else if (proc_flag) atom->nlocal--;
      if (force->kspace) force->kspace->qsum_qsq();
      if (force->pair->tail_flag) force->pair->reinit();
      energy_stored = energy_before;
    }
  }

  // update counts

  update_gas_atoms_list();
}

/* ----------------------------------------------------------------------
   copied directly from fix_gcmc.cpp
------------------------------------------------------------------------- */

void FixGEMC::attempt_atomic_translation_full()
{
  ntranslation_attempts += 1.0;

  if (natom_total == 0) return;

  double energy_before = energy_stored;

  int i = pick_random_gas_atom();

  double **x = atom->x;
  double xtmp[3];

  xtmp[0] = xtmp[1] = xtmp[2] = 0.0;

  tagint tmptag = -1;
  imageint tmpimage = 0;

  if (i >= 0) {

    double rsq = 1.1;
    double rx, ry, rz;
    rx = ry = rz = 0.0;
    double coord[3];
    while (rsq > 1.0) {
      rx = 2 * random_proc->uniform() - 1.0;
      ry = 2 * random_proc->uniform() - 1.0;
      rz = 2 * random_proc->uniform() - 1.0;
      rsq = rx * rx + ry * ry + rz * rz;
    }
    coord[0] = x[i][0] + displace * rx;
    coord[1] = x[i][1] + displace * ry;
    coord[2] = x[i][2] + displace * rz;

    if (!domain->inside_nonperiodic(coord)) error->one(FLERR, "Fix gcmc put atom outside box");
    xtmp[0] = x[i][0];
    xtmp[1] = x[i][1];
    xtmp[2] = x[i][2];
    x[i][0] = coord[0];
    x[i][1] = coord[1];
    x[i][2] = coord[2];

    tmptag = atom->tag[i];
    tmpimage = atom->image[i];
  }

  double energy_after = energy_full();

  if (energy_after < MAXENERGYTEST &&
      random_world->uniform() < exp(beta * (energy_before - energy_after))) {

    energy_stored = energy_after;
    ntranslation_successes += 1.0;
  } else {

    tagint tmptag_all;
    MPI_Allreduce(&tmptag, &tmptag_all, 1, MPI_LMP_TAGINT, MPI_MAX, world);

    double xtmp_all[3];
    MPI_Allreduce(&xtmp, &xtmp_all, 3, MPI_DOUBLE, MPI_SUM, world);
    imageint tmpimage_all;
    MPI_Allreduce(&tmpimage, &tmpimage_all, 1, MPI_LMP_IMAGEINT, MPI_SUM, world);

    // energy_full() may have reallocated the per-atom arrays,
    // and may have wrapped the atom around a periodic boundary

    x = atom->x;
    for (int i = 0; i < atom->nlocal; i++) {
      if (tmptag_all == atom->tag[i]) {
        x[i][0] = xtmp_all[0];
        x[i][1] = xtmp_all[1];
        x[i][2] = xtmp_all[2];
        atom->image[i] = tmpimage_all;
      }
    }

    energy_stored = energy_before;

    // this remapping is necessary, but is not clear why
    // also, not clear why it is *not* necessary in fix gcmc

    if (triclinic_flag) domain->x2lamda(atom->nlocal);
    domain->pbc();
    comm->exchange();
    atom->nghost = 0;
    comm->borders();
    if (triclinic_flag) domain->lamda2x(atom->nlocal + atom->nghost);
    if (modify->n_pre_neighbor) modify->pre_neighbor();
    neighbor->build(1);
  }
  update_gas_atoms_list();
}

/* ----------------------------------------------------------------------
------------------------------------------------------------------------- */

int FixGEMC::pick_random_gas_atom()
{
  int i = -1;
  int iwhichglobal = static_cast<int>(natom_total * random_world->uniform());
  if ((iwhichglobal >= natom_lower) && (iwhichglobal < natom_lower + natom_local)) {
    i = local_gas_list[iwhichglobal - natom_lower];
  }

  return i;
}

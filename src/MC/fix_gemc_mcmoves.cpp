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
#include "fix.h"
#include "force.h"
#include "math_const.h"
#include "math_extra.h"
#include "modify.h"
#include "molecule.h"
#include "neighbor.h"
#include "pair.h"
#include "random_park.h"

#include <algorithm>
#include <cmath>

using namespace LAMMPS_NS;
using MathConst::MY_2PI;

// trial energies above this value (or not a number) are always rejected

static constexpr double MAXENERGYTEST = 1.0e50;

// values per atom in gather_molecule() and in the exchanged molecule data

static constexpr int NMOLDATA = 10;
static constexpr int NMOLINFO = 8;

/* ----------------------------------------------------------------------
  Shrink one box and expand the other one by the same volume
------------------------------------------------------------------------- */

void FixGEMC::attempt_volume_change_full()
{
  nvolume_attempts += 1.0;

  // all owned atoms must be in the atom arrays

  flush_pending();

  // random walk in logvolratio = log(V1/V2), identical in both boxes
  // V1 = Vtotal/(1+exp(-logvolratio)), V2 = Vtotal/(1+exp(logvolratio))
  // box volumes are always positive and Vtotal is conserved

  double dlogvolratio = max_dlogvolratio * (2.0 * random_universe->uniform() - 1.0);

  // fvolume = Vnew/Vold of my box

  double fvolume;
  if (myworld == 0)
    fvolume = (1.0 + exp(-logvolratio)) / (1.0 + exp(-(logvolratio + dlogvolratio)));
  else
    fvolume = (1.0 + exp(logvolratio)) / (1.0 + exp(logvolratio + dlogvolratio));

  double s = cbrt(fvolume);

  // scale the box uniformly about its lower corner (keeping its shape)
  // atoms are scaled individually, molecules with their center of mass

  // moves that make a box narrower than twice the size of a molecule are
  // rejected in both boxes

  bigint nunits = scale_positions(s, 1);
  if (nunits < 0) return;
  set_box(xlo + (xhi - xlo) * s, ylo + (yhi - ylo) * s, zlo + (zhi - zlo) * s, boxxy * s, boxxz * s,
          boxyz * s);

  // Frenkel & Smit, 3rd Ed. (2023), Eq. (6.6.10) for a random walk in log(V1/V2):
  // acc = (V1new/V1old)^(N1+1) * (V2new/V2old)^(N2+1) * exp(-beta*(dU1+dU2))
  // N = number of independent units (molecules or atoms) in a box
  // each box contributes dU - (N+1)*kT*log(Vnew/Vold)

  double energy_after = energy_full();
  double dU = energy_after - energy_stored - (nunits + 1) * force->boltz * box_temp * log(fvolume);
  int overflow = !(energy_after < MAXENERGYTEST);

  if (accept_both(dU, overflow)) {
    nvolume_successes += 1.0;
    logvolratio += dlogvolratio;
    energy_stored = energy_after;
    xhi = domain->boxhi[0];
    yhi = domain->boxhi[1];
    zhi = domain->boxhi[2];
    boxxy = domain->xy;
    boxxz = domain->xz;
    boxyz = domain->yz;

  } else {

    // rejected: scale positions back while the trial box is still set, restore box

    scale_positions(1.0 / s, 0);
    set_box(xhi, yhi, zhi, boxxy, boxxz, boxyz);
    ghosts_stale = 1;
  }

  // atoms may have migrated between ranks, so the local list must be regenerated

  update_gas_atoms_list();
}

/* ----------------------------------------------------------------------
  Move a randomly chosen atom from one box to a random position in the other box
------------------------------------------------------------------------- */

void FixGEMC::attempt_atomic_exchange_full()
{
  nexchange_attempts += 1.0;

  // choose donor box with equal probability, identical in both boxes

  int donor = (random_universe->uniform() < 0.5) ? 0 : 1;
  int sender = (myworld == donor) ? 1 : 0;

  // each box independently decides whether it can use single-atom energies

  int local = use_local();

  double energy_before = energy_stored;
  double volume = box_volume();
  bigint nold = natom_total;
  double denergy = 0.0;

  // donor box: pick an atom
  // with full energy: remove it and keep its complete state for a possible restore
  // with local energy: only compute its energy, remove it after acceptance
  // donor_info = {box is empty, atom type, atom mask, atom charge, velocity}

  double donor_info[7] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  int iremove = -1;
  tagint itag = 0;
  int owner = -1;
  if (sender) {
    if (natom_total == 0) {
      donor_info[0] = 1.0;
    } else {
      iremove = pick_random_gas_atom();
      if (local) {
        itag = (iremove >= 0) ? atom->tag[iremove] : 0;
        int bad = (iremove >= 0) && (subdomain_excess(atom->x[iremove]) > ghost_skin);
        if (need_refresh(bad, itag)) iremove = local_index(itag);
      }

      // info = {type, mask, charge, velocity, owning rank + 1}, only from the owning rank

      double info[7] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
      double eone = 0.0;
      if (iremove >= 0) {
        info[0] = atom->type[iremove];
        info[1] = atom->mask[iremove];
        info[2] = atom->q_flag ? atom->q[iremove] : 0.0;
        info[3] = atom->v[iremove][0];
        info[4] = atom->v[iremove][1];
        info[5] = atom->v[iremove][2];
        info[6] = me + 1;
        if (local) {
          eone = energy_local(iremove, atom->type[iremove], atom->tag[iremove], atom->x[iremove]);
        } else {
          size_t nbuf = atom->avec->maxexchange + 1024;
          for (const auto &ifix : modify->get_fix_list()) nbuf += ifix->maxexchange;
          if (exchange_buf.size() < nbuf) exchange_buf.resize(nbuf);
          atom->avec->pack_exchange(iremove, exchange_buf.data());
          atom->avec->copy(atom->nlocal - 1, iremove, 1);
          atom->nlocal--;
        }
      }
      // exactly one rank owns the atom, all others contribute zeros
      double info_all[7];
      MPI_Allreduce(info, info_all, 7, MPI_DOUBLE, MPI_SUM, world);
      for (int k = 0; k < 6; k++) donor_info[k + 1] = info_all[k];
      owner = static_cast<int>(info_all[6]) - 1;
      if (local) {
        MPI_Allreduce(&eone, &denergy, 1, MPI_DOUBLE, MPI_SUM, world);
        denergy = -denergy;
      } else {
        atom->natoms--;
        changed_atoms();
      }
    }
  }

  // receiver box needs the donor info

  if (me == 0) {
    double other[7];
    MPI_Sendrecv(donor_info, 7, MPI_DOUBLE, 1 - myworld, 0, other, 7, MPI_DOUBLE, 1 - myworld, 0,
                 comm_replica, MPI_STATUS_IGNORE);
    if (!sender)
      for (int k = 0; k < 7; k++) donor_info[k] = other[k];
  }
  MPI_Bcast(donor_info, 7, MPI_DOUBLE, 0, world);

  // donor box is empty: reject without doing anything

  if (donor_info[0] != 0.0) return;

  // receiver box: insert an atom of the same type, group membership, and charge
  // at a random position
  // with local energy: compute its energy in a scratch slot and insert it after acceptance

  int itype = static_cast<int>(donor_info[1]);
  int imask = static_cast<int>(donor_info[2]);
  double iq = donor_info[3];
  double *iv = &donor_info[4];
  double coord[3];
  int proc_flag = 0;
  tagint newtag = 0;
  if (!sender) {
    random_point(coord);
    proc_flag = owns(coord);

    if (local) {

      // {energy of the inserted atom, owning rank + 1}, only from the owning rank

      double eone[2] = {0.0, 0.0};
      if (proc_flag) {
        int ii = nstored;
        grow_stored(ii + 1);
        atom->type[ii] = itype;
        atom->mask[ii] = imask;
        atom->tag[ii] = 0;
        if (atom->q_flag) atom->q[ii] = iq;
        atom->x[ii][0] = coord[0];
        atom->x[ii][1] = coord[1];
        atom->x[ii][2] = coord[2];
        eone[0] = energy_local(ii, itype, 0, coord);
        eone[1] = me + 1;
      }
      double eone_all[2];
      MPI_Allreduce(eone, eone_all, 2, MPI_DOUBLE, MPI_SUM, world);
      denergy = eone_all[0];
      owner = static_cast<int>(eone_all[1]) - 1;
    } else {
      newtag = insert_atom(itype, imask, iq, iv, coord, proc_flag);
      changed_atoms();
    }
  }

  // Frenkel & Smit, 3rd Ed. (2023), Eq. (6.6.11)
  // acc = (Ndonor,old * Vreceiver) / (Vdonor * (Nreceiver,old+1)) * exp(-beta*(dU1+dU2))
  // each box contributes dU - kT*log(its factor)

  double energy_after = local ? energy_before + denergy : energy_full();
  double logfactor;
  if (sender)
    logfactor = log(nold / volume);
  else
    logfactor = log(volume / (nold + 1));
  double dU = energy_after - energy_before - force->boltz * box_temp * logfactor;
  int overflow = !(energy_after < MAXENERGYTEST);

  if (accept_both(dU, overflow)) {
    nexchange_successes += 1.0;
    energy_stored = energy_after;

    // with local energy the move has not been applied yet: update the grid
    // on all ranks, the atom arrays are updated later by flush_pending()

    if (local) {
      if (sender) {
        remove_images(itag);
        if (iremove >= 0) {
          gas_remove(iremove);
          if (iremove < atom->nlocal) {
            removed.push_back(iremove);
          } else {
            for (auto &k : pending)
              if (k == iremove) k = -1;
          }
        }
        natom_total--;
        if (me > owner) natom_lower--;
        atom->natoms--;
      } else {
        if (maxtag_box >= MAXTAGINT)
          error->all(FLERR, "Fix gemc ran out of atom IDs; use reset_atoms id between runs");
        tagint newid = ++maxtag_box;
        insert_pending(newid, itype, imask, iq, coord, iv, proc_flag);
        natom_total++;
        if (me > owner) natom_lower++;
        atom->natoms++;
      }
      pending_changes = 1;
    }

  } else if (!local) {

    // rejected: put the removed atom back or delete the inserted atom

    if (sender) {
      if (iremove >= 0) atom->avec->unpack_exchange(exchange_buf.data());
      atom->natoms++;
    } else {
      int i = local_index(newtag);
      if (i >= 0) {
        atom->avec->copy(atom->nlocal - 1, i, 1);
        atom->nlocal--;
      }
      atom->natoms--;
    }
    changed_atoms();
    energy_stored = energy_before;
    ghosts_stale = 1;
  }

  if (!local) update_gas_atoms_list();
}

/* ----------------------------------------------------------------------
   insert an atom of type itype with group mask imask, charge iq, and velocity iv
   at coord on the owning rank (proc_flag = 1) and assign a new atom ID.
   the velocity of the removed atom is kept, so that time integration
   continues consistently in hybrid MD/MC simulations.
   must be called by all ranks of a box, returns the new atom ID.
------------------------------------------------------------------------- */

tagint FixGEMC::insert_atom(int itype, int imask, double iq, double *iv, double *coord,
                            int proc_flag)
{
  tagint mytag = 0;
  if (proc_flag) {
    atom->avec->create_atom(itype, coord);
    int m = atom->nlocal - 1;
    atom->mask[m] = imask;
    if (atom->q_flag) atom->q[m] = iq;
    atom->v[m][0] = iv[0];
    atom->v[m][1] = iv[1];
    atom->v[m][2] = iv[2];
    modify->create_attribute(m);
    atom->tag_extend();
    mytag = atom->tag[m];
  } else {
    atom->tag_extend();
  }
  tagint newtag;
  MPI_Allreduce(&mytag, &newtag, 1, MPI_LMP_TAGINT, MPI_MAX, world);
  if (newtag == 0) error->all(FLERR, "Fix gemc failed to insert atom");
  atom->natoms++;
  return newtag;
}

/* ----------------------------------------------------------------------
   Displace a randomly chosen atom within a box
------------------------------------------------------------------------- */

void FixGEMC::attempt_atomic_translation_full()
{
  if (natom_total == 0) return;

  ntranslation_attempts += 1.0;

  int local = use_local();
  double energy_before = energy_stored;
  double dmax = local ? MIN(displace, move_limit) : max_translation();

  int i = pick_random_gas_atom();

  // displacement drawn by the rank that owns the atom

  double step[3] = {0.0, 0.0, 0.0};
  if (i >= 0) {
    double rsq = 1.1;
    double rx, ry, rz;
    rx = ry = rz = 0.0;
    while (rsq > 1.0) {
      rx = 2.0 * random_proc->uniform() - 1.0;
      ry = 2.0 * random_proc->uniform() - 1.0;
      rz = 2.0 * random_proc->uniform() - 1.0;
      rsq = rx * rx + ry * ry + rz * rz;
    }
    step[0] = dmax * rx;
    step[1] = dmax * ry;
    step[2] = dmax * rz;
  }

  if (local) {

    // the energies at the old and new position are covered by the stored atoms,
    // if both are close enough to the subdomain of the owning rank, and the atom
    // must stay within reach of the neighboring subdomains. otherwise rebuild the
    // ghost atoms, the atom may then be owned by another rank.

    tagint itag = (i >= 0) ? atom->tag[i] : 0;
    int bad = 0;
    if (i >= 0) {
      double xnew[3], dnew[3];
      MathExtra::add3(atom->x[i], step, xnew);
      MathExtra::add3(&dacc[3 * i], step, dnew);
      bad = (subdomain_excess(atom->x[i]) > ghost_skin) || (subdomain_excess(xnew) > ghost_skin) ||
          (MathExtra::len3(dnew) > move_limit);
    }
    if (need_refresh(bad, itag)) {
      double step_all[3];
      MPI_Allreduce(step, step_all, 3, MPI_DOUBLE, MPI_SUM, world);
      i = local_index(itag);
      for (int d = 0; d < 3; d++) step[d] = (i >= 0) ? step_all[d] : 0.0;
    }

    // buf = {energy change, displacement, new position, type, mask, charge}
    // only the owning rank contributes

    double buf[10] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    if (i >= 0) {
      int itype = atom->type[i];
      double *xi = atom->x[i];
      double xnew[3];
      MathExtra::add3(xi, step, xnew);
      buf[0] = energy_local(i, itype, itag, xnew) - energy_local(i, itype, itag, xi);
      buf[1] = step[0];
      buf[2] = step[1];
      buf[3] = step[2];
      buf[4] = xnew[0];
      buf[5] = xnew[1];
      buf[6] = xnew[2];
      buf[7] = itype;
      buf[8] = atom->mask[i];
      buf[9] = atom->q_flag ? atom->q[i] : 0.0;
    }
    double buf_all[10];
    MPI_Allreduce(buf, buf_all, 10, MPI_DOUBLE, MPI_SUM, world);
    double energy_after = energy_before + buf_all[0];

    if ((energy_after < MAXENERGYTEST) &&
        (random_world->uniform() < exp(beta * (energy_before - energy_after)))) {
      energy_stored = energy_after;
      ntranslation_successes += 1.0;

      // move the atom and all its images on all ranks, and store the images
      // that are now within the ghost cutoff of a rank

      move_images(itag, &buf_all[1]);
      double vzero[3] = {0.0, 0.0, 0.0};
      add_images(itag, static_cast<int>(buf_all[7]), static_cast<int>(buf_all[8]), buf_all[9],
                 &buf_all[4], vzero);
      if (i >= 0) MathExtra::add3(&dacc[3 * i], &buf_all[1], &dacc[3 * i]);
    }
    return;
  }

  double xold[3] = {0.0, 0.0, 0.0};
  imageint imageold = 0;
  tagint tagold = 0;

  if (i >= 0) {
    double **x = atom->x;
    xold[0] = x[i][0];
    xold[1] = x[i][1];
    xold[2] = x[i][2];
    imageold = atom->image[i];
    tagold = atom->tag[i];
    x[i][0] += step[0];
    x[i][1] += step[1];
    x[i][2] += step[2];
  }

  double energy_after = energy_full();

  if ((energy_after < MAXENERGYTEST) &&
      (random_world->uniform() < exp(beta * (energy_before - energy_after)))) {
    energy_stored = energy_after;
    ntranslation_successes += 1.0;

  } else {

    // rejected: restore position and image flags of the atom,
    // which may have moved to a different rank or across a periodic boundary

    tagint tagold_all;
    MPI_Allreduce(&tagold, &tagold_all, 1, MPI_LMP_TAGINT, MPI_MAX, world);
    double xold_all[3];
    MPI_Allreduce(xold, xold_all, 3, MPI_DOUBLE, MPI_SUM, world);
    imageint imageold_all;
    MPI_Allreduce(&imageold, &imageold_all, 1, MPI_LMP_IMAGEINT, MPI_SUM, world);

    int j = local_index(tagold_all);
    if (j >= 0) {
      atom->x[j][0] = xold_all[0];
      atom->x[j][1] = xold_all[1];
      atom->x[j][2] = xold_all[2];
      atom->image[j] = imageold_all;
    }
    energy_stored = energy_before;
    ghosts_stale = 1;
  }
  update_gas_atoms_list();
}

/* ----------------------------------------------------------------------
   return molecule ID of a randomly chosen molecule in my box
   all molecules have the same number of atoms, so picking a random group
   atom picks each molecule with equal probability
------------------------------------------------------------------------- */

tagint FixGEMC::pick_random_molecule()
{
  int i = pick_random_gas_atom();
  tagint molid = (i >= 0) ? atom->molecule[i] : 0;
  tagint molid_all;
  MPI_Allreduce(&molid, &molid_all, 1, MPI_LMP_TAGINT, MPI_MAX, world);
  return molid_all;
}

/* ----------------------------------------------------------------------
   smallest atom ID of molecule molid, on all ranks of the box
   molecules of the fix group have consecutive atom IDs in template order,
   so atom ID - first atom ID is the index of an atom in its molecule
------------------------------------------------------------------------- */

tagint FixGEMC::first_atom(tagint molid)
{
  tagint first = MAXTAGINT;
  for (int i = 0; i < atom->nlocal; i++)
    if (atom->molecule[i] == molid) first = MIN(first, atom->tag[i]);
  tagint first_all;
  MPI_Allreduce(&first, &first_all, 1, MPI_LMP_TAGINT, MPI_MIN, world);
  return first_all;
}

/* ----------------------------------------------------------------------
   gather data of all atoms of molecule molid on all ranks, in template order
   NMOLDATA values per atom: 1 (count), unwrapped x,y,z, charge, mask, mass, velocity
   returns the first atom ID of the molecule
------------------------------------------------------------------------- */

tagint FixGEMC::gather_molecule(tagint molid, std::vector<double> &data)
{
  const int n = natoms_per_molecule;
  tagint first = first_atom(molid);

  // each atom fills its own slot, all other ranks contribute zeros

  std::vector<double> mine(NMOLDATA * n, 0.0);
  int bad = 0;
  double xu[3];
  for (int i = 0; i < atom->nlocal; i++) {
    if (atom->molecule[i] != molid) continue;
    bigint k = atom->tag[i] - first;
    if ((k < 0) || (k >= n)) {
      bad = 1;
      continue;
    }
    double *slot = &mine[NMOLDATA * k];
    domain->unmap(atom->x[i], atom->image[i], xu);
    slot[0] += 1.0;
    slot[1] = xu[0];
    slot[2] = xu[1];
    slot[3] = xu[2];
    slot[4] = atom->q_flag ? atom->q[i] : 0.0;
    slot[5] = atom->mask[i];
    slot[6] = atom->mass[atom->type[i]];
    slot[7] = atom->v[i][0];
    slot[8] = atom->v[i][1];
    slot[9] = atom->v[i][2];
  }
  data.resize(NMOLDATA * n);
  MPI_Allreduce(mine.data(), data.data(), NMOLDATA * n, MPI_DOUBLE, MPI_SUM, world);
  for (int k = 0; k < n; k++)
    if (data[NMOLDATA * k] != 1.0) bad = 1;
  int bad_all;
  MPI_Allreduce(&bad, &bad_all, 1, MPI_INT, MPI_MAX, world);
  if (bad_all)
    error->all(FLERR,
               "Fix gemc found molecule {} whose atoms do not match molecule template {} with "
               "consecutive atom IDs",
               molid, idmol);
  return first;
}

/* ----------------------------------------------------------------------
   Displace a randomly chosen molecule within a box
------------------------------------------------------------------------- */

void FixGEMC::attempt_molecule_translation_full()
{
  // all owned atoms must be in the atom arrays

  flush_pending();

  if (natom_total == 0) return;

  ntranslation_attempts += 1.0;

  double energy_before = energy_stored;
  tagint molid = pick_random_molecule();

  double d[3];
  double rsq = 1.1;
  while (rsq > 1.0) {
    d[0] = 2.0 * random_world->uniform() - 1.0;
    d[1] = 2.0 * random_world->uniform() - 1.0;
    d[2] = 2.0 * random_world->uniform() - 1.0;
    rsq = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
  }
  double dmax = max_translation();
  d[0] *= dmax;
  d[1] *= dmax;
  d[2] *= dmax;

  for (int i = 0; i < atom->nlocal; i++) {
    if (atom->molecule[i] == molid) {
      atom->x[i][0] += d[0];
      atom->x[i][1] += d[1];
      atom->x[i][2] += d[2];
    }
  }

  double energy_after = energy_full();

  if ((energy_after < MAXENERGYTEST) &&
      (random_world->uniform() < exp(beta * (energy_before - energy_after)))) {
    energy_stored = energy_after;
    ntranslation_successes += 1.0;
  } else {

    // rejected: shifting back is consistent with the image flags,
    // even if energy_full() wrapped atoms around periodic boundaries

    double **x = atom->x;
    for (int i = 0; i < atom->nlocal; i++) {
      if (atom->molecule[i] == molid) {
        x[i][0] -= d[0];
        x[i][1] -= d[1];
        x[i][2] -= d[2];
      }
    }
    energy_stored = energy_before;
  }
  update_gas_atoms_list();
}

/* ----------------------------------------------------------------------
   Rotate a randomly chosen molecule about its center of mass
------------------------------------------------------------------------- */

void FixGEMC::attempt_molecule_rotation_full()
{
  // all owned atoms must be in the atom arrays

  flush_pending();

  if (natom_total == 0) return;

  nrotation_attempts += 1.0;

  double energy_before = energy_stored;
  tagint molid = pick_random_molecule();

  // center of mass from unwrapped coordinates

  std::vector<double> data;
  tagint first = gather_molecule(molid, data);
  const int n = natoms_per_molecule;
  double com[3] = {0.0, 0.0, 0.0};
  double mtotal = 0.0;
  for (int k = 0; k < n; k++) {
    double m = data[NMOLDATA * k + 6];
    com[0] += m * data[NMOLDATA * k + 1];
    com[1] += m * data[NMOLDATA * k + 2];
    com[2] += m * data[NMOLDATA * k + 3];
    mtotal += m;
  }
  com[0] /= mtotal;
  com[1] /= mtotal;
  com[2] /= mtotal;

  // random axis and angle in [-maxangle,maxangle], a symmetric proposal

  double r[3], quat[4], rotmat[3][3];
  double rsq = 1.1;
  while ((rsq > 1.0) || (rsq < 1.0e-6)) {
    r[0] = 2.0 * random_world->uniform() - 1.0;
    r[1] = 2.0 * random_world->uniform() - 1.0;
    r[2] = 2.0 * random_world->uniform() - 1.0;
    rsq = MathExtra::dot3(r, r);
  }
  MathExtra::norm3(r);
  double theta = maxangle * (2.0 * random_world->uniform() - 1.0);
  MathExtra::axisangle_to_quat(r, theta, quat);
  MathExtra::quat_to_mat(quat, rotmat);

  // reject rotations that move an atom farther than allowed by the subdomain
  // size. the reverse rotation moves each atom by the same distance, so the
  // rejection is symmetric. all ranks hold all atoms of the molecule.

  double dmove = max_move();
  for (int k = 0; k < n; k++) {
    double dx[3], dxr[3];
    dx[0] = data[NMOLDATA * k + 1] - com[0];
    dx[1] = data[NMOLDATA * k + 2] - com[1];
    dx[2] = data[NMOLDATA * k + 3] - com[2];
    MathExtra::matvec(rotmat, dx, dxr);
    MathExtra::sub3(dxr, dx, dxr);
    if (MathExtra::lensq3(dxr) >= dmove * dmove) return;
  }

  // save old positions and image flags, apply rotation to unwrapped coordinates

  const imageint imagezero =
      ((imageint) IMGMAX << IMG2BITS) | ((imageint) IMGMAX << IMGBITS) | IMGMAX;
  std::vector<double> xsaved(3 * n, 0.0);
  std::vector<imageint> imagesaved(n, 0);
  double **x = atom->x;
  imageint *image = atom->image;
  for (int i = 0; i < atom->nlocal; i++) {
    if (atom->molecule[i] != molid) continue;
    bigint k = atom->tag[i] - first;
    xsaved[3 * k] = x[i][0];
    xsaved[3 * k + 1] = x[i][1];
    xsaved[3 * k + 2] = x[i][2];
    imagesaved[k] = image[i];
    double xu[3], dx[3];
    domain->unmap(x[i], image[i], xu);
    dx[0] = xu[0] - com[0];
    dx[1] = xu[1] - com[1];
    dx[2] = xu[2] - com[2];
    MathExtra::matvec(rotmat, dx, x[i]);
    x[i][0] += com[0];
    x[i][1] += com[1];
    x[i][2] += com[2];
    image[i] = imagezero;
    domain->remap(x[i], image[i]);
  }

  double energy_after = energy_full();

  if ((energy_after < MAXENERGYTEST) &&
      (random_world->uniform() < exp(beta * (energy_before - energy_after)))) {
    energy_stored = energy_after;
    nrotation_successes += 1.0;
  } else {

    // rejected: energy_full() may have moved atoms of the molecule to other
    // ranks, so restore saved positions and image flags by atom ID.
    // each atom was saved by exactly one rank, all others contributed zeros

    std::vector<double> xsaved_all(3 * n);
    std::vector<imageint> imagesaved_all(n);
    MPI_Allreduce(xsaved.data(), xsaved_all.data(), 3 * n, MPI_DOUBLE, MPI_SUM, world);
    MPI_Allreduce(imagesaved.data(), imagesaved_all.data(), n, MPI_LMP_IMAGEINT, MPI_SUM, world);
    x = atom->x;
    image = atom->image;
    for (int k = 0; k < n; k++) {
      int i = local_index(first + k);
      if (i >= 0) {
        x[i][0] = xsaved_all[3 * k];
        x[i][1] = xsaved_all[3 * k + 1];
        x[i][2] = xsaved_all[3 * k + 2];
        image[i] = imagesaved_all[k];
      }
    }
    energy_stored = energy_before;
  }
  update_gas_atoms_list();
}

/* ----------------------------------------------------------------------
  Move a randomly chosen molecule from one box to the other box.
  the molecule keeps its conformation and is inserted at a random position
  with a random orientation, so the move is symmetric for flexible molecules
------------------------------------------------------------------------- */

void FixGEMC::attempt_molecule_exchange_full()
{
  nexchange_attempts += 1.0;

  // all owned atoms must be in the atom arrays

  flush_pending();

  const int n = natoms_per_molecule;

  // choose donor box with equal probability, identical in both boxes

  int donor = (random_universe->uniform() < 0.5) ? 0 : 1;
  int sender = (myworld == donor) ? 1 : 0;

  double energy_before = energy_stored;
  double volume = box_volume();
  bigint nold = natom_total / n;

  // donor box: pick a molecule, record its conformation relative to its
  // center of mass, then remove its atoms and keep their state for a restore
  // molinfo = {box is empty, n x (dx, dy, dz, charge, mask, vx, vy, vz)}

  std::vector<double> molinfo(1 + NMOLINFO * n, 0.0);
  int nremoved = 0;
  if (sender) {
    if (nold == 0) {
      molinfo[0] = 1.0;
    } else {
      tagint molid = pick_random_molecule();
      std::vector<double> data;
      gather_molecule(molid, data);
      double com[3] = {0.0, 0.0, 0.0};
      double mtotal = 0.0;
      for (int k = 0; k < n; k++) {
        double m = data[NMOLDATA * k + 6];
        com[0] += m * data[NMOLDATA * k + 1];
        com[1] += m * data[NMOLDATA * k + 2];
        com[2] += m * data[NMOLDATA * k + 3];
        mtotal += m;
      }
      for (int j = 0; j < 3; j++) com[j] /= mtotal;
      for (int k = 0; k < n; k++) {
        double *info = &molinfo[1 + NMOLINFO * k];
        info[0] = data[NMOLDATA * k + 1] - com[0];
        info[1] = data[NMOLDATA * k + 2] - com[1];
        info[2] = data[NMOLDATA * k + 3] - com[2];
        info[3] = data[NMOLDATA * k + 4];
        info[4] = data[NMOLDATA * k + 5];
        info[5] = data[NMOLDATA * k + 7];
        info[6] = data[NMOLDATA * k + 8];
        info[7] = data[NMOLDATA * k + 9];
      }

      // remove atoms of the molecule, save their complete state

      size_t nbuf = atom->avec->maxexchange + 1024 + 2 * atom->bond_per_atom +
          4 * atom->angle_per_atom + 5 * atom->dihedral_per_atom + 5 * atom->improper_per_atom +
          atom->maxspecial;
      for (const auto &ifix : modify->get_fix_list()) nbuf += ifix->maxexchange;
      size_t used = 0;
      int i = 0;
      while (i < atom->nlocal) {
        if (atom->molecule[i] == molid) {
          if (exchange_buf.size() < used + nbuf) exchange_buf.resize(used + nbuf);
          used += atom->avec->pack_exchange(i, &exchange_buf[used]);
          atom->avec->copy(atom->nlocal - 1, i, 1);
          atom->nlocal--;
          nremoved++;
        } else {
          i++;
        }
      }
      atom->natoms -= n;
      atom->nbonds -= onemol->nbonds;
      atom->nangles -= onemol->nangles;
      atom->ndihedrals -= onemol->ndihedrals;
      atom->nimpropers -= onemol->nimpropers;
      changed_atoms();
    }
  }

  // receiver box needs the donor info

  if (me == 0) {
    std::vector<double> other(1 + NMOLINFO * n);
    MPI_Sendrecv(molinfo.data(), 1 + NMOLINFO * n, MPI_DOUBLE, 1 - myworld, 0, other.data(),
                 1 + NMOLINFO * n, MPI_DOUBLE, 1 - myworld, 0, comm_replica, MPI_STATUS_IGNORE);
    if (!sender) molinfo = other;
  }
  MPI_Bcast(molinfo.data(), 1 + NMOLINFO * n, MPI_DOUBLE, 0, world);

  // donor box is empty: reject without doing anything

  if (molinfo[0] != 0.0) return;

  // put the removed atoms back into the donor box

  auto restore_donor = [&]() {
    size_t pos = 0;
    for (int k = 0; k < nremoved; k++) pos += atom->avec->unpack_exchange(&exchange_buf[pos]);
    atom->natoms += n;
    atom->nbonds += onemol->nbonds;
    atom->nangles += onemol->nangles;
    atom->ndihedrals += onemol->ndihedrals;
    atom->nimpropers += onemol->nimpropers;
    changed_atoms();
  };

  // reject if the receiving box is not wider than twice the molecule
  // (the same condition as for volume moves, see scale_positions()),
  // so that both moves sample the same set of allowed states

  double rmaxsq = 0.0;
  for (int k = 0; k < n; k++) rmaxsq = MAX(rmaxsq, MathExtra::lensq3(&molinfo[1 + NMOLINFO * k]));
  if (any_box(!sender && (min_box_width() <= 4.0 * sqrt(rmaxsq)))) {
    if (sender) restore_donor();
    update_gas_atoms_list();
    return;
  }

  // receiver box: insert the molecule at a random position with a uniformly
  // distributed random orientation (random unit quaternion).
  // velocities are rotated with the molecule, so that time integration
  // continues consistently in hybrid MD/MC simulations

  tagint newmol = 0;
  if (!sender) {
    tagint maxtag = 0, maxmol = 0;
    for (int i = 0; i < atom->nlocal; i++) {
      maxtag = MAX(maxtag, atom->tag[i]);
      maxmol = MAX(maxmol, atom->molecule[i]);
    }
    tagint maxtag_all, maxmol_all;
    MPI_Allreduce(&maxtag, &maxtag_all, 1, MPI_LMP_TAGINT, MPI_MAX, world);
    MPI_Allreduce(&maxmol, &maxmol_all, 1, MPI_LMP_TAGINT, MPI_MAX, world);
    if ((maxtag_all + n >= MAXTAGINT) || (maxmol_all + 1 >= MAXTAGINT))
      error->all(FLERR, "Fix gemc ran out of atom or molecule IDs");
    newmol = maxmol_all + 1;

    double com[3];
    random_point(com);
    double u1 = random_world->uniform();
    double u2 = random_world->uniform();
    double u3 = random_world->uniform();
    double quat[4] = {sqrt(1.0 - u1) * sin(MY_2PI * u2), sqrt(1.0 - u1) * cos(MY_2PI * u2),
                      sqrt(u1) * sin(MY_2PI * u3), sqrt(u1) * cos(MY_2PI * u3)};
    double rotmat[3][3];
    MathExtra::quat_to_mat(quat, rotmat);

    const imageint imagezero =
        ((imageint) IMGMAX << IMG2BITS) | ((imageint) IMGMAX << IMGBITS) | IMGMAX;
    int ncreated = 0;
    for (int k = 0; k < n; k++) {
      double xnew[3];
      double *info = &molinfo[1 + NMOLINFO * k];
      MathExtra::matvec(rotmat, info, xnew);
      xnew[0] += com[0];
      xnew[1] += com[1];
      xnew[2] += com[2];
      imageint imagenew = imagezero;
      domain->remap(xnew, imagenew);
      if (!owns(xnew)) continue;
      ncreated++;

      int itype = onemol->type[k];
      atom->avec->create_atom(itype, xnew);
      int m = atom->nlocal - 1;
      atom->mask[m] = static_cast<int>(info[4]);
      atom->image[m] = imagenew;
      atom->molecule[m] = newmol;
      atom->tag[m] = maxtag_all + k + 1;
      MathExtra::matvec(rotmat, &info[5], atom->v[m]);
      atom->add_molecule_atom(onemol, k, m, maxtag_all);
      if (atom->q_flag) atom->q[m] = info[3];
      modify->create_attribute(m);
    }
    int ncreated_all = 0;
    MPI_Allreduce(&ncreated, &ncreated_all, 1, MPI_INT, MPI_SUM, world);
    if (ncreated_all != n)
      error->all(FLERR, "Fix gemc inserted {} instead of {} atoms of a molecule", ncreated_all, n);
    atom->natoms += n;
    atom->nbonds += onemol->nbonds;
    atom->nangles += onemol->nangles;
    atom->ndihedrals += onemol->ndihedrals;
    atom->nimpropers += onemol->nimpropers;
    changed_atoms();
  }

  // same acceptance rule as for atoms, with N = number of molecules

  double energy_after = energy_full();
  double logfactor;
  if (sender)
    logfactor = log(nold / volume);
  else
    logfactor = log(volume / (nold + 1));
  double dU = energy_after - energy_before - force->boltz * box_temp * logfactor;
  int overflow = !(energy_after < MAXENERGYTEST);

  if (accept_both(dU, overflow)) {
    nexchange_successes += 1.0;
    energy_stored = energy_after;

  } else {

    // rejected: put the removed atoms back or delete the inserted molecule

    if (sender) {
      restore_donor();
    } else {
      int i = 0;
      while (i < atom->nlocal) {
        if (atom->molecule[i] == newmol) {
          atom->avec->copy(atom->nlocal - 1, i, 1);
          atom->nlocal--;
        } else {
          i++;
        }
      }
      atom->natoms -= n;
      atom->nbonds -= onemol->nbonds;
      atom->nangles -= onemol->nangles;
      atom->ndihedrals -= onemol->ndihedrals;
      atom->nimpropers -= onemol->nimpropers;
      changed_atoms();
    }
    energy_stored = energy_before;
  }

  update_gas_atoms_list();
}

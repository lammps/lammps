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
   fix msevb eigensystem:
   Hamiltonian construction, coupling, eigensolve, weight-based
   Hellmann-Feynman force mixing, Fermi-Dirac occupancies.

   Weight-based mixing replaces the dH tensor approach:
   - Grimme coupling → per-state scalar weights
   - Raiteri/Vuilleumier → sparse 3-atom corrections
   Communication: O(natoms*3) allreduce instead of O(ns*natoms*3)
------------------------------------------------------------------------- */

#include "fix_msevb.h"

#include "angle.h"
#include "atom.h"
#include "bond.h"
#include "comm.h"
#include "coupling_grimme2015.h"
#include "coupling_raiteri2011.h"
#include "coupling_vuilleumier1998.h"
#include "dihedral.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "improper.h"
#include "jacobi_eigen.h"
#include "kspace.h"
#include "memory.h"
#include "modify.h"
#include "neighbor.h"
#include "pair.h"
#include "special.h"
#include "universe.h"
#include "update.h"

#include <cmath>
#include <cstring>
#include <vector>

using namespace LAMMPS_NS;

/* ----------------------------------------------------------------------
   Grow dynamically-sized eigensystem arrays to nstates if needed.
   Called at the top of post_force() after nstates is known for this step.
   All arrays that scale with nstates x nstates are managed here.
---------------------------------------------------------------------- */

void FixMSEVB::grow_eigensystem_arrays()
{
  if (nstates <= eigensys_nmax) return;

  memory->destroy(hamiltonian);
  memory->destroy(H_work);
  memory->destroy(eigenvalues);
  memory->destroy(eigenvectors);
  memory->destroy(amplitudes);
  memory->destroy(fd_occ);
  memory->destroy(weights);

  eigensys_nmax = nstates;

  memory->create(hamiltonian, eigensys_nmax * eigensys_nmax, "msevb:hamiltonian");
  memory->create(H_work, eigensys_nmax * eigensys_nmax, "msevb:H_work");
  memory->create(eigenvalues, eigensys_nmax, "msevb:eigenvalues");
  memory->create(eigenvectors, eigensys_nmax * eigensys_nmax, "msevb:eigenvectors");
  memory->create(amplitudes, eigensys_nmax, "msevb:amplitudes");
  memory->create(weights, eigensys_nmax, "msevb:weights");
  if (fermi_dirac_enabled) memory->create(fd_occ, eigensys_nmax, "msevb:fd_occ");

  weights_nmax = eigensys_nmax;
}

/* ----------------------------------------------------------------------
   Grow epot[] to hold nstates entries when nstates has increased.
   my_epot stays npartitions-sized — it is only used for the parallel
   MPI_Allreduce in gather_potential_energies().
---------------------------------------------------------------------- */

void FixMSEVB::grow_epot_arrays()
{
  if (nstates <= epot_nmax) return;

  // Preserve existing values (parallel states were just gathered).
  double *old_epot = epot;
  memory->create(epot, nstates, "msevb:epot");
  for (int i = 0; i < epot_nmax; i++) epot[i] = old_epot[i];
  for (int i = epot_nmax; i < nstates; i++) epot[i] = 0.0;
  memory->destroy(old_epot);
  epot_nmax = nstates;
}

/* ----------------------------------------------------------------------
   Grow excess_forces buffer for excess state force storage (Approach A).
   Layout: excess_forces[serial_idx * ef_nmax * 3 + atom * 3 + xyz]
   Size: nsites_serial * ef_nmax * 3  (vs old dH_all: ns² * nmax * 3)
---------------------------------------------------------------------- */

void FixMSEVB::grow_excess_forces(int nlocal_max)
{
  if (nsites_serial <= excess_forces_nmax_serial && nlocal_max <= excess_forces_nmax) return;

  if (nlocal_max > excess_forces_nmax) excess_forces_nmax = nlocal_max;
  if (nsites_serial > excess_forces_nmax_serial) excess_forces_nmax_serial = nsites_serial;

  delete[] excess_forces;
  const bigint sz = (bigint) excess_forces_nmax_serial * excess_forces_nmax * 3;
  excess_forces = new double[sz]();
}

/* ----------------------------------------------------------------------
   Approach A: Evaluate excess EVB states with save/restore of atom->f.
   Stores per-excess-state forces in excess_forces buffer (diagonal only).
   Uses compute_scalar for Grimme coupling (no off-diagonal dH needed).

   Batched parallel: distribute excess states across all partitions in
   batches of npartitions, so each batch does one parallel force eval.
---------------------------------------------------------------------- */

void FixMSEVB::compute_excess_states()
{
  if (nsites_serial == 0) return;

  const int ns = nstates;
  const int np = npartitions;
  const int eflag = 1, vflag = 0;

  std::vector<int> active_serial;
  active_serial.reserve(nsites_serial);
  for (int b = 0; b < nsites_serial; b++) active_serial.push_back(b);

  const int n_active = static_cast<int>(active_serial.size());
  if (n_active == 0) return;

  const int nbatches = (n_active + np - 1) / np;

  // Save force energy and virial caches — excess evaluations overwrite these,
  // but thermo reads them later via compute_pe / compute_pressure.  Restore
  // after all batches so that thermo sees the partition's own state values.
  double saved_pair_vdwl = 0, saved_pair_coul = 0;
  double saved_bond_energy = 0, saved_angle_energy = 0;
  double saved_dihedral_energy = 0, saved_improper_energy = 0;
  double saved_kspace_energy = 0;
  double saved_pair_virial[6] = {0}, saved_bond_virial[6] = {0};
  double saved_angle_virial[6] = {0}, saved_dihedral_virial[6] = {0};
  double saved_improper_virial[6] = {0}, saved_kspace_virial[6] = {0};
  if (force->pair) {
    saved_pair_vdwl = force->pair->eng_vdwl;
    saved_pair_coul = force->pair->eng_coul;
    memcpy(saved_pair_virial, force->pair->virial, 6 * sizeof(double));
  }
  if (force->bond) {
    saved_bond_energy = force->bond->energy;
    memcpy(saved_bond_virial, force->bond->virial, 6 * sizeof(double));
  }
  if (force->angle) {
    saved_angle_energy = force->angle->energy;
    memcpy(saved_angle_virial, force->angle->virial, 6 * sizeof(double));
  }
  if (force->dihedral) {
    saved_dihedral_energy = force->dihedral->energy;
    memcpy(saved_dihedral_virial, force->dihedral->virial, 6 * sizeof(double));
  }
  if (force->improper) {
    saved_improper_energy = force->improper->energy;
    memcpy(saved_improper_virial, force->improper->virial, 6 * sizeof(double));
  }
  if (force->kspace) {
    saved_kspace_energy = force->kspace->energy;
    memcpy(saved_kspace_virial, force->kspace->virial, 6 * sizeof(double));
  }

  // Save parallel state forces (restored at exit)
  const int save_nlocal = atom->nlocal;
  memory->grow(saved_forces, atom->nmax * 3, "msevb:saved_forces");
  memcpy(saved_forces, &atom->f[0][0], save_nlocal * 3 * sizeof(double));

  // Grow excess_forces buffer
  {
    int my_nmax_atoms = atom->nmax;
    int global_nmax_atoms;
    MPI_Allreduce(&my_nmax_atoms, &global_nmax_atoms, 1, MPI_INT, MPI_MAX, universe->uworld);
    grow_excess_forces(global_nmax_atoms);
  }

  // Zero buffer entries for skipped states
  const int ef_stride = excess_forces_nmax * 3;
  for (int b = 0; b < nsites_serial; b++) {
    bool is_active = false;
    for (int a = 0; a < n_active; a++) {
      if (active_serial[a] == b) {
        is_active = true;
        break;
      }
    }
    if (!is_active) memset(&excess_forces[b * ef_stride], 0, ef_stride * sizeof(double));
  }

  for (int batch = 0; batch < nbatches; batch++) {
    const int batch_start = batch * np;
    const int batch_count = (n_active - batch_start < np) ? n_active - batch_start : np;

    // Which excess state does this partition evaluate?
    const bool active = (ipartition < batch_count);
    const int my_active_idx = active ? batch_start + ipartition : -1;
    const int my_serial_idx = active ? active_serial[my_active_idx] : -1;
    const int sk = active ? nsites_parallel + my_serial_idx : -1;
    const int state_k = active ? sk + 1 : -1;

    // ---- 1. Restore reference topology (all partitions) ---------------
    restore_reference_topology();
    comm->forward_comm(this);

    // ---- 2. Apply excess state topology (active partitions only) ------
    if (active) apply_site_changes(sk);

    // ---- 3. Rebuild specials + neighbors (all partitions) -------------
    {
      Special special_obj(lmp);
      special_obj.build(true);
      if (force->kspace) force->kspace->qsum_qsq(0);
      modified_topology_on_host();
      if (domain->triclinic) domain->x2lamda(atom->nlocal);
      domain->pbc();
      comm->exchange();
      comm->borders();
      if (domain->triclinic) domain->lamda2x(atom->nlocal + atom->nghost);
      sync_before_neighbor_build();
      neighbor->build(1);
    }

    // Grow excess_forces if atom count changed after exchange.
    if (atom->nlocal > excess_forces_nmax) {
      int global_nmax;
      MPI_Allreduce(&atom->nlocal, &global_nmax, 1, MPI_INT, MPI_MAX, world);
      grow_excess_forces(global_nmax);
    }
    const int cur_nlocal = atom->nlocal;

    // ---- 4. Forward comm + zero forces (all partitions) ---------------
    comm->forward_comm(this);
    double **f = atom->f;
    for (int i = 0; i < cur_nlocal; i++) f[i][0] = f[i][1] = f[i][2] = 0.0;
    modified_forces_on_host();

    // ---- 5. Evaluate forces + energy (all partitions in parallel) -----
    sync_before_force_compute();
    if (force->pair && force->pair->compute_flag) force->pair->compute(eflag, vflag);

    if (atom->molecular != Atom::ATOMIC) {
      if (force->bond) force->bond->compute(eflag, vflag);
      if (force->angle) force->angle->compute(eflag, vflag);
      if (force->dihedral) force->dihedral->compute(eflag, vflag);
      if (force->improper) force->improper->compute(eflag, vflag);
    }

    if (force->kspace && force->kspace->compute_flag) force->kspace->compute(eflag, vflag);

    if (force->newton) comm->reverse_comm();
    modified_after_force_compute();

    // ---- 6. Reduce energy within each partition -----------------------
    double eng = 0.0;
    if (force->pair) eng += force->pair->eng_vdwl + force->pair->eng_coul;
    if (atom->molecular != Atom::ATOMIC) {
      if (force->bond) eng += force->bond->energy;
      if (force->angle) eng += force->angle->energy;
      if (force->dihedral) eng += force->dihedral->energy;
      if (force->improper) eng += force->improper->energy;
    }
    MPI_Allreduce(MPI_IN_PLACE, &eng, 1, MPI_DOUBLE, MPI_SUM, world);
    if (force->kspace) eng += force->kspace->energy;

    // ---- 7. Share energies + excess forces across partitions ----------
    for (int b = 0; b < batch_count; b++) {
      const int b_active_idx = batch_start + b;
      const int b_serial_idx = active_serial[b_active_idx];
      const int b_sk = nsites_parallel + b_serial_idx;
      const int b_state = b_sk + 1;

      double b_eng = (ipartition == b) ? eng : 0.0;
      MPI_Bcast(&b_eng, 1, MPI_DOUBLE, b, samerank);
      epot[b_state] = b_eng;
      // Serial-state diagonal = that state's total potential energy (any species
      // offset is already in epot via pair_style template/offset).
      hamiltonian[b_state * ns + b_state] = b_eng;

      // Store excess forces in compact buffer (replaces dH_all diagonal)
      sync_forces_to_host();
      const int ef_off = b_serial_idx * ef_stride;
      if (ipartition == b) {
        for (int a = 0; a < cur_nlocal; a++) {
          excess_forces[ef_off + a * 3 + 0] = f[a][0];
          excess_forces[ef_off + a * 3 + 1] = f[a][1];
          excess_forces[ef_off + a * 3 + 2] = f[a][2];
        }
      }
      MPI_Bcast(&excess_forces[ef_off], ef_stride, MPI_DOUBLE, b, samerank);
    }

    // ---- 8. Compute coupling for each state in this batch -------------
    for (int b = 0; b < batch_count; b++) {
      const int b_active_idx = batch_start + b;
      const int b_sk = nsites_parallel + active_serial[b_active_idx];
      const int b_state = b_sk + 1;

      if (sites[b_sk].n_components > 0) {
        for (int ci = 0; ci < sites[b_sk].n_components; ci++) {
          // Couple this product to the state with component ci removed
          // (neighbor_state[ci]) via component ci's transfer.
          const int par_i = sites[b_sk].neighbor_state[ci];
          const int comp_j = sites[b_sk].components[ci];
          const ReactionDef &rxn_cpl = rxndefs[sites[comp_j].rxn_idx];

          if (rxn_cpl.coupling_type == COUPLING_NONE) continue;

          if (rxn_cpl.coupling_type == COUPLING_GRIMME2015) {
            double V, g;
            if (rxn_cpl.coupling_taper > 0.0) {
              int gidx_H = atom->map(sites[comp_j].tag_H);
              int gidx_Y = atom->map(sites[comp_j].tag_Y);
              double tfH[3], tfY[3];
              if (gidx_H >= 0 && gidx_Y >= 0)
                Grimme2015::compute_scalar_tapered(
                    rxn_cpl.coupling_a, rxn_cpl.coupling_b, hamiltonian[par_i * ns + par_i],
                    hamiltonian[b_state * ns + b_state], atom->x[gidx_H], atom->x[gidx_Y], domain,
                    rxn_cpl.coupling_taper, std::sqrt(rxn_cpl.cutoff_sq), V, g, tfH, tfY);
              else
                Grimme2015::compute_scalar(rxn_cpl.coupling_a, rxn_cpl.coupling_b,
                                           hamiltonian[par_i * ns + par_i],
                                           hamiltonian[b_state * ns + b_state], V, g);
            } else {
              Grimme2015::compute_scalar(rxn_cpl.coupling_a, rxn_cpl.coupling_b,
                                         hamiltonian[par_i * ns + par_i],
                                         hamiltonian[b_state * ns + b_state], V, g);
            }
            hamiltonian[par_i * ns + b_state] = V;
            hamiltonian[b_state * ns + par_i] = V;
          } else {
            int idx_H = atom->map(sites[comp_j].tag_H);
            int idx_X = atom->map(sites[comp_j].tag_X);
            int idx_Y = atom->map(sites[comp_j].tag_Y);

            if (idx_H >= 0 && idx_X >= 0 && idx_Y >= 0) {
              double **x = atom->x;
              double V, fH[3], fX[3], fY[3];
              double cutoff = std::sqrt(rxn_cpl.cutoff_sq);

              if (rxn_cpl.coupling_type == COUPLING_RAITERI2011) {
                if (rxn_cpl.coupling_taper > 0.0)
                  Raiteri2011::compute_tapered(rxn_cpl.coupling_lambda, rxn_cpl.coupling_zeta,
                                               x[idx_H], x[idx_X], x[idx_Y], domain,
                                               rxn_cpl.coupling_taper, cutoff, V, fH, fX, fY);
                else
                  Raiteri2011::compute(rxn_cpl.coupling_lambda, rxn_cpl.coupling_zeta, x[idx_H],
                                       x[idx_X], x[idx_Y], domain, V, fH, fX, fY);
              } else {
                if (rxn_cpl.coupling_taper > 0.0)
                  Vuilleumier1998::compute_tapered(rxn_cpl.coupling_v12, rxn_cpl.coupling_alpha,
                                                   rxn_cpl.coupling_gamma_v, x[idx_H], x[idx_X],
                                                   x[idx_Y], domain, rxn_cpl.coupling_taper, cutoff,
                                                   V, fH, fX, fY);
                else
                  Vuilleumier1998::compute(rxn_cpl.coupling_v12, rxn_cpl.coupling_alpha,
                                           rxn_cpl.coupling_gamma_v, x[idx_H], x[idx_X], x[idx_Y],
                                           domain, V, fH, fX, fY);
              }

              if (comm->me == 0) {
                hamiltonian[par_i * ns + b_state] = V;
                hamiltonian[b_state * ns + par_i] = V;
              }
              MPI_Allreduce(MPI_IN_PLACE, &hamiltonian[par_i * ns + b_state], 1, MPI_DOUBLE,
                            MPI_SUM, world);
              hamiltonian[b_state * ns + par_i] = hamiltonian[par_i * ns + b_state];
            }
          }
        }
      } else {
        const int parent = sites[b_sk].parent_state;
        const ReactionDef &rxn_cpl = rxndefs[sites[b_sk].rxn_idx];

        if (rxn_cpl.coupling_type == COUPLING_NONE) continue;

        if (rxn_cpl.coupling_type == COUPLING_GRIMME2015) {
          double V, g;
          if (rxn_cpl.coupling_taper > 0.0) {
            int gidx_H = atom->map(sites[b_sk].tag_H);
            int gidx_Y = atom->map(sites[b_sk].tag_Y);
            double tfH[3], tfY[3];
            if (gidx_H >= 0 && gidx_Y >= 0)
              Grimme2015::compute_scalar_tapered(
                  rxn_cpl.coupling_a, rxn_cpl.coupling_b, hamiltonian[parent * ns + parent],
                  hamiltonian[b_state * ns + b_state], atom->x[gidx_H], atom->x[gidx_Y], domain,
                  rxn_cpl.coupling_taper, std::sqrt(rxn_cpl.cutoff_sq), V, g, tfH, tfY);
            else
              Grimme2015::compute_scalar(rxn_cpl.coupling_a, rxn_cpl.coupling_b,
                                         hamiltonian[parent * ns + parent],
                                         hamiltonian[b_state * ns + b_state], V, g);
          } else {
            Grimme2015::compute_scalar(rxn_cpl.coupling_a, rxn_cpl.coupling_b,
                                       hamiltonian[parent * ns + parent],
                                       hamiltonian[b_state * ns + b_state], V, g);
          }
          hamiltonian[parent * ns + b_state] = V;
          hamiltonian[b_state * ns + parent] = V;
        } else {
          int idx_H = atom->map(sites[b_sk].tag_H);
          int idx_X = atom->map(sites[b_sk].tag_X);
          int idx_Y = atom->map(sites[b_sk].tag_Y);

          if (idx_H >= 0 && idx_X >= 0 && idx_Y >= 0) {
            double **x = atom->x;
            double V, fH[3], fX[3], fY[3];
            double cutoff = std::sqrt(rxn_cpl.cutoff_sq);

            if (rxn_cpl.coupling_type == COUPLING_RAITERI2011) {
              if (rxn_cpl.coupling_taper > 0.0)
                Raiteri2011::compute_tapered(rxn_cpl.coupling_lambda, rxn_cpl.coupling_zeta,
                                             x[idx_H], x[idx_X], x[idx_Y], domain,
                                             rxn_cpl.coupling_taper, cutoff, V, fH, fX, fY);
              else
                Raiteri2011::compute(rxn_cpl.coupling_lambda, rxn_cpl.coupling_zeta, x[idx_H],
                                     x[idx_X], x[idx_Y], domain, V, fH, fX, fY);
            } else {
              if (rxn_cpl.coupling_taper > 0.0)
                Vuilleumier1998::compute_tapered(rxn_cpl.coupling_v12, rxn_cpl.coupling_alpha,
                                                 rxn_cpl.coupling_gamma_v, x[idx_H], x[idx_X],
                                                 x[idx_Y], domain, rxn_cpl.coupling_taper, cutoff,
                                                 V, fH, fX, fY);
              else
                Vuilleumier1998::compute(rxn_cpl.coupling_v12, rxn_cpl.coupling_alpha,
                                         rxn_cpl.coupling_gamma_v, x[idx_H], x[idx_X], x[idx_Y],
                                         domain, V, fH, fX, fY);
            }

            if (comm->me == 0) {
              hamiltonian[parent * ns + b_state] = V;
              hamiltonian[b_state * ns + parent] = V;
            }
            MPI_Allreduce(MPI_IN_PLACE, &hamiltonian[parent * ns + b_state], 1, MPI_DOUBLE, MPI_SUM,
                          world);
            hamiltonian[b_state * ns + parent] = hamiltonian[parent * ns + b_state];
          }
        }
      }
    }
  }    // end batch loop

  // Restore reference topology on all partitions after all batches.
  restore_reference_topology();
  {
    Special special_obj(lmp);
    special_obj.build(true);
    // Excess-state evaluation changed the charges (and thus qsqsum) via the
    // per-batch qsum_qsq() above; restoring the reference topology must also
    // restore the reference qsqsum, otherwise partition 0 keeps the last excess
    // state's qsqsum and the reference kspace self-energy is wrong on the next
    // step.
    if (force->kspace) force->kspace->qsum_qsq(0);
    if (domain->triclinic) domain->x2lamda(atom->nlocal);
    domain->pbc();
    comm->exchange();
    comm->borders();
    if (domain->triclinic) domain->lamda2x(atom->nlocal + atom->nghost);
    neighbor->build(1);
  }

  // Restore force energy and virial caches so thermo sees partition's own values
  if (force->pair) {
    force->pair->eng_vdwl = saved_pair_vdwl;
    force->pair->eng_coul = saved_pair_coul;
    memcpy(force->pair->virial, saved_pair_virial, 6 * sizeof(double));
  }
  if (force->bond) {
    force->bond->energy = saved_bond_energy;
    memcpy(force->bond->virial, saved_bond_virial, 6 * sizeof(double));
  }
  if (force->angle) {
    force->angle->energy = saved_angle_energy;
    memcpy(force->angle->virial, saved_angle_virial, 6 * sizeof(double));
  }
  if (force->dihedral) {
    force->dihedral->energy = saved_dihedral_energy;
    memcpy(force->dihedral->virial, saved_dihedral_virial, 6 * sizeof(double));
  }
  if (force->improper) {
    force->improper->energy = saved_improper_energy;
    memcpy(force->improper->virial, saved_improper_virial, 6 * sizeof(double));
  }
  if (force->kspace) {
    force->kspace->energy = saved_kspace_energy;
    memcpy(force->kspace->virial, saved_kspace_virial, 6 * sizeof(double));
  }

  // Restore parallel state forces (saved at entry)
  memcpy(&atom->f[0][0], saved_forces, save_nlocal * 3 * sizeof(double));
}

/* ----------------------------------------------------------------------
   Approach B, Phase 1: Evaluate excess state energies and coupling values
   only.  Does NOT store forces — atom->f is left in an undefined state.
   Called before eigensystem solve (only needs Hamiltonian entries).
---------------------------------------------------------------------- */

void FixMSEVB::compute_excess_energies()
{
  if (nsites_serial == 0) return;

  const int ns = nstates;
  const int np = npartitions;
  const int eflag = 1, vflag = 0;

  std::vector<int> active_serial;
  active_serial.reserve(nsites_serial);
  for (int b = 0; b < nsites_serial; b++) active_serial.push_back(b);

  const int n_active = static_cast<int>(active_serial.size());
  if (n_active == 0) return;

  const int nbatches = (n_active + np - 1) / np;

  // Save force energy and virial caches (same rationale as compute_excess_states)
  double saved_pair_vdwl = 0, saved_pair_coul = 0;
  double saved_bond_energy = 0, saved_angle_energy = 0;
  double saved_dihedral_energy = 0, saved_improper_energy = 0;
  double saved_kspace_energy = 0;
  double saved_pair_virial[6] = {0}, saved_bond_virial[6] = {0};
  double saved_angle_virial[6] = {0}, saved_dihedral_virial[6] = {0};
  double saved_improper_virial[6] = {0}, saved_kspace_virial[6] = {0};
  if (force->pair) {
    saved_pair_vdwl = force->pair->eng_vdwl;
    saved_pair_coul = force->pair->eng_coul;
    memcpy(saved_pair_virial, force->pair->virial, 6 * sizeof(double));
  }
  if (force->bond) {
    saved_bond_energy = force->bond->energy;
    memcpy(saved_bond_virial, force->bond->virial, 6 * sizeof(double));
  }
  if (force->angle) {
    saved_angle_energy = force->angle->energy;
    memcpy(saved_angle_virial, force->angle->virial, 6 * sizeof(double));
  }
  if (force->dihedral) {
    saved_dihedral_energy = force->dihedral->energy;
    memcpy(saved_dihedral_virial, force->dihedral->virial, 6 * sizeof(double));
  }
  if (force->improper) {
    saved_improper_energy = force->improper->energy;
    memcpy(saved_improper_virial, force->improper->virial, 6 * sizeof(double));
  }
  if (force->kspace) {
    saved_kspace_energy = force->kspace->energy;
    memcpy(saved_kspace_virial, force->kspace->virial, 6 * sizeof(double));
  }

  // Save parallel state forces (excess eval destroys atom->f)
  const int save_nlocal = atom->nlocal;
  memory->grow(saved_forces, atom->nmax * 3, "msevb:saved_forces");
  memcpy(saved_forces, &atom->f[0][0], save_nlocal * 3 * sizeof(double));

  for (int batch = 0; batch < nbatches; batch++) {
    const int batch_start = batch * np;
    const int batch_count = (n_active - batch_start < np) ? n_active - batch_start : np;

    const bool active = (ipartition < batch_count);
    const int my_active_idx = active ? batch_start + ipartition : -1;
    const int my_serial_idx = active ? active_serial[my_active_idx] : -1;
    const int sk = active ? nsites_parallel + my_serial_idx : -1;
    const int state_k = active ? sk + 1 : -1;

    // 1. Restore reference topology
    restore_reference_topology();
    comm->forward_comm(this);

    // 2. Apply excess state topology
    if (active) apply_site_changes(sk);

    // 3. Rebuild specials + neighbors
    {
      Special special_obj(lmp);
      special_obj.build(true);
      if (force->kspace) force->kspace->qsum_qsq(0);
      modified_topology_on_host();
      if (domain->triclinic) domain->x2lamda(atom->nlocal);
      domain->pbc();
      comm->exchange();
      comm->borders();
      if (domain->triclinic) domain->lamda2x(atom->nlocal + atom->nghost);
      sync_before_neighbor_build();
      neighbor->build(1);
    }

    const int cur_nlocal = atom->nlocal;

    // 4. Forward comm + zero forces
    comm->forward_comm(this);
    double **f = atom->f;
    for (int i = 0; i < cur_nlocal; i++) f[i][0] = f[i][1] = f[i][2] = 0.0;
    modified_forces_on_host();

    // 5. Evaluate forces + energy
    sync_before_force_compute();
    if (force->pair && force->pair->compute_flag) force->pair->compute(eflag, vflag);

    if (atom->molecular != Atom::ATOMIC) {
      if (force->bond) force->bond->compute(eflag, vflag);
      if (force->angle) force->angle->compute(eflag, vflag);
      if (force->dihedral) force->dihedral->compute(eflag, vflag);
      if (force->improper) force->improper->compute(eflag, vflag);
    }

    if (force->kspace && force->kspace->compute_flag) force->kspace->compute(eflag, vflag);

    // No reverse_comm needed — forces not stored
    modified_after_force_compute();

    // 6. Reduce energy
    double eng = 0.0;
    if (force->pair) eng += force->pair->eng_vdwl + force->pair->eng_coul;
    if (atom->molecular != Atom::ATOMIC) {
      if (force->bond) eng += force->bond->energy;
      if (force->angle) eng += force->angle->energy;
      if (force->dihedral) eng += force->dihedral->energy;
      if (force->improper) eng += force->improper->energy;
    }
    MPI_Allreduce(MPI_IN_PLACE, &eng, 1, MPI_DOUBLE, MPI_SUM, world);
    if (force->kspace) eng += force->kspace->energy;

    // 7. Share energies across partitions (no force storage)
    for (int b = 0; b < batch_count; b++) {
      const int b_active_idx = batch_start + b;
      const int b_serial_idx = active_serial[b_active_idx];
      const int b_sk = nsites_parallel + b_serial_idx;
      const int b_state = b_sk + 1;

      double b_eng = (ipartition == b) ? eng : 0.0;
      MPI_Bcast(&b_eng, 1, MPI_DOUBLE, b, samerank);
      epot[b_state] = b_eng;
      // Serial-state diagonal = that state's total potential energy (any species
      // offset is already in epot via pair_style template/offset).
      hamiltonian[b_state * ns + b_state] = b_eng;
    }

    // 8. Compute coupling (scalar only — no dH needed)
    for (int b = 0; b < batch_count; b++) {
      const int b_active_idx = batch_start + b;
      const int b_sk = nsites_parallel + active_serial[b_active_idx];
      const int b_state = b_sk + 1;

      if (sites[b_sk].n_components > 0) {
        for (int ci = 0; ci < sites[b_sk].n_components; ci++) {
          // Couple this product to the state with component ci removed
          // (neighbor_state[ci]) via component ci's transfer.
          const int par_i = sites[b_sk].neighbor_state[ci];
          const int comp_j = sites[b_sk].components[ci];
          const ReactionDef &rxn_cpl = rxndefs[sites[comp_j].rxn_idx];

          if (rxn_cpl.coupling_type == COUPLING_NONE) continue;

          if (rxn_cpl.coupling_type == COUPLING_GRIMME2015) {
            double V, g;
            if (rxn_cpl.coupling_taper > 0.0) {
              int gidx_H = atom->map(sites[comp_j].tag_H);
              int gidx_Y = atom->map(sites[comp_j].tag_Y);
              double tfH[3], tfY[3];
              if (gidx_H >= 0 && gidx_Y >= 0)
                Grimme2015::compute_scalar_tapered(
                    rxn_cpl.coupling_a, rxn_cpl.coupling_b, hamiltonian[par_i * ns + par_i],
                    hamiltonian[b_state * ns + b_state], atom->x[gidx_H], atom->x[gidx_Y], domain,
                    rxn_cpl.coupling_taper, std::sqrt(rxn_cpl.cutoff_sq), V, g, tfH, tfY);
              else
                Grimme2015::compute_scalar(rxn_cpl.coupling_a, rxn_cpl.coupling_b,
                                           hamiltonian[par_i * ns + par_i],
                                           hamiltonian[b_state * ns + b_state], V, g);
            } else {
              Grimme2015::compute_scalar(rxn_cpl.coupling_a, rxn_cpl.coupling_b,
                                         hamiltonian[par_i * ns + par_i],
                                         hamiltonian[b_state * ns + b_state], V, g);
            }
            hamiltonian[par_i * ns + b_state] = V;
            hamiltonian[b_state * ns + par_i] = V;
          } else {
            int idx_H = atom->map(sites[comp_j].tag_H);
            int idx_X = atom->map(sites[comp_j].tag_X);
            int idx_Y = atom->map(sites[comp_j].tag_Y);

            if (idx_H >= 0 && idx_X >= 0 && idx_Y >= 0) {
              double **x = atom->x;
              double V, fH[3], fX[3], fY[3];
              double cutoff = std::sqrt(rxn_cpl.cutoff_sq);

              if (rxn_cpl.coupling_type == COUPLING_RAITERI2011) {
                if (rxn_cpl.coupling_taper > 0.0)
                  Raiteri2011::compute_tapered(rxn_cpl.coupling_lambda, rxn_cpl.coupling_zeta,
                                               x[idx_H], x[idx_X], x[idx_Y], domain,
                                               rxn_cpl.coupling_taper, cutoff, V, fH, fX, fY);
                else
                  Raiteri2011::compute(rxn_cpl.coupling_lambda, rxn_cpl.coupling_zeta, x[idx_H],
                                       x[idx_X], x[idx_Y], domain, V, fH, fX, fY);
              } else {
                if (rxn_cpl.coupling_taper > 0.0)
                  Vuilleumier1998::compute_tapered(rxn_cpl.coupling_v12, rxn_cpl.coupling_alpha,
                                                   rxn_cpl.coupling_gamma_v, x[idx_H], x[idx_X],
                                                   x[idx_Y], domain, rxn_cpl.coupling_taper, cutoff,
                                                   V, fH, fX, fY);
                else
                  Vuilleumier1998::compute(rxn_cpl.coupling_v12, rxn_cpl.coupling_alpha,
                                           rxn_cpl.coupling_gamma_v, x[idx_H], x[idx_X], x[idx_Y],
                                           domain, V, fH, fX, fY);
              }

              if (comm->me == 0) {
                hamiltonian[par_i * ns + b_state] = V;
                hamiltonian[b_state * ns + par_i] = V;
              }
              MPI_Allreduce(MPI_IN_PLACE, &hamiltonian[par_i * ns + b_state], 1, MPI_DOUBLE,
                            MPI_SUM, world);
              hamiltonian[b_state * ns + par_i] = hamiltonian[par_i * ns + b_state];
            }
          }
        }
      } else {
        const int parent = sites[b_sk].parent_state;
        const ReactionDef &rxn_cpl = rxndefs[sites[b_sk].rxn_idx];

        if (rxn_cpl.coupling_type == COUPLING_NONE) continue;

        if (rxn_cpl.coupling_type == COUPLING_GRIMME2015) {
          double V, g;
          if (rxn_cpl.coupling_taper > 0.0) {
            int gidx_H = atom->map(sites[b_sk].tag_H);
            int gidx_Y = atom->map(sites[b_sk].tag_Y);
            double tfH[3], tfY[3];
            if (gidx_H >= 0 && gidx_Y >= 0)
              Grimme2015::compute_scalar_tapered(
                  rxn_cpl.coupling_a, rxn_cpl.coupling_b, hamiltonian[parent * ns + parent],
                  hamiltonian[b_state * ns + b_state], atom->x[gidx_H], atom->x[gidx_Y], domain,
                  rxn_cpl.coupling_taper, std::sqrt(rxn_cpl.cutoff_sq), V, g, tfH, tfY);
            else
              Grimme2015::compute_scalar(rxn_cpl.coupling_a, rxn_cpl.coupling_b,
                                         hamiltonian[parent * ns + parent],
                                         hamiltonian[b_state * ns + b_state], V, g);
          } else {
            Grimme2015::compute_scalar(rxn_cpl.coupling_a, rxn_cpl.coupling_b,
                                       hamiltonian[parent * ns + parent],
                                       hamiltonian[b_state * ns + b_state], V, g);
          }
          hamiltonian[parent * ns + b_state] = V;
          hamiltonian[b_state * ns + parent] = V;
        } else {
          int idx_H = atom->map(sites[b_sk].tag_H);
          int idx_X = atom->map(sites[b_sk].tag_X);
          int idx_Y = atom->map(sites[b_sk].tag_Y);

          if (idx_H >= 0 && idx_X >= 0 && idx_Y >= 0) {
            double **x = atom->x;
            double V, fH[3], fX[3], fY[3];
            double cutoff = std::sqrt(rxn_cpl.cutoff_sq);

            if (rxn_cpl.coupling_type == COUPLING_RAITERI2011) {
              if (rxn_cpl.coupling_taper > 0.0)
                Raiteri2011::compute_tapered(rxn_cpl.coupling_lambda, rxn_cpl.coupling_zeta,
                                             x[idx_H], x[idx_X], x[idx_Y], domain,
                                             rxn_cpl.coupling_taper, cutoff, V, fH, fX, fY);
              else
                Raiteri2011::compute(rxn_cpl.coupling_lambda, rxn_cpl.coupling_zeta, x[idx_H],
                                     x[idx_X], x[idx_Y], domain, V, fH, fX, fY);
            } else {
              if (rxn_cpl.coupling_taper > 0.0)
                Vuilleumier1998::compute_tapered(rxn_cpl.coupling_v12, rxn_cpl.coupling_alpha,
                                                 rxn_cpl.coupling_gamma_v, x[idx_H], x[idx_X],
                                                 x[idx_Y], domain, rxn_cpl.coupling_taper, cutoff,
                                                 V, fH, fX, fY);
              else
                Vuilleumier1998::compute(rxn_cpl.coupling_v12, rxn_cpl.coupling_alpha,
                                         rxn_cpl.coupling_gamma_v, x[idx_H], x[idx_X], x[idx_Y],
                                         domain, V, fH, fX, fY);
            }

            if (comm->me == 0) {
              hamiltonian[parent * ns + b_state] = V;
              hamiltonian[b_state * ns + parent] = V;
            }
            MPI_Allreduce(MPI_IN_PLACE, &hamiltonian[parent * ns + b_state], 1, MPI_DOUBLE, MPI_SUM,
                          world);
            hamiltonian[b_state * ns + parent] = hamiltonian[parent * ns + b_state];
          }
        }
      }
    }
  }    // end batch loop

  // Restore reference topology
  restore_reference_topology();
  {
    Special special_obj(lmp);
    special_obj.build(true);
    // Restore the reference qsqsum too (excess-state evaluation changed the
    // charges); otherwise partition 0 keeps the last excess state's qsqsum and
    // the reference kspace self-energy is wrong on the next step.
    if (force->kspace) force->kspace->qsum_qsq(0);
    if (domain->triclinic) domain->x2lamda(atom->nlocal);
    domain->pbc();
    comm->exchange();
    comm->borders();
    if (domain->triclinic) domain->lamda2x(atom->nlocal + atom->nghost);
    neighbor->build(1);
  }

  // Restore force energy and virial caches so thermo sees partition's own values
  if (force->pair) {
    force->pair->eng_vdwl = saved_pair_vdwl;
    force->pair->eng_coul = saved_pair_coul;
    memcpy(force->pair->virial, saved_pair_virial, 6 * sizeof(double));
  }
  if (force->bond) {
    force->bond->energy = saved_bond_energy;
    memcpy(force->bond->virial, saved_bond_virial, 6 * sizeof(double));
  }
  if (force->angle) {
    force->angle->energy = saved_angle_energy;
    memcpy(force->angle->virial, saved_angle_virial, 6 * sizeof(double));
  }
  if (force->dihedral) {
    force->dihedral->energy = saved_dihedral_energy;
    memcpy(force->dihedral->virial, saved_dihedral_virial, 6 * sizeof(double));
  }
  if (force->improper) {
    force->improper->energy = saved_improper_energy;
    memcpy(force->improper->virial, saved_improper_virial, 6 * sizeof(double));
  }
  if (force->kspace) {
    force->kspace->energy = saved_kspace_energy;
    memcpy(force->kspace->virial, saved_kspace_virial, 6 * sizeof(double));
  }

  // Restore parallel state forces
  memcpy(&atom->f[0][0], saved_forces, save_nlocal * 3 * sizeof(double));
}

/* ----------------------------------------------------------------------
   Approach B, Phase 2: Re-evaluate excess state forces with known weights
   and accumulate directly into atom->f.  No force buffer needed — forces
   are applied immediately with their weights.

   Called after weight_based_hellmann_feynman_forces() has mixed the
   parallel state forces.  atom->f contains the mixed parallel forces.
---------------------------------------------------------------------- */

void FixMSEVB::apply_excess_forces()
{
  if (nsites_serial == 0) return;

  const int ns = nstates;
  const int np = npartitions;
  const int eflag = 0, vflag = 0;    // energy already known, only need forces

  std::vector<int> active_serial;
  active_serial.reserve(nsites_serial);
  for (int b = 0; b < nsites_serial; b++) active_serial.push_back(b);

  const int n_active = static_cast<int>(active_serial.size());
  if (n_active == 0) return;

  const int nbatches = (n_active + np - 1) / np;

  // Save mixed parallel forces (will accumulate excess contributions)
  const int save_nlocal = atom->nlocal;
  memory->grow(saved_forces, atom->nmax * 3, "msevb:saved_forces");
  memcpy(saved_forces, &atom->f[0][0], save_nlocal * 3 * sizeof(double));

  for (int batch = 0; batch < nbatches; batch++) {
    const int batch_start = batch * np;
    const int batch_count = (n_active - batch_start < np) ? n_active - batch_start : np;

    const bool active = (ipartition < batch_count);
    const int my_active_idx = active ? batch_start + ipartition : -1;
    const int my_serial_idx = active ? active_serial[my_active_idx] : -1;
    const int sk = active ? nsites_parallel + my_serial_idx : -1;
    const int state_k = active ? sk + 1 : -1;

    // 1. Restore reference topology
    restore_reference_topology();
    comm->forward_comm(this);

    // 2. Apply excess state topology
    if (active) apply_site_changes(sk);

    // 3. Rebuild specials + neighbors
    {
      Special special_obj(lmp);
      special_obj.build(true);
      if (force->kspace) force->kspace->qsum_qsq(0);
      modified_topology_on_host();
      if (domain->triclinic) domain->x2lamda(atom->nlocal);
      domain->pbc();
      comm->exchange();
      comm->borders();
      if (domain->triclinic) domain->lamda2x(atom->nlocal + atom->nghost);
      sync_before_neighbor_build();
      neighbor->build(1);
    }

    const int cur_nlocal = atom->nlocal;

    // 4. Forward comm + zero forces
    comm->forward_comm(this);
    double **f = atom->f;
    for (int i = 0; i < cur_nlocal; i++) f[i][0] = f[i][1] = f[i][2] = 0.0;
    modified_forces_on_host();

    // 5. Evaluate forces only (no energy needed)
    sync_before_force_compute();
    if (force->pair && force->pair->compute_flag) force->pair->compute(eflag, vflag);

    if (atom->molecular != Atom::ATOMIC) {
      if (force->bond) force->bond->compute(eflag, vflag);
      if (force->angle) force->angle->compute(eflag, vflag);
      if (force->dihedral) force->dihedral->compute(eflag, vflag);
      if (force->improper) force->improper->compute(eflag, vflag);
    }

    if (force->kspace && force->kspace->compute_flag) force->kspace->compute(eflag, vflag);

    if (force->newton) comm->reverse_comm();
    modified_after_force_compute();
    sync_forces_to_host();

    // 6. Scale forces by weight, allreduce across partitions, accumulate.
    //    Active partitions scale by their excess state weight; inactive zero.
    if (active) {
      double wk = weights[state_k];
      for (int a = 0; a < cur_nlocal; a++) {
        f[a][0] *= wk;
        f[a][1] *= wk;
        f[a][2] *= wk;
      }
    } else {
      for (int a = 0; a < cur_nlocal; a++) f[a][0] = f[a][1] = f[a][2] = 0.0;
    }

    MPI_Allreduce(MPI_IN_PLACE, &f[0][0], cur_nlocal * 3, MPI_DOUBLE, MPI_SUM, samerank);

    for (int a = 0; a < cur_nlocal; a++) {
      saved_forces[a * 3 + 0] += f[a][0];
      saved_forces[a * 3 + 1] += f[a][1];
      saved_forces[a * 3 + 2] += f[a][2];
    }
  }    // end batch loop

  // Restore reference topology
  restore_reference_topology();
  {
    Special special_obj(lmp);
    special_obj.build(true);
    // Restore the reference qsqsum too (excess-state evaluation changed the
    // charges); otherwise partition 0 keeps the last excess state's qsqsum and
    // the reference kspace self-energy is wrong on the next step.
    if (force->kspace) force->kspace->qsum_qsq(0);
    if (domain->triclinic) domain->x2lamda(atom->nlocal);
    domain->pbc();
    comm->exchange();
    comm->borders();
    if (domain->triclinic) domain->lamda2x(atom->nlocal + atom->nghost);
    neighbor->build(1);
  }

  // Restore mixed forces (parallel + excess contributions)
  memcpy(&atom->f[0][0], saved_forces, save_nlocal * 3 * sizeof(double));
}

/* ----------------------------------------------------------------------
   Build Hamiltonian matrix from per-partition potential energies.
   Diagonal only — off-diagonal coupling added in compute_coupling_values().
   Uses nstates x nstates layout (not npartitions x npartitions).
---------------------------------------------------------------------- */

void FixMSEVB::build_hamiltonian()
{
  const int ns = nstates;

  for (int i = 0; i < ns * ns; i++) hamiltonian[i] = 0.0;

  // Each diagonal is the total potential energy of that state's topology.  Any
  // per-species diabatic energy shift is supplied externally by adding a
  // pair_style template/offset to the force field, so it is already included in
  // epot[i] for every state -- the fix itself carries no offset bookkeeping.
  for (int i = 0; i < ns; i++) hamiltonian[i * ns + i] = epot[i];
}

/* ----------------------------------------------------------------------
   Evaluate coupling for a single parent-child edge.
   Dispatches to the correct coupling model based on rxn.coupling_type.
   Returns CouplingResult with V, g (Grimme), forces (position-based or
   Grimme taper), and validity flags.
---------------------------------------------------------------------- */

FixMSEVB::CouplingResult FixMSEVB::evaluate_coupling(const ReactionDef &rxn, double E_parent,
                                                     double E_child, tagint tag_H, tagint tag_X,
                                                     tagint tag_Y)
{
  CouplingResult cr = {};

  switch (rxn.coupling_type) {
    case COUPLING_GRIMME2015: {
      cr.is_grimme = true;
      cr.valid = true;
      if (rxn.coupling_taper > 0.0) {
        int idx_H = atom->map(tag_H);
        int idx_Y = atom->map(tag_Y);
        if (idx_H >= 0 && idx_Y >= 0) {
          cr.has_forces = true;
          Grimme2015::compute_scalar_tapered(
              rxn.coupling_a, rxn.coupling_b, E_parent, E_child, atom->x[idx_H], atom->x[idx_Y],
              domain, rxn.coupling_taper, std::sqrt(rxn.cutoff_sq), cr.V, cr.g, cr.fH, cr.fY);
        } else {
          Grimme2015::compute_scalar(rxn.coupling_a, rxn.coupling_b, E_parent, E_child, cr.V, cr.g);
        }
      } else {
        Grimme2015::compute_scalar(rxn.coupling_a, rxn.coupling_b, E_parent, E_child, cr.V, cr.g);
      }
      break;
    }

    case COUPLING_RAITERI2011: {
      int idx_H = atom->map(tag_H);
      int idx_X = atom->map(tag_X);
      int idx_Y = atom->map(tag_Y);
      if (idx_H >= 0 && idx_X >= 0 && idx_Y >= 0) {
        cr.valid = true;
        cr.has_forces = true;
        double **x = atom->x;
        if (rxn.coupling_taper > 0.0)
          Raiteri2011::compute_tapered(rxn.coupling_lambda, rxn.coupling_zeta, x[idx_H], x[idx_X],
                                       x[idx_Y], domain, rxn.coupling_taper,
                                       std::sqrt(rxn.cutoff_sq), cr.V, cr.fH, cr.fX, cr.fY);
        else
          Raiteri2011::compute(rxn.coupling_lambda, rxn.coupling_zeta, x[idx_H], x[idx_X], x[idx_Y],
                               domain, cr.V, cr.fH, cr.fX, cr.fY);
      }
      break;
    }

    case COUPLING_VUILLEUMIER1998: {
      int idx_H = atom->map(tag_H);
      int idx_X = atom->map(tag_X);
      int idx_Y = atom->map(tag_Y);
      if (idx_H >= 0 && idx_X >= 0 && idx_Y >= 0) {
        cr.valid = true;
        cr.has_forces = true;
        double **x = atom->x;
        if (rxn.coupling_taper > 0.0)
          Vuilleumier1998::compute_tapered(rxn.coupling_v12, rxn.coupling_alpha,
                                           rxn.coupling_gamma_v, x[idx_H], x[idx_X], x[idx_Y],
                                           domain, rxn.coupling_taper, std::sqrt(rxn.cutoff_sq),
                                           cr.V, cr.fH, cr.fX, cr.fY);
        else
          Vuilleumier1998::compute(rxn.coupling_v12, rxn.coupling_alpha, rxn.coupling_gamma_v,
                                   x[idx_H], x[idx_X], x[idx_Y], domain, cr.V, cr.fH, cr.fX, cr.fY);
      }
      break;
    }

    default:
      break;
  }

  return cr;
}

/* ----------------------------------------------------------------------
   Compute coupling values (V) for all parallel sites and fill
   hamiltonian off-diagonal.  All partitions compute independently
   (identical positions and energies after sync).
---------------------------------------------------------------------- */

void FixMSEVB::compute_coupling_values()
{
  if (nsites == 0) return;

  const int ns = nstates;

  // Helper lambda: process one coupling edge, fill hamiltonian
  auto fill_hamiltonian_edge = [&](int parent, int child, const ReactiveSite &site) {
    const ReactionDef &rxn = rxndefs[site.rxn_idx];
    if (rxn.coupling_type == COUPLING_NONE) return;

    auto cr =
        evaluate_coupling(rxn, hamiltonian[parent * ns + parent], hamiltonian[child * ns + child],
                          site.tag_H, site.tag_X, site.tag_Y);

    if (cr.is_grimme) {
      // Grimme: all ranks compute identically from synced energies
      hamiltonian[parent * ns + child] = cr.V;
      hamiltonian[child * ns + parent] = cr.V;
    } else if (cr.valid) {
      // Position-based: rank 0 contributes, Allreduce sums
      if (comm->me == 0) {
        hamiltonian[parent * ns + child] = cr.V;
        hamiltonian[child * ns + parent] = cr.V;
      }
      MPI_Allreduce(MPI_IN_PLACE, &hamiltonian[parent * ns + child], 1, MPI_DOUBLE, MPI_SUM, world);
      hamiltonian[child * ns + parent] = hamiltonian[parent * ns + child];
    }
  };

  for (int k = 0; k < nsites_parallel; k++) {
    int state_k = k + 1;
    if (sites[k].n_components > 0) {
      // Couple this product to each state that differs from it by one
      // component's transfer: neighbor_state[j] is that (n-1)-component state,
      // and the edge carries component j's transfer coupling.  For a 2-way
      // product this is exactly the two single-component edges.
      for (int j = 0; j < sites[k].n_components; j++) {
        fill_hamiltonian_edge(sites[k].neighbor_state[j], state_k, sites[sites[k].components[j]]);
      }
    } else {
      fill_hamiltonian_edge(sites[k].parent_state, state_k, sites[k]);
    }
  }

  // Broadcast position-based coupling values across samerank
  // (multi-rank partitions may have different contributing ranks)
  auto bcast_edge = [&](int parent, int child, const ReactiveSite &site) {
    const ReactionDef &rxn = rxndefs[site.rxn_idx];
    if (rxn.coupling_type != COUPLING_NONE && rxn.coupling_type != COUPLING_GRIMME2015) {
      MPI_Bcast(&hamiltonian[parent * ns + child], 1, MPI_DOUBLE, 0, samerank);
      hamiltonian[child * ns + parent] = hamiltonian[parent * ns + child];
    }
  };

  for (int k = 0; k < nsites_parallel; k++) {
    int state_k = k + 1;
    if (sites[k].n_components > 0) {
      for (int j = 0; j < sites[k].n_components; j++) {
        bcast_edge(sites[k].neighbor_state[j], state_k, sites[sites[k].components[j]]);
      }
    } else {
      bcast_edge(sites[k].parent_state, state_k, sites[k]);
    }
  }
}

/* ----------------------------------------------------------------------
   Solve the eigensystem of the nstates x nstates Hamiltonian.
---------------------------------------------------------------------- */

void FixMSEVB::solve_eigensystem()
{
  const int ns = nstates;

  for (int i = 0; i < ns; i++)
    for (int j = 0; j < ns; j++) H_work[i * ns + j] = hamiltonian[i * ns + j];

  for (int i = 0; i < ns; i++) {
    eigenvalues[i] = 0.0;
    amplitudes[i] = 0.0;
  }
  for (int i = 0; i < ns * ns; i++) eigenvectors[i] = 0.0;

  int info = JacobiEigen::solve(ns, H_work, eigenvalues, eigenvectors);

  if (info != 0) error->universe_all(FLERR, "Fix msevb: Jacobi eigensolver did not converge");

  epot_ground = eigenvalues[0];

  for (int i = 0; i < ns; i++) amplitudes[i] = eigenvectors[i * ns + 0] * eigenvectors[i * ns + 0];
}

/* ----------------------------------------------------------------------
   Compute per-state mixing weights and sparse coupling corrections.

   For the density matrix rho[i][j]:
     F[a] = sum_ij rho[i][j] * dH[i][j][a]

   Diagonal blocks: dH[k][k][a] = f_k[a]  (per-partition forces)
     → weights[k] = rho[k][k]

   Grimme off-diagonal: dH[p][c][a] = -g*(f_p[a] - f_c[a])
     → weights[p] -= 2*rho[p][c]*g,  weights[c] += 2*rho[p][c]*g

   Position-based off-diagonal: dH[p][c] is sparse (3 atoms only)
     → stored as CouplingCorrection for direct application
---------------------------------------------------------------------- */

void FixMSEVB::compute_mixing_weights()
{
  const int ns = nstates;

  // Build density matrix (unchanged)
  std::vector<double> rho(ns * ns, 0.0);

  if (fermi_dirac_enabled && ns > 1) {
    fermi_dirac_occupancies(eigenvalues, ns, fd_occ);

    epot_ground = 0.0;
    for (int k = 0; k < ns; k++) epot_ground += fd_occ[k] * eigenvalues[k];

    for (int k = 0; k < ns; k++) {
      if (fd_occ[k] < 1.0e-15) continue;
      for (int i = 0; i < ns; i++)
        for (int j = 0; j < ns; j++)
          rho[i * ns + j] += fd_occ[k] * eigenvectors[i * ns + k] * eigenvectors[j * ns + k];
    }

    for (int i = 0; i < ns; i++) amplitudes[i] = rho[i * ns + i];
  } else {
    for (int i = 0; i < ns; i++)
      for (int j = 0; j < ns; j++)
        rho[i * ns + j] = eigenvectors[i * ns + 0] * eigenvectors[j * ns + 0];
  }

  // Initialize weights from diagonal of density matrix
  for (int k = 0; k < ns; k++) weights[k] = rho[k * ns + k];

  // Coupling adjustments and corrections for all sites (parallel + serial).
  // Single pass replaces separate Grimme weight adjustment + position-based
  // correction loops.  evaluate_coupling() handles all dispatch.
  coupling_corrections.clear();

  auto process_edge = [&](int parent, int child, const ReactiveSite &site) {
    const ReactionDef &rxn = rxndefs[site.rxn_idx];
    if (rxn.coupling_type == COUPLING_NONE) return;

    double rho_pc = rho[parent * ns + child];

    auto cr =
        evaluate_coupling(rxn, hamiltonian[parent * ns + parent], hamiltonian[child * ns + child],
                          site.tag_H, site.tag_X, site.tag_Y);

    // Grimme weight adjustment: -2*rho_pc*g to parent, +2*rho_pc*g to child
    if (cr.is_grimme) {
      double w_adj = 2.0 * rho_pc * cr.g;
      weights[parent] -= w_adj;
      weights[child] += w_adj;
    }

    // Sparse force corrections (position-based coupling or Grimme taper)
    if (cr.has_forces && std::fabs(rho_pc) >= 1.0e-30) {
      CouplingCorrection cc;
      cc.tag_H = site.tag_H;
      cc.tag_X = site.tag_X;
      cc.tag_Y = site.tag_Y;
      cc.rho_factor = 2.0 * rho_pc;
      for (int d = 0; d < 3; d++) {
        cc.fH[d] = cr.fH[d];
        cc.fX[d] = cr.fX[d];
        cc.fY[d] = cr.fY[d];
      }
      coupling_corrections.push_back(cc);
    }
  };

  for (int k = 0; k < nsites_parallel + nsites_serial; k++) {
    int state_k = k + 1;
    if (sites[k].n_components > 0) {
      for (int ci = 0; ci < sites[k].n_components; ci++) {
        int par_i = sites[k].components[ci] + 1;
        int comp_j = sites[k].components[1 - ci];
        process_edge(par_i, state_k, sites[comp_j]);
      }
    } else {
      process_edge(sites[k].parent_state, state_k, sites[k]);
    }
  }
}

/* ----------------------------------------------------------------------
   Apply sparse position-based coupling corrections to atom->f and
   accumulate coupling_virial.  Called from weight_based_hellmann_feynman_forces().
---------------------------------------------------------------------- */

void FixMSEVB::apply_coupling_corrections()
{
  const int nlocal = atom->nlocal;
  double **f = atom->f;
  double **x = atom->x;

  for (int k = 0; k < 6; k++) coupling_virial[k] = 0.0;

  for (const auto &cc : coupling_corrections) {
    int idx_H = atom->map(cc.tag_H);
    int idx_X = atom->map(cc.tag_X);
    int idx_Y = atom->map(cc.tag_Y);

    if (idx_H >= 0 && idx_H < nlocal) {
      f[idx_H][0] += cc.rho_factor * cc.fH[0];
      f[idx_H][1] += cc.rho_factor * cc.fH[1];
      f[idx_H][2] += cc.rho_factor * cc.fH[2];
    }
    if (idx_X >= 0 && idx_X < nlocal) {
      f[idx_X][0] += cc.rho_factor * cc.fX[0];
      f[idx_X][1] += cc.rho_factor * cc.fX[1];
      f[idx_X][2] += cc.rho_factor * cc.fX[2];
    }
    if (idx_Y >= 0 && idx_Y < nlocal) {
      f[idx_Y][0] += cc.rho_factor * cc.fY[0];
      f[idx_Y][1] += cc.rho_factor * cc.fY[1];
      f[idx_Y][2] += cc.rho_factor * cc.fY[2];
    }

    // Accumulate virial: vir[ab] = dHX[a]*FH[b] + dYX[a]*FY[b]
    if (idx_H >= 0 && idx_H < nlocal && idx_X >= 0 && idx_Y >= 0) {
      double dHX[3] = {x[idx_H][0] - x[idx_X][0], x[idx_H][1] - x[idx_X][1],
                       x[idx_H][2] - x[idx_X][2]};
      double dYX[3] = {x[idx_Y][0] - x[idx_X][0], x[idx_Y][1] - x[idx_X][1],
                       x[idx_Y][2] - x[idx_X][2]};
      domain->minimum_image(FLERR, dHX[0], dHX[1], dHX[2]);
      domain->minimum_image(FLERR, dYX[0], dYX[1], dYX[2]);

      double FH[3] = {cc.rho_factor * cc.fH[0], cc.rho_factor * cc.fH[1], cc.rho_factor * cc.fH[2]};
      double FY[3] = {cc.rho_factor * cc.fY[0], cc.rho_factor * cc.fY[1], cc.rho_factor * cc.fY[2]};

      coupling_virial[0] += dHX[0] * FH[0] + dYX[0] * FY[0];
      coupling_virial[1] += dHX[1] * FH[1] + dYX[1] * FY[1];
      coupling_virial[2] += dHX[2] * FH[2] + dYX[2] * FY[2];
      coupling_virial[3] += dHX[0] * FH[1] + dYX[0] * FY[1];
      coupling_virial[4] += dHX[0] * FH[2] + dYX[0] * FY[2];
      coupling_virial[5] += dHX[1] * FH[2] + dYX[1] * FY[2];
    }
  }

  MPI_Allreduce(MPI_IN_PLACE, coupling_virial, 6, MPI_DOUBLE, MPI_SUM, world);
}

/* ----------------------------------------------------------------------
   Weight-based Hellmann-Feynman force mixing.

   1. Scale local forces by this partition's weight (or zero if excess).
   2. MPI_Allreduce forces across samerank communicator.
   3. Apply sparse position-based coupling corrections.

   Communication: single Allreduce of natoms*3 doubles (vs ns*natoms*3).
---------------------------------------------------------------------- */

void FixMSEVB::weight_based_hellmann_feynman_forces()
{
  sync_forces_to_host();
  const int nlocal = atom->nlocal;
  double **f = atom->f;

  // Scale this partition's forces by its weight
  double w = (ipartition < nstates) ? weights[ipartition] : 0.0;

  for (int a = 0; a < nlocal; a++) {
    f[a][0] *= w;
    f[a][1] *= w;
    f[a][2] *= w;
  }

  // Allreduce weighted forces across partitions (samerank comm)
  MPI_Allreduce(MPI_IN_PLACE, &f[0][0], 3 * nlocal, MPI_DOUBLE, MPI_SUM, samerank);

  // Add excess state diagonal force contributions from excess_forces buffer.
  // Only used by Approach A (compute_excess_states fills excess_forces).
  // Approach B applies excess forces directly in apply_excess_forces().
  if (nsites_serial > 0 && excess_forces) {
    const int ef_stride = excess_forces_nmax * 3;
    for (int b = 0; b < nsites_serial; b++) {
      int sk = nsites_parallel + b;
      int state_k = sk + 1;
      double wk = weights[state_k];
      if (std::fabs(wk) < 1.0e-30) continue;
      const int ef_off = b * ef_stride;
      for (int a = 0; a < nlocal; a++) {
        f[a][0] += wk * excess_forces[ef_off + a * 3 + 0];
        f[a][1] += wk * excess_forces[ef_off + a * 3 + 1];
        f[a][2] += wk * excess_forces[ef_off + a * 3 + 2];
      }
    }
  }

  apply_coupling_corrections();
  modified_forces_on_host();
}

/* ----------------------------------------------------------------------
   Fermi-Dirac occupancy solver.
---------------------------------------------------------------------- */

void FixMSEVB::fermi_dirac_occupancies(double *evals, int ns, double *occ)
{
  // Shift eigenvalues by evals[0] so all values are small (0 to range).
  // Occupancies depend only on differences (evals[i] - mu), so this shift
  // is exact but eliminates catastrophic cancellation when |evals| >> fd_RT.
  const double e0 = evals[0];
  double mu = (ns > 1) ? 0.5 * (evals[1] - e0) : 0.0;

  for (int iter = 0; iter < 150; iter++) {
    double sumq = 0.0, dsumq = 0.0;
    for (int i = 0; i < ns; i++) {
      double x = ((evals[i] - e0) - mu) / fd_RT;
      double f, df;
      if (x < -100.0) {
        f = 1.0;
        df = 0.0;
      } else if (x > 100.0) {
        f = 0.0;
        df = 0.0;
      } else {
        double ex = exp(x);
        f = 1.0 / (1.0 + ex);
        df = ex / (fd_RT * (ex + 1.0) * (ex + 1.0));
      }
      occ[i] = f;
      sumq += f;
      dsumq += df;
    }
    if (fabs(sumq - 1.0) < 1.0e-10) return;
    if (fabs(dsumq) > 1.0e-20) mu += (1.0 - sumq) / dsumq;
  }

  error->universe_all(FLERR, "Fix msevb: Fermi-Dirac occupancy solver did not converge");
}

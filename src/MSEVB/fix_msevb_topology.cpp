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
   fix msevb topology management:
   Reference snapshot/restore, reactive site detection, per-partition
   state changes, bond/angle/type modifications.
------------------------------------------------------------------------- */

#include "fix_msevb.h"
#include "fix_msevb_superimpose.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "kspace.h"
#include "memory.h"
#include "molecule.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "special.h"
#include "universe.h"
#include "update.h"

#include <algorithm>
#include <cstring>
#include <functional>
#include <map>

using namespace LAMMPS_NS;

/* ----------------------------------------------------------------------
   Ensure sites[] has capacity for at least n entries.
   Uses memory->destroy/create to stay consistent with LAMMPS allocator
   conventions. Content is NOT preserved across a reallocation — callers
   must refill from scratch, which is what detect_reactive_sites() does
   every timestep anyway.
---------------------------------------------------------------------- */

void FixMSEVB::grow_sites(int n)
{
  if (n <= sites_nmax) return;

  int newmax = MAX(n, MAX(2 * sites_nmax, 4));

  ReactiveSite *old_sites = sites;
  int old_nmax = sites_nmax;

  memory->create(sites, newmax, "msevb:sites");
  sites_nmax = newmax;

  if (old_sites && old_nmax > 0) {
    int ncopy = MIN(old_nmax, newmax);
    memcpy(sites, old_sites, ncopy * sizeof(ReactiveSite));
  }
  memory->destroy(old_sites);

  // Co-grow flat chain arrays (size = sites_nmax * max_shells).
  const int chain_size = newmax * max_shells;
  tagint *old_H = chain_H_flat, *old_X = chain_X_flat, *old_Y = chain_Y_flat;
  int *old_rxn = chain_rxn_flat;
  memory->create(chain_H_flat, chain_size, "msevb:chain_H_flat");
  memory->create(chain_X_flat, chain_size, "msevb:chain_X_flat");
  memory->create(chain_Y_flat, chain_size, "msevb:chain_Y_flat");
  memory->create(chain_rxn_flat, chain_size, "msevb:chain_rxn_flat");

  if (old_H && old_nmax > 0) {
    int ncopy = MIN(old_nmax, newmax) * max_shells;
    memcpy(chain_H_flat, old_H, ncopy * sizeof(tagint));
    memcpy(chain_X_flat, old_X, ncopy * sizeof(tagint));
    memcpy(chain_Y_flat, old_Y, ncopy * sizeof(tagint));
    memcpy(chain_rxn_flat, old_rxn, ncopy * sizeof(int));
  }
  memory->destroy(old_H);
  memory->destroy(old_X);
  memory->destroy(old_Y);
  memory->destroy(old_rxn);

  // Co-grow glove_flat and chain_glove_flat (size = sites_nmax * glove_nmax).
  // Only allocated when glove_nmax > 0.
  if (glove_nmax > 0) {
    tagint *old_glove = glove_flat;
    tagint *old_cglove = chain_glove_flat;
    const int gsize = newmax * glove_nmax;
    const int cgsize = newmax * max_shells * glove_nmax;
    memory->create(glove_flat, gsize, "msevb:glove_flat");
    memory->create(chain_glove_flat, cgsize, "msevb:chain_glove_flat");
    memset(glove_flat, 0, gsize * sizeof(tagint));
    memset(chain_glove_flat, 0, cgsize * sizeof(tagint));
    if (old_glove && old_nmax > 0) {
      int ncopy = MIN(old_nmax, newmax) * glove_nmax;
      memcpy(glove_flat, old_glove, ncopy * sizeof(tagint));
    }
    if (old_cglove && old_nmax > 0) {
      int ncopy = MIN(old_nmax, newmax) * max_shells * glove_nmax;
      memcpy(chain_glove_flat, old_cglove, ncopy * sizeof(tagint));
    }
    memory->destroy(old_glove);
    memory->destroy(old_cglove);
  }
}

/* ----------------------------------------------------------------------
   Initialise a site entry as a shell-1 site.
   Caller must set sites[site_idx].rxn_idx and
   chain_rxn_flat[site_idx * max_shells + 0] after this call.
---------------------------------------------------------------------- */

void FixMSEVB::init_site(int site_idx, tagint tag_H, tagint tag_X, tagint tag_Y, double dist_sq)
{
  ReactiveSite *s = &sites[site_idx];
  s->tag_H = tag_H;
  s->tag_X = tag_X;
  s->tag_Y = tag_Y;
  s->dist_sq = dist_sq;
  s->rxn_idx = 0;
  s->parent_state = 0;
  s->shell = 1;
  s->chain_len = 1;
  s->n_components = 0;
  for (int c = 0; c < MAX_COMPONENTS; c++) s->components[c] = -1;
  tagint *cH = chain_H_flat + site_idx * max_shells;
  tagint *cX = chain_X_flat + site_idx * max_shells;
  tagint *cY = chain_Y_flat + site_idx * max_shells;
  int *cR = chain_rxn_flat + site_idx * max_shells;
  cH[0] = tag_H;
  cX[0] = tag_X;
  cY[0] = tag_Y;
  cR[0] = 0;    // caller must overwrite with actual rxn_idx
  for (int d = 1; d < max_shells; d++) {
    cH[d] = 0;
    cX[d] = 0;
    cY[d] = 0;
    cR[d] = 0;
  }
  // Zero chain_glove_flat for this site (all depths).
  if (glove_nmax > 0)
    memset(chain_glove_flat + site_idx * max_shells * glove_nmax, 0,
           max_shells * glove_nmax * sizeof(tagint));
}

/* ----------------------------------------------------------------------
   check consistency of owned atom count and ordering across partitions
---------------------------------------------------------------------- */

void FixMSEVB::check_consistency_atoms()
{
  const int nlocal = atom->nlocal;
  int np = npartitions;

  for (int i = 0; i < np; i++) {
    my_nlocal[i] = 0;
    all_nlocal[i] = 0;
  }
  my_nlocal[ipartition] = nlocal;
  MPI_Allreduce(my_nlocal, all_nlocal, npartitions, MPI_INT, MPI_SUM, samerank);

  int fail = 0;
  for (int i = 1; i < npartitions; i++)
    if (all_nlocal[i] != all_nlocal[0]) fail = 1;

  int allfail = 0;
  MPI_Allreduce(&fail, &allfail, 1, MPI_INT, MPI_MAX, universe->uworld);
  if (allfail)
    error->universe_all(FLERR, "Fix msevb: local atom count is inconsistent across partitions");

  tagint *tagbuf = (tagint *) commbuf;
  tagint *tag = atom->tag;
  if (nlocal > 0) {
    if (ipartition == 0)
      for (int i = 0; i < nlocal; i++) tagbuf[i] = tag[i];
    MPI_Bcast(tagbuf, nlocal, MPI_LMP_TAGINT, 0, samerank);
  }

  fail = allfail = 0;
  if (ipartition > 0)
    for (int i = 0; i < nlocal; i++)
      if (tag[i] != tagbuf[i]) fail = 1;
  MPI_Allreduce(&fail, &allfail, 1, MPI_INT, MPI_MAX, universe->uworld);
  if (allfail)
    error->universe_all(FLERR, "Fix msevb: local atom ordering is inconsistent across partitions");
}

/* ----------------------------------------------------------------------
   Destroy all ref_* flat arrays (safe to call when pointers are null).
---------------------------------------------------------------------- */

void FixMSEVB::destroy_ref_topology()
{
  memory->destroy(ref_type);
  memory->destroy(ref_charge);
  memory->destroy(ref_molecule);
  memory->destroy(ref_num_bond);
  memory->destroy(ref_bond_type_flat);
  memory->destroy(ref_bond_atom_flat);
  memory->destroy(ref_num_angle);
  memory->destroy(ref_angle_type_flat);
  memory->destroy(ref_angle_atom1_flat);
  memory->destroy(ref_angle_atom2_flat);
  memory->destroy(ref_angle_atom3_flat);
  memory->destroy(ref_num_dihedral);
  memory->destroy(ref_dihedral_type_flat);
  memory->destroy(ref_dihedral_atom1_flat);
  memory->destroy(ref_dihedral_atom2_flat);
  memory->destroy(ref_dihedral_atom3_flat);
  memory->destroy(ref_dihedral_atom4_flat);
  memory->destroy(ref_num_improper);
  memory->destroy(ref_improper_type_flat);
  memory->destroy(ref_improper_atom1_flat);
  memory->destroy(ref_improper_atom2_flat);
  memory->destroy(ref_improper_atom3_flat);
  memory->destroy(ref_improper_atom4_flat);
  memory->destroy(ref_nspecial_flat);
  memory->destroy(ref_special_flat);
}

/* ----------------------------------------------------------------------
   Allocate all ref_* flat arrays for the given sizes.
   The caller must update ref_maxtag, ref_*_per_atom, ref_maxspecial
   after returning.
---------------------------------------------------------------------- */

void FixMSEVB::allocate_ref_topology(tagint sz, int bpa, int apa, int dpa, int ipa, int mspc)
{
  memory->create(ref_type, sz, "msevb:ref_type");
  memory->create(ref_charge, sz, "msevb:ref_charge");
  memory->create(ref_molecule, sz, "msevb:ref_molecule");
  memory->create(ref_num_bond, sz, "msevb:ref_num_bond");
  memory->create(ref_bond_type_flat, sz * bpa, "msevb:ref_bond_type_flat");
  memory->create(ref_bond_atom_flat, sz * bpa, "msevb:ref_bond_atom_flat");
  memory->create(ref_num_angle, sz, "msevb:ref_num_angle");
  memory->create(ref_angle_type_flat, sz * apa, "msevb:ref_angle_type_flat");
  memory->create(ref_angle_atom1_flat, sz * apa, "msevb:ref_angle_atom1_flat");
  memory->create(ref_angle_atom2_flat, sz * apa, "msevb:ref_angle_atom2_flat");
  memory->create(ref_angle_atom3_flat, sz * apa, "msevb:ref_angle_atom3_flat");
  memory->create(ref_num_dihedral, sz, "msevb:ref_num_dihedral");
  memory->create(ref_dihedral_type_flat, sz * dpa, "msevb:ref_dihedral_type_flat");
  memory->create(ref_dihedral_atom1_flat, sz * dpa, "msevb:ref_dihedral_atom1_flat");
  memory->create(ref_dihedral_atom2_flat, sz * dpa, "msevb:ref_dihedral_atom2_flat");
  memory->create(ref_dihedral_atom3_flat, sz * dpa, "msevb:ref_dihedral_atom3_flat");
  memory->create(ref_dihedral_atom4_flat, sz * dpa, "msevb:ref_dihedral_atom4_flat");
  memory->create(ref_num_improper, sz, "msevb:ref_num_improper");
  memory->create(ref_improper_type_flat, sz * ipa, "msevb:ref_improper_type_flat");
  memory->create(ref_improper_atom1_flat, sz * ipa, "msevb:ref_improper_atom1_flat");
  memory->create(ref_improper_atom2_flat, sz * ipa, "msevb:ref_improper_atom2_flat");
  memory->create(ref_improper_atom3_flat, sz * ipa, "msevb:ref_improper_atom3_flat");
  memory->create(ref_improper_atom4_flat, sz * ipa, "msevb:ref_improper_atom4_flat");
  memory->create(ref_nspecial_flat, sz * 3, "msevb:ref_nspecial_flat");
  memory->create(ref_special_flat, sz * mspc, "msevb:ref_special_flat");
}

/* ----------------------------------------------------------------------
   Snapshot the current topology into ref_* flat arrays.
---------------------------------------------------------------------- */

void FixMSEVB::snapshot_reference_topology()
{
  const int nlocal = atom->nlocal;
  const int bpa = atom->bond_per_atom;
  const int apa = atom->angle_per_atom;
  const int dpa = atom->dihedral_per_atom;
  const int ipa = atom->improper_per_atom;
  const int mspc = atom->maxspecial;
  const tagint maxtag = atom->map_tag_max;

  if (maxtag != ref_maxtag || bpa != ref_bond_per_atom || apa != ref_angle_per_atom ||
      dpa != ref_dihedral_per_atom || ipa != ref_improper_per_atom || mspc != ref_maxspecial) {
    destroy_ref_topology();
    allocate_ref_topology(maxtag + 1, bpa, apa, dpa, ipa, mspc);
    ref_maxtag = maxtag;
    ref_bond_per_atom = bpa;
    ref_angle_per_atom = apa;
    ref_dihedral_per_atom = dpa;
    ref_improper_per_atom = ipa;
    ref_maxspecial = mspc;
  }

  const tagint sz = maxtag + 1;
  memset(ref_type, 0, sizeof(int) * sz);
  memset(ref_charge, 0, sizeof(double) * sz);
  memset(ref_molecule, 0, sizeof(tagint) * sz);
  memset(ref_num_bond, 0, sizeof(int) * sz);
  memset(ref_bond_type_flat, 0, sizeof(int) * sz * bpa);
  memset(ref_bond_atom_flat, 0, sizeof(tagint) * sz * bpa);
  memset(ref_num_angle, 0, sizeof(int) * sz);
  memset(ref_angle_type_flat, 0, sizeof(int) * sz * apa);
  memset(ref_angle_atom1_flat, 0, sizeof(tagint) * sz * apa);
  memset(ref_angle_atom2_flat, 0, sizeof(tagint) * sz * apa);
  memset(ref_angle_atom3_flat, 0, sizeof(tagint) * sz * apa);
  if (dpa > 0) {
    memset(ref_num_dihedral, 0, sizeof(int) * sz);
    memset(ref_dihedral_type_flat, 0, sizeof(int) * sz * dpa);
    memset(ref_dihedral_atom1_flat, 0, sizeof(tagint) * sz * dpa);
    memset(ref_dihedral_atom2_flat, 0, sizeof(tagint) * sz * dpa);
    memset(ref_dihedral_atom3_flat, 0, sizeof(tagint) * sz * dpa);
    memset(ref_dihedral_atom4_flat, 0, sizeof(tagint) * sz * dpa);
  }
  if (ipa > 0) {
    memset(ref_num_improper, 0, sizeof(int) * sz);
    memset(ref_improper_type_flat, 0, sizeof(int) * sz * ipa);
    memset(ref_improper_atom1_flat, 0, sizeof(tagint) * sz * ipa);
    memset(ref_improper_atom2_flat, 0, sizeof(tagint) * sz * ipa);
    memset(ref_improper_atom3_flat, 0, sizeof(tagint) * sz * ipa);
    memset(ref_improper_atom4_flat, 0, sizeof(tagint) * sz * ipa);
  }
  memset(ref_nspecial_flat, 0, sizeof(int) * sz * 3);
  memset(ref_special_flat, 0, sizeof(tagint) * sz * mspc);

  tagint *tag = atom->tag;
  int *type = atom->type;
  double *q = atom->q;
  tagint *molecule = atom->molecule;
  int *num_bond = atom->num_bond;
  int **bond_type = atom->bond_type;
  tagint **bond_atom = atom->bond_atom;
  int *num_angle = atom->num_angle;
  int **angle_type = atom->angle_type;
  tagint **angle_atom1 = atom->angle_atom1;
  tagint **angle_atom2 = atom->angle_atom2;
  tagint **angle_atom3 = atom->angle_atom3;
  int *num_dihedral = atom->num_dihedral;
  int **dihedral_type = atom->dihedral_type;
  tagint **dihedral_atom1 = atom->dihedral_atom1;
  tagint **dihedral_atom2 = atom->dihedral_atom2;
  tagint **dihedral_atom3 = atom->dihedral_atom3;
  tagint **dihedral_atom4 = atom->dihedral_atom4;
  int *num_improper = atom->num_improper;
  int **improper_type = atom->improper_type;
  tagint **improper_atom1 = atom->improper_atom1;
  tagint **improper_atom2 = atom->improper_atom2;
  tagint **improper_atom3 = atom->improper_atom3;
  tagint **improper_atom4 = atom->improper_atom4;
  int **nspecial = atom->nspecial;
  tagint **special = atom->special;

  for (int i = 0; i < nlocal; i++) {
    const tagint t = tag[i];
    ref_type[t] = type[i];
    ref_charge[t] = q[i];
    ref_molecule[t] = molecule[i];
    ref_num_bond[t] = num_bond[i];
    for (int k = 0; k < bpa; k++) {
      ref_bond_type_flat[t * bpa + k] = bond_type[i][k];
      ref_bond_atom_flat[t * bpa + k] = bond_atom[i][k];
    }
    ref_num_angle[t] = num_angle[i];
    for (int k = 0; k < apa; k++) {
      ref_angle_type_flat[t * apa + k] = angle_type[i][k];
      ref_angle_atom1_flat[t * apa + k] = angle_atom1[i][k];
      ref_angle_atom2_flat[t * apa + k] = angle_atom2[i][k];
      ref_angle_atom3_flat[t * apa + k] = angle_atom3[i][k];
    }
    if (dpa > 0 && num_dihedral) {
      ref_num_dihedral[t] = num_dihedral[i];
      for (int k = 0; k < dpa; k++) {
        ref_dihedral_type_flat[t * dpa + k] = dihedral_type[i][k];
        ref_dihedral_atom1_flat[t * dpa + k] = dihedral_atom1[i][k];
        ref_dihedral_atom2_flat[t * dpa + k] = dihedral_atom2[i][k];
        ref_dihedral_atom3_flat[t * dpa + k] = dihedral_atom3[i][k];
        ref_dihedral_atom4_flat[t * dpa + k] = dihedral_atom4[i][k];
      }
    }
    if (ipa > 0 && num_improper) {
      ref_num_improper[t] = num_improper[i];
      for (int k = 0; k < ipa; k++) {
        ref_improper_type_flat[t * ipa + k] = improper_type[i][k];
        ref_improper_atom1_flat[t * ipa + k] = improper_atom1[i][k];
        ref_improper_atom2_flat[t * ipa + k] = improper_atom2[i][k];
        ref_improper_atom3_flat[t * ipa + k] = improper_atom3[i][k];
        ref_improper_atom4_flat[t * ipa + k] = improper_atom4[i][k];
      }
    }
    for (int k = 0; k < 3; k++) ref_nspecial_flat[t * 3 + k] = nspecial[i][k];
    int ns = nspecial[i][2];
    for (int k = 0; k < ns; k++) ref_special_flat[t * mspc + k] = special[i][k];
  }

  if (comm->nprocs > 1) {
    MPI_Allreduce(MPI_IN_PLACE, ref_type, sz, MPI_INT, MPI_SUM, world);
    MPI_Allreduce(MPI_IN_PLACE, ref_charge, sz, MPI_DOUBLE, MPI_SUM, world);
    MPI_Allreduce(MPI_IN_PLACE, ref_molecule, sz, MPI_LMP_TAGINT, MPI_SUM, world);
    MPI_Allreduce(MPI_IN_PLACE, ref_num_bond, sz, MPI_INT, MPI_SUM, world);
    MPI_Allreduce(MPI_IN_PLACE, ref_bond_type_flat, sz * bpa, MPI_INT, MPI_SUM, world);
    MPI_Allreduce(MPI_IN_PLACE, ref_bond_atom_flat, sz * bpa, MPI_LMP_TAGINT, MPI_SUM, world);
    MPI_Allreduce(MPI_IN_PLACE, ref_num_angle, sz, MPI_INT, MPI_SUM, world);
    MPI_Allreduce(MPI_IN_PLACE, ref_angle_type_flat, sz * apa, MPI_INT, MPI_SUM, world);
    MPI_Allreduce(MPI_IN_PLACE, ref_angle_atom1_flat, sz * apa, MPI_LMP_TAGINT, MPI_SUM, world);
    MPI_Allreduce(MPI_IN_PLACE, ref_angle_atom2_flat, sz * apa, MPI_LMP_TAGINT, MPI_SUM, world);
    MPI_Allreduce(MPI_IN_PLACE, ref_angle_atom3_flat, sz * apa, MPI_LMP_TAGINT, MPI_SUM, world);
    if (dpa > 0) {
      MPI_Allreduce(MPI_IN_PLACE, ref_num_dihedral, sz, MPI_INT, MPI_SUM, world);
      MPI_Allreduce(MPI_IN_PLACE, ref_dihedral_type_flat, sz * dpa, MPI_INT, MPI_SUM, world);
      MPI_Allreduce(MPI_IN_PLACE, ref_dihedral_atom1_flat, sz * dpa, MPI_LMP_TAGINT, MPI_SUM,
                    world);
      MPI_Allreduce(MPI_IN_PLACE, ref_dihedral_atom2_flat, sz * dpa, MPI_LMP_TAGINT, MPI_SUM,
                    world);
      MPI_Allreduce(MPI_IN_PLACE, ref_dihedral_atom3_flat, sz * dpa, MPI_LMP_TAGINT, MPI_SUM,
                    world);
      MPI_Allreduce(MPI_IN_PLACE, ref_dihedral_atom4_flat, sz * dpa, MPI_LMP_TAGINT, MPI_SUM,
                    world);
    }
    if (ipa > 0) {
      MPI_Allreduce(MPI_IN_PLACE, ref_num_improper, sz, MPI_INT, MPI_SUM, world);
      MPI_Allreduce(MPI_IN_PLACE, ref_improper_type_flat, sz * ipa, MPI_INT, MPI_SUM, world);
      MPI_Allreduce(MPI_IN_PLACE, ref_improper_atom1_flat, sz * ipa, MPI_LMP_TAGINT, MPI_SUM,
                    world);
      MPI_Allreduce(MPI_IN_PLACE, ref_improper_atom2_flat, sz * ipa, MPI_LMP_TAGINT, MPI_SUM,
                    world);
      MPI_Allreduce(MPI_IN_PLACE, ref_improper_atom3_flat, sz * ipa, MPI_LMP_TAGINT, MPI_SUM,
                    world);
      MPI_Allreduce(MPI_IN_PLACE, ref_improper_atom4_flat, sz * ipa, MPI_LMP_TAGINT, MPI_SUM,
                    world);
    }
    MPI_Allreduce(MPI_IN_PLACE, ref_nspecial_flat, sz * 3, MPI_INT, MPI_SUM, world);
    MPI_Allreduce(MPI_IN_PLACE, ref_special_flat, sz * mspc, MPI_LMP_TAGINT, MPI_SUM, world);
  }

  ref_snapshot_valid = 1;
}

/* ----------------------------------------------------------------------
   Restore topology from ref_* snapshot.
---------------------------------------------------------------------- */

void FixMSEVB::restore_reference_topology()
{
  const int nlocal = atom->nlocal;
  const int bpa = ref_bond_per_atom;
  const int apa = ref_angle_per_atom;
  const int dpa = ref_dihedral_per_atom;
  const int ipa = ref_improper_per_atom;
  const int mspc = ref_maxspecial;

  tagint *tag = atom->tag;
  int *type = atom->type;
  double *q = atom->q;
  tagint *molecule = atom->molecule;
  int *num_bond = atom->num_bond;
  int **bond_type = atom->bond_type;
  tagint **bond_atom = atom->bond_atom;
  int *num_angle = atom->num_angle;
  int **angle_type = atom->angle_type;
  tagint **angle_atom1 = atom->angle_atom1;
  tagint **angle_atom2 = atom->angle_atom2;
  tagint **angle_atom3 = atom->angle_atom3;
  int **nspecial = atom->nspecial;
  tagint **special_arr = atom->special;

  for (int i = 0; i < nlocal; i++) {
    const tagint t = tag[i];
    if (t < 0 || t > ref_maxtag) continue;
    type[i] = ref_type[t];
    q[i] = ref_charge[t];
    molecule[i] = ref_molecule[t];
    num_bond[i] = ref_num_bond[t];
    for (int k = 0; k < bpa; k++) {
      bond_type[i][k] = ref_bond_type_flat[t * bpa + k];
      bond_atom[i][k] = ref_bond_atom_flat[t * bpa + k];
    }
    num_angle[i] = ref_num_angle[t];
    for (int k = 0; k < apa; k++) {
      angle_type[i][k] = ref_angle_type_flat[t * apa + k];
      angle_atom1[i][k] = ref_angle_atom1_flat[t * apa + k];
      angle_atom2[i][k] = ref_angle_atom2_flat[t * apa + k];
      angle_atom3[i][k] = ref_angle_atom3_flat[t * apa + k];
    }
    if (dpa > 0 && atom->num_dihedral) {
      atom->num_dihedral[i] = ref_num_dihedral[t];
      for (int k = 0; k < dpa; k++) {
        atom->dihedral_type[i][k] = ref_dihedral_type_flat[t * dpa + k];
        atom->dihedral_atom1[i][k] = ref_dihedral_atom1_flat[t * dpa + k];
        atom->dihedral_atom2[i][k] = ref_dihedral_atom2_flat[t * dpa + k];
        atom->dihedral_atom3[i][k] = ref_dihedral_atom3_flat[t * dpa + k];
        atom->dihedral_atom4[i][k] = ref_dihedral_atom4_flat[t * dpa + k];
      }
    }
    if (ipa > 0 && atom->num_improper) {
      atom->num_improper[i] = ref_num_improper[t];
      for (int k = 0; k < ipa; k++) {
        atom->improper_type[i][k] = ref_improper_type_flat[t * ipa + k];
        atom->improper_atom1[i][k] = ref_improper_atom1_flat[t * ipa + k];
        atom->improper_atom2[i][k] = ref_improper_atom2_flat[t * ipa + k];
        atom->improper_atom3[i][k] = ref_improper_atom3_flat[t * ipa + k];
        atom->improper_atom4[i][k] = ref_improper_atom4_flat[t * ipa + k];
      }
    }
    for (int k = 0; k < 3; k++) nspecial[i][k] = ref_nspecial_flat[t * 3 + k];
    int ns_tot = nspecial[i][2];
    for (int k = 0; k < ns_tot; k++) special_arr[i][k] = ref_special_flat[t * mspc + k];
  }
}

/* ----------------------------------------------------------------------
   Detect reactive sites using occasional neighbor list.
   Run on partition 0 only, broadcast results to all partitions.

   All reactions are template-based (use msevb_superimpose BFS matching).
   sites[] is grown via grow_sites(); nsites = total sites detected.
---------------------------------------------------------------------- */

void FixMSEVB::detect_reactive_sites()
{
  nsites = 0;
  nstates = 1;

  if (ipartition == 0) {
    if (!list) goto broadcast;

    const int nlocal = atom->nlocal;
    int *type = atom->type;
    double **x = atom->x;
    tagint *tag = atom->tag;
    int *ilist = list->ilist;
    int *numneigh = list->numneigh;
    int **firstneigh = list->firstneigh;
    int inum = list->inum;

    // ---- Template site detection ------------------------------------
    // Stride per candidate in the flat double buffer:
    //   [0]          dist_sq
    //   [1]          tag_H  (as double)
    //   [2]          tag_X  (as double; 0 if ix_bonding < 0)
    //   [3]          tag_Y  (as double)
    //   [4]          rxn_idx (as double)
    //   [5..4+GN]   glove[0..glove_nmax-1] (as doubles)
    const int GN = glove_nmax;    // 0 if no reactions (shouldn't happen)
    const int TSTRIDE = 5 + GN;

    tagint *glove_tmp = new tagint[MAX(GN, 1)]();

    static const int MAX_TMPL_CAND = 256;
    double *tmpl_local = new double[MAX_TMPL_CAND * TSTRIDE]();
    int n_tmpl_local = 0;

    for (int ri = 0; ri < (int) rxndefs.size(); ri++) {
      const ReactionDef &rxn = rxndefs[ri];
      if (rxn.glove_n == 0) continue;

      for (int ii = 0; ii < inum; ii++) {
        int ih = ilist[ii];
        if (ih >= nlocal) continue;
        if (type[ih] != rxn.type_H) continue;

        int *jlistH = firstneigh[ih];
        int jnumH = numneigh[ih];
        for (int jj = 0; jj < jnumH; jj++) {
          int j = jlistH[jj] & NEIGHMASK;
          if (type[j] != rxn.type_Y) continue;

          double dx = x[ih][0] - x[j][0];
          double dy = x[ih][1] - x[j][1];
          double dz = x[ih][2] - x[j][2];
          domain->minimum_image(FLERR, dx, dy, dz);
          double rsq = dx * dx + dy * dy + dz * dz;
          if (rsq >= rxn.cutoff_sq) continue;

          // Run superimpose: BFS from (ih→ibonding, j→jbonding).
          // Pass ref topology (special list = symmetric 1-2 neighbors) so
          // ghost/off-rank atoms can be matched in multi-rank partitions.
          memset(glove_tmp, 0, GN * sizeof(tagint));
          RefTopo ref_topo{ref_type, ref_nspecial_flat, ref_special_flat, ref_maxspecial,
                           ref_maxtag};
          bool ok = msevb_superimpose(lmp, rxn.pre_mol, rxn.is_edge.data(), rxn.ibonding,
                                      rxn.jbonding, tag[ih], tag[j], glove_tmp, nullptr, &ref_topo);
          if (!ok) continue;

          if (n_tmpl_local >= MAX_TMPL_CAND) continue;

          tagint tX = (rxn.ix_bonding >= 0) ? glove_tmp[rxn.ix_bonding] : 0;

          int b = n_tmpl_local * TSTRIDE;
          tmpl_local[b + 0] = rsq;
          tmpl_local[b + 1] = static_cast<double>(tag[ih]);
          tmpl_local[b + 2] = static_cast<double>(tX);
          tmpl_local[b + 3] = static_cast<double>(tag[j]);
          tmpl_local[b + 4] = static_cast<double>(ri);
          for (int g = 0; g < GN; g++)
            tmpl_local[b + 5 + g] = static_cast<double>((g < rxn.glove_n) ? glove_tmp[g] : 0);
          n_tmpl_local++;
        }
      }
    }

    delete[] glove_tmp;

    // ---- Gather candidates from all ranks in partition 0 ------------
    int nranks_p0;
    MPI_Comm_size(world, &nranks_p0);
    int *tl_nc = new int[nranks_p0];
    MPI_Allgather(&n_tmpl_local, 1, MPI_INT, tl_nc, 1, MPI_INT, world);

    int total_tmpl = 0;
    int *tl_sc = new int[nranks_p0];
    int *tl_rd = new int[nranks_p0];
    for (int r = 0; r < nranks_p0; r++) {
      tl_rd[r] = total_tmpl * TSTRIDE;
      tl_sc[r] = tl_nc[r] * TSTRIDE;
      total_tmpl += tl_nc[r];
    }

    double *tmpl_all = new double[MAX(total_tmpl, 1) * TSTRIDE]();
    MPI_Allgatherv(tmpl_local, n_tmpl_local * TSTRIDE, MPI_DOUBLE, tmpl_all, tl_sc, tl_rd,
                   MPI_DOUBLE, world);

    delete[] tmpl_local;
    delete[] tl_nc;
    delete[] tl_sc;
    delete[] tl_rd;

    // ---- Insertion-sort by dist_sq -----------------------------------
    // Maximum TSTRIDE is 5+GN; GN typically ~7 for H3O+, bounded by 128
    // for safety in the sort buffer.
    const int SORT_STRIDE = TSTRIDE;
    std::vector<double> tmp_e(SORT_STRIDE);
    for (int i = 1; i < total_tmpl; i++) {
      for (int k = 0; k < SORT_STRIDE; k++) tmp_e[k] = tmpl_all[i * TSTRIDE + k];
      int j = i - 1;
      while (j >= 0 && tmpl_all[j * TSTRIDE] > tmp_e[0]) {
        for (int k = 0; k < SORT_STRIDE; k++)
          tmpl_all[(j + 1) * TSTRIDE + k] = tmpl_all[j * TSTRIDE + k];
        j--;
      }
      for (int k = 0; k < SORT_STRIDE; k++) tmpl_all[(j + 1) * TSTRIDE + k] = tmp_e[k];
    }

    // ---- Pick sites: unique (donor H, acceptor Y) transfers, sorted by
    // distance.  A state is defined by which proton moves to which acceptor, so
    // the identity of a transfer is the (tag_H, tag_Y) pair -- NOT tag_Y alone.
    // Deduping on tag_Y only would collapse distinct states in which several
    // donors hand a proton to the same acceptor (e.g. two waters donating to
    // one hydroxide, or a bicarbonate proton and a water proton both targeting
    // the same OH), keeping just the closest and dropping the rest.
    grow_sites(MAX(total_tmpl, 1));
    nsites = 0;
    for (int i = 0; i < total_tmpl; i++) {
      int b = i * TSTRIDE;
      tagint tY = static_cast<tagint>(tmpl_all[b + 3]);
      tagint tH = static_cast<tagint>(tmpl_all[b + 1]);
      int dup = 0;
      for (int k = 0; k < nsites; k++)
        if (sites[k].tag_Y == tY && sites[k].tag_H == tH) {
          dup = 1;
          break;
        }
      if (dup) continue;

      tagint tX = static_cast<tagint>(tmpl_all[b + 2]);
      double rsq = tmpl_all[b + 0];
      int ri = static_cast<int>(tmpl_all[b + 4]);

      grow_sites(nsites + 1);
      init_site(nsites, tH, tX, tY, rsq);
      sites[nsites].rxn_idx = ri;
      chain_rxn_flat[nsites * max_shells + 0] = ri;
      // Copy glove into glove_flat and chain_glove_flat depth 0.
      tagint *gptr = glove_flat + nsites * GN;
      tagint *cgptr = chain_glove_flat + nsites * max_shells * GN;
      for (int g = 0; g < GN; g++) {
        gptr[g] = static_cast<tagint>(tmpl_all[b + 5 + g]);
        cgptr[g] = gptr[g];
      }
      nsites++;
    }

    delete[] tmpl_all;

    // ---- Multi-shell + multi-site notice ----------------------------
    // A multi-shell reaction acting on several distinct donor molecules is
    // allowed: each donor's chain is detected independently and enters the
    // Hamiltonian as its own state, and product-state combination is now gated
    // by an exact atom-set disjointness test (chains_disjoint()), so disjoint
    // chains combine safely and overlapping ones are simply not combined.  We
    // therefore no longer error on this combination.  We do emit a one-time
    // note when product_states is enabled, because combining multi-shell chains
    // into product states has not yet been validated against a reference
    // implementation.
    //
    // "Distinct donor molecule" is determined by glove membership: two shell-1
    // sites belong to the same donor molecule if one site's tag_X appears in
    // the other's glove.  This correctly handles molecules with multiple
    // donatable atoms on different heavy atoms (e.g. imidazolium has two N-H
    // bonds, giving different tag_X values that still belong to the same
    // molecule), while still catching genuinely separate donor molecules.
    if (max_shells > 1 && nsites > 1 && enumerate_product_states &&
        !product_states_multishell_warning) {
      for (size_t ri = 0; ri < rxndefs.size(); ri++) {
        if (rxndefs[ri].shells <= 1) continue;
        // Count distinct donor molecules among this reaction's shell-1 sites.
        // Grown without bound so the count is exact regardless of how many
        // donors are present.
        std::vector<tagint> donor_X;
        std::vector<const tagint *> donor_glove;
        for (int s = 0; s < nsites; s++) {
          if (sites[s].rxn_idx != (int) ri || sites[s].shell != 1) continue;
          tagint tX = sites[s].tag_X;
          const tagint *sg = glove_flat + s * GN;
          int same_mol = 0;
          for (size_t d = 0; d < donor_X.size() && !same_mol; d++) {
            if (donor_X[d] == tX) {
              same_mol = 1;
              break;
            }
            // tX in the glove of donor d -> same molecule
            for (int gi = 0; gi < GN && !same_mol; gi++)
              if (donor_glove[d][gi] == tX) same_mol = 1;
            // donor_X[d] in this site's glove -> same molecule
            for (int gi = 0; gi < GN && !same_mol; gi++)
              if (sg[gi] == donor_X[d]) same_mol = 1;
          }
          if (!same_mol) {
            donor_X.push_back(tX);
            donor_glove.push_back(sg);
          }
        }
        int n_distinct_donors = (int) donor_X.size();
        if (n_distinct_donors > 1 && universe->me == 0) {
          utils::logmesg(lmp,
                         fmt::format("WARNING: Fix msevb: reaction {} has shells={} and "
                                     "{} distinct donor molecules.  Combining multi-shell "
                                     "chains with product_states is experimental and has "
                                     "not been validated against a reference implementation."
                                     "  Proceed with caution.\n",
                                     ri, rxndefs[ri].shells, n_distinct_donors));
          product_states_multishell_warning = true;
          break;    // one notice is enough
        }
      }
    }

    // ---- Multi-shell detection (shell 2+) ---------------------------
    // Uses the type changes from the reaction template to build a VirtualTopo
    // overlay that represents the post-transfer bonding state, then searches
    // for further transfer candidates from the new acceptor.
    for (int shell = 2; shell <= max_shells; shell++) {
      int prev_nsites = nsites;
      for (int s = 0; s < prev_nsites; s++) {
        if (sites[s].shell != shell - 1) continue;
        if (sites[s].chain_len >= rxndefs[sites[s].rxn_idx].shells) continue;

        int s_rxn_idx = sites[s].rxn_idx;
        const ReactionDef &prxn = rxndefs[s_rxn_idx];
        const tagint *s_glove = glove_flat + s * GN;

        // Build VirtualTopo for the post-transfer state of site s.
        // This allows superimpose to see the correct types/bonds after
        // transfer.  For shell >= 3, we must include type changes from
        // ALL ancestor transfers in the chain (earliest first) so that
        // superimpose sees the cumulative post-state types.
        VirtualTopo vtopo;

        // Walk the ancestor chain and apply type changes from each
        // ancestor's transfer (oldest first).  This is needed because
        // atoms outside the current site's glove may have been retyped
        // by earlier transfers in the chain.
        {
          // Collect ancestor sites from root to parent.
          std::vector<int> ancestors;
          int cur = s;
          while (cur >= 0 && sites[cur].shell > 1) {
            // Find the parent site index.
            int parent_state = sites[cur].parent_state;
            int parent_site = parent_state - 1;    // state N -> site N-1
            ancestors.push_back(parent_site);
            cur = parent_site;
          }
          // Apply from oldest ancestor to most recent (reverse order).
          for (int ai = (int) ancestors.size() - 1; ai >= 0; ai--) {
            int anc = ancestors[ai];
            const ReactionDef &arxn = rxndefs[sites[anc].rxn_idx];
            const tagint *a_glove = glove_flat + anc * GN;
            for (const auto &tc : arxn.type_changes) {
              if (tc.pre_idx < 0 || tc.pre_idx >= arxn.glove_n) continue;
              tagint r = a_glove[tc.pre_idx];
              if (r == 0) continue;
              vtopo.type_override[r] = tc.new_type;
            }
          }
        }

        // Apply type changes from this site's transfer (overrides
        // any ancestor changes on the same atoms).
        for (const auto &tc : prxn.type_changes) {
          if (tc.pre_idx < 0 || tc.pre_idx >= prxn.glove_n) continue;
          tagint r = s_glove[tc.pre_idx];
          if (r == 0) continue;
          vtopo.type_override[r] = tc.new_type;
        }

        // Helper: lazily initialise bond_override from ref topology.
        // Uses the symmetric special list (ref_nspecial_flat /
        // ref_special_flat) rather than the one-sided bond_atom list, because
        // LAMMPS stores each bond on only one atom — the acceptor Y may not
        // store its own bonds.
        auto get_vtopo_bonds = [&](tagint t_tag) -> std::vector<tagint> & {
          if (!vtopo.bond_override.count(t_tag)) {
            std::vector<tagint> bonds;
            if (t_tag > 0 && t_tag <= ref_maxtag) {
              int nb = ref_nspecial_flat[t_tag * 3 + 0];    // 1-2 neighbor count
              for (int b = 0; b < nb; b++) {
                tagint ba = ref_special_flat[(tagint) t_tag * ref_maxspecial + b];
                if (ba != 0) bonds.push_back(ba);
              }
            }
            vtopo.bond_override[t_tag] = bonds;
          }
          return vtopo.bond_override[t_tag];
        };

        // Apply bond breaks and creates from the template diff.
        for (const auto &bb : prxn.bond_breaks) {
          if (bb.pre_idx1 < 0 || bb.pre_idx1 >= prxn.glove_n) continue;
          if (bb.pre_idx2 < 0 || bb.pre_idx2 >= prxn.glove_n) continue;
          tagint r1 = s_glove[bb.pre_idx1];
          tagint r2 = s_glove[bb.pre_idx2];
          if (r1 == 0 || r2 == 0) continue;
          auto &b1 = get_vtopo_bonds(r1);
          b1.erase(std::remove(b1.begin(), b1.end(), r2), b1.end());
          auto &b2 = get_vtopo_bonds(r2);
          b2.erase(std::remove(b2.begin(), b2.end(), r1), b2.end());
        }
        for (const auto &bc : prxn.bond_creates) {
          if (bc.pre_idx1 < 0 || bc.pre_idx1 >= prxn.glove_n) continue;
          if (bc.pre_idx2 < 0 || bc.pre_idx2 >= prxn.glove_n) continue;
          tagint r1 = s_glove[bc.pre_idx1];
          tagint r2 = s_glove[bc.pre_idx2];
          if (r1 == 0 || r2 == 0) continue;
          auto &b1 = get_vtopo_bonds(r1);
          if (std::find(b1.begin(), b1.end(), r2) == b1.end()) b1.push_back(r2);
          auto &b2 = get_vtopo_bonds(r2);
          if (std::find(b2.begin(), b2.end(), r1) == b2.end()) b2.push_back(r1);
        }

        // The new donor X is the acceptor Y of the parent site.
        tagint new_X_tag = sites[s].tag_Y;

        // Gather bonded neighbors of new_X in vtopo.
        const std::vector<tagint> &new_X_bonds = get_vtopo_bonds(new_X_tag);

        for (tagint bp_tag : new_X_bonds) {
          // Determine type of bp_tag in vtopo.
          int bp_type;
          {
            auto vit = vtopo.type_override.find(bp_tag);
            if (vit != vtopo.type_override.end())
              bp_type = vit->second;
            else {
              int bp_local = atom->map(bp_tag);
              if (bp_local < 0) continue;
              bp_type = type[bp_local];
            }
          }
          // bp_tag must not be the parent H (already transferred).
          if (bp_tag == sites[s].tag_H) continue;

          // Look for acceptors Y near bp_tag in neighbor list.
          int bp_local = atom->map(bp_tag);
          if (bp_local < 0 || bp_local >= nlocal) continue;

          int *jlist_bp = firstneigh[bp_local];
          int jnum_bp = numneigh[bp_local];
          for (int jj = 0; jj < jnum_bp; jj++) {
            int j = jlist_bp[jj] & NEIGHMASK;

            // Determine type in vtopo.
            int j_type;
            {
              auto vit2 = vtopo.type_override.find(tag[j]);
              if (vit2 != vtopo.type_override.end())
                j_type = vit2->second;
              else
                j_type = type[j];
            }

            tagint j_tag = tag[j];

            // Not an ancestor.
            int is_ancestor = 0;
            for (int d = 0; d < sites[s].chain_len; d++)
              if (j_tag == chain_X_flat[s * max_shells + d] ||
                  j_tag == chain_Y_flat[s * max_shells + d]) {
                is_ancestor = 1;
                break;
              }
            if (is_ancestor) continue;

            // Try all reactions for this (H, Y) candidate.
            for (int ri = 0; ri < (int) rxndefs.size(); ri++) {
              const ReactionDef &crxn = rxndefs[ri];
              if (bp_type != crxn.type_H) continue;
              if (j_type != crxn.type_Y) continue;

              double dx = x[bp_local][0] - x[j][0];
              double dy = x[bp_local][1] - x[j][1];
              double dz = x[bp_local][2] - x[j][2];
              domain->minimum_image(FLERR, dx, dy, dz);
              double rsq = dx * dx + dy * dy + dz * dz;
              if (rsq >= crxn.cutoff_sq) continue;

              // Skip duplicate (donor H, acceptor Y) transfers within this
              // shell.  As in the shell-1 selection, the transfer identity is
              // the (H, Y) pair, so distinct donors to the same acceptor are
              // kept as separate states.
              int dup = 0;
              for (int k = prev_nsites; k < nsites; k++)
                if (sites[k].tag_Y == j_tag && sites[k].tag_H == bp_tag) {
                  dup = 1;
                  break;
                }
              if (dup) continue;

              // Run superimpose with vtopo (virtual bonds) and ref (fallback
              // for ghost/off-rank atoms in multi-rank partitions).
              std::vector<tagint> new_glove(GN, 0);
              RefTopo ref_topo2{ref_type, ref_nspecial_flat, ref_special_flat, ref_maxspecial,
                                ref_maxtag};
              bool ok = msevb_superimpose(lmp, crxn.pre_mol, crxn.is_edge.data(), crxn.ibonding,
                                          crxn.jbonding, bp_tag, j_tag, new_glove.data(), &vtopo,
                                          &ref_topo2);
              if (!ok) continue;

              // Snapshot parent chain before grow_sites() may reallocate.
              int s_chain_len = sites[s].chain_len;
              std::vector<tagint> s_chain_H(s_chain_len), s_chain_X(s_chain_len),
                  s_chain_Y(s_chain_len);
              std::vector<int> s_chain_rxn(s_chain_len);
              std::vector<tagint> s_chain_glove(s_chain_len * GN);
              for (int d = 0; d < s_chain_len; d++) {
                s_chain_H[d] = chain_H_flat[s * max_shells + d];
                s_chain_X[d] = chain_X_flat[s * max_shells + d];
                s_chain_Y[d] = chain_Y_flat[s * max_shells + d];
                s_chain_rxn[d] = chain_rxn_flat[s * max_shells + d];
                if (GN > 0)
                  memcpy(&s_chain_glove[d * GN],
                         chain_glove_flat + (s * max_shells + d) * glove_nmax, GN * sizeof(tagint));
              }

              grow_sites(nsites + 1);
              ReactiveSite *ns_ptr = &sites[nsites];
              tagint tX_new = (crxn.ix_bonding >= 0) ? new_glove[crxn.ix_bonding] : 0;
              ns_ptr->tag_H = bp_tag;
              ns_ptr->tag_X = tX_new;
              ns_ptr->tag_Y = j_tag;
              ns_ptr->dist_sq = rsq;
              ns_ptr->rxn_idx = ri;
              ns_ptr->parent_state = s + 1;
              ns_ptr->shell = shell;
              ns_ptr->chain_len = s_chain_len + 1;
              ns_ptr->n_components = 0;
              for (int c = 0; c < MAX_COMPONENTS; c++) ns_ptr->components[c] = -1;

              tagint *nH = chain_H_flat + nsites * max_shells;
              tagint *nX = chain_X_flat + nsites * max_shells;
              tagint *nY = chain_Y_flat + nsites * max_shells;
              int *nR = chain_rxn_flat + nsites * max_shells;
              for (int d = 0; d < s_chain_len; d++) {
                nH[d] = s_chain_H[d];
                nX[d] = s_chain_X[d];
                nY[d] = s_chain_Y[d];
                nR[d] = s_chain_rxn[d];
              }
              nH[s_chain_len] = bp_tag;
              nX[s_chain_len] = tX_new;
              nY[s_chain_len] = j_tag;
              nR[s_chain_len] = ri;
              for (int d = s_chain_len + 1; d < max_shells; d++) {
                nH[d] = 0;
                nX[d] = 0;
                nY[d] = 0;
                nR[d] = 0;
              }

              // Store glove for this shell-N site + chain_glove for all depths.
              if (GN > 0) {
                tagint *gptr = glove_flat + nsites * GN;
                for (int g = 0; g < GN; g++) gptr[g] = new_glove[g];
                // Copy parent chain glove entries for depths 0..chain_len-2
                for (int d = 0; d < s_chain_len; d++)
                  memcpy(chain_glove_flat + (nsites * max_shells + d) * glove_nmax,
                         &s_chain_glove[d * GN], GN * sizeof(tagint));
                // Store new glove at depth chain_len-1 (= s_chain_len)
                memcpy(chain_glove_flat + (nsites * max_shells + s_chain_len) * glove_nmax,
                       new_glove.data(), GN * sizeof(tagint));
                // Zero remaining depths
                for (int d = s_chain_len + 1; d < max_shells; d++)
                  memset(chain_glove_flat + (nsites * max_shells + d) * glove_nmax, 0,
                         glove_nmax * sizeof(tagint));
              }

              nsites++;
            }
          }
        }
      }

      // ---- Gather shell-N sites from all ranks within partition 0 -----
      if (comm->nprocs > 1) {
        // SFIELD: 8 header + 3*MSD chain_H/X/Y + MSD chain_rxn + GN glove
        //         + MSD*GN chain_glove
        const int SFIELD = 8 + 4 * max_shells + GN + max_shells * GN;
        int my_new = nsites - prev_nsites;

        int nr = comm->nprocs;
        int *all_new = new int[nr];
        MPI_Allgather(&my_new, 1, MPI_INT, all_new, 1, MPI_INT, world);

        int total_new = 0;
        int *scounts = new int[nr];
        int *rdispls = new int[nr];
        for (int r = 0; r < nr; r++) {
          rdispls[r] = total_new * SFIELD;
          scounts[r] = all_new[r] * SFIELD;
          total_new += all_new[r];
        }

        double *sendbuf = new double[MAX(my_new, 1) * SFIELD]();
        for (int k = 0; k < my_new; k++) {
          const ReactiveSite &sv = sites[prev_nsites + k];
          int b = k * SFIELD;
          sendbuf[b + 0] = sv.dist_sq;
          sendbuf[b + 1] = (double) sv.tag_H;
          sendbuf[b + 2] = (double) sv.tag_X;
          sendbuf[b + 3] = (double) sv.tag_Y;
          sendbuf[b + 4] = (double) sv.parent_state;
          sendbuf[b + 5] = (double) sv.shell;
          sendbuf[b + 6] = (double) sv.chain_len;
          sendbuf[b + 7] = (double) sv.rxn_idx;
          const tagint *svH = chain_H_flat + (prev_nsites + k) * max_shells;
          const tagint *svX = chain_X_flat + (prev_nsites + k) * max_shells;
          const tagint *svY = chain_Y_flat + (prev_nsites + k) * max_shells;
          const int *svR = chain_rxn_flat + (prev_nsites + k) * max_shells;
          for (int d = 0; d < max_shells; d++) {
            sendbuf[b + 8 + d] = (double) svH[d];
            sendbuf[b + 8 + max_shells + d] = (double) svX[d];
            sendbuf[b + 8 + 2 * max_shells + d] = (double) svY[d];
            sendbuf[b + 8 + 3 * max_shells + d] = (double) svR[d];
          }
          if (GN > 0) {
            const tagint *gv = glove_flat + (prev_nsites + k) * GN;
            for (int g = 0; g < GN; g++) sendbuf[b + 8 + 4 * max_shells + g] = (double) gv[g];
          }
          if (GN > 0) {
            for (int d = 0; d < max_shells; d++) {
              const tagint *cg =
                  chain_glove_flat + ((prev_nsites + k) * max_shells + d) * glove_nmax;
              for (int g = 0; g < GN; g++)
                sendbuf[b + 8 + 4 * max_shells + GN + d * GN + g] = (double) cg[g];
            }
          }
        }

        double *recvbuf = new double[MAX(total_new, 1) * SFIELD]();
        MPI_Allgatherv(sendbuf, my_new * SFIELD, MPI_DOUBLE, recvbuf, scounts, rdispls, MPI_DOUBLE,
                       world);

        nsites = prev_nsites;
        for (int i = 0; i < total_new; i++) {
          int b = i * SFIELD;
          tagint tY = (tagint) recvbuf[b + 3];
          tagint tH = (tagint) recvbuf[b + 1];
          // Dedup on the (donor H, acceptor Y) transfer identity, matching the
          // shell-1 selection: distinct donors to the same acceptor are distinct
          // states.  (Keying on tag_Y alone would re-collapse them here in
          // multi-rank runs, undoing the shell-1/within-rank dedup fix.)
          int dup = 0;
          for (int k = prev_nsites; k < nsites; k++)
            if (sites[k].tag_Y == tY && sites[k].tag_H == tH) {
              dup = 1;
              break;
            }
          if (dup) continue;
          grow_sites(nsites + 1);
          ReactiveSite &ns = sites[nsites];
          ns.dist_sq = recvbuf[b + 0];
          ns.tag_H = (tagint) recvbuf[b + 1];
          ns.tag_X = (tagint) recvbuf[b + 2];
          ns.tag_Y = tY;
          ns.parent_state = (int) recvbuf[b + 4];
          ns.shell = (int) recvbuf[b + 5];
          ns.chain_len = (int) recvbuf[b + 6];
          ns.rxn_idx = (int) recvbuf[b + 7];
          ns.n_components = 0;
          for (int c = 0; c < MAX_COMPONENTS; c++) ns.components[c] = -1;
          tagint *nH = chain_H_flat + nsites * max_shells;
          tagint *nX = chain_X_flat + nsites * max_shells;
          tagint *nY = chain_Y_flat + nsites * max_shells;
          int *nR = chain_rxn_flat + nsites * max_shells;
          for (int d = 0; d < max_shells; d++) {
            nH[d] = (tagint) recvbuf[b + 8 + d];
            nX[d] = (tagint) recvbuf[b + 8 + max_shells + d];
            nY[d] = (tagint) recvbuf[b + 8 + 2 * max_shells + d];
            nR[d] = (int) recvbuf[b + 8 + 3 * max_shells + d];
          }
          if (GN > 0) {
            tagint *gptr = glove_flat + nsites * GN;
            for (int g = 0; g < GN; g++) gptr[g] = (tagint) recvbuf[b + 8 + 4 * max_shells + g];
          }
          if (GN > 0) {
            for (int d = 0; d < max_shells; d++) {
              tagint *cg = chain_glove_flat + (nsites * max_shells + d) * glove_nmax;
              for (int g = 0; g < GN; g++)
                cg[g] = (tagint) recvbuf[b + 8 + 4 * max_shells + GN + d * GN + g];
            }
          }
          nsites++;
        }

        delete[] all_new;
        delete[] scounts;
        delete[] rdispls;
        delete[] sendbuf;
        delete[] recvbuf;
      }
    }
    // ---- Product state enumeration ----------------------------------
    // Enumerate all compatible pairs of chain-tip sites and create
    // combination states that apply both transfer chains simultaneously.
    // Works for any shell depth: components[ci] is the index of a site
    // whose full chain lives in chain_*_flat[comp * max_shells + d].
    // Compatibility: the atom sets of the two chains must be fully
    // disjoint — no shared H, X, or Y tag across any depth of either chain.
    if (enumerate_product_states) {
      // Enumerate every set of >=2 chain-tip sites whose atom sets are pairwise
      // disjoint (chains_disjoint) and turn each into a product state that
      // applies all of those chains simultaneously.  Two sites of the same
      // reactive complex overlap, so at most one per complex is ever chosen;
      // sites of different complexes are disjoint and combine.  This yields the
      // full n-way product space (the Cartesian product across complexes, minus
      // overlapping combinations), with the order emerging from the data rather
      // than a fixed keyword.  Order is capped only by MAX_COMPONENTS (storage).
      const int n_tip = nsites;

      // Collect the disjoint combinations first: creating the product sites
      // below may reallocate the flat arrays that chains_disjoint() reads.
      std::vector<std::vector<int>> combos;
      std::vector<int> cur;
      std::function<void(int)> enumerate = [&](int start) {
        for (int j = start; j < n_tip; j++) {
          if (sites[j].n_components != 0) continue;
          bool ok = true;
          for (int c : cur)
            if (!chains_disjoint(c, j)) {
              ok = false;
              break;
            }
          if (!ok) continue;
          cur.push_back(j);
          if ((int) cur.size() >= 2) combos.push_back(cur);
          if ((int) cur.size() < MAX_COMPONENTS) enumerate(j + 1);
          cur.pop_back();
        }
      };
      enumerate(0);

      for (const auto &combo : combos) {
        const int a = combo[0];
        grow_sites(nsites + 1);
        init_site(nsites, chain_H_flat[a * max_shells + 0], chain_X_flat[a * max_shells + 0],
                  chain_Y_flat[a * max_shells + 0], 0.0);
        chain_rxn_flat[nsites * max_shells + 0] = chain_rxn_flat[a * max_shells + 0];
        ReactiveSite &ps = sites[nsites];
        ps.rxn_idx = sites[a].rxn_idx;
        ps.shell = sites[a].shell;
        ps.chain_len = sites[a].chain_len;
        ps.n_components = (int) combo.size();
        for (int ci = 0; ci < (int) combo.size(); ci++) ps.components[ci] = combo[ci];
        for (int c = (int) combo.size(); c < MAX_COMPONENTS; c++) ps.components[c] = -1;
        if (GN > 0) {
          memcpy(glove_flat + nsites * GN, glove_flat + a * GN, GN * sizeof(tagint));
          memcpy(chain_glove_flat + nsites * max_shells * glove_nmax,
                 chain_glove_flat + a * max_shells * glove_nmax, glove_nmax * sizeof(tagint));
        }
        nsites++;
      }
    }

    // Prune to the max_states budget (ranked by collective shells, then
    // distance) before the sites are broadcast.
    prune_states_to_max();

  }    // end ipartition == 0

broadcast:
  MPI_Bcast(&nsites, 1, MPI_INT, 0, samerank);
  MPI_Bcast(&nsites, 1, MPI_INT, 0, world);

  nstates = 1 + nsites;

  // Expose one global-vector entry per EVB state this step (reference + all
  // reactive states), so f_evb has length nstates rather than npartitions.
  size_vector = nstates;

  nsites_parallel = MIN(nsites, npartitions - 1);
  nsites_serial = nsites - nsites_parallel;

  if (nsites_serial > 0 && universe->me == 0) {
    if (!partition_warning) {
      auto msg = fmt::format("WARNING: Fix msevb detected {} reactive states but only {} "
                             "partitions are available; \nWARNING: {} excess state(s) will be "
                             "evaluated "
                             "in {} batched round(s). \nWARNING: This will only affect simulation "
                             "speed "
                             "and will not impact results. \nWARNING: This warning will not "
                             "appear again. "
                             "\n",
                             nstates, npartitions, nsites_serial,
                             (nsites_serial + npartitions - 1) / npartitions);
      utils::logmesg(lmp, msg);
      partition_warning = true;
    }
  }

  // Ensure all partitions have sites[] large enough, then broadcast site data.
  grow_sites(MAX(nsites, 1));

  if (nsites > 0) {
    // Per site in broadcast buffer:
    //   tag_H, tag_X, tag_Y, parent_state, shell, chain_len, rxn_idx,
    //   chain_H[MSD], chain_X[MSD], chain_Y[MSD], chain_rxn[MSD],
    //   glove[GN], chain_glove[MSD*GN],
    //   n_components, components[MPC]
    const int MSD = max_shells;
    const int GN = glove_nmax;
    const int MPC = MAX_COMPONENTS;
    const int PER_SITE = 7 + 4 * MSD + GN + MSD * GN + 1 + MPC;

    int buf_needed = nsites * PER_SITE;
    tagint *tagbuf;
    bool used_tmp = false;
    if ((int) (nmax * sizeof(double) / sizeof(tagint)) >= buf_needed) {
      tagbuf = reinterpret_cast<tagint *>(commbuf);
    } else {
      tagbuf = new tagint[buf_needed];
      used_tmp = true;
    }

    if (ipartition == 0) {
      for (int k = 0; k < nsites; k++) {
        int base = k * PER_SITE;
        tagbuf[base + 0] = sites[k].tag_H;
        tagbuf[base + 1] = sites[k].tag_X;
        tagbuf[base + 2] = sites[k].tag_Y;
        tagbuf[base + 3] = static_cast<tagint>(sites[k].parent_state);
        tagbuf[base + 4] = static_cast<tagint>(sites[k].shell);
        tagbuf[base + 5] = static_cast<tagint>(sites[k].chain_len);
        tagbuf[base + 6] = static_cast<tagint>(sites[k].rxn_idx);
        for (int d = 0; d < MSD; d++) {
          tagbuf[base + 7 + d] = chain_H_flat[k * max_shells + d];
          tagbuf[base + 7 + MSD + d] = chain_X_flat[k * max_shells + d];
          tagbuf[base + 7 + 2 * MSD + d] = chain_Y_flat[k * max_shells + d];
          tagbuf[base + 7 + 3 * MSD + d] = static_cast<tagint>(chain_rxn_flat[k * max_shells + d]);
        }
        if (GN > 0) {
          const tagint *gptr = glove_flat + k * GN;
          for (int g = 0; g < GN; g++) tagbuf[base + 7 + 4 * MSD + g] = gptr[g];
          for (int d = 0; d < MSD; d++) {
            const tagint *cg = chain_glove_flat + (k * max_shells + d) * glove_nmax;
            for (int g = 0; g < GN; g++) tagbuf[base + 7 + 4 * MSD + GN + d * GN + g] = cg[g];
          }
        }
        tagbuf[base + 7 + 4 * MSD + GN + MSD * GN] = static_cast<tagint>(sites[k].n_components);
        for (int c = 0; c < MPC; c++)
          tagbuf[base + 7 + 4 * MSD + GN + MSD * GN + 1 + c] =
              static_cast<tagint>((c < sites[k].n_components) ? sites[k].components[c] : -1);
      }
    }

    MPI_Bcast(tagbuf, buf_needed, MPI_LMP_TAGINT, 0, samerank);
    MPI_Bcast(tagbuf, buf_needed, MPI_LMP_TAGINT, 0, world);

    // On non-p0 partitions, ensure glove_flat is sized.
    if (GN > 0) grow_sites(MAX(nsites, 1));

    for (int k = 0; k < nsites; k++) {
      int base = k * PER_SITE;
      sites[k].tag_H = tagbuf[base + 0];
      sites[k].tag_X = tagbuf[base + 1];
      sites[k].tag_Y = tagbuf[base + 2];
      sites[k].parent_state = static_cast<int>(tagbuf[base + 3]);
      sites[k].shell = static_cast<int>(tagbuf[base + 4]);
      sites[k].chain_len = static_cast<int>(tagbuf[base + 5]);
      sites[k].rxn_idx = static_cast<int>(tagbuf[base + 6]);
      for (int d = 0; d < MSD; d++) {
        chain_H_flat[k * max_shells + d] = tagbuf[base + 7 + d];
        chain_X_flat[k * max_shells + d] = tagbuf[base + 7 + MSD + d];
        chain_Y_flat[k * max_shells + d] = tagbuf[base + 7 + 2 * MSD + d];
        chain_rxn_flat[k * max_shells + d] = static_cast<int>(tagbuf[base + 7 + 3 * MSD + d]);
      }
      if (GN > 0) {
        tagint *gptr = glove_flat + k * GN;
        for (int g = 0; g < GN; g++) gptr[g] = tagbuf[base + 7 + 4 * MSD + g];
        for (int d = 0; d < MSD; d++) {
          tagint *cg = chain_glove_flat + (k * max_shells + d) * glove_nmax;
          for (int g = 0; g < GN; g++) cg[g] = tagbuf[base + 7 + 4 * MSD + GN + d * GN + g];
        }
      }
      sites[k].n_components = static_cast<int>(tagbuf[base + 7 + 4 * MSD + GN + MSD * GN]);
      for (int c = 0; c < MPC; c++)
        sites[k].components[c] =
            static_cast<int>(tagbuf[base + 7 + 4 * MSD + GN + MSD * GN + 1 + c]);
    }

    if (used_tmp) delete[] tagbuf;
  }

  // Resolve global tags to local indices on this rank (-1 if not owned here).
  for (int k = 0; k < nsites; k++) {
    sites[k].idx_H = atom->map(sites[k].tag_H);
    sites[k].idx_X = atom->map(sites[k].tag_X);
    sites[k].idx_Y = atom->map(sites[k].tag_Y);
  }

  // Resolve product-state coupling neighbours (needs the finalized, broadcast
  // sites[]; deterministic, so every partition computes the same result).
  build_product_neighbors();
}

/* ----------------------------------------------------------------------
   Atom-set of a reactive site's transfer chain: the real atom tags the chain
   modifies, i.e. the non-edge glove atoms at every chain depth.  Edge atoms are
   matched by the template but not rewritten, so they are excluded.  This is the
   Chain::atom_set() primitive used for product-state (and future multi-site)
   compatibility decisions.
---------------------------------------------------------------------- */

std::vector<tagint> FixMSEVB::chain_atom_set(int site) const
{
  std::vector<tagint> tags;
  const int clen = sites[site].chain_len;
  for (int d = 0; d < clen; d++) {
    const ReactionDef &rx = rxndefs[chain_rxn_flat[site * max_shells + d]];
    const tagint *g = chain_glove_flat + (size_t) (site * max_shells + d) * glove_nmax;
    for (int pi = 0; pi < rx.glove_n; pi++) {
      if (rx.is_edge[pi]) continue;
      if (g[pi] != 0) tags.push_back(g[pi]);
    }
  }
  return tags;
}

/* ---------------------------------------------------------------------- */

bool FixMSEVB::chains_disjoint(int site_a, int site_b) const
{
  const std::vector<tagint> ta = chain_atom_set(site_a);
  const std::vector<tagint> tb = chain_atom_set(site_b);
  for (tagint x : ta)
    for (tagint y : tb)
      if (x == y) return false;
  return true;
}

/* ----------------------------------------------------------------------
   Ranking keys used to prune to max_states: collective shell depth (sum of the
   chain lengths of a state's components; a single-site state has one component
   = itself) and aggregate transfer distance (sum of the components' seed
   distances).  Lower of each is kept preferentially.
---------------------------------------------------------------------- */

int FixMSEVB::site_collective_shells(int site) const
{
  const ReactiveSite &s = sites[site];
  if (s.n_components == 0) return s.chain_len;
  int tot = 0;
  for (int ci = 0; ci < s.n_components; ci++) tot += sites[s.components[ci]].chain_len;
  return tot;
}

double FixMSEVB::site_aggregate_dist(int site) const
{
  const ReactiveSite &s = sites[site];
  if (s.n_components == 0) return s.dist_sq;
  double tot = 0.0;
  for (int ci = 0; ci < s.n_components; ci++) tot += sites[s.components[ci]].dist_sq;
  return tot;
}

/* ----------------------------------------------------------------------
   For each product state, resolve neighbor_state[j] = the state index of the
   configuration obtained by removing component j (the state that differs from
   this product by exactly component j's transfer).  Runs on every partition
   after sites[] is finalized and broadcast, so the coupling routines can place
   an n-way product's off-diagonals onto its (n-1)-component neighbours.

   Config key = the sorted set of component-site indices: {i} for a single
   (state i+1), {i,j,...} for a product.  Removing one component of a 2-way
   product leaves a single; of a 3-way leaves a 2-way; etc.  Every such subset
   is itself an enumerated state, so the lookup always resolves (the pruning
   ordering keeps the kept set closed under component removal).
---------------------------------------------------------------------- */

void FixMSEVB::build_product_neighbors()
{
  // Map each state's canonical component-set to its site index.
  std::map<std::vector<int>, int> config_to_site;
  for (int k = 0; k < nsites; k++) {
    std::vector<int> key;
    if (sites[k].n_components == 0) {
      key.push_back(k);
    } else {
      for (int ci = 0; ci < sites[k].n_components; ci++) key.push_back(sites[k].components[ci]);
      std::sort(key.begin(), key.end());
    }
    config_to_site[key] = k;
  }

  for (int k = 0; k < nsites; k++) {
    ReactiveSite &s = sites[k];
    if (s.n_components == 0) continue;
    for (int j = 0; j < s.n_components; j++) {
      std::vector<int> key;
      for (int ci = 0; ci < s.n_components; ci++)
        if (ci != j) key.push_back(s.components[ci]);
      std::sort(key.begin(), key.end());
      auto it = config_to_site.find(key);
      s.neighbor_state[j] = (it != config_to_site.end()) ? it->second + 1 : 0;
    }
    for (int j = s.n_components; j < MAX_COMPONENTS; j++) s.neighbor_state[j] = 0;
  }
}

/* ----------------------------------------------------------------------
   Prune the reactive states to at most max_states-1 (reserving one slot for
   the reference), keeping the lowest collective shell depth first and, within
   a tier, the smallest aggregate distance.  Survivors are re-indexed in that
   ranked order (so the most important states occupy the parallel partitions).

   The ranking is closed under removing a component or shortening a chain: a
   product's collective depth is strictly greater than any of its components',
   and a chain state's is strictly greater than its parent's, so every kept
   state's coupling neighbours (which have strictly smaller depth) are kept too.
   Runs on partition 0 only, before the sites are broadcast.
---------------------------------------------------------------------- */

void FixMSEVB::prune_states_to_max()
{
  if (max_states <= 0 || nsites + 1 <= max_states) return;
  const int keep = max_states - 1;    // reactive states (reference is state 0)
  const int GN = glove_nmax;
  const int ms = max_shells;

  // Rank all sites by (collective shells, aggregate distance, index).
  std::vector<int> order(nsites);
  for (int i = 0; i < nsites; i++) order[i] = i;
  std::stable_sort(order.begin(), order.end(), [&](int a, int b) {
    const int ca = site_collective_shells(a), cb = site_collective_shells(b);
    if (ca != cb) return ca < cb;
    const double da = site_aggregate_dist(a), db = site_aggregate_dist(b);
    if (da != db) return da < db;
    return a < b;
  });

  std::vector<int> old2new(nsites, -1);
  for (int n = 0; n < keep; n++) old2new[order[n]] = n;

  // Rebuild the kept sites (and their chain/glove data) in ranked order,
  // remapping component and parent references through old2new.
  std::vector<ReactiveSite> new_sites(keep);
  std::vector<tagint> nH((size_t) keep * ms), nX((size_t) keep * ms), nY((size_t) keep * ms);
  std::vector<int> nR((size_t) keep * ms);
  std::vector<tagint> ng((size_t) keep * MAX(GN, 1));
  std::vector<tagint> ncg((size_t) keep * ms * MAX(glove_nmax, 1));

  for (int n = 0; n < keep; n++) {
    const int old = order[n];
    ReactiveSite s = sites[old];
    for (int ci = 0; ci < s.n_components; ci++) {
      int nc = old2new[s.components[ci]];
      s.components[ci] = (nc >= 0) ? nc : 0;
    }
    if (s.parent_state > 0) {
      int np = old2new[s.parent_state - 1];
      s.parent_state = (np >= 0) ? np + 1 : 0;
    }
    new_sites[n] = s;
    for (int d = 0; d < ms; d++) {
      nH[(size_t) n * ms + d] = chain_H_flat[(size_t) old * ms + d];
      nX[(size_t) n * ms + d] = chain_X_flat[(size_t) old * ms + d];
      nY[(size_t) n * ms + d] = chain_Y_flat[(size_t) old * ms + d];
      nR[(size_t) n * ms + d] = chain_rxn_flat[(size_t) old * ms + d];
    }
    if (GN > 0) {
      for (int g = 0; g < GN; g++) ng[(size_t) n * GN + g] = glove_flat[(size_t) old * GN + g];
      for (int d = 0; d < ms; d++)
        for (int g = 0; g < GN; g++)
          ncg[((size_t) n * ms + d) * glove_nmax + g] =
              chain_glove_flat[((size_t) old * ms + d) * glove_nmax + g];
    }
  }

  // Copy back into the live arrays.
  for (int n = 0; n < keep; n++) {
    sites[n] = new_sites[n];
    for (int d = 0; d < ms; d++) {
      chain_H_flat[(size_t) n * ms + d] = nH[(size_t) n * ms + d];
      chain_X_flat[(size_t) n * ms + d] = nX[(size_t) n * ms + d];
      chain_Y_flat[(size_t) n * ms + d] = nY[(size_t) n * ms + d];
      chain_rxn_flat[(size_t) n * ms + d] = nR[(size_t) n * ms + d];
    }
    if (GN > 0) {
      for (int g = 0; g < GN; g++) glove_flat[(size_t) n * GN + g] = ng[(size_t) n * GN + g];
      for (int d = 0; d < ms; d++)
        for (int g = 0; g < GN; g++)
          chain_glove_flat[((size_t) n * ms + d) * glove_nmax + g] =
              ncg[((size_t) n * ms + d) * glove_nmax + g];
    }
  }

  nsites = keep;
}

/* ----------------------------------------------------------------------
   Apply per-partition state changes: types, charges, bonds, angles,
   then rebuild special bonds via Special::build().
---------------------------------------------------------------------- */

void FixMSEVB::apply_per_partition_state_changes()
{
  if (!reaction_enabled || nsites == 0) return;

  int istate = (ipartition > 0 && ipartition <= nsites_parallel) ? ipartition : 0;
  if (istate == 0) return;

  int sk = istate - 1;

  apply_site_changes(sk);

  Special special(lmp);
  special.build(true);

  // Charges changed — kspace (PPPM) caches qsqsum for its self-energy
  // correction and only updates it when atom count changes.  Force a
  // recalculation so the kspace energy reflects the new charges.
  if (force->kspace) force->kspace->qsum_qsq(0);
  modified_topology_on_host();
}

/* ----------------------------------------------------------------------
   Apply state changes for all chain depths of site sk.
   For product states (n_components > 0): iterates every component and
   every depth within each component, with forward_comm between depths
   and between components.
   For single sites: iterates all chain depths, with forward_comm between
   consecutive depths.
---------------------------------------------------------------------- */

void FixMSEVB::apply_site_changes(int sk)
{
  if (sites[sk].n_components > 0) {
    for (int ci = 0; ci < sites[sk].n_components; ci++) {
      int comp = sites[sk].components[ci];
      const int comp_len = sites[comp].chain_len;
      for (int d = 0; d < comp_len; d++) {
        const ReactionDef &rxn_c = rxndefs[chain_rxn_flat[comp * max_shells + d]];
        const tagint *glove_c = chain_glove_flat + (comp * max_shells + d) * glove_nmax;
        tagint tX = chain_X_flat[comp * max_shells + d];
        tagint tH = chain_H_flat[comp * max_shells + d];
        tagint tY = chain_Y_flat[comp * max_shells + d];
        apply_state_change(glove_c, atom->map(tX), atom->map(tH), atom->map(tY), tX, tH, tY, rxn_c);
        if (d < comp_len - 1) comm->forward_comm(this);
      }
      if (ci < sites[sk].n_components - 1) comm->forward_comm(this);
    }
  } else {
    const int chain_len = sites[sk].chain_len;
    for (int d = 0; d < chain_len; d++) {
      const ReactionDef &rxn_d = rxndefs[chain_rxn_flat[sk * max_shells + d]];
      const tagint *glove_d = chain_glove_flat + (sk * max_shells + d) * glove_nmax;
      tagint tX = chain_X_flat[sk * max_shells + d];
      tagint tH = chain_H_flat[sk * max_shells + d];
      tagint tY = chain_Y_flat[sk * max_shells + d];
      apply_state_change(glove_d, atom->map(tX), atom->map(tH), atom->map(tY), tX, tH, tY, rxn_d);
      if (d < chain_len - 1) comm->forward_comm(this);
    }
  }
}

/* ----------------------------------------------------------------------
   Apply state change for one template transfer step.
   Uses the glove stored in glove_flat[site_idx * glove_nmax + *] to map
   pre-template atom indices to real global tags, then applies:
     1. Type/charge changes from rd.type_changes
     2. Bond breaks from rd.bond_breaks
     3. Bond creates from rd.bond_creates
     4. Rebuild angles for all non-edge, local atoms in the glove.
     5. Update molecule ID for H (transferred atom).
---------------------------------------------------------------------- */

void FixMSEVB::apply_state_change(const tagint *glove, int idx_X, int idx_H, int idx_Y,
                                  tagint tag_X, tagint tag_H, tagint tag_Y, const ReactionDef &rxn)
{
  int *type = atom->type;
  double *q = atom->q;
  tagint *molecule = atom->molecule;
  int *num_bond_arr = atom->num_bond;
  int **bond_type_arr = atom->bond_type;
  tagint **bond_atom_arr = atom->bond_atom;
  const int nlocal = atom->nlocal;

  if (!glove) return;

  const int GN = rxn.glove_n;

  // Helper: remove a bond stored on atom idx pointing to partner_tag.
  auto remove_bond_on = [&](int idx, tagint partner_tag) {
    if (idx < 0 || idx >= nlocal) return;
    for (int k = 0; k < num_bond_arr[idx]; k++) {
      if (bond_atom_arr[idx][k] == partner_tag) {
        int last = num_bond_arr[idx] - 1;
        if (k != last) {
          bond_atom_arr[idx][k] = bond_atom_arr[idx][last];
          bond_type_arr[idx][k] = bond_type_arr[idx][last];
        }
        num_bond_arr[idx]--;
        return;
      }
    }
  };

  // Step 1: Apply type/charge changes.
  for (const auto &tc : rxn.type_changes) {
    if (tc.pre_idx < 0 || tc.pre_idx >= GN) continue;
    tagint r = glove[tc.pre_idx];
    if (r == 0) continue;
    int local_r = atom->map(r);
    if (local_r < 0 || local_r >= nlocal) continue;
    type[local_r] = tc.new_type;
    // Set charge from the post-template; look it up from post_mol if available.
    if (rxn.post_mol && rxn.pre_to_post[tc.pre_idx] >= 0) {
      int post_i = rxn.pre_to_post[tc.pre_idx];
      if (rxn.post_mol->qflag) q[local_r] = rxn.post_mol->q[post_i];
    }
  }

  // Step 2: Bond breaks.
  for (const auto &bb : rxn.bond_breaks) {
    if (bb.pre_idx1 < 0 || bb.pre_idx1 >= GN) continue;
    if (bb.pre_idx2 < 0 || bb.pre_idx2 >= GN) continue;
    tagint r1 = glove[bb.pre_idx1];
    tagint r2 = glove[bb.pre_idx2];
    if (r1 == 0 || r2 == 0) continue;
    int local_r1 = atom->map(r1);
    int local_r2 = atom->map(r2);
    remove_bond_on(local_r1, r2);
    remove_bond_on(local_r2, r1);
  }

  // Step 3: Bond creates.
  // Store each bond on exactly ONE atom: the one with the lower global tag,
  // but only if that atom is local to THIS rank.  This avoids double-storing
  // the bond when both endpoints are owned by different ranks in a multi-rank
  // partition (without this rule each rank would independently add the bond
  // to its local endpoint, giving a double-counted bond energy).
  for (const auto &bc : rxn.bond_creates) {
    if (bc.pre_idx1 < 0 || bc.pre_idx1 >= GN) continue;
    if (bc.pre_idx2 < 0 || bc.pre_idx2 >= GN) continue;
    tagint r1 = glove[bc.pre_idx1];
    tagint r2 = glove[bc.pre_idx2];
    if (r1 == 0 || r2 == 0) continue;

    // Determine the canonical "owner": atom with the lower global tag.
    tagint owner_tag = (r1 < r2) ? r1 : r2;
    tagint partner_tag = (r1 < r2) ? r2 : r1;
    int local_owner = atom->map(owner_tag);
    if (local_owner < 0 || local_owner >= nlocal)
      continue;    // owner is on another rank; that rank will store the bond

    if (num_bond_arr[local_owner] >= atom->bond_per_atom)
      error->one(FLERR,
                 "Fix msevb: bond count on an atom exceeds bond_per_atom during a "
                 "template transfer; increase it via the data file 'extra/bond/per/atom' "
                 "keyword (or create_box's bond/per/atom setting)");

    bond_atom_arr[local_owner][num_bond_arr[local_owner]] = partner_tag;
    bond_type_arr[local_owner][num_bond_arr[local_owner]] = bc.bond_type;
    num_bond_arr[local_owner]++;
  }

  // Step 3.5: Bond retypes — bonds that persist between the same pair of atoms
  // but change type (e.g., O-H bonds changing from H3O+ type to H2O type).
  // Scan bond lists of both endpoints; update whichever side stores the bond.
  for (const auto &br : rxn.bond_retypes) {
    if (br.pre_idx1 < 0 || br.pre_idx1 >= GN) continue;
    if (br.pre_idx2 < 0 || br.pre_idx2 >= GN) continue;
    tagint r1 = glove[br.pre_idx1];
    tagint r2 = glove[br.pre_idx2];
    if (r1 == 0 || r2 == 0) continue;
    // Update on whichever local atom stores the bond.
    for (int side = 0; side < 2; side++) {
      tagint owner_tag = (side == 0) ? r1 : r2;
      tagint partner_tag = (side == 0) ? r2 : r1;
      int local_owner = atom->map(owner_tag);
      if (local_owner < 0 || local_owner >= nlocal) continue;
      for (int k = 0; k < num_bond_arr[local_owner]; k++) {
        if (bond_atom_arr[local_owner][k] == partner_tag) {
          bond_type_arr[local_owner][k] = br.new_bond_type;
          break;
        }
      }
    }
  }

  // Step 4: Rebuild angles from the post-template for all non-edge glove atoms.
  //
  // Strategy:
  //   a. Zero num_angle for every non-edge glove atom that is local (so stale
  //      pre-reaction angles — e.g. the extra H3O+ angles on the donor O after
  //      it becomes H2O — are discarded).
  //   b. Walk the post_mol angle list.  For each post angle, map the three
  //      atom indices through post_to_pre → glove to obtain real global tags.
  //      Store the angle on the central atom (atom2) if it is local.
  //
  // This mirrors exactly what MuStaRD does: angle types come from the post-mol
  // template and are fully determined by the topology diff.
  // Build post_to_pre map once (inverse of rxn.pre_to_post).
  // Shared by angle, dihedral, and improper rebuild steps below.
  const int natoms_post = rxn.post_mol ? rxn.post_mol->natoms : 0;
  std::vector<int> post_to_pre(natoms_post, -1);
  if (rxn.post_mol) {
    for (int i = 0; i < GN; i++)
      if (rxn.pre_to_post[i] >= 0) post_to_pre[rxn.pre_to_post[i]] = i;
  }

  // NOTE: the guard must NOT require post_mol->nangles > 0.  Step (a) below
  // zeroes stale angle counts on the reactive atoms; if a product template has
  // zero angles (e.g. reacting to a species with none) but the reference does
  // not, skipping this block would leave the reference's angles on the daughter
  // atoms and corrupt its energy.  The add-loop is a no-op when nangles == 0.
  if (atom->nangles > 0 && rxn.post_mol) {
    int *num_angle_arr = atom->num_angle;
    int **angle_type_arr = atom->angle_type;
    tagint **angle_atom1_arr = atom->angle_atom1;
    tagint **angle_atom2_arr = atom->angle_atom2;
    tagint **angle_atom3_arr = atom->angle_atom3;

    // a. Zero angle counts for non-edge local glove atoms.
    for (int pi = 0; pi < GN; pi++) {
      if (rxn.is_edge[pi]) continue;
      tagint r = glove[pi];
      if (r == 0) continue;
      int local_r = atom->map(r);
      if (local_r >= 0 && local_r < nlocal) num_angle_arr[local_r] = 0;
    }

    // b. Add angles from post_mol, storing on the central atom (atom2).
    // In LAMMPS Molecule, angles are stored per-atom (central atom = atom2).
    // Iterate per post-mol atom m; only process entries where atom2-1 == m
    // to visit each unique angle exactly once.
    for (int m = 0; m < natoms_post; m++) {
      int pre_m = post_to_pre[m];
      if (pre_m < 0) continue;    // no pre-atom maps here
      for (int k = 0; k < rxn.post_mol->num_angle[m]; k++) {
        int pa1 = rxn.post_mol->angle_atom1[m][k] - 1;    // 0-based post indices
        int pa2 = rxn.post_mol->angle_atom2[m][k] - 1;    // central atom
        int pa3 = rxn.post_mol->angle_atom3[m][k] - 1;
        if (pa2 != m) continue;    // only process when m is the central atom
        int pre1 = (pa1 >= 0 && pa1 < natoms_post) ? post_to_pre[pa1] : -1;
        int pre2 = pre_m;    // pa2 == m
        int pre3 = (pa3 >= 0 && pa3 < natoms_post) ? post_to_pre[pa3] : -1;
        if (pre1 < 0 || pre3 < 0) continue;
        tagint r1 = glove[pre1], r2 = glove[pre2], r3 = glove[pre3];
        if (r1 == 0 || r2 == 0 || r3 == 0) continue;
        int local_r2 = atom->map(r2);
        if (local_r2 < 0 || local_r2 >= nlocal) continue;
        int na = num_angle_arr[local_r2];
        if (na >= atom->angle_per_atom)
          error->one(FLERR,
                     "Fix msevb: angle count on an atom exceeds angle_per_atom during a "
                     "template transfer; increase it via the data file 'extra/angle/per/atom' "
                     "keyword (or create_box's angle/per/atom setting)");
        angle_type_arr[local_r2][na] = rxn.post_mol->angle_type[m][k];
        angle_atom1_arr[local_r2][na] = r1;
        angle_atom2_arr[local_r2][na] = r2;
        angle_atom3_arr[local_r2][na] = r3;
        num_angle_arr[local_r2]++;
      }
    }
  }

  // Step 4.5: Rebuild dihedrals from the post-template for all non-edge
  // glove atoms.  Same pattern as angles: zero counts, walk post_mol list,
  // filter by storage atom (atom2 in LAMMPS Molecule dihedral convention),
  // map through post_to_pre -> glove to real tags.
  // See the angle note above: the guard must not require post_mol->ndihedrals
  // > 0, or a product with no dihedrals would leave the reference's dihedrals
  // (e.g. O-C-O-H torsions) on the daughter atoms and inflate its energy.
  if (force->dihedral && rxn.post_mol) {
    int *num_dihedral_arr = atom->num_dihedral;
    int **dihedral_type_arr = atom->dihedral_type;
    tagint **dihedral_atom1_arr = atom->dihedral_atom1;
    tagint **dihedral_atom2_arr = atom->dihedral_atom2;
    tagint **dihedral_atom3_arr = atom->dihedral_atom3;
    tagint **dihedral_atom4_arr = atom->dihedral_atom4;

    // a. Zero dihedral counts for non-edge local glove atoms.
    for (int pi = 0; pi < GN; pi++) {
      if (rxn.is_edge[pi]) continue;
      tagint r = glove[pi];
      if (r == 0) continue;
      int local_r = atom->map(r);
      if (local_r >= 0 && local_r < nlocal) num_dihedral_arr[local_r] = 0;
    }

    // b. Add dihedrals from post_mol, storing on atom2 (LAMMPS convention).
    for (int m = 0; m < natoms_post; m++) {
      int pre_m = post_to_pre[m];
      if (pre_m < 0) continue;
      for (int k = 0; k < rxn.post_mol->num_dihedral[m]; k++) {
        int pa1 = rxn.post_mol->dihedral_atom1[m][k] - 1;
        int pa2 = rxn.post_mol->dihedral_atom2[m][k] - 1;
        int pa3 = rxn.post_mol->dihedral_atom3[m][k] - 1;
        int pa4 = rxn.post_mol->dihedral_atom4[m][k] - 1;
        if (pa2 != m) continue;    // only process when m is the storage atom
        int pre1 = (pa1 >= 0 && pa1 < natoms_post) ? post_to_pre[pa1] : -1;
        int pre2 = pre_m;
        int pre3 = (pa3 >= 0 && pa3 < natoms_post) ? post_to_pre[pa3] : -1;
        int pre4 = (pa4 >= 0 && pa4 < natoms_post) ? post_to_pre[pa4] : -1;
        if (pre1 < 0 || pre3 < 0 || pre4 < 0) continue;
        tagint r1 = glove[pre1], r2 = glove[pre2];
        tagint r3 = glove[pre3], r4 = glove[pre4];
        if (r1 == 0 || r2 == 0 || r3 == 0 || r4 == 0) continue;
        int local_r2 = atom->map(r2);
        if (local_r2 < 0 || local_r2 >= nlocal) continue;
        int nd = num_dihedral_arr[local_r2];
        if (nd >= atom->dihedral_per_atom)
          error->one(FLERR,
                     "Fix msevb: dihedral count on an atom exceeds dihedral_per_atom during a "
                     "template transfer; increase it via the data file 'extra/dihedral/per/atom' "
                     "keyword (or create_box's dihedral/per/atom setting)");
        dihedral_type_arr[local_r2][nd] = rxn.post_mol->dihedral_type[m][k];
        dihedral_atom1_arr[local_r2][nd] = r1;
        dihedral_atom2_arr[local_r2][nd] = r2;
        dihedral_atom3_arr[local_r2][nd] = r3;
        dihedral_atom4_arr[local_r2][nd] = r4;
        num_dihedral_arr[local_r2]++;
      }
    }
  }

  // Step 4.6: Rebuild impropers from the post-template for all non-edge
  // glove atoms.  Same pattern as dihedrals: impropers are stored on atom2.
  // See the angle note above: the guard must not require post_mol->nimpropers
  // > 0, or a product with no impropers would leave the reference's impropers on
  // the daughter atoms and inflate its energy.
  if (force->improper && rxn.post_mol) {
    int *num_improper_arr = atom->num_improper;
    int **improper_type_arr = atom->improper_type;
    tagint **improper_atom1_arr = atom->improper_atom1;
    tagint **improper_atom2_arr = atom->improper_atom2;
    tagint **improper_atom3_arr = atom->improper_atom3;
    tagint **improper_atom4_arr = atom->improper_atom4;

    // a. Zero improper counts for non-edge local glove atoms.
    for (int pi = 0; pi < GN; pi++) {
      if (rxn.is_edge[pi]) continue;
      tagint r = glove[pi];
      if (r == 0) continue;
      int local_r = atom->map(r);
      if (local_r >= 0 && local_r < nlocal) num_improper_arr[local_r] = 0;
    }

    // b. Add impropers from post_mol, storing on atom2.
    for (int m = 0; m < natoms_post; m++) {
      int pre_m = post_to_pre[m];
      if (pre_m < 0) continue;
      for (int k = 0; k < rxn.post_mol->num_improper[m]; k++) {
        int pa1 = rxn.post_mol->improper_atom1[m][k] - 1;
        int pa2 = rxn.post_mol->improper_atom2[m][k] - 1;
        int pa3 = rxn.post_mol->improper_atom3[m][k] - 1;
        int pa4 = rxn.post_mol->improper_atom4[m][k] - 1;
        if (pa2 != m) continue;
        int pre1 = (pa1 >= 0 && pa1 < natoms_post) ? post_to_pre[pa1] : -1;
        int pre2 = pre_m;
        int pre3 = (pa3 >= 0 && pa3 < natoms_post) ? post_to_pre[pa3] : -1;
        int pre4 = (pa4 >= 0 && pa4 < natoms_post) ? post_to_pre[pa4] : -1;
        if (pre1 < 0 || pre3 < 0 || pre4 < 0) continue;
        tagint r1 = glove[pre1], r2 = glove[pre2];
        tagint r3 = glove[pre3], r4 = glove[pre4];
        if (r1 == 0 || r2 == 0 || r3 == 0 || r4 == 0) continue;
        int local_r2 = atom->map(r2);
        if (local_r2 < 0 || local_r2 >= nlocal) continue;
        int ni = num_improper_arr[local_r2];
        if (ni >= atom->improper_per_atom)
          error->one(FLERR,
                     "Fix msevb: improper count on an atom exceeds improper_per_atom during a "
                     "template transfer; increase it via the data file 'extra/improper/per/atom' "
                     "keyword (or create_box's improper/per/atom setting)");
        improper_type_arr[local_r2][ni] = rxn.post_mol->improper_type[m][k];
        improper_atom1_arr[local_r2][ni] = r1;
        improper_atom2_arr[local_r2][ni] = r2;
        improper_atom3_arr[local_r2][ni] = r3;
        improper_atom4_arr[local_r2][ni] = r4;
        num_improper_arr[local_r2]++;
      }
    }
  }

  // Step 5: Update molecule ID for the transferred H atom.
  // H (ibonding) moves from X's molecule to Y's molecule.
  tagint tag_H_used = (tag_H != 0) ? tag_H
                                   : ((glove_nmax > rxn.ibonding) ? glove[rxn.ibonding] : 0);
  tagint tag_Y_used = (tag_Y != 0) ? tag_Y
                                   : ((glove_nmax > rxn.jbonding) ? glove[rxn.jbonding] : 0);
  if (tag_H_used != 0) {
    int local_H = atom->map(tag_H_used);
    if (local_H >= 0 && local_H < nlocal) {
      tagint mol_Y = -1;
      if (tag_Y_used > 0 && tag_Y_used <= ref_maxtag)
        mol_Y = ref_molecule[tag_Y_used];
      else {
        int local_Y = atom->map(tag_Y_used);
        if (local_Y >= 0) mol_Y = molecule[local_Y];
      }
      if (mol_Y > 0) molecule[local_H] = mol_Y;
    }
  }
}

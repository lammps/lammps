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

#include "bond_bpm.h"

#include "atom.h"
#include "citeme.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "fix_bond_history.h"
#include "fix_store_local.h"
#include "fix_update_special_bonds.h"
#include "force.h"
#include "memory.h"
#include "modify.h"
#include "update.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;

static const char cite_bpm[] =
  "BPM bond style: doi:10.1039/D3SM01373A\n\n"
  "@Article{Clemmer2024,\n"
  " author =  {Clemmer, Joel T. and Monti, Joseph M. and Lechman, Jeremy B.},\n"
  " title =   {A soft departure from jamming: the compaction of deformable\n"
  "            granular matter under high pressures},\n"
  " journal = {Soft Matter},\n"
  " year =    2024,\n"
  " volume =  20,\n"
  " number =  8,\n"
  " pages =   {1702--1718}\n"
  "}\n\n";

/* ---------------------------------------------------------------------- */

BondBPM::BondBPM(LAMMPS *_lmp) :
    Bond(_lmp), id_fix_dummy(nullptr), id_fix_dummy2(nullptr), id_fix_update(nullptr),
    id_fix_bond_history(nullptr), id_fix_store_local(nullptr), id_fix_prop_atom(nullptr),
    fix_store_local(nullptr), fix_bond_history(nullptr), fix_update_special_bonds(nullptr),
    pack_choice(nullptr), output_data(nullptr)
{
  overlay_flag = 0;
  prop_atom_flag = 0;
  break_flag = 1;
  nvalues = 0;

  r0_max_estimate = 0.0;
  max_stretch = 1.0;

  // create dummy fix as placeholder for FixUpdateSpecialBonds & BondHistory
  // this is so final order of Modify:fix will conform to input script
  // BondHistory technically only needs this if updateflag = 1

  id_fix_dummy = utils::strdup("BPM_DUMMY");
  modify->add_fix(fmt::format("{} all DUMMY ", id_fix_dummy));

  id_fix_dummy2 = utils::strdup("BPM_DUMMY2");
  modify->add_fix(fmt::format("{} all DUMMY ", id_fix_dummy2));

  if (lmp->citeme) lmp->citeme->add(cite_bpm);
}

/* ---------------------------------------------------------------------- */

BondBPM::~BondBPM()
{
  delete[] pack_choice;

  if (id_fix_dummy) modify->delete_fix(id_fix_dummy);
  if (id_fix_dummy2) modify->delete_fix(id_fix_dummy2);
  if (id_fix_update) modify->delete_fix(id_fix_update);
  if (id_fix_bond_history) modify->delete_fix(id_fix_bond_history);
  if (id_fix_store_local) modify->delete_fix(id_fix_store_local);
  if (id_fix_prop_atom) modify->delete_fix(id_fix_prop_atom);

  delete[] id_fix_dummy;
  delete[] id_fix_dummy2;
  delete[] id_fix_update;
  delete[] id_fix_bond_history;
  delete[] id_fix_store_local;
  delete[] id_fix_prop_atom;

  memory->destroy(output_data);
}

/* ---------------------------------------------------------------------- */

void BondBPM::init_style()
{
  if (id_fix_store_local) {
    auto ifix = modify->get_fix_by_id(id_fix_store_local);
    if (!ifix) error->all(FLERR, "Cannot find fix STORE/LOCAL id {}", id_fix_store_local);
    if (strcmp(ifix->style, "STORE/LOCAL") != 0)
      error->all(FLERR, "Incorrect fix style matched, not STORE/LOCAL: {}", ifix->style);
    fix_store_local = dynamic_cast<FixStoreLocal *>(ifix);
    fix_store_local->nvalues = nvalues;
  }

  if (overlay_flag) {
    if (force->special_lj[1] != 1.0 || force->special_lj[2] != 1.0 || force->special_lj[3] != 1.0 ||
        force->special_coul[1] != 1.0 || force->special_coul[2] != 1.0 || force->special_coul[3] != 1.0)
      error->all(FLERR,
                 "With overlay/pair yes, BPM bond styles require a value of 1.0 for all special_bonds weights");
    if (id_fix_update) {
      modify->delete_fix(id_fix_update);
      delete[] id_fix_update;
      id_fix_update = nullptr;
    }
  } else {
    // Require atoms know about all of their bonds and if they break
    if (force->newton_bond && break_flag)
      error->all(FLERR, "With overlay/pair no, or break yes, BPM bond styles require Newton bond off");

    // special lj must be 0 1 1 to censor pair forces between bonded particles
    if (force->special_lj[1] != 0.0 || force->special_lj[2] != 1.0 || force->special_lj[3] != 1.0)
      error->all(FLERR,
                 "With overlay/pair no, BPM bond styles require special LJ weights = 0,1,1");
    // if bonds can break, special coulomb must be 1 1 1 to ensure all pairs are included in the
    //    neighbor list and 1-3 and 1-4 special bond lists are skipped
    if (break_flag && (force->special_coul[1] != 1.0 || force->special_coul[2] != 1.0 ||
        force->special_coul[3] != 1.0))
      error->all(FLERR,
                 "With overlay/pair no, and break yes, BPM bond styles requires special Coulomb weights = 1,1,1");

    if (id_fix_dummy && break_flag) {
      id_fix_update = utils::strdup("BPM_UPDATE_SPECIAL_BONDS");
      fix_update_special_bonds = dynamic_cast<FixUpdateSpecialBonds *>(modify->replace_fix(
          id_fix_dummy, fmt::format("{} all UPDATE_SPECIAL_BONDS", id_fix_update), 1));
      delete[] id_fix_dummy;
      id_fix_dummy = nullptr;
    }
  }

  if (force->angle || force->dihedral || force->improper)
    error->all(FLERR, "Bond style bpm cannot be used with 3,4-body interactions");
  if (atom->molecular == 2)
    error->all(FLERR, "Bond style bpm cannot be used with atom style template");

  // special 1-3 and 1-4 weights must be 1 to prevent building 1-3 and 1-4 special bond lists
  if (force->special_lj[2] != 1.0 || force->special_lj[3] != 1.0 || force->special_coul[2] != 1.0 ||
      force->special_coul[3] != 1.0)
    error->all(FLERR, "Bond style bpm requires 1-3 and 1-4 special weights of 1.0");
}

/* ----------------------------------------------------------------------
   global settings
   All args before store/local command are saved for potential args
     for specific bond BPM substyles
   All args after optional stode/local command are variables stored
     in the compute store/local
------------------------------------------------------------------------- */

void BondBPM::settings(int narg, char **arg)
{
  leftover_iarg.clear();

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "store/local") == 0) {
      nvalues = 0;
      id_fix_store_local = utils::strdup(arg[iarg + 1]);
      store_local_freq = utils::inumeric(FLERR, arg[iarg + 2], false, lmp);
      pack_choice = new FnPtrPack[narg - iarg - 1];
      iarg += 3;
      while (iarg < narg) {
        if (strcmp(arg[iarg], "id1") == 0) {
          pack_choice[nvalues++] = &BondBPM::pack_id1;
        } else if (strcmp(arg[iarg], "id2") == 0) {
          pack_choice[nvalues++] = &BondBPM::pack_id2;
        } else if (strcmp(arg[iarg], "time") == 0) {
          pack_choice[nvalues++] = &BondBPM::pack_time;
        } else if (strcmp(arg[iarg], "x") == 0) {
          pack_choice[nvalues++] = &BondBPM::pack_x;
        } else if (strcmp(arg[iarg], "y") == 0) {
          pack_choice[nvalues++] = &BondBPM::pack_y;
        } else if (strcmp(arg[iarg], "z") == 0) {
          pack_choice[nvalues++] = &BondBPM::pack_z;
        } else if (strcmp(arg[iarg], "x/ref") == 0) {
          pack_choice[nvalues++] = &BondBPM::pack_x_ref;
          prop_atom_flag = 1;
        } else if (strcmp(arg[iarg], "y/ref") == 0) {
          pack_choice[nvalues++] = &BondBPM::pack_y_ref;
          prop_atom_flag = 1;
        } else if (strcmp(arg[iarg], "z/ref") == 0) {
          pack_choice[nvalues++] = &BondBPM::pack_z_ref;
          prop_atom_flag = 1;
        } else {
          break;
        }
        iarg++;
      }
    } else if (strcmp(arg[iarg], "overlay/pair") == 0) {
      if (iarg + 1 > narg) error->all(FLERR, "Illegal bond bpm command, missing option for overlay/pair");
      overlay_flag = utils::logical(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "break") == 0) {
      if (iarg + 1 > narg) error->all(FLERR, "Illegal bond bpm command, missing option for break");
      break_flag = utils::logical(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else {
      leftover_iarg.push_back(iarg);
      iarg++;
    }
  }

  if (id_fix_store_local) {

    if (nvalues == 0)
      error->all(FLERR, "Storing local data must include at least one value to output");
    memory->create(output_data, nvalues, "bond/bpm:output_data");

    auto ifix = modify->get_fix_by_id(id_fix_store_local);
    if (!ifix)
      ifix = modify->add_fix(
          fmt::format("{} all STORE/LOCAL {} {}", id_fix_store_local, store_local_freq, nvalues));
    fix_store_local = dynamic_cast<FixStoreLocal *>(ifix);

    // Use property/atom to save reference positions as it can transfer to ghost atoms
    // This won't work for instances where bonds are added (e.g. fix pour) but in those cases
    // a reference state isn't well defined
    if (prop_atom_flag == 1) {

      id_fix_prop_atom = utils::strdup("BPM_property_atom");
      char *x_ref_id = utils::strdup("BPM_X_REF");
      char *y_ref_id = utils::strdup("BPM_Y_REF");
      char *z_ref_id = utils::strdup("BPM_Z_REF");

      ifix = modify->get_fix_by_id(id_fix_prop_atom);
      if (!ifix)
        ifix = modify->add_fix(fmt::format("{} all property/atom d_{} d_{} d_{} ghost yes",
                                           id_fix_prop_atom, x_ref_id, y_ref_id, z_ref_id));

      int type_flag;
      int col_flag;
      index_x_ref = atom->find_custom(x_ref_id, type_flag, col_flag);
      index_y_ref = atom->find_custom(y_ref_id, type_flag, col_flag);
      index_z_ref = atom->find_custom(z_ref_id, type_flag, col_flag);

      delete[] x_ref_id;
      delete[] y_ref_id;
      delete[] z_ref_id;

      if (ifix->restart_reset) {
        ifix->restart_reset = 0;
      } else {
        double *x_ref = atom->dvector[index_x_ref];
        double *y_ref = atom->dvector[index_y_ref];
        double *z_ref = atom->dvector[index_z_ref];

        double **x = atom->x;
        for (int i = 0; i < atom->nlocal; i++) {
          x_ref[i] = x[i][0];
          y_ref[i] = x[i][1];
          z_ref[i] = x[i][2];
        }
      }
    }
  }
}

/* ----------------------------------------------------------------------
   used to check bond communiction cutoff - not perfect, estimates based on local-local only
------------------------------------------------------------------------- */

double BondBPM::equilibrium_distance(int /*i*/)
{
  // Ghost atoms may not yet be communicated, this may only be an estimate
  if (r0_max_estimate == 0) {
    if (!fix_bond_history->restart_reset) {
      int type, j;
      double delx, dely, delz, r;
      double **x = atom->x;
      for (int i = 0; i < atom->nlocal; i++) {
        for (int m = 0; m < atom->num_bond[i]; m++) {
          type = atom->bond_type[i][m];
          if (type == 0) continue;

          j = atom->map(atom->bond_atom[i][m]);
          if (j == -1) continue;

          delx = x[i][0] - x[j][0];
          dely = x[i][1] - x[j][1];
          delz = x[i][2] - x[j][2];
          domain->minimum_image(delx, dely, delz);

          r = sqrt(delx * delx + dely * dely + delz * delz);
          if (r > r0_max_estimate) r0_max_estimate = r;
        }
      }
    } else {
      int type, j;
      double r;
      for (int i = 0; i < atom->nlocal; i++) {
        for (int m = 0; m < atom->num_bond[i]; m++) {
          type = atom->bond_type[i][m];
          if (type == 0) continue;

          j = atom->map(atom->bond_atom[i][m]);
          if (j == -1) continue;

          // First value must always be reference length
          r = fix_bond_history->get_atom_value(i, m, 0);
          if (r > r0_max_estimate) r0_max_estimate = r;
        }
      }
    }

    double temp;
    MPI_Allreduce(&r0_max_estimate, &temp, 1, MPI_DOUBLE, MPI_MAX, world);
    r0_max_estimate = temp;
  }

  // Divide out heuristic prefactor added in comm class
  return max_stretch * r0_max_estimate / 1.5;
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
 ------------------------------------------------------------------------- */

void BondBPM::write_restart(FILE *fp)
{
  fwrite(&overlay_flag, sizeof(int), 1, fp);
}

/* ----------------------------------------------------------------------
    proc 0 reads from restart file, bcasts
 ------------------------------------------------------------------------- */

void BondBPM::read_restart(FILE *fp)
{
  if (comm->me == 0) utils::sfread(FLERR, &overlay_flag, sizeof(int), 1, fp, nullptr, error);
  MPI_Bcast(&overlay_flag, 1, MPI_INT, 0, world);
}

/* ---------------------------------------------------------------------- */

void BondBPM::process_broken(int i, int j)
{
  if (!break_flag)
    error->one(FLERR, "BPM bond broke with break no option");
  if (fix_store_local) {
    for (int n = 0; n < nvalues; n++) (this->*pack_choice[n])(n, i, j);

    fix_store_local->add_data(output_data, i, j);
  }

  if (fix_update_special_bonds) fix_update_special_bonds->add_broken_bond(i, j);

  // Manually search and remove from atom arrays
  // need to remove in case special bonds arrays rebuilt
  int m, n;
  int nlocal = atom->nlocal;

  tagint *tag = atom->tag;
  tagint **bond_atom = atom->bond_atom;
  int **bond_type = atom->bond_type;
  int *num_bond = atom->num_bond;

  if (i < nlocal) {
    for (m = 0; m < num_bond[i]; m++) {
      if (bond_atom[i][m] == tag[j]) {
        n = num_bond[i];
        bond_type[i][m] = bond_type[i][n - 1];
        bond_atom[i][m] = bond_atom[i][n - 1];
        fix_bond_history->shift_history(i, m, n - 1);
        fix_bond_history->delete_history(i, n - 1);
        num_bond[i]--;
        break;
      }
    }
  }

  if (j < nlocal) {
    for (m = 0; m < num_bond[j]; m++) {
      if (bond_atom[j][m] == tag[i]) {
        n = num_bond[j];
        bond_type[j][m] = bond_type[j][n - 1];
        bond_atom[j][m] = bond_atom[j][n - 1];
        fix_bond_history->shift_history(j, m, n - 1);
        fix_bond_history->delete_history(j, n - 1);
        num_bond[j]--;
        break;
      }
    }
  }
}

/* ----------------------------------------------------------------------
   one method for every keyword bond bpm can output
   the atom property is packed into array or vector
------------------------------------------------------------------------- */

void BondBPM::pack_id1(int n, int i, int /*j*/)
{
  tagint *tag = atom->tag;
  output_data[n] = tag[i];
}

/* ---------------------------------------------------------------------- */

void BondBPM::pack_id2(int n, int /*i*/, int j)
{
  tagint *tag = atom->tag;
  output_data[n] = tag[j];
}

/* ---------------------------------------------------------------------- */

void BondBPM::pack_time(int n, int /*i*/, int /*j*/)
{
  bigint time = update->ntimestep;
  output_data[n] = time;
}

/* ---------------------------------------------------------------------- */

void BondBPM::pack_x(int n, int i, int j)
{
  double **x = atom->x;
  output_data[n] = (x[i][0] + x[j][0]) * 0.5;
}

/* ---------------------------------------------------------------------- */

void BondBPM::pack_y(int n, int i, int j)
{
  double **x = atom->x;
  output_data[n] = (x[i][1] + x[j][1]) * 0.5;
}

/* ---------------------------------------------------------------------- */

void BondBPM::pack_z(int n, int i, int j)
{
  double **x = atom->x;
  output_data[n] = (x[i][2] + x[j][2]) * 0.5;
}

/* ---------------------------------------------------------------------- */

void BondBPM::pack_x_ref(int n, int i, int j)
{
  double *x = atom->dvector[index_x_ref];
  output_data[n] = (x[i] + x[j]) * 0.5;
}

/* ---------------------------------------------------------------------- */

void BondBPM::pack_y_ref(int n, int i, int j)
{
  double *y = atom->dvector[index_y_ref];
  output_data[n] = (y[i] + y[j]) * 0.5;
}

/* ---------------------------------------------------------------------- */

void BondBPM::pack_z_ref(int n, int i, int j)
{
  double *z = atom->dvector[index_z_ref];
  output_data[n] = (z[i] + z[j]) * 0.5;
}

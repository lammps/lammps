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

#include "atom_vec_hybrid.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "tokenizer.h"

#include <algorithm>
#include <cstring>

using namespace LAMMPS_NS;

enum { ELLIPSOID, LINE, TRIANGLE, BODY };    // also in WriteData

/* ---------------------------------------------------------------------- */

AtomVecHybrid::AtomVecHybrid(LAMMPS *lmp) : AtomVec(lmp)
{
  nstyles = 0;
  styles = nullptr;
  keywords = nullptr;

  bonus_flag = 0;
  nstyles_bonus = 0;
  styles_bonus = nullptr;

  // field strings will be concatenated from sub-style strings
  // fields_data_atom & fields_data_vel start with fields common to all styles

  fields_data_atom = {"id", "type", "x"};
  fields_data_vel = {"id", "v"};
}

/* ---------------------------------------------------------------------- */

AtomVecHybrid::~AtomVecHybrid()
{
  for (int k = 0; k < nstyles; k++) delete styles[k];
  delete[] styles;
  for (int k = 0; k < nstyles; k++) delete[] keywords[k];
  delete[] keywords;
  delete[] styles_bonus;
}

/* ----------------------------------------------------------------------
   process sub-style args
------------------------------------------------------------------------- */

void AtomVecHybrid::process_args(int narg, char **arg)
{
  // allocate list of sub-styles as big as possibly needed if no extra args

  styles = new AtomVec *[narg];
  keywords = new char *[narg];

  // allocate each sub-style
  // call process_args() with set of args that are not atom style names
  // use known_style() to determine which args these are

  int dummy;
  int iarg = 0;
  nstyles = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "hybrid") == 0)
      error->all(FLERR, "Atom style hybrid cannot have hybrid as an argument");
    for (int i = 0; i < nstyles; i++)
      if (strcmp(arg[iarg], keywords[i]) == 0)
        error->all(FLERR, "Atom style hybrid cannot use same atom style twice");
    styles[nstyles] = atom->new_avec(arg[iarg], 1, dummy);
    keywords[nstyles] = utils::strdup(arg[iarg]);

    // determine list of arguments for atom style settings
    // by looking for the next known atom style name.

    int jarg = iarg + 1;
    while ((jarg < narg) && !atom->avec_map->count(arg[jarg]) &&
           !lmp->match_style("atom", arg[jarg]))
      jarg++;

    styles[nstyles]->process_args(jarg - iarg - 1, &arg[iarg + 1]);
    iarg = jarg;
    nstyles++;
  }

  // hybrid settings are MAX or MIN of sub-style settings
  // check for both mass_type = 0 and 1, so can warn

  molecular = Atom::ATOMIC;
  maxexchange = 0;

  for (int k = 0; k < nstyles; k++) {
    if ((styles[k]->molecular == Atom::MOLECULAR && molecular == Atom::TEMPLATE) ||
        (styles[k]->molecular == Atom::TEMPLATE && molecular == Atom::MOLECULAR))
      error->all(FLERR, "Cannot mix molecular and molecule template atom styles");
    molecular = MAX(molecular, styles[k]->molecular);

    bonds_allow = MAX(bonds_allow, styles[k]->bonds_allow);
    angles_allow = MAX(angles_allow, styles[k]->angles_allow);
    dihedrals_allow = MAX(dihedrals_allow, styles[k]->dihedrals_allow);
    impropers_allow = MAX(impropers_allow, styles[k]->impropers_allow);
    mass_type = MAX(mass_type, styles[k]->mass_type);
    dipole_type = MAX(dipole_type, styles[k]->dipole_type);
    forceclearflag = MAX(forceclearflag, styles[k]->forceclearflag);
    maxexchange += styles[k]->maxexchange;

    if (styles[k]->molecular == Atom::TEMPLATE) onemols = styles[k]->onemols;
  }

  // issue a warning if both per-type mass and per-atom rmass are defined

  int mass_pertype = 0;
  int mass_peratom = 0;

  for (int k = 0; k < nstyles; k++) {
    if (styles[k]->mass_type == 0) mass_peratom = 1;
    if (styles[k]->mass_type == 1) mass_pertype = 1;
  }

  if (mass_pertype && mass_peratom && comm->me == 0)
    error->warning(FLERR,
                   "Atom style hybrid defines both, per-type "
                   "and per-atom masses; both must be set, but only "
                   "per-atom masses will be used");

  // merge field strings from all sub-styles
  // save concat_grow to check for duplicates of special-case fields

  std::vector<std::string> concat_grow;
  std::vector<std::string> concat_dummy;

  for (int k = 0; k < nstyles; k++) {
    merge_fields(fields_grow, styles[k]->fields_grow, 1, concat_grow);
    merge_fields(fields_copy, styles[k]->fields_copy, 0, concat_dummy);
    merge_fields(fields_comm, styles[k]->fields_comm, 0, concat_dummy);
    merge_fields(fields_comm_vel, styles[k]->fields_comm_vel, 0, concat_dummy);
    merge_fields(fields_reverse, styles[k]->fields_reverse, 0, concat_dummy);
    merge_fields(fields_border, styles[k]->fields_border, 0, concat_dummy);
    merge_fields(fields_border_vel, styles[k]->fields_border_vel, 0, concat_dummy);
    merge_fields(fields_exchange, styles[k]->fields_exchange, 0, concat_dummy);
    merge_fields(fields_restart, styles[k]->fields_restart, 0, concat_dummy);
    merge_fields(fields_create, styles[k]->fields_create, 0, concat_dummy);
    merge_fields(fields_data_atom, styles[k]->fields_data_atom, 0, concat_dummy);
    merge_fields(fields_data_vel, styles[k]->fields_data_vel, 0, concat_dummy);
  }

  // check concat_grow for multiple special-case fields
  // may cause issues with style-specific create_atom() and data_atom() methods
  // issue warnings if appear in multiple sub-styles

  std::vector<std::string> dupfield = {"radius", "rmass"};

  for (const auto &idup : dupfield) {
    if ((comm->me == 0) && (std::count(concat_grow.begin(), concat_grow.end(), idup) > 1))
      error->warning(FLERR,
                     "Per-atom field {} is used in multiple sub-styles; must be used consistently",
                     idup);
  }

  // set bonus_flag if any substyle has bonus data
  // set nstyles_bonus & styles_bonus
  // sum two sizes over contributions from each substyle with bonus data.

  nstyles_bonus = 0;
  for (int k = 0; k < nstyles; k++)
    if (styles[k]->bonus_flag) nstyles_bonus++;

  if (nstyles_bonus) {
    bonus_flag = 1;
    styles_bonus = new AtomVec *[nstyles_bonus];
    nstyles_bonus = 0;
    size_forward_bonus = 0;
    size_border_bonus = 0;
    for (int k = 0; k < nstyles; k++) {
      if (styles[k]->bonus_flag) {
        styles_bonus[nstyles_bonus++] = styles[k];
        size_forward_bonus += styles[k]->size_forward_bonus;
        size_border_bonus += styles[k]->size_border_bonus;
      }
    }
  }

  // parent AtomVec can now operate on merged fields

  setup_fields();
}

/* ---------------------------------------------------------------------- */

void AtomVecHybrid::init()
{
  AtomVec::init();
  for (int k = 0; k < nstyles; k++) styles[k]->init();
}

/* ---------------------------------------------------------------------- */

void AtomVecHybrid::grow_pointers()
{
  for (int k = 0; k < nstyles; k++) styles[k]->grow_pointers();
}

/* ---------------------------------------------------------------------- */

void AtomVecHybrid::force_clear(int n, size_t nbytes)
{
  for (int k = 0; k < nstyles; k++)
    if (styles[k]->forceclearflag) styles[k]->force_clear(n, nbytes);
}

/* ---------------------------------------------------------------------- */

void AtomVecHybrid::copy_bonus(int i, int j, int delflag)
{
  for (int k = 0; k < nstyles_bonus; k++) styles_bonus[k]->copy_bonus(i, j, delflag);
}

/* ---------------------------------------------------------------------- */

void AtomVecHybrid::clear_bonus()
{
  for (int k = 0; k < nstyles_bonus; k++) styles_bonus[k]->clear_bonus();
}

/* ---------------------------------------------------------------------- */

int AtomVecHybrid::pack_comm_bonus(int n, int *list, double *buf)
{
  int m = 0;
  for (int k = 0; k < nstyles_bonus; k++) m += styles_bonus[k]->pack_comm_bonus(n, list, buf);
  return m;
}

/* ---------------------------------------------------------------------- */

void AtomVecHybrid::unpack_comm_bonus(int n, int first, double *buf)
{
  for (int k = 0; k < nstyles_bonus; k++) styles_bonus[k]->unpack_comm_bonus(n, first, buf);
}

/* ---------------------------------------------------------------------- */

int AtomVecHybrid::pack_border_bonus(int n, int *list, double *buf)
{
  int m = 0;
  for (int k = 0; k < nstyles_bonus; k++) m += styles_bonus[k]->pack_border_bonus(n, list, buf);
  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecHybrid::unpack_border_bonus(int n, int first, double *buf)
{
  int m = 0;
  for (int k = 0; k < nstyles_bonus; k++) m += styles_bonus[k]->unpack_border_bonus(n, first, buf);
  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecHybrid::pack_exchange_bonus(int i, double *buf)
{
  int m = 0;
  for (int k = 0; k < nstyles_bonus; k++) m += styles_bonus[k]->pack_exchange_bonus(i, buf);
  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecHybrid::unpack_exchange_bonus(int ilocal, double *buf)
{
  int m = 0;
  for (int k = 0; k < nstyles_bonus; k++) m += styles_bonus[k]->unpack_exchange_bonus(ilocal, buf);
  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecHybrid::size_restart_bonus()
{
  int n = 0;
  for (int k = 0; k < nstyles_bonus; k++) n += styles_bonus[k]->size_restart_bonus();
  return n;
}

/* ---------------------------------------------------------------------- */

int AtomVecHybrid::pack_restart_bonus(int i, double *buf)
{
  int m = 0;
  for (int k = 0; k < nstyles_bonus; k++) m += styles_bonus[k]->pack_restart_bonus(i, buf);
  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecHybrid::unpack_restart_bonus(int ilocal, double *buf)
{
  int m = 0;
  for (int k = 0; k < nstyles_bonus; k++) m += styles_bonus[k]->unpack_restart_bonus(ilocal, buf);
  return m;
}

/* ---------------------------------------------------------------------- */

double AtomVecHybrid::memory_usage_bonus()
{
  double bytes = 0;
  for (int k = 0; k < nstyles_bonus; k++) bytes += styles_bonus[k]->memory_usage_bonus();
  return bytes;
}

/* ----------------------------------------------------------------------
   modify values for AtomVec::pack_restart() to pack
------------------------------------------------------------------------- */

void AtomVecHybrid::pack_restart_pre(int ilocal)
{
  for (int k = 0; k < nstyles; k++) styles[k]->pack_restart_pre(ilocal);
}

/* ----------------------------------------------------------------------
   unmodify values packed by AtomVec::pack_restart()
------------------------------------------------------------------------- */

void AtomVecHybrid::pack_restart_post(int ilocal)
{
  for (int k = 0; k < nstyles; k++) styles[k]->pack_restart_post(ilocal);
}

/* ----------------------------------------------------------------------
   initialize other atom quantities after AtomVec::unpack_restart()
------------------------------------------------------------------------- */

void AtomVecHybrid::unpack_restart_init(int ilocal)
{
  for (int k = 0; k < nstyles; k++) styles[k]->unpack_restart_init(ilocal);
}

/* ----------------------------------------------------------------------
   initialize non-zero atom quantities
------------------------------------------------------------------------- */

void AtomVecHybrid::create_atom_post(int ilocal)
{
  for (int k = 0; k < nstyles; k++) styles[k]->create_atom_post(ilocal);
}

/* ----------------------------------------------------------------------
   modify what AtomVec::data_atom() just unpacked
   or initialize other atom quantities
------------------------------------------------------------------------- */

void AtomVecHybrid::data_atom_post(int ilocal)
{
  for (int k = 0; k < nstyles; k++) styles[k]->data_atom_post(ilocal);
}

/* ----------------------------------------------------------------------
   modify what AtomVec::data_bonds() just unpacked
   or initialize other bond quantities
------------------------------------------------------------------------- */
void AtomVecHybrid::data_bonds_post(int m, int num_bond, tagint atom1, tagint atom2,
                                    tagint id_offset)
{
  for (int k = 0; k < nstyles; k++)
    styles[k]->data_bonds_post(m, num_bond, atom1, atom2, id_offset);
}

/* ----------------------------------------------------------------------
   modify values for AtomVec::pack_data() to pack
------------------------------------------------------------------------- */

void AtomVecHybrid::pack_data_pre(int ilocal)
{
  for (int k = 0; k < nstyles; k++) styles[k]->pack_data_pre(ilocal);
}

/* ----------------------------------------------------------------------
   unmodify values packed by AtomVec::pack_data()
------------------------------------------------------------------------- */

void AtomVecHybrid::pack_data_post(int ilocal)
{
  for (int k = 0; k < nstyles; k++) styles[k]->pack_data_post(ilocal);
}

/* ----------------------------------------------------------------------
   pack bonus info for writing to data file, match flag to sub-style
------------------------------------------------------------------------- */

int AtomVecHybrid::pack_data_bonus(double *buf, int flag)
{
  for (int k = 0; k < nstyles; k++) {
    if (flag == ELLIPSOID && strcmp(keywords[k], "ellipsoid") != 0) continue;
    if (flag == LINE && strcmp(keywords[k], "line") != 0) continue;
    if (flag == TRIANGLE && strcmp(keywords[k], "tri") != 0) continue;
    if (flag == BODY && strcmp(keywords[k], "body") != 0) continue;

    return styles[k]->pack_data_bonus(buf, flag);
  }
  return 0;
}

/* ----------------------------------------------------------------------
   write bonus info to data file, match flag to sub-style
------------------------------------------------------------------------- */

void AtomVecHybrid::write_data_bonus(FILE *fp, int n, double *buf, int flag)
{
  for (int k = 0; k < nstyles; k++) {
    if (flag == ELLIPSOID && strcmp(keywords[k], "ellipsoid") != 0) continue;
    if (flag == LINE && strcmp(keywords[k], "line") != 0) continue;
    if (flag == TRIANGLE && strcmp(keywords[k], "tri") != 0) continue;
    if (flag == BODY && strcmp(keywords[k], "body") != 0) continue;

    styles[k]->write_data_bonus(fp, n, buf, flag);
  }
}

/* ----------------------------------------------------------------------
   assign an index to named atom property and return index
   returned value encodes which sub-style and index returned by sub-style
   return -1 if name is unknown to any sub-styles
------------------------------------------------------------------------- */

int AtomVecHybrid::property_atom(const std::string &name)
{
  for (int k = 0; k < nstyles; k++) {
    int index = styles[k]->property_atom(name);
    if (index >= 0) return index * nstyles + k;
  }
  return -1;
}

/* ----------------------------------------------------------------------
   pack per-atom data into buf for ComputePropertyAtom
   index maps to data specific to this atom style
------------------------------------------------------------------------- */

void AtomVecHybrid::pack_property_atom(int multiindex, double *buf, int nvalues, int groupbit)
{
  int k = multiindex % nstyles;
  int index = multiindex / nstyles;
  styles[k]->pack_property_atom(index, buf, nvalues, groupbit);
}

// ----------------------------------------------------------------------
// internal methods
// ----------------------------------------------------------------------

/* ----------------------------------------------------------------------
   merge fields into root vector and remove duplicate fields
   if concat_flag set, also return concat (w/ duplicates)
   so caller can check for problematic fields
------------------------------------------------------------------------- */

void AtomVecHybrid::merge_fields(std::vector<std::string> &root,
                                 const std::vector<std::string> &fields, int concat_flag,
                                 std::vector<std::string> &concat)
{
  // grow vector with all words combined with dedup and

  for (const auto &field : fields) {
    if (concat_flag) concat.push_back(field);
    if (std::find(root.begin(), root.end(), field) == root.end()) root.push_back(field);
  }
}

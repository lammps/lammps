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

#include "group.h"

#include "atom.h"
#include "comm.h"
#include "compute.h"
#include "domain.h"
#include "dump.h"
#include "error.h"
#include "exceptions.h"
#include "fix.h"
#include "force.h"
#include "input.h"
#include "math_eigen.h"
#include "math_extra.h"
#include "memory.h"
#include "modify.h"
#include "output.h"
#include "region.h"
#include "tokenizer.h"
#include "variable.h"

#include <cmath>
#include <cstring>
#include <map>
#include <utility>

using namespace LAMMPS_NS;

static constexpr int MAX_GROUP = 32;
static constexpr double EPSILON = 1.0e-6;

enum { NONE, TYPE, MOLECULE, ID };
enum { LT, LE, GT, GE, EQ, NEQ, BETWEEN };

static constexpr double BIG = 1.0e20;

/* ----------------------------------------------------------------------
   initialize group memory
------------------------------------------------------------------------- */

Group::Group(LAMMPS *lmp) : Pointers(lmp)
{
  MPI_Comm_rank(world, &me);

  names = new char *[MAX_GROUP];
  bitmask = new int[MAX_GROUP];
  inversemask = new int[MAX_GROUP];
  dynamic = new int[MAX_GROUP];

  for (int i = 0; i < MAX_GROUP; i++) names[i] = nullptr;
  for (int i = 0; i < MAX_GROUP; i++) bitmask[i] = 1 << i;
  for (int i = 0; i < MAX_GROUP; i++) inversemask[i] = bitmask[i] ^ ~0;
  for (int i = 0; i < MAX_GROUP; i++) dynamic[i] = 0;

  // create "all" group

  names[0] = utils::strdup("all");
  ngroup = 1;
}

/* ----------------------------------------------------------------------
   free all memory
------------------------------------------------------------------------- */

Group::~Group()
{
  for (int i = 0; i < MAX_GROUP; i++) delete[] names[i];
  delete[] names;
  delete[] bitmask;
  delete[] inversemask;
  delete[] dynamic;
}

/* ----------------------------------------------------------------------
   assign atoms to a new or existing group
------------------------------------------------------------------------- */

void Group::assign(int narg, char **arg)
{
  int i;

  if (domain->box_exist == 0) error->all(FLERR, "Group command before simulation box is defined");
  if (narg < 2) utils::missing_cmd_args(FLERR, "group", error);

  // delete the group if not being used elsewhere
  // clear mask of each atom assigned to this group

  if (strcmp(arg[1], "delete") == 0) {
    if (narg != 2) error->all(FLERR, "Illegal group delete command: too many arguments");
    int igroup = find(arg[0]);
    if (igroup == -1) error->all(FLERR, "Could not find group delete group ID {}", arg[0]);
    if (igroup == 0) error->all(FLERR, "Cannot delete group all");
    for (const auto &i : modify->get_fix_list())
      if (i->igroup == igroup)
        error->all(FLERR, "Cannot delete group {} currently used by fix ID {}", arg[0], i->id);
    for (const auto &i : modify->get_compute_list())
      if (i->igroup == igroup)
        error->all(FLERR, "Cannot delete group {} currently used by compute ID {}", arg[0], i->id);
    for (const auto &i : output->get_dump_list())
      if (i->igroup == igroup)
        error->all(FLERR, "Cannot delete group {} currently used by dump ID {}", arg[0], i->id);
    if (atom->firstgroupname && strcmp(arg[0], atom->firstgroupname) == 0)
      error->all(FLERR, "Cannot delete group {} currently used by atom_modify first", arg[0]);

    int *mask = atom->mask;
    int nlocal = atom->nlocal;
    int bits = inversemask[igroup];
    for (i = 0; i < nlocal; i++) mask[i] &= bits;

    if (dynamic[igroup]) modify->delete_fix(std::string("GROUP_") + names[igroup]);

    delete[] names[igroup];
    names[igroup] = nullptr;
    dynamic[igroup] = 0;
    ngroup--;

    return;
  }

  // clear the group

  if (strcmp(arg[1], "clear") == 0) {
    int igroup = find(arg[0]);
    if (igroup == -1) error->all(FLERR, "Could not find group clear group ID {}", arg[0]);
    if (igroup == 0) error->all(FLERR, "Cannot clear group all");

    int *mask = atom->mask;
    int nlocal = atom->nlocal;
    int bits = inversemask[igroup];
    for (i = 0; i < nlocal; i++) mask[i] &= bits;

    return;
  }

  // find group in existing list
  // add a new group if igroup = -1

  int igroup = find(arg[0]);
  bool created = false;

  if (igroup == -1) {
    if (ngroup == MAX_GROUP) error->all(FLERR, "Too many groups (max {})", MAX_GROUP);
    igroup = find_unused();
    names[igroup] = utils::strdup(arg[0]);
    ngroup++;
    created = true;
  }

  double **x = atom->x;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int bit = bitmask[igroup];

  try {
    // style = region
    // add to group if atom is in region

    if (strcmp(arg[1], "region") == 0) {

      if (narg != 3) error->all(FLERR, "Illegal group region command");

      auto region = domain->get_region_by_id(arg[2]);
      if (!region) error->all(FLERR, "Region {} for group region does not exist", arg[2]);
      region->init();
      region->prematch();

      for (i = 0; i < nlocal; i++)
        if (region->match(x[i][0], x[i][1], x[i][2])) mask[i] |= bit;

      // create an empty group

    } else if (strcmp(arg[1], "empty") == 0) {

      if (narg != 2) error->all(FLERR, "Illegal group empty command");
      // nothing else to do here

      // style = type, molecule, id
      // add to group if atom matches type/molecule/id or condition

    } else if (strcmp(arg[1], "type") == 0 || strcmp(arg[1], "molecule") == 0 ||
               strcmp(arg[1], "id") == 0) {

      if (narg < 3) utils::missing_cmd_args(FLERR, std::string("group ") + arg[1], error);

      int category = NONE;
      if (strcmp(arg[1], "type") == 0)
        category = TYPE;
      else if (strcmp(arg[1], "molecule") == 0)
        category = MOLECULE;
      else if (strcmp(arg[1], "id") == 0)
        category = ID;

      if ((category == MOLECULE) && (!atom->molecule_flag))
        error->all(FLERR, "Group molecule command requires atom attribute molecule");

      if ((category == ID) && (!atom->tag_enable))
        error->all(FLERR, "Group id command requires atom IDs");

      // args = logical condition

      if (narg > 3 &&
          (strcmp(arg[2], "<") == 0 || strcmp(arg[2], ">") == 0 || strcmp(arg[2], "<=") == 0 ||
           strcmp(arg[2], ">=") == 0 || strcmp(arg[2], "==") == 0 || strcmp(arg[2], "!=") == 0 ||
           strcmp(arg[2], "<>") == 0)) {

        int condition = -1;
        if (strcmp(arg[2], "<") == 0)
          condition = LT;
        else if (strcmp(arg[2], "<=") == 0)
          condition = LE;
        else if (strcmp(arg[2], ">") == 0)
          condition = GT;
        else if (strcmp(arg[2], ">=") == 0)
          condition = GE;
        else if (strcmp(arg[2], "==") == 0)
          condition = EQ;
        else if (strcmp(arg[2], "!=") == 0)
          condition = NEQ;
        else if (strcmp(arg[2], "<>") == 0)
          condition = BETWEEN;
        else
          error->all(FLERR, "Illegal group command");

        tagint bound1, bound2;
        if (category == TYPE)
          bound1 = (tagint) utils::expand_type_int(FLERR, arg[3], Atom::ATOM, lmp);
        else
          bound1 = utils::tnumeric(FLERR, arg[3], false, lmp);
        bound2 = -1;

        if (condition == BETWEEN) {
          if (narg != 5) error->all(FLERR, "Illegal group command");
          if (category == TYPE)
            bound2 = (tagint) utils::expand_type_int(FLERR, arg[4], Atom::ATOM, lmp);
          else
            bound2 = utils::tnumeric(FLERR, arg[4], false, lmp);
        } else if (narg != 4)
          error->all(FLERR, "Illegal group command");

        int *attribute = nullptr;
        tagint *tattribute = nullptr;
        if (category == TYPE)
          attribute = atom->type;
        else if (category == MOLECULE)
          tattribute = atom->molecule;
        else if (category == ID)
          tattribute = atom->tag;

        // add to group if meets condition

        if (attribute) {
          if (condition == LT) {
            for (i = 0; i < nlocal; i++)
              if (attribute[i] < bound1) mask[i] |= bit;
          } else if (condition == LE) {
            for (i = 0; i < nlocal; i++)
              if (attribute[i] <= bound1) mask[i] |= bit;
          } else if (condition == GT) {
            for (i = 0; i < nlocal; i++)
              if (attribute[i] > bound1) mask[i] |= bit;
          } else if (condition == GE) {
            for (i = 0; i < nlocal; i++)
              if (attribute[i] >= bound1) mask[i] |= bit;
          } else if (condition == EQ) {
            for (i = 0; i < nlocal; i++)
              if (attribute[i] == bound1) mask[i] |= bit;
          } else if (condition == NEQ) {
            for (i = 0; i < nlocal; i++)
              if (attribute[i] != bound1) mask[i] |= bit;
          } else if (condition == BETWEEN) {
            for (i = 0; i < nlocal; i++)
              if (attribute[i] >= bound1 && attribute[i] <= bound2) mask[i] |= bit;
          }
        } else {
          if (condition == LT) {
            for (i = 0; i < nlocal; i++)
              if (tattribute[i] < bound1) mask[i] |= bit;
          } else if (condition == LE) {
            for (i = 0; i < nlocal; i++)
              if (tattribute[i] <= bound1) mask[i] |= bit;
          } else if (condition == GT) {
            for (i = 0; i < nlocal; i++)
              if (tattribute[i] > bound1) mask[i] |= bit;
          } else if (condition == GE) {
            for (i = 0; i < nlocal; i++)
              if (tattribute[i] >= bound1) mask[i] |= bit;
          } else if (condition == EQ) {
            for (i = 0; i < nlocal; i++)
              if (tattribute[i] == bound1) mask[i] |= bit;
          } else if (condition == NEQ) {
            for (i = 0; i < nlocal; i++)
              if (tattribute[i] != bound1) mask[i] |= bit;
          } else if (condition == BETWEEN) {
            for (i = 0; i < nlocal; i++)
              if (tattribute[i] >= bound1 && tattribute[i] <= bound2) mask[i] |= bit;
          }
        }

        // args = list of values

      } else {
        int *attribute = nullptr;
        tagint *tattribute = nullptr;
        if (category == TYPE)
          attribute = atom->type;
        else if (category == MOLECULE)
          tattribute = atom->molecule;
        else if (category == ID)
          tattribute = atom->tag;

        char *typestr = nullptr;
        tagint start, stop, delta;

        for (int iarg = 2; iarg < narg; iarg++) {
          delta = 1;
          if (category == TYPE) {
            delete[] typestr;
            typestr = utils::expand_type(FLERR, arg[iarg], Atom::ATOM, lmp);
            if (typestr) stop = start = utils::tnumeric(FLERR, typestr, false, lmp);
          }
          if (typestr == nullptr) {
            try {
              ValueTokenizer values(arg[iarg], ":");
              start = values.next_tagint();
              if (utils::strmatch(arg[iarg], "^-?\\d+$")) {
                stop = start;
              } else if (utils::strmatch(arg[iarg], "^-?\\d+:-?\\d+$")) {
                stop = values.next_tagint();
              } else if (utils::strmatch(arg[iarg], "^-?\\d+:-?\\d+:\\d+$")) {
                stop = values.next_tagint();
                delta = values.next_tagint();
              } else
                throw TokenizerException("Syntax error", "");
            } catch (TokenizerException &e) {
              error->all(FLERR, "Incorrect range string '{}': {}", arg[iarg], e.what());
            }
            if (delta < 1) error->all(FLERR, "Illegal range increment value");
          }

          // add to group if attribute matches value or sequence

          if (attribute) {
            for (i = 0; i < nlocal; i++)
              if (attribute[i] >= start && attribute[i] <= stop &&
                  (attribute[i] - start) % delta == 0)
                mask[i] |= bit;
          } else {
            for (i = 0; i < nlocal; i++)
              if (tattribute[i] >= start && tattribute[i] <= stop &&
                  (tattribute[i] - start) % delta == 0)
                mask[i] |= bit;
          }
        }

        delete[] typestr;
      }

      // style = variable
      // add to group if atom-atyle variable is non-zero

    } else if (strcmp(arg[1], "variable") == 0) {

      int ivar = input->variable->find(arg[2]);
      if (ivar < 0) error->all(FLERR, "Variable name {} for group does not exist", arg[2]);
      if (!input->variable->atomstyle(ivar))
        error->all(FLERR, "Variable {} for group is invalid style", arg[2]);

      double *aflag;

      // aflag = evaluation of per-atom variable

      memory->create(aflag, nlocal, "group:aflag");
      input->variable->compute_atom(ivar, 0, aflag, 1, 0);

      // add to group if per-atom variable evaluated to non-zero

      for (i = 0; i < nlocal; i++)
        if (aflag[i] != 0.0) mask[i] |= bit;

      memory->destroy(aflag);

      // style = include

    } else if (strcmp(arg[1], "include") == 0) {

      if (narg != 3) error->all(FLERR, "Illegal group include command");
      if (strcmp(arg[2], "molecule") == 0) {
        if (!atom->molecule_flag)
          error->all(FLERR, "Group include molecule command requires atom attribute molecule");

        add_molecules(igroup, bit);

      } else
        error->all(FLERR, "Unknown group include keyword {}", arg[2]);

      // style = subtract

    } else if (strcmp(arg[1], "subtract") == 0) {

      if (narg < 4) utils::missing_cmd_args(FLERR, "group subtract", error);

      int length = narg - 2;
      std::vector<int> list(length);

      int jgroup;
      for (int iarg = 2; iarg < narg; iarg++) {
        jgroup = find(arg[iarg]);
        if (jgroup == -1) error->all(FLERR, "Group ID {} does not exist", arg[iarg]);
        if (dynamic[jgroup]) error->all(FLERR, "Cannot subtract dynamic groups");
        list[iarg - 2] = jgroup;
      }

      // add to group if in 1st group in list

      int otherbit = bitmask[list[0]];

      for (i = 0; i < nlocal; i++)
        if (mask[i] & otherbit) mask[i] |= bit;

      // remove atoms if they are in any of the other groups
      // AND with inverse mask removes the atom from group

      int inverse = inversemask[igroup];

      for (int ilist = 1; ilist < length; ilist++) {
        otherbit = bitmask[list[ilist]];
        for (i = 0; i < nlocal; i++)
          if (mask[i] & otherbit) mask[i] &= inverse;
      }

      // style = union

    } else if (strcmp(arg[1], "union") == 0) {

      if (narg < 3) utils::missing_cmd_args(FLERR, "group union", error);

      int length = narg - 2;
      std::vector<int> list(length);

      int jgroup;
      for (int iarg = 2; iarg < narg; iarg++) {
        jgroup = find(arg[iarg]);
        if (jgroup == -1) error->all(FLERR, "Group ID {} does not exist", arg[iarg]);
        if (dynamic[jgroup]) error->all(FLERR, "Cannot union groups from a dynamic group");
        list[iarg - 2] = jgroup;
      }

      // add to group if in any other group in list

      int otherbit;

      for (int ilist = 0; ilist < length; ilist++) {
        otherbit = bitmask[list[ilist]];
        for (i = 0; i < nlocal; i++)
          if (mask[i] & otherbit) mask[i] |= bit;
      }

      // style = intersect

    } else if (strcmp(arg[1], "intersect") == 0) {

      if (narg < 4) utils::missing_cmd_args(FLERR, "group intersect", error);

      int length = narg - 2;
      std::vector<int> list(length);

      int jgroup;
      for (int iarg = 2; iarg < narg; iarg++) {
        jgroup = find(arg[iarg]);
        if (jgroup == -1) error->all(FLERR, "Group ID {} does not exist", arg[iarg]);
        if (dynamic[jgroup]) error->all(FLERR, "Cannot intersect groups using a dynamic group");
        list[iarg - 2] = jgroup;
      }

      // add to group if in all groups in list

      int otherbit, ok, ilist;

      for (i = 0; i < nlocal; i++) {
        ok = 1;
        for (ilist = 0; ilist < length; ilist++) {
          otherbit = bitmask[list[ilist]];
          if ((mask[i] & otherbit) == 0) ok = 0;
        }
        if (ok) mask[i] |= bit;
      }

      // style = dynamic
      // create a new FixGroup to dynamically determine atoms in group

    } else if (strcmp(arg[1], "dynamic") == 0) {

      if (narg < 4) error->all(FLERR, "Illegal group command");
      if (strcmp(arg[0], arg[2]) == 0) error->all(FLERR, "Group dynamic cannot reference itself");
      if (find(arg[2]) < 0)
        error->all(FLERR, "Group dynamic parent group {} does not exist", arg[2]);
      if (igroup == 0) error->all(FLERR, "Group all cannot be made dynamic");

      // if group is already dynamic, delete existing FixGroup

      if (dynamic[igroup]) modify->delete_fix(std::string("GROUP_") + names[igroup]);

      dynamic[igroup] = 1;

      std::string fixcmd = "GROUP_";
      fixcmd += fmt::format("{} {} GROUP", names[igroup], arg[2]);
      for (i = 3; i < narg; i++) fixcmd += std::string(" ") + arg[i];
      modify->add_fix(fixcmd);

      // style = static
      // remove dynamic FixGroup if necessary

    } else if (strcmp(arg[1], "static") == 0) {

      if (narg != 2) error->all(FLERR, "Illegal group static command");

      if (dynamic[igroup]) modify->delete_fix(std::string("GROUP_") + names[igroup]);

      dynamic[igroup] = 0;

      // not a valid group style

    } else
      error->all(FLERR, "Unknown group command keyword: {}", arg[1]);

  } catch (LAMMPSException &e) {
    // undo created group in case of an error
    if (created) {
      delete[] names[igroup];
      names[igroup] = nullptr;
      ngroup--;
    }
    throw e;
  }

  // print stats for changed group

  int n;
  n = 0;
  for (i = 0; i < nlocal; i++)
    if (mask[i] & bit) n++;

  double rlocal = n;
  double all;
  MPI_Allreduce(&rlocal, &all, 1, MPI_DOUBLE, MPI_SUM, world);

  if (me == 0) {
    if (dynamic[igroup])
      utils::logmesg(lmp, "dynamic group {} defined\n", names[igroup]);
    else
      utils::logmesg(lmp, "{:.15g} atoms in group {}\n", all, names[igroup]);
  }
}

/* ----------------------------------------------------------------------
   convenience function to allow assigning to groups from a single string
------------------------------------------------------------------------- */

void Group::assign(const std::string &groupcmd)
{
  auto args = utils::split_words(groupcmd);
  std::vector<char *> newarg(args.size());
  int i = 0;
  for (const auto &arg : args) { newarg[i++] = (char *) arg.c_str(); }
  assign(args.size(), newarg.data());
}

/* ----------------------------------------------------------------------
   add flagged atoms to a new or existing group
------------------------------------------------------------------------- */

void Group::create(const std::string &name, int *flag)
{
  int i;

  // find group in existing list
  // add a new group if igroup = -1

  int igroup = find(name);

  if (igroup == -1) {
    if (ngroup == MAX_GROUP) error->all(FLERR, "Too many groups (max {})", MAX_GROUP);
    igroup = find_unused();
    names[igroup] = utils::strdup(name);
    ngroup++;
  }

  // add atoms to group whose flags are set

  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int bit = bitmask[igroup];

  for (i = 0; i < nlocal; i++)
    if (flag[i]) mask[i] |= bit;
}

/* ----------------------------------------------------------------------
   return group index if name matches existing group, -1 if no such group
------------------------------------------------------------------------- */

int Group::find(const std::string &name)
{
  for (int igroup = 0; igroup < MAX_GROUP; igroup++)
    if (names[igroup] && (name == names[igroup])) return igroup;
  return -1;
}

/* ----------------------------------------------------------------------
   find group with name or create group if it doesn't exist
   return group index
------------------------------------------------------------------------- */

int Group::find_or_create(const char *name)
{
  int igroup = find(name);
  if (igroup >= 0) return igroup;

  if (ngroup == MAX_GROUP) error->all(FLERR, "Too many groups (max {})", MAX_GROUP);
  igroup = find_unused();
  names[igroup] = utils::strdup(name);
  ngroup++;

  return igroup;
}

/* ----------------------------------------------------------------------
   return index of first available group
   should never be called when group limit has been reached
------------------------------------------------------------------------- */

int Group::find_unused()
{
  for (int igroup = 0; igroup < MAX_GROUP; igroup++)
    if (names[igroup] == nullptr) return igroup;
  return -1;
}

/* ----------------------------------------------------------------------
   return group bitmask for given group id. Error out if group is not found.
------------------------------------------------------------------------- */

int Group::get_bitmask_by_id(const std::string &file, int line, const std::string &name,
                             const std::string &caller)
{
  int igroup = 0;
  for (; igroup < MAX_GROUP; ++igroup) {
    if (names[igroup] && (name == names[igroup])) break;
  }
  if (igroup == MAX_GROUP)
    error->all(file, line, "Group ID {} requested by {} does not exist", name, caller);
  return bitmask[igroup];
}

/* ----------------------------------------------------------------------
   add atoms to group that are in same molecules as atoms already in group
   do not include molID = 0
------------------------------------------------------------------------- */

void Group::add_molecules(int /*igroup*/, int bit)
{
  // hash = unique molecule IDs of atoms already in group

  hash = new std::map<tagint, int>();

  tagint *molecule = atom->molecule;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & bit) {
      if (molecule[i] == 0) continue;
      if (hash->find(molecule[i]) == hash->end()) (*hash)[molecule[i]] = 1;
    }

  // list = set of unique molecule IDs for atoms to add
  // pass list to all other procs via comm->ring()

  int n = hash->size();
  tagint *list;
  memory->create(list, n, "group:list");

  n = 0;
  std::map<tagint, int>::iterator pos;
  for (pos = hash->begin(); pos != hash->end(); ++pos) list[n++] = pos->first;

  molbit = bit;
  comm->ring(n, sizeof(tagint), list, 1, molring, nullptr, (void *) this);

  delete hash;
  memory->destroy(list);
}

/* ----------------------------------------------------------------------
   callback from comm->ring()
   cbuf = list of N molecule IDs, put them in hash
   loop over my atoms, if matches molecule ID in hash,
     add atom to group flagged by molbit
------------------------------------------------------------------------- */

void Group::molring(int n, char *cbuf, void *ptr)
{
  auto gptr = (Group *) ptr;
  auto list = (tagint *) cbuf;
  std::map<tagint, int> *hash = gptr->hash;
  int nlocal = gptr->atom->nlocal;
  tagint *molecule = gptr->atom->molecule;
  int *mask = gptr->atom->mask;
  int molbit = gptr->molbit;

  hash->clear();
  for (int i = 0; i < n; i++) (*hash)[list[i]] = 1;

  for (int i = 0; i < nlocal; i++)
    if (hash->find(molecule[i]) != hash->end()) mask[i] |= molbit;
}

/* ----------------------------------------------------------------------
   write group info to a restart file
   only called by proc 0
------------------------------------------------------------------------- */

void Group::write_restart(FILE *fp)
{
  fwrite(&ngroup, sizeof(int), 1, fp);

  // use count to not change restart format with deleted groups
  // remove this on next major release

  int n;
  int count = 0;
  for (int i = 0; i < MAX_GROUP; i++) {
    if (names[i])
      n = strlen(names[i]) + 1;
    else
      n = 0;
    fwrite(&n, sizeof(int), 1, fp);
    if (n) {
      fwrite(names[i], sizeof(char), n, fp);
      count++;
    }
    if (count == ngroup) break;
  }
}

/* ----------------------------------------------------------------------
   read group info from a restart file
   proc 0 reads, bcast to all procs
------------------------------------------------------------------------- */

void Group::read_restart(FILE *fp)
{
  int i, n;

  // delete existing group names
  // atom masks will be overwritten by reading of restart file

  for (i = 0; i < MAX_GROUP; i++) delete[] names[i];

  if (me == 0) utils::sfread(FLERR, &ngroup, sizeof(int), 1, fp, nullptr, error);
  MPI_Bcast(&ngroup, 1, MPI_INT, 0, world);

  // use count to not change restart format with deleted groups
  // remove this on next major release

  int count = 0;
  for (i = 0; i < MAX_GROUP; i++) {
    if (count == ngroup) {
      names[i] = nullptr;
      continue;
    }
    if (me == 0) utils::sfread(FLERR, &n, sizeof(int), 1, fp, nullptr, error);
    MPI_Bcast(&n, 1, MPI_INT, 0, world);
    if (n) {
      names[i] = new char[n];
      if (me == 0) utils::sfread(FLERR, names[i], sizeof(char), n, fp, nullptr, error);
      MPI_Bcast(names[i], n, MPI_CHAR, 0, world);
      count++;
    } else
      names[i] = nullptr;
  }
}

// ----------------------------------------------------------------------
// computations on a group of atoms
// ----------------------------------------------------------------------

/* ----------------------------------------------------------------------
   count atoms in group all
------------------------------------------------------------------------- */

bigint Group::count_all()
{
  bigint nme = atom->nlocal;
  bigint nall;
  MPI_Allreduce(&nme, &nall, 1, MPI_LMP_BIGINT, MPI_SUM, world);
  return nall;
}

/* ----------------------------------------------------------------------
   count atoms in group
------------------------------------------------------------------------- */

bigint Group::count(int igroup)
{
  int groupbit = bitmask[igroup];

  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  int n = 0;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) n++;

  bigint nsingle = n;
  bigint nall;
  MPI_Allreduce(&nsingle, &nall, 1, MPI_LMP_BIGINT, MPI_SUM, world);
  return nall;
}

/* ----------------------------------------------------------------------
   count atoms in group and region
------------------------------------------------------------------------- */

bigint Group::count(int igroup, Region *region)
{
  region->prematch();

  const int groupbit = bitmask[igroup];
  double **x = atom->x;
  int *mask = atom->mask;
  const int nlocal = atom->nlocal;

  int n = 0;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit && region->match(x[i][0], x[i][1], x[i][2])) n++;

  bigint nsingle = n;
  bigint nall;
  MPI_Allreduce(&nsingle, &nall, 1, MPI_LMP_BIGINT, MPI_SUM, world);
  return nall;
}

/* ----------------------------------------------------------------------
   compute the total mass of group of atoms
   use either per-type mass or per-atom rmass
------------------------------------------------------------------------- */

double Group::mass(int igroup)
{
  int groupbit = bitmask[igroup];

  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int *mask = atom->mask;
  int *type = atom->type;
  int nlocal = atom->nlocal;

  double one = 0.0;

  if (rmass) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) one += rmass[i];
  } else {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) one += mass[type[i]];
  }

  double all;
  MPI_Allreduce(&one, &all, 1, MPI_DOUBLE, MPI_SUM, world);
  return all;
}

/* ----------------------------------------------------------------------
   compute the total mass of group of atoms in region
   use either per-type mass or per-atom rmass
------------------------------------------------------------------------- */

double Group::mass(int igroup, Region *region)
{
  int groupbit = bitmask[igroup];
  region->prematch();

  double **x = atom->x;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int *mask = atom->mask;
  int *type = atom->type;
  int nlocal = atom->nlocal;

  double one = 0.0;

  if (rmass) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit && region->match(x[i][0], x[i][1], x[i][2])) one += rmass[i];
  } else {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit && region->match(x[i][0], x[i][1], x[i][2])) one += mass[type[i]];
  }

  double all;
  MPI_Allreduce(&one, &all, 1, MPI_DOUBLE, MPI_SUM, world);
  return all;
}

/* ----------------------------------------------------------------------
   compute the total charge of group of atoms
------------------------------------------------------------------------- */

double Group::charge(int igroup)
{
  int groupbit = bitmask[igroup];

  double *q = atom->q;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double qone = 0.0;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) qone += q[i];

  double qall;
  MPI_Allreduce(&qone, &qall, 1, MPI_DOUBLE, MPI_SUM, world);
  return qall;
}

/* ----------------------------------------------------------------------
   compute the total charge of group of atoms in region
------------------------------------------------------------------------- */

double Group::charge(int igroup, Region *region)
{
  int groupbit = bitmask[igroup];
  region->prematch();

  double **x = atom->x;
  double *q = atom->q;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double qone = 0.0;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit && region->match(x[i][0], x[i][1], x[i][2])) qone += q[i];

  double qall;
  MPI_Allreduce(&qone, &qall, 1, MPI_DOUBLE, MPI_SUM, world);
  return qall;
}

/* ----------------------------------------------------------------------
   compute the coordinate bounds of the group of atoms
   periodic images are not considered, so atoms are NOT unwrapped
------------------------------------------------------------------------- */

void Group::bounds(int igroup, double *minmax)
{
  int groupbit = bitmask[igroup];

  double extent[6];
  extent[0] = extent[2] = extent[4] = BIG;
  extent[1] = extent[3] = extent[5] = -BIG;

  double **x = atom->x;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      extent[0] = MIN(extent[0], x[i][0]);
      extent[1] = MAX(extent[1], x[i][0]);
      extent[2] = MIN(extent[2], x[i][1]);
      extent[3] = MAX(extent[3], x[i][1]);
      extent[4] = MIN(extent[4], x[i][2]);
      extent[5] = MAX(extent[5], x[i][2]);
    }
  }

  // compute extent across all procs
  // flip sign of MIN to do it in one Allreduce MAX
  // set box by extent in shrink-wrapped dims

  extent[0] = -extent[0];
  extent[2] = -extent[2];
  extent[4] = -extent[4];

  MPI_Allreduce(extent, minmax, 6, MPI_DOUBLE, MPI_MAX, world);

  minmax[0] = -minmax[0];
  minmax[2] = -minmax[2];
  minmax[4] = -minmax[4];
}

/* ----------------------------------------------------------------------
   compute the coordinate bounds of the group of atoms in region
   periodic images are not considered, so atoms are NOT unwrapped
------------------------------------------------------------------------- */

void Group::bounds(int igroup, double *minmax, Region *region)
{
  int groupbit = bitmask[igroup];
  region->prematch();

  double extent[6];
  extent[0] = extent[2] = extent[4] = BIG;
  extent[1] = extent[3] = extent[5] = -BIG;

  double **x = atom->x;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit && region->match(x[i][0], x[i][1], x[i][2])) {
      extent[0] = MIN(extent[0], x[i][0]);
      extent[1] = MAX(extent[1], x[i][0]);
      extent[2] = MIN(extent[2], x[i][1]);
      extent[3] = MAX(extent[3], x[i][1]);
      extent[4] = MIN(extent[4], x[i][2]);
      extent[5] = MAX(extent[5], x[i][2]);
    }
  }

  // compute extent across all procs
  // flip sign of MIN to do it in one Allreduce MAX
  // set box by extent in shrink-wrapped dims

  extent[0] = -extent[0];
  extent[2] = -extent[2];
  extent[4] = -extent[4];

  MPI_Allreduce(extent, minmax, 6, MPI_DOUBLE, MPI_MAX, world);

  minmax[0] = -minmax[0];
  minmax[2] = -minmax[2];
  minmax[4] = -minmax[4];
}

/* ----------------------------------------------------------------------
   compute the center-of-mass coords of group of atoms
   masstotal = total mass
   return center-of-mass coords in cm[]
   must unwrap atoms to compute center-of-mass correctly
------------------------------------------------------------------------- */

void Group::xcm(int igroup, double masstotal, double *cm)
{
  int groupbit = bitmask[igroup];

  double **x = atom->x;
  int *mask = atom->mask;
  int *type = atom->type;
  imageint *image = atom->image;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int nlocal = atom->nlocal;

  double cmone[3];
  cmone[0] = cmone[1] = cmone[2] = 0.0;

  double massone;
  double unwrap[3];

  if (rmass) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        massone = rmass[i];
        domain->unmap(x[i], image[i], unwrap);
        cmone[0] += unwrap[0] * massone;
        cmone[1] += unwrap[1] * massone;
        cmone[2] += unwrap[2] * massone;
      }
  } else {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        massone = mass[type[i]];
        domain->unmap(x[i], image[i], unwrap);
        cmone[0] += unwrap[0] * massone;
        cmone[1] += unwrap[1] * massone;
        cmone[2] += unwrap[2] * massone;
      }
  }

  MPI_Allreduce(cmone, cm, 3, MPI_DOUBLE, MPI_SUM, world);
  if (masstotal > 0.0) {
    cm[0] /= masstotal;
    cm[1] /= masstotal;
    cm[2] /= masstotal;
  }
}

/* ----------------------------------------------------------------------
   compute the center-of-mass coords of group of atoms in region
   mastotal = total mass
   return center-of-mass coords in cm[]
   must unwrap atoms to compute center-of-mass correctly
------------------------------------------------------------------------- */

void Group::xcm(int igroup, double masstotal, double *cm, Region *region)
{
  int groupbit = bitmask[igroup];
  region->prematch();

  double **x = atom->x;
  int *mask = atom->mask;
  int *type = atom->type;
  imageint *image = atom->image;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int nlocal = atom->nlocal;

  double cmone[3];
  cmone[0] = cmone[1] = cmone[2] = 0.0;

  double massone;
  double unwrap[3];

  if (rmass) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit && region->match(x[i][0], x[i][1], x[i][2])) {
        massone = rmass[i];
        domain->unmap(x[i], image[i], unwrap);
        cmone[0] += unwrap[0] * massone;
        cmone[1] += unwrap[1] * massone;
        cmone[2] += unwrap[2] * massone;
      }
  } else {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit && region->match(x[i][0], x[i][1], x[i][2])) {
        massone = mass[type[i]];
        domain->unmap(x[i], image[i], unwrap);
        cmone[0] += unwrap[0] * massone;
        cmone[1] += unwrap[1] * massone;
        cmone[2] += unwrap[2] * massone;
      }
  }

  MPI_Allreduce(cmone, cm, 3, MPI_DOUBLE, MPI_SUM, world);
  if (masstotal > 0.0) {
    cm[0] /= masstotal;
    cm[1] /= masstotal;
    cm[2] /= masstotal;
  }
}

/* ----------------------------------------------------------------------
   compute the center-of-mass velocity of group of atoms
   masstotal = total mass
   return center-of-mass velocity in cm[]
------------------------------------------------------------------------- */

void Group::vcm(int igroup, double masstotal, double *cm)
{
  int groupbit = bitmask[igroup];

  double **v = atom->v;
  int *mask = atom->mask;
  int *type = atom->type;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int nlocal = atom->nlocal;

  double p[3], massone;
  p[0] = p[1] = p[2] = 0.0;

  if (rmass) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        massone = rmass[i];
        p[0] += v[i][0] * massone;
        p[1] += v[i][1] * massone;
        p[2] += v[i][2] * massone;
      }
  } else {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        massone = mass[type[i]];
        p[0] += v[i][0] * massone;
        p[1] += v[i][1] * massone;
        p[2] += v[i][2] * massone;
      }
  }

  MPI_Allreduce(p, cm, 3, MPI_DOUBLE, MPI_SUM, world);
  if (masstotal > 0.0) {
    cm[0] /= masstotal;
    cm[1] /= masstotal;
    cm[2] /= masstotal;
  }
}

/* ----------------------------------------------------------------------
   compute the center-of-mass velocity of group of atoms in region
   masstotal = total mass
   return center-of-mass velocity in cm[]
------------------------------------------------------------------------- */

void Group::vcm(int igroup, double masstotal, double *cm, Region *region)
{
  int groupbit = bitmask[igroup];
  region->prematch();

  double **x = atom->x;
  double **v = atom->v;
  int *mask = atom->mask;
  int *type = atom->type;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int nlocal = atom->nlocal;

  double p[3], massone;
  p[0] = p[1] = p[2] = 0.0;

  if (rmass) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit && region->match(x[i][0], x[i][1], x[i][2])) {
        massone = rmass[i];
        p[0] += v[i][0] * massone;
        p[1] += v[i][1] * massone;
        p[2] += v[i][2] * massone;
      }
  } else {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit && region->match(x[i][0], x[i][1], x[i][2])) {
        massone = mass[type[i]];
        p[0] += v[i][0] * massone;
        p[1] += v[i][1] * massone;
        p[2] += v[i][2] * massone;
      }
  }

  MPI_Allreduce(p, cm, 3, MPI_DOUBLE, MPI_SUM, world);
  if (masstotal > 0.0) {
    cm[0] /= masstotal;
    cm[1] /= masstotal;
    cm[2] /= masstotal;
  }
}

/* ----------------------------------------------------------------------
   compute the total force on group of atoms
------------------------------------------------------------------------- */

void Group::fcm(int igroup, double *cm)
{
  int groupbit = bitmask[igroup];

  double **f = atom->f;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double flocal[3];
  flocal[0] = flocal[1] = flocal[2] = 0.0;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      flocal[0] += f[i][0];
      flocal[1] += f[i][1];
      flocal[2] += f[i][2];
    }

  MPI_Allreduce(flocal, cm, 3, MPI_DOUBLE, MPI_SUM, world);
}

/* ----------------------------------------------------------------------
   compute the total force on group of atoms in region
------------------------------------------------------------------------- */

void Group::fcm(int igroup, double *cm, Region *region)
{
  int groupbit = bitmask[igroup];
  region->prematch();

  double **x = atom->x;
  double **f = atom->f;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double flocal[3];
  flocal[0] = flocal[1] = flocal[2] = 0.0;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit && region->match(x[i][0], x[i][1], x[i][2])) {
      flocal[0] += f[i][0];
      flocal[1] += f[i][1];
      flocal[2] += f[i][2];
    }

  MPI_Allreduce(flocal, cm, 3, MPI_DOUBLE, MPI_SUM, world);
}

/* ----------------------------------------------------------------------
   compute the total kinetic energy of group of atoms and return it
------------------------------------------------------------------------- */

double Group::ke(int igroup)
{
  int groupbit = bitmask[igroup];

  double **v = atom->v;
  int *mask = atom->mask;
  int *type = atom->type;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int nlocal = atom->nlocal;

  double one = 0.0;

  if (rmass) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit)
        one += (v[i][0] * v[i][0] + v[i][1] * v[i][1] + v[i][2] * v[i][2]) * rmass[i];
  } else {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit)
        one += (v[i][0] * v[i][0] + v[i][1] * v[i][1] + v[i][2] * v[i][2]) * mass[type[i]];
  }

  double all;
  MPI_Allreduce(&one, &all, 1, MPI_DOUBLE, MPI_SUM, world);
  all *= 0.5 * force->mvv2e;
  return all;
}

/* ----------------------------------------------------------------------
   compute the total kinetic energy of group of atoms in region and return it
------------------------------------------------------------------------- */

double Group::ke(int igroup, Region *region)
{
  int groupbit = bitmask[igroup];
  region->prematch();

  double **x = atom->x;
  double **v = atom->v;
  int *mask = atom->mask;
  int *type = atom->type;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int nlocal = atom->nlocal;

  double one = 0.0;

  if (rmass) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit && region->match(x[i][0], x[i][1], x[i][2]))
        one += (v[i][0] * v[i][0] + v[i][1] * v[i][1] + v[i][2] * v[i][2]) * rmass[i];
  } else {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit && region->match(x[i][0], x[i][1], x[i][2]))
        one += (v[i][0] * v[i][0] + v[i][1] * v[i][1] + v[i][2] * v[i][2]) * mass[type[i]];
  }

  double all;
  MPI_Allreduce(&one, &all, 1, MPI_DOUBLE, MPI_SUM, world);
  all *= 0.5 * force->mvv2e;
  return all;
}

/* ----------------------------------------------------------------------
   compute the radius-of-gyration of group of atoms
   around center-of-mass cm
   must unwrap atoms to compute Rg correctly
------------------------------------------------------------------------- */

double Group::gyration(int igroup, double masstotal, double *cm)
{
  int groupbit = bitmask[igroup];

  double **x = atom->x;
  int *mask = atom->mask;
  int *type = atom->type;
  imageint *image = atom->image;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int nlocal = atom->nlocal;

  double dx, dy, dz, massone;
  double unwrap[3];
  double rg = 0.0;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      domain->unmap(x[i], image[i], unwrap);
      dx = unwrap[0] - cm[0];
      dy = unwrap[1] - cm[1];
      dz = unwrap[2] - cm[2];
      if (rmass)
        massone = rmass[i];
      else
        massone = mass[type[i]];
      rg += (dx * dx + dy * dy + dz * dz) * massone;
    }
  double rg_all;
  MPI_Allreduce(&rg, &rg_all, 1, MPI_DOUBLE, MPI_SUM, world);

  if (masstotal > 0.0) return sqrt(rg_all / masstotal);
  return 0.0;
}

/* ----------------------------------------------------------------------
   compute the radius-of-gyration of group of atoms in region
   around center-of-mass cm
   must unwrap atoms to compute Rg correctly
------------------------------------------------------------------------- */

double Group::gyration(int igroup, double masstotal, double *cm, Region *region)
{
  int groupbit = bitmask[igroup];
  region->prematch();

  double **x = atom->x;
  int *mask = atom->mask;
  int *type = atom->type;
  imageint *image = atom->image;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int nlocal = atom->nlocal;

  double dx, dy, dz, massone;
  double unwrap[3];
  double rg = 0.0;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit && region->match(x[i][0], x[i][1], x[i][2])) {
      domain->unmap(x[i], image[i], unwrap);
      dx = unwrap[0] - cm[0];
      dy = unwrap[1] - cm[1];
      dz = unwrap[2] - cm[2];
      if (rmass)
        massone = rmass[i];
      else
        massone = mass[type[i]];
      rg += (dx * dx + dy * dy + dz * dz) * massone;
    }
  double rg_all;
  MPI_Allreduce(&rg, &rg_all, 1, MPI_DOUBLE, MPI_SUM, world);

  if (masstotal > 0.0) return sqrt(rg_all / masstotal);
  return 0.0;
}

/* ----------------------------------------------------------------------
   compute the angular momentum L (lmom) of group
   around center-of-mass cm
   must unwrap atoms to compute L correctly
------------------------------------------------------------------------- */

void Group::angmom(int igroup, double *cm, double *lmom)
{
  int groupbit = bitmask[igroup];

  double **x = atom->x;
  double **v = atom->v;
  int *mask = atom->mask;
  int *type = atom->type;
  imageint *image = atom->image;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int nlocal = atom->nlocal;

  double dx, dy, dz, massone;
  double unwrap[3];

  double p[3];
  p[0] = p[1] = p[2] = 0.0;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      domain->unmap(x[i], image[i], unwrap);
      dx = unwrap[0] - cm[0];
      dy = unwrap[1] - cm[1];
      dz = unwrap[2] - cm[2];
      if (rmass)
        massone = rmass[i];
      else
        massone = mass[type[i]];
      p[0] += massone * (dy * v[i][2] - dz * v[i][1]);
      p[1] += massone * (dz * v[i][0] - dx * v[i][2]);
      p[2] += massone * (dx * v[i][1] - dy * v[i][0]);
    }

  MPI_Allreduce(p, lmom, 3, MPI_DOUBLE, MPI_SUM, world);
}

/* ----------------------------------------------------------------------
   compute the angular momentum L (lmom) of group of atoms in region
   around center-of-mass cm
   must unwrap atoms to compute L correctly
------------------------------------------------------------------------- */

void Group::angmom(int igroup, double *cm, double *lmom, Region *region)
{
  int groupbit = bitmask[igroup];
  region->prematch();

  double **x = atom->x;
  double **v = atom->v;
  int *mask = atom->mask;
  int *type = atom->type;
  imageint *image = atom->image;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int nlocal = atom->nlocal;

  double dx, dy, dz, massone;
  double unwrap[3];

  double p[3];
  p[0] = p[1] = p[2] = 0.0;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit && region->match(x[i][0], x[i][1], x[i][2])) {
      domain->unmap(x[i], image[i], unwrap);
      dx = unwrap[0] - cm[0];
      dy = unwrap[1] - cm[1];
      dz = unwrap[2] - cm[2];
      if (rmass)
        massone = rmass[i];
      else
        massone = mass[type[i]];
      p[0] += massone * (dy * v[i][2] - dz * v[i][1]);
      p[1] += massone * (dz * v[i][0] - dx * v[i][2]);
      p[2] += massone * (dx * v[i][1] - dy * v[i][0]);
    }

  MPI_Allreduce(p, lmom, 3, MPI_DOUBLE, MPI_SUM, world);
}

/* ----------------------------------------------------------------------
   compute the torque T (tq) on group
   around center-of-mass cm
   must unwrap atoms to compute T correctly
------------------------------------------------------------------------- */

void Group::torque(int igroup, double *cm, double *tq)
{
  int groupbit = bitmask[igroup];

  double **x = atom->x;
  double **f = atom->f;
  int *mask = atom->mask;
  imageint *image = atom->image;
  int nlocal = atom->nlocal;

  double dx, dy, dz;
  double unwrap[3];

  double tlocal[3];
  tlocal[0] = tlocal[1] = tlocal[2] = 0.0;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      domain->unmap(x[i], image[i], unwrap);
      dx = unwrap[0] - cm[0];
      dy = unwrap[1] - cm[1];
      dz = unwrap[2] - cm[2];
      tlocal[0] += dy * f[i][2] - dz * f[i][1];
      tlocal[1] += dz * f[i][0] - dx * f[i][2];
      tlocal[2] += dx * f[i][1] - dy * f[i][0];
    }

  MPI_Allreduce(tlocal, tq, 3, MPI_DOUBLE, MPI_SUM, world);
}

/* ----------------------------------------------------------------------
   compute the torque T (tq) on group of atoms in region
   around center-of-mass cm
   must unwrap atoms to compute T correctly
------------------------------------------------------------------------- */

void Group::torque(int igroup, double *cm, double *tq, Region *region)
{
  int groupbit = bitmask[igroup];
  region->prematch();

  double **x = atom->x;
  double **f = atom->f;
  int *mask = atom->mask;
  imageint *image = atom->image;
  int nlocal = atom->nlocal;

  double dx, dy, dz;
  double unwrap[3];

  double tlocal[3];
  tlocal[0] = tlocal[1] = tlocal[2] = 0.0;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit && region->match(x[i][0], x[i][1], x[i][2])) {
      domain->unmap(x[i], image[i], unwrap);
      dx = unwrap[0] - cm[0];
      dy = unwrap[1] - cm[1];
      dz = unwrap[2] - cm[2];
      tlocal[0] += dy * f[i][2] - dz * f[i][1];
      tlocal[1] += dz * f[i][0] - dx * f[i][2];
      tlocal[2] += dx * f[i][1] - dy * f[i][0];
    }

  MPI_Allreduce(tlocal, tq, 3, MPI_DOUBLE, MPI_SUM, world);
}

/* ----------------------------------------------------------------------
   compute moment of inertia tensor around center-of-mass cm of group
   must unwrap atoms to compute itensor correctly
------------------------------------------------------------------------- */

void Group::inertia(int igroup, double *cm, double itensor[3][3])
{
  int i, j;

  int groupbit = bitmask[igroup];

  double **x = atom->x;
  int *mask = atom->mask;
  int *type = atom->type;
  imageint *image = atom->image;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int nlocal = atom->nlocal;

  double dx, dy, dz, massone;
  double unwrap[3];

  double ione[3][3];
  for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++) ione[i][j] = 0.0;

  for (i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      domain->unmap(x[i], image[i], unwrap);
      dx = unwrap[0] - cm[0];
      dy = unwrap[1] - cm[1];
      dz = unwrap[2] - cm[2];
      if (rmass)
        massone = rmass[i];
      else
        massone = mass[type[i]];
      ione[0][0] += massone * (dy * dy + dz * dz);
      ione[1][1] += massone * (dx * dx + dz * dz);
      ione[2][2] += massone * (dx * dx + dy * dy);
      ione[0][1] -= massone * dx * dy;
      ione[1][2] -= massone * dy * dz;
      ione[0][2] -= massone * dx * dz;
    }
  ione[1][0] = ione[0][1];
  ione[2][1] = ione[1][2];
  ione[2][0] = ione[0][2];

  MPI_Allreduce(&ione[0][0], &itensor[0][0], 9, MPI_DOUBLE, MPI_SUM, world);
}

/* ----------------------------------------------------------------------
   compute moment of inertia tensor around cm of group of atoms in region
   must unwrap atoms to compute itensor correctly
------------------------------------------------------------------------- */

void Group::inertia(int igroup, double *cm, double itensor[3][3], Region *region)
{
  int i, j;

  int groupbit = bitmask[igroup];
  region->prematch();

  double **x = atom->x;
  int *mask = atom->mask;
  int *type = atom->type;
  imageint *image = atom->image;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int nlocal = atom->nlocal;

  double dx, dy, dz, massone;
  double unwrap[3];

  double ione[3][3];
  for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++) ione[i][j] = 0.0;

  for (i = 0; i < nlocal; i++)
    if (mask[i] & groupbit && region->match(x[i][0], x[i][1], x[i][2])) {
      domain->unmap(x[i], image[i], unwrap);
      dx = unwrap[0] - cm[0];
      dy = unwrap[1] - cm[1];
      dz = unwrap[2] - cm[2];
      if (rmass)
        massone = rmass[i];
      else
        massone = mass[type[i]];
      ione[0][0] += massone * (dy * dy + dz * dz);
      ione[1][1] += massone * (dx * dx + dz * dz);
      ione[2][2] += massone * (dx * dx + dy * dy);
      ione[0][1] -= massone * dx * dy;
      ione[1][2] -= massone * dy * dz;
      ione[0][2] -= massone * dx * dz;
    }
  ione[1][0] = ione[0][1];
  ione[2][1] = ione[1][2];
  ione[2][0] = ione[0][2];

  MPI_Allreduce(&ione[0][0], &itensor[0][0], 9, MPI_DOUBLE, MPI_SUM, world);
}

/* ----------------------------------------------------------------------
   compute angular velocity omega from L and I
------------------------------------------------------------------------- */

void Group::omega(double *angmom, double inertia[3][3], double *w)
{
  double idiag[3], ex[3], ey[3], ez[3], cross[3];
  double evectors[3][3], inverse[3][3];

  // determinant = triple product of rows of inertia matrix

  double determinant = inertia[0][0] * inertia[1][1] * inertia[2][2] +
      inertia[0][1] * inertia[1][2] * inertia[2][0] +
      inertia[0][2] * inertia[1][0] * inertia[2][1] -
      inertia[0][0] * inertia[1][2] * inertia[2][1] -
      inertia[0][1] * inertia[1][0] * inertia[2][2] - inertia[2][0] * inertia[1][1] * inertia[0][2];

  // non-singular I matrix
  // use L = Iw, inverting I to solve for w
  // this should give exact zeroing of angular momentum by velocity command

  if (determinant > EPSILON) {

    inverse[0][0] = inertia[1][1] * inertia[2][2] - inertia[1][2] * inertia[2][1];
    inverse[0][1] = -(inertia[0][1] * inertia[2][2] - inertia[0][2] * inertia[2][1]);
    inverse[0][2] = inertia[0][1] * inertia[1][2] - inertia[0][2] * inertia[1][1];

    inverse[1][0] = -(inertia[1][0] * inertia[2][2] - inertia[1][2] * inertia[2][0]);
    inverse[1][1] = inertia[0][0] * inertia[2][2] - inertia[0][2] * inertia[2][0];
    inverse[1][2] = -(inertia[0][0] * inertia[1][2] - inertia[0][2] * inertia[1][0]);

    inverse[2][0] = inertia[1][0] * inertia[2][1] - inertia[1][1] * inertia[2][0];
    inverse[2][1] = -(inertia[0][0] * inertia[2][1] - inertia[0][1] * inertia[2][0]);
    inverse[2][2] = inertia[0][0] * inertia[1][1] - inertia[0][1] * inertia[1][0];

    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++) inverse[i][j] /= determinant;

    w[0] = inverse[0][0] * angmom[0] + inverse[0][1] * angmom[1] + inverse[0][2] * angmom[2];
    w[1] = inverse[1][0] * angmom[0] + inverse[1][1] * angmom[1] + inverse[1][2] * angmom[2];
    w[2] = inverse[2][0] * angmom[0] + inverse[2][1] * angmom[1] + inverse[2][2] * angmom[2];

    // handle (nearly) singular I matrix
    // typically due to 2-atom group or linear molecule
    // use jacobi3() and angmom_to_omega() to calculate valid omega
    // less exact answer than matrix inversion, due to iterative Jacobi method

  } else {
    int ierror = MathEigen::jacobi3(inertia, idiag, evectors);
    if (ierror) error->all(FLERR, "Insufficient Jacobi rotations for group::omega");

    ex[0] = evectors[0][0];
    ex[1] = evectors[1][0];
    ex[2] = evectors[2][0];
    ey[0] = evectors[0][1];
    ey[1] = evectors[1][1];
    ey[2] = evectors[2][1];
    ez[0] = evectors[0][2];
    ez[1] = evectors[1][2];
    ez[2] = evectors[2][2];

    // enforce 3 evectors as a right-handed coordinate system
    // flip 3rd vector if needed

    MathExtra::cross3(ex, ey, cross);
    if (MathExtra::dot3(cross, ez) < 0.0) MathExtra::negate3(ez);

    // if any principal moment < scaled EPSILON, set to 0.0

    double max;
    max = MAX(idiag[0], idiag[1]);
    max = MAX(max, idiag[2]);

    if (idiag[0] < EPSILON * max) idiag[0] = 0.0;
    if (idiag[1] < EPSILON * max) idiag[1] = 0.0;
    if (idiag[2] < EPSILON * max) idiag[2] = 0.0;

    // calculate omega using diagonalized inertia matrix

    MathExtra::angmom_to_omega(angmom, ex, ey, ez, idiag, w);
  }
}

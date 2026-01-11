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

#include "delete_atoms.h"

#include "atom.h"
#include "atom_vec.h"
#include "atom_vec_body.h"
#include "atom_vec_ellipsoid.h"
#include "atom_vec_line.h"
#include "atom_vec_tri.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "fix_bond_history.h"
#include "force.h"
#include "group.h"
#include "input.h"
#include "memory.h"
#include "modify.h"
#include "molecule.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "random_mars.h"
#include "region.h"
#include "variable.h"

#include <algorithm>
#include <cstring>
#include <numeric>
#include <set>
#include <unordered_set>
#include <utility>
#include <vector>

using namespace LAMMPS_NS;

enum { UNKNOWN, FRACTION, COUNT };
enum { TAGADD = 100, TAGSET };

// compare two tags but consider tag == 0 very large so zeroes move to the end
static bool tagcomp(tagint a, tagint b)
{
  if (a == 0) return false;
  if (b == 0) return true;
  return (a < b);
}

/* ---------------------------------------------------------------------- */

DeleteAtoms::DeleteAtoms(LAMMPS *lmp) :
    Command(lmp), dlist(nullptr), tagproc(nullptr), newtags(nullptr)
{
}

/* ---------------------------------------------------------------------- */

void DeleteAtoms::command(int narg, char **arg)
{
  if (domain->box_exist == 0)
    error->all(FLERR, Error::NOLASTLINE,
               "Delete_atoms command before simulation box is defined" + utils::errorurl(33));
  if (narg < 1) utils::missing_cmd_args(FLERR, "delete_atoms", error);
  if (atom->tag_enable == 0)
    error->all(FLERR, Error::COMMAND, "Cannot use delete_atoms unless atoms have IDs");

  // store state before delete

  bigint natoms_previous = atom->natoms;
  bigint nbonds_previous = atom->nbonds;
  bigint nangles_previous = atom->nangles;
  bigint ndihedrals_previous = atom->ndihedrals;
  bigint nimpropers_previous = atom->nimpropers;

  // flag atoms for deletion

  allflag = 0;

  if (strcmp(arg[0], "group") == 0)
    delete_group(narg, arg);
  else if (strcmp(arg[0], "region") == 0)
    delete_region(narg, arg);
  else if (strcmp(arg[0], "overlap") == 0)
    delete_overlap(narg, arg);
  else if (strcmp(arg[0], "random") == 0)
    delete_random(narg, arg);
  // deprecated porosity option, now included in new partial option
  else if (strcmp(arg[0], "porosity") == 0) {
    error->all(FLERR, Error::ARGZERO,
               "The delete_atoms 'porosity' keyword has been removed.\n"
               "Please use: delete_atoms random fraction frac exact group-ID region-ID seed\n");
  } else if (strcmp(arg[0], "variable") == 0)
    delete_variable(narg, arg);
  else
    error->all(FLERR, Error::ARGZERO, "Unknown delete_atoms sub-command: {}", arg[0]);

  if (allflag) {
    int igroupbit = group->get_bitmask_by_id(FLERR, "all", "delete_atoms");
    if (modify->check_rigid_group_overlap(igroupbit))
      error->warning(FLERR, "Attempting to delete atoms in rigid bodies");
  } else {
    if (modify->check_rigid_list_overlap(dlist))
      if (comm->me == 0) error->warning(FLERR, "Attempting to delete atoms in rigid bodies");
  }

  // if allflag = 1, just reset atom->nlocal
  // else delete atoms one by one

  if (allflag)
    atom->nlocal = 0;
  else {

    // optionally delete additional bonds or atoms in molecules

    if (bond_flag) delete_bond();
    if (mol_flag) delete_molecule();

    // delete local atoms flagged in dlist
    // reset nlocal

    AtomVec *avec = atom->avec;
    int nlocal = atom->nlocal;

    int i = 0;
    while (i < nlocal) {
      if (dlist[i]) {
        avec->copy(nlocal - 1, i, 1);
        dlist[i] = dlist[nlocal - 1];
        nlocal--;
      } else
        i++;
    }

    atom->nlocal = nlocal;
    memory->destroy(dlist);
  }

  // if non-molecular system and compress flag set:
  // reset atom tags to be contiguous
  // set all atom IDs to 0, call tag_extend()
  // for molecular system call reset_atoms id, unless there is a fix that stores atom IDs.

  // if condense flag set, collect atom IDs for all atoms in array distributed across all ranks
  // sort complete array according to tag value, tags that are zero and determine new tag.
  // build new array with tags for all local atoms and then apply.

  if (compress_flag) {

    if (atom->molecular == Atom::ATOMIC) {
      tagint *tag = atom->tag;
      int nlocal = atom->nlocal;
      for (int i = 0; i < nlocal; i++) tag[i] = 0;
      atom->tag_extend();
    } else {
      bool can_compress = true;
      for (const auto &ifix : modify->get_fix_list())
        if (ifix->stores_ids) can_compress = false;

      if (can_compress) {
        input->one("reset_atoms id");
      } else {
        if (comm->me == 0)
          error->warning(FLERR, "Ignoring 'compress yes' because of a fix storing atom IDs");
      }
    }
  } else if (condense_flag)
    condense_tags();

  // reset atom->natoms and also topology counts

  bigint nblocal = atom->nlocal;
  MPI_Allreduce(&nblocal, &atom->natoms, 1, MPI_LMP_BIGINT, MPI_SUM, world);

  // reset bonus data counts

  auto *avec_ellipsoid = dynamic_cast<AtomVecEllipsoid *>(atom->style_match("ellipsoid"));
  auto *avec_line = dynamic_cast<AtomVecLine *>(atom->style_match("line"));
  auto *avec_tri = dynamic_cast<AtomVecTri *>(atom->style_match("tri"));
  auto *avec_body = dynamic_cast<AtomVecBody *>(atom->style_match("body"));
  bigint nlocal_bonus;

  if (atom->nellipsoids > 0) {
    nlocal_bonus = avec_ellipsoid->nlocal_bonus;
    MPI_Allreduce(&nlocal_bonus, &atom->nellipsoids, 1, MPI_LMP_BIGINT, MPI_SUM, world);
  }
  if (atom->nlines > 0) {
    nlocal_bonus = avec_line->nlocal_bonus;
    MPI_Allreduce(&nlocal_bonus, &atom->nlines, 1, MPI_LMP_BIGINT, MPI_SUM, world);
  }
  if (atom->ntris > 0) {
    nlocal_bonus = avec_tri->nlocal_bonus;
    MPI_Allreduce(&nlocal_bonus, &atom->ntris, 1, MPI_LMP_BIGINT, MPI_SUM, world);
  }
  if (atom->nbodies > 0) {
    nlocal_bonus = avec_body->nlocal_bonus;
    MPI_Allreduce(&nlocal_bonus, &atom->nbodies, 1, MPI_LMP_BIGINT, MPI_SUM, world);
  }

  // reset atom->map if it exists
  // set nghost to 0 so old ghosts of deleted atoms won't be mapped

  if (atom->map_style != Atom::MAP_NONE) {
    atom->nghost = 0;
    atom->map_init();
    atom->map_set();
  }

  recount_topology();

  // print before and after atom and topology counts

  bigint ndelete = natoms_previous - atom->natoms;
  bigint ndelete_bonds = nbonds_previous - atom->nbonds;
  bigint ndelete_angles = nangles_previous - atom->nangles;
  bigint ndelete_dihedrals = ndihedrals_previous - atom->ndihedrals;
  bigint ndelete_impropers = nimpropers_previous - atom->nimpropers;

  if (comm->me == 0) {
    std::string mesg = fmt::format("Deleted {} atoms, new total = {}\n", ndelete, atom->natoms);
    if (bond_flag || mol_flag) {
      if (nbonds_previous)
        mesg += fmt::format("Deleted {} bonds, new total = {}\n", ndelete_bonds, atom->nbonds);
      if (nangles_previous)
        mesg += fmt::format("Deleted {} angles, new total = {}\n", ndelete_angles, atom->nangles);
      if (ndihedrals_previous)
        mesg += fmt::format("Deleted {} dihedrals, new total = {}\n", ndelete_dihedrals,
                            atom->ndihedrals);
      if (nimpropers_previous)
        mesg += fmt::format("Deleted {} impropers, new total = {}\n", ndelete_impropers,
                            atom->nimpropers);
    }
    utils::logmesg(lmp, mesg);
  }
}

/* ----------------------------------------------------------------------
   delete all atoms in group, group will still exist
------------------------------------------------------------------------- */

void DeleteAtoms::delete_group(int narg, char **arg)
{
  if (narg < 2) utils::missing_cmd_args(FLERR, "delete_atoms group", error);

  int groupbit = group->get_bitmask_by_id(FLERR, arg[1], "delete_atoms");
  options(narg - 2, &arg[2]);

  // check for special case of group = all

  if (strcmp(arg[1], "all") == 0) {
    allflag = 1;
    return;
  }

  // allocate and initialize deletion list

  int nlocal = atom->nlocal;
  memory->create(dlist, nlocal, "delete_atoms:dlist");
  for (int i = 0; i < nlocal; i++) dlist[i] = 0;

  int *mask = atom->mask;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) dlist[i] = 1;
}

/* ----------------------------------------------------------------------
   delete all atoms in region
------------------------------------------------------------------------- */

void DeleteAtoms::delete_region(int narg, char **arg)
{
  if (narg < 2) utils::missing_cmd_args(FLERR, "delete_atoms region", error);

  auto *iregion = domain->get_region_by_id(arg[1]);
  if (!iregion) error->all(FLERR, "Could not find delete_atoms region ID {}", arg[1]);
  iregion->prematch();

  options(narg - 2, &arg[2]);

  // allocate and initialize deletion list

  int nlocal = atom->nlocal;
  memory->create(dlist, nlocal, "delete_atoms:dlist");
  for (int i = 0; i < nlocal; i++) dlist[i] = 0;

  double **x = atom->x;

  for (int i = 0; i < nlocal; i++)
    if (iregion->match(x[i][0], x[i][1], x[i][2])) dlist[i] = 1;
}

/* ----------------------------------------------------------------------
   delete atoms so there are no pairs within cutoff
   which atoms are deleted depends on ordering of atoms within proc
   deletions can vary with processor count
   no guarantee that minimium number of atoms will be deleted
------------------------------------------------------------------------- */

void DeleteAtoms::delete_overlap(int narg, char **arg)
{
  if (narg < 4) utils::missing_cmd_args(FLERR, "delete_atoms overlap", error);

  // read args

  const double cut = utils::numeric(FLERR, arg[1], false, lmp);
  const double cutsq = cut * cut;
  const int group1bit = group->get_bitmask_by_id(FLERR, arg[2], "delete_atoms");
  const int group2bit = group->get_bitmask_by_id(FLERR, arg[3], "delete_atoms");

  options(narg - 4, &arg[4]);

  if (comm->me == 0) utils::logmesg(lmp, "System init for delete_atoms ...\n");

  // request a full neighbor list for use by this command

  neighbor->add_request(this, "delete_atoms", NeighConst::REQ_FULL);

  // init entire system since comm->borders and neighbor->build is done
  // comm::init needs neighbor::init needs pair::init needs kspace::init, etc

  lmp->init();

  // error check on cutoff
  // if no pair style, neighbor list will be empty

  if (force->pair == nullptr) error->all(FLERR, "Delete_atoms requires a pair style be defined");
  if (cut > neighbor->cutneighmax) error->all(FLERR, "Delete_atoms cutoff > max neighbor cutoff");
  if ((cut > neighbor->cutneighmin) && comm->me == 0)
    error->warning(FLERR, "Delete_atoms cutoff > minimum neighbor cutoff");

  // setup domain, communication and neighboring
  // acquire ghosts and build standard neighbor lists

  if (domain->triclinic) domain->x2lamda(atom->nlocal);
  domain->pbc();
  domain->reset_box();
  comm->setup();
  if (neighbor->style) neighbor->setup_bins();
  comm->exchange();
  comm->borders();
  if (domain->triclinic) domain->lamda2x(atom->nlocal + atom->nghost);
  neighbor->build(1);

  // build neighbor list this command needs based on the earlier request

  auto *list = neighbor->find_list(this);
  neighbor->build_one(list);

  // allocate and initialize deletion list
  // must be after exchange potentially changes nlocal

  int nlocal = atom->nlocal;
  memory->create(dlist, nlocal, "delete_atoms:dlist");
  for (int i = 0; i < nlocal; i++) dlist[i] = 0;

  // double loop over owned atoms and their full neighbor list
  // at end of loop, there are no more overlaps
  // only ever delete owned atom I in I loop iteration, never J even if owned

  tagint *tag = atom->tag;
  int *mask = atom->mask;
  double **x = atom->x;
  double *special_coul = force->special_coul;
  double *special_lj = force->special_lj;

  int i, j, ii, jj, inum, jnum;
  double xtmp, ytmp, ztmp, delx, dely, delz, rsq;
  int *ilist, *jlist, *numneigh, **firstneigh;
  double factor_lj, factor_coul;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    if (!(mask[i] & (group1bit | group2bit))) continue;
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      factor_coul = special_coul[sbmask(j)];
      j &= NEIGHMASK;
      if (!(mask[j] & (group1bit | group2bit))) continue;

      // if both weighting factors are 0, skip this pair
      // could be 0 and still be in neigh list for long-range Coulombics
      // want consistency with non-charged pairs which wouldn't be in list

      if (factor_lj == 0.0 && factor_coul == 0.0) continue;

      // only consider deletion if I,J distance < cutoff
      // compute rsq identically on both I,J loop iterations
      // ignoring possibility that I,J tags are equal

      if (tag[i] < tag[j]) {
        delx = xtmp - x[j][0];
        dely = ytmp - x[j][1];
        delz = ztmp - x[j][2];
      } else {
        delx = x[j][0] - xtmp;
        dely = x[j][1] - ytmp;
        delz = x[j][2] - ztmp;
      }
      rsq = delx * delx + dely * dely + delz * delz;
      if (rsq >= cutsq) continue;

      // only consider deletion if I,J are in groups 1,2 respectively
      // true whether J is owned or ghost atom

      if (!(mask[i] & group1bit)) continue;
      if (!(mask[j] & group2bit)) continue;

      // J is owned atom:
      //   delete atom I if atom J has not already been deleted
      // J is ghost atom:
      //   delete atom I if J,I is not a candidate deletion pair
      //     due to being in groups 1,2 respectively
      //   if they are candidate pair, then either:
      //      another proc owns J and could delete J
      //      J is a ghost of another of my owned atoms, and I could delete J
      //   test on tags of I,J ensures that only I or J is deleted

      if (j < nlocal) {
        if (dlist[j]) continue;
      } else if ((mask[i] & group2bit) && (mask[j] & group1bit)) {
        if (tag[i] > tag[j]) continue;
      }

      dlist[i] = 1;
      break;
    }
  }
  neighbor->init();
}

/* ----------------------------------------------------------------------
   delete specified portion of atoms within group and/or region
------------------------------------------------------------------------- */

void DeleteAtoms::delete_random(int narg, char **arg)
{
  if (narg < 7) utils::missing_cmd_args(FLERR, "delete_atoms random", error);

  int random_style = UNKNOWN;
  bool exactflag = false;
  bool errorflag = false;
  bigint ncount = 0;
  double fraction = 0.0;

  if (strcmp(arg[1], "fraction") == 0) {
    random_style = FRACTION;
    fraction = utils::numeric(FLERR, arg[2], false, lmp);
    exactflag = utils::logical(FLERR, arg[3], false, lmp);
    if (fraction < 0.0 || fraction > 1.0)
      error->all(FLERR, "Delete_atoms random fraction has invalid value: {}", fraction);
  } else if (strcmp(arg[1], "count") == 0) {
    random_style = COUNT;
    ncount = utils::bnumeric(FLERR, arg[2], false, lmp);
    errorflag = utils::logical(FLERR, arg[3], false, lmp);
    if (ncount < 0) error->all(FLERR, "Delete_atoms random count has invalid value: {}", ncount);
    exactflag = true;
  } else {
    error->all(FLERR, "Unknown delete_atoms random style: {}", arg[1]);
  }

  int groupbit = group->get_bitmask_by_id(FLERR, arg[4], "delete_atoms");
  auto *region = domain->get_region_by_id(arg[5]);
  if (!region && (strcmp(arg[5], "NULL") != 0))
    error->all(FLERR, "Could not find delete_atoms random region ID {}", arg[5]);

  int seed = utils::inumeric(FLERR, arg[6], false, lmp);
  options(narg - 7, &arg[7]);

  auto *ranmars = new RanMars(lmp, seed + comm->me);

  // allocate and initialize deletion list

  int nlocal = atom->nlocal;
  memory->create(dlist, nlocal, "delete_atoms:dlist");
  for (int i = 0; i < nlocal; i++) dlist[i] = 0;

  // setup

  double **x = atom->x;
  int *mask = atom->mask;

  if (region) region->prematch();

  // delete approximate fraction of atoms in both group and region

  if (random_style == FRACTION && !exactflag) {
    for (int i = 0; i < nlocal; i++) {
      if (!(mask[i] & groupbit)) continue;
      if (region && !region->match(x[i][0], x[i][1], x[i][2])) continue;
      if (ranmars->uniform() <= fraction) dlist[i] = 1;
    }

    // delete exact fraction or count of atoms in both group and region

  } else {
    double **x = atom->x;
    int *mask = atom->mask;

    // count = number of atoms this proc owns in both group and region

    int count = 0;
    for (int i = 0; i < nlocal; i++) {
      if (!(mask[i] & groupbit)) continue;
      if (region && !region->match(x[i][0], x[i][1], x[i][2])) continue;
      count++;
    }

    // convert specified fraction to ncount

    bigint bcount = count;
    bigint allcount;
    MPI_Allreduce(&bcount, &allcount, 1, MPI_LMP_BIGINT, MPI_SUM, world);

    if (random_style == FRACTION) {
      ncount = static_cast<bigint>(fraction * allcount);
    } else if (random_style == COUNT) {
      if (ncount > allcount) {
        auto mesg = fmt::format("Delete_atoms count of {} exceeds number of eligible atoms {}",
                                ncount, allcount);
        ncount = allcount;
        if (errorflag) {
          error->all(FLERR, mesg);
        } else {
          if (comm->me == 0) error->warning(FLERR, mesg);
        }
      }
    }

    // make selection

    int *flag = memory->create(flag, count, "delete_atoms:flag");
    int *work = memory->create(work, count, "delete_atoms:work");

    ranmars->select_subset(ncount, count, flag, work);

    // set dlist for atom indices in flag
    // flag vector from select_subset() is only for eligible atoms

    int j = 0;
    for (int i = 0; i < nlocal; i++) {
      if (!(mask[i] & groupbit)) continue;
      if (region && !region->match(x[i][0], x[i][1], x[i][2])) continue;
      if (flag[j]) dlist[i] = 1;
      j++;
    }

    memory->destroy(flag);
    memory->destroy(work);
  }

  // delete RN generator

  delete ranmars;
}

/* ----------------------------------------------------------------------
   delete all as flagged by non-zero atom style variable
------------------------------------------------------------------------- */

void DeleteAtoms::delete_variable(int narg, char **arg)
{
  if (narg < 2) utils::missing_cmd_args(FLERR, "delete_atoms variable", error);

  int ivar = input->variable->find(arg[1]);
  if (ivar < 0) error->all(FLERR, "Variable name {} for delete_atoms does not exist", arg[1]);
  if (!input->variable->atomstyle(ivar))
    error->all(FLERR, "Variable {} for delete_atoms is invalid style", arg[1]);

  // consume remaining options

  options(narg - 2, &arg[2]);

  // aflag = evaluation of per-atom variable

  const int nlocal = atom->nlocal;
  double *aflag;
  memory->create(dlist, nlocal, "delete_atoms:dlist");
  memory->create(aflag, nlocal, "group:aflag");
  input->variable->compute_atom(ivar, 0, aflag, 1, 0);

  // delete if per-atom variable evaluated to non-zero

  for (int i = 0; i < nlocal; i++) dlist[i] = (aflag[i] == 0.0) ? 0 : 1;

  memory->destroy(aflag);
}

/* ----------------------------------------------------------------------
   delete all topology interactions that include deleted atoms
------------------------------------------------------------------------- */

void DeleteAtoms::delete_bond()
{
  // list = set of unique molecule IDs from which I deleted atoms
  // pass list to all other procs via comm->ring()

  tagint *tag = atom->tag;
  int nlocal = atom->nlocal;

  int n = 0;
  for (int i = 0; i < nlocal; i++)
    if (dlist[i]) n++;
  tagint *list;
  memory->create(list, n, "delete_atoms:list");

  n = 0;
  for (int i = 0; i < nlocal; i++)
    if (dlist[i]) list[n++] = tag[i];

  comm->ring(n, sizeof(tagint), list, 1, bondring, nullptr, (void *) this);

  memory->destroy(list);
}

/* ----------------------------------------------------------------------
   delete atoms in molecules with any deletions
   use dlist marked with atom deletions, and mark additional atoms
   do not include molID = 0
------------------------------------------------------------------------- */

void DeleteAtoms::delete_molecule()
{
  // hash = unique molecule IDs from which I deleted atoms

  std::set<tagint> hash;

  tagint *molecule = atom->molecule;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (dlist[i] && molecule[i] != 0) hash.insert(molecule[i]);
  }

  // list = set of unique molecule IDs from which I deleted atoms
  // pass list to all other procs via comm->ring()

  auto n = hash.size();
  tagint *list;
  memory->create(list, n, "delete_atoms:list");

  n = 0;
  for (const auto pos : hash) list[n++] = pos;

  comm->ring(n, sizeof(tagint), list, 1, molring, nullptr, (void *) this);

  memory->destroy(list);
}

/* ----------------------------------------------------------------------
   update bond,angle,etc counts
   different for atom->molecular = Atom::MOLECULAR or Atom::TEMPLATE
------------------------------------------------------------------------- */

void DeleteAtoms::recount_topology()
{
  bigint nbonds = 0;
  bigint nangles = 0;
  bigint ndihedrals = 0;
  bigint nimpropers = 0;

  if (atom->molecular == Atom::MOLECULAR) {
    int *num_bond = atom->num_bond;
    int *num_angle = atom->num_angle;
    int *num_dihedral = atom->num_dihedral;
    int *num_improper = atom->num_improper;
    int nlocal = atom->nlocal;

    for (int i = 0; i < nlocal; i++) {
      if (num_bond) nbonds += num_bond[i];
      if (num_angle) nangles += num_angle[i];
      if (num_dihedral) ndihedrals += num_dihedral[i];
      if (num_improper) nimpropers += num_improper[i];
    }

  } else if (atom->molecular == Atom::TEMPLATE) {
    Molecule **onemols = atom->avec->onemols;
    int *molindex = atom->molindex;
    int *molatom = atom->molatom;
    int nlocal = atom->nlocal;

    int imol, iatom;

    for (int i = 0; i < nlocal; i++) {
      imol = molindex[i];
      iatom = molatom[i];
      if (imol < 0) continue;
      nbonds += onemols[imol]->num_bond[iatom];
      nangles += onemols[imol]->num_angle[iatom];
      ndihedrals += onemols[imol]->num_dihedral[iatom];
      nimpropers += onemols[imol]->num_improper[iatom];
    }
  }

  if (atom->avec->bonds_allow) {
    MPI_Allreduce(&nbonds, &atom->nbonds, 1, MPI_LMP_BIGINT, MPI_SUM, world);
    if (!force->newton_bond) atom->nbonds /= 2;
  }
  if (atom->avec->angles_allow) {
    MPI_Allreduce(&nangles, &atom->nangles, 1, MPI_LMP_BIGINT, MPI_SUM, world);
    if (!force->newton_bond) atom->nangles /= 3;
  }
  if (atom->avec->dihedrals_allow) {
    MPI_Allreduce(&ndihedrals, &atom->ndihedrals, 1, MPI_LMP_BIGINT, MPI_SUM, world);
    if (!force->newton_bond) atom->ndihedrals /= 4;
  }
  if (atom->avec->impropers_allow) {
    MPI_Allreduce(&nimpropers, &atom->nimpropers, 1, MPI_LMP_BIGINT, MPI_SUM, world);
    if (!force->newton_bond) atom->nimpropers /= 4;
  }
}

/* ----------------------------------------------------------------------
   callback from comm->ring() in delete_bond()
------------------------------------------------------------------------- */

void DeleteAtoms::bondring(int nbuf, char *cbuf, void *ptr)
{
  auto *daptr = (DeleteAtoms *) ptr;
  auto *list = (tagint *) cbuf;

  int *num_bond = daptr->atom->num_bond;
  int *num_angle = daptr->atom->num_angle;
  int *num_dihedral = daptr->atom->num_dihedral;
  int *num_improper = daptr->atom->num_improper;

  int **bond_type = daptr->atom->bond_type;
  tagint **bond_atom = daptr->atom->bond_atom;

  int **angle_type = daptr->atom->angle_type;
  tagint **angle_atom1 = daptr->atom->angle_atom1;
  tagint **angle_atom2 = daptr->atom->angle_atom2;
  tagint **angle_atom3 = daptr->atom->angle_atom3;

  int **dihedral_type = daptr->atom->dihedral_type;
  tagint **dihedral_atom1 = daptr->atom->dihedral_atom1;
  tagint **dihedral_atom2 = daptr->atom->dihedral_atom2;
  tagint **dihedral_atom3 = daptr->atom->dihedral_atom3;
  tagint **dihedral_atom4 = daptr->atom->dihedral_atom4;

  int **improper_type = daptr->atom->improper_type;
  tagint **improper_atom1 = daptr->atom->improper_atom1;
  tagint **improper_atom2 = daptr->atom->improper_atom2;
  tagint **improper_atom3 = daptr->atom->improper_atom3;
  tagint **improper_atom4 = daptr->atom->improper_atom4;

  int nlocal = daptr->atom->nlocal;

  // find instances of bond history to delete data
  auto histories = daptr->modify->get_fix_by_style("BOND_HISTORY");
  int n_histories = histories.size();

  // cbuf = list of N deleted atom IDs from other proc, put them in hash
  std::unordered_set<tagint> hash(list, list + nbuf);

  // loop over my atoms and their bond topology lists
  // if any atom in an interaction matches atom ID in hash, delete interaction

  int m, n;
  for (int i = 0; i < nlocal; i++) {
    if (num_bond) {
      m = 0;
      n = num_bond[i];
      while (m < n) {
        if (hash.find(bond_atom[i][m]) != hash.end()) {
          bond_type[i][m] = bond_type[i][n - 1];
          bond_atom[i][m] = bond_atom[i][n - 1];
          if (n_histories > 0)
            for (auto &ihistory : histories) {
              dynamic_cast<FixBondHistory *>(ihistory)->shift_history(i, m, n - 1);
              dynamic_cast<FixBondHistory *>(ihistory)->delete_history(i, n - 1);
            }
          n--;
        } else
          m++;
      }
      num_bond[i] = n;
    }

    if (num_angle) {
      m = 0;
      n = num_angle[i];
      while (m < n) {
        if (hash.find(angle_atom1[i][m]) != hash.end() ||
            hash.find(angle_atom2[i][m]) != hash.end() ||
            hash.find(angle_atom3[i][m]) != hash.end()) {
          angle_type[i][m] = angle_type[i][n - 1];
          angle_atom1[i][m] = angle_atom1[i][n - 1];
          angle_atom2[i][m] = angle_atom2[i][n - 1];
          angle_atom3[i][m] = angle_atom3[i][n - 1];
          n--;
        } else
          m++;
      }
      num_angle[i] = n;
    }

    if (num_dihedral) {
      m = 0;
      n = num_dihedral[i];
      while (m < n) {
        if (hash.find(dihedral_atom1[i][m]) != hash.end() ||
            hash.find(dihedral_atom2[i][m]) != hash.end() ||
            hash.find(dihedral_atom3[i][m]) != hash.end() ||
            hash.find(dihedral_atom4[i][m]) != hash.end()) {
          dihedral_type[i][m] = dihedral_type[i][n - 1];
          dihedral_atom1[i][m] = dihedral_atom1[i][n - 1];
          dihedral_atom2[i][m] = dihedral_atom2[i][n - 1];
          dihedral_atom3[i][m] = dihedral_atom3[i][n - 1];
          dihedral_atom4[i][m] = dihedral_atom4[i][n - 1];
          n--;
        } else
          m++;
      }
      num_dihedral[i] = n;
    }

    if (num_improper) {
      m = 0;
      n = num_improper[i];
      while (m < n) {
        if (hash.find(improper_atom1[i][m]) != hash.end() ||
            hash.find(improper_atom2[i][m]) != hash.end() ||
            hash.find(improper_atom3[i][m]) != hash.end() ||
            hash.find(improper_atom4[i][m]) != hash.end()) {
          improper_type[i][m] = improper_type[i][n - 1];
          improper_atom1[i][m] = improper_atom1[i][n - 1];
          improper_atom2[i][m] = improper_atom2[i][n - 1];
          improper_atom3[i][m] = improper_atom3[i][n - 1];
          improper_atom4[i][m] = improper_atom4[i][n - 1];
          n--;
        } else
          m++;
      }
      num_improper[i] = n;
    }
  }
}

/* ----------------------------------------------------------------------
   callback from comm->ring() in delete_molecule()
------------------------------------------------------------------------- */

void DeleteAtoms::molring(int n, char *cbuf, void *ptr)
{
  auto *daptr = (DeleteAtoms *) ptr;
  auto *list = (tagint *) cbuf;
  int *dlist = daptr->dlist;
  int nlocal = daptr->atom->nlocal;
  tagint *molecule = daptr->atom->molecule;

  // cbuf = list of N molecule IDs from other proc, put them in hash

  std::unordered_set<tagint> hash(list, list + n);

  // loop over my atoms, if matches molecule ID in hash, delete that atom

  for (int i = 0; i < nlocal; i++)
    if (hash.find(molecule[i]) != hash.end()) dlist[i] = 1;
}

/* ----------------------------------------------------------------------
   callback from comm->ring() for adding tags from remaining atoms to buffer
------------------------------------------------------------------------- */

void DeleteAtoms::addtags(int nbuf, char *cbuf, void *ptr)
{
  auto *daptr = (DeleteAtoms *) ptr;
  auto *taglist = (tagint *) cbuf;

  int nlocal = daptr->atom->nlocal;
  tagint *tag = daptr->atom->tag;
  int buf_rank = (daptr->ringrank + daptr->comm->nprocs) % daptr->comm->nprocs;
  int idx = 0;
  tagint mytagmin = (tagint) buf_rank * nbuf;
  tagint mytagmax = mytagmin + nbuf;

  // skip over list entries that are already set
  while ((taglist[idx] != 0) && (idx < nbuf)) ++idx;

  // add remaining tags that fit the range into buffer
  for (int i = 0; i < nlocal; ++i) {
    if ((tag[i] > mytagmin) && (tag[i] <= mytagmax)) taglist[idx++] = tag[i];
  }

  daptr->ringrank--;    // next buffer will be from previous rank.
}

/* ----------------------------------------------------------------------
   callback from comm->ring() for setting new tags
------------------------------------------------------------------------- */

void DeleteAtoms::settags(int nbuf, char *cbuf, void *ptr)
{
  auto *daptr = (DeleteAtoms *) ptr;
  auto *taglist = (tagint *) cbuf;

  int nlocal = daptr->atom->nlocal;
  int buf_rank = (daptr->ringrank + daptr->comm->nprocs) % daptr->comm->nprocs;
  tagint *newtags = daptr->newtags;

  // loop over tags in this buffer and determine their new tag
  // if the atom with the old tag exists as a local atom, assign it in list
  int idx = 0;
  tagint tagval = daptr->tagproc[buf_rank];
  for (int i = 0; i < nbuf; ++i) {
    if (taglist[i] == 0) break;    // end of tags; only zeros beyond this.
    ++tagval;
    idx = daptr->atom->map(taglist[i]);
    if ((idx < 0) || (idx >= nlocal)) continue;
    newtags[idx] = tagval;
  }

  daptr->ringrank--;    // next buffer will be from previous rank
}

/* ----------------------------------------------------------------------
   process command options
------------------------------------------------------------------------- */

void DeleteAtoms::options(int narg, char **arg)
{
  // default to "compress yes" for atomic systems and "compress no" for molecular systems
  if (atom->molecular == Atom::ATOMIC)
    compress_flag = 1;
  else
    compress_flag = 0;

  bond_flag = mol_flag = 0;
  condense_flag = 0;

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "compress") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "delete_atoms compress", error);
      compress_flag = utils::logical(FLERR, arg[iarg + 1], false, lmp);
      if (condense_flag && compress_flag) condense_flag = 0;
      iarg += 2;
    } else if (strcmp(arg[iarg], "condense") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "delete_atoms condense", error);
      condense_flag = utils::logical(FLERR, arg[iarg + 1], false, lmp);
      if (compress_flag && condense_flag) compress_flag = 0;
      iarg += 2;
    } else if (strcmp(arg[iarg], "bond") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "delete_atoms bond", error);
      bond_flag = utils::logical(FLERR, arg[iarg + 1], false, lmp);
      if (bond_flag && (atom->molecular == Atom::ATOMIC))
        error->all(FLERR, iarg, "Cannot use delete_atoms bond yes for non-molecular systems");
      if (bond_flag && (atom->molecular == Atom::TEMPLATE))
        error->all(FLERR, iarg, "Cannot use delete_atoms bond yes with atom_style template");
      iarg += 2;
    } else if (strcmp(arg[iarg], "mol") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "delete_atoms mol", error);
      mol_flag = utils::logical(FLERR, arg[iarg + 1], false, lmp);
      if (mol_flag && (atom->molecule_flag == 0))
        error->all(FLERR, iarg, "Delete_atoms mol yes requires atom attribute molecule");
      iarg += 2;
    } else
      error->all(FLERR, iarg, "Unknown delete_atoms option: {}", arg[iarg]);
  }
}

/* ----------------------------------------------------------------------
   condense tags for ATOMIC and MOLECULAR
   Algorithm: Parallel Sample Sort (Scalable O(N/P))
   Preserves global ID ordering: OldID_A < OldID_B ==> NewID_A < NewID_B
------------------------------------------------------------------------- */
void DeleteAtoms::condense_tags() {

    if (atom->map_style == Atom::MAP_NONE)
        error->all(FLERR, "Using 'condense yes' option requires an atom map");

    // ---------------------------------------------------------
    // 1. REFRESH COMM LISTS (Prevents Segfault/BusError)
    // ---------------------------------------------------------
    if (domain->triclinic) domain->x2lamda(atom->nlocal);
    domain->pbc();
    domain->reset_box();

    // Force re-binning and list rebuild
    if (neighbor->style) neighbor->setup_bins();
    comm->setup();

    if (domain->triclinic) domain->lamda2x(atom->nlocal + atom->nghost);

    atom->map_init();
    atom->map_set();

    int nlocal = atom->nlocal;
    int nmax = atom->nmax;
    int nprocs, me;
    MPI_Comm_size(world, &nprocs);
    MPI_Comm_rank(world, &me);

    // ---------------------------------------------------------
    // 2. PARALLEL SAMPLE SORT
    // ---------------------------------------------------------

    struct SortItem {
        tagint tag;
        int origin_index;
        int origin_rank;
    };

    // A. Local Sort & Sampling
    std::vector<SortItem> my_items(nlocal);
    for(int i=0; i<nlocal; i++) {
        my_items[i] = {atom->tag[i], i, me};
    }
    // Sort local atoms by Tag
    std::sort(my_items.begin(), my_items.end(),
              [](const SortItem &a, const SortItem &b) { return a.tag < b.tag; });

    // Select P samples per proc
    int num_samples = nprocs;
    std::vector<tagint> local_samples(num_samples, MAXTAGINT);
    if (nlocal > 0) {
        for(int i=0; i<num_samples; i++) {
            int idx = (i * nlocal) / num_samples;
            local_samples[i] = my_items[idx].tag;
        }
    }

    // B. Gather Samples & Pick Splitters (Root only)
    std::vector<tagint> all_samples(nprocs * num_samples);
    MPI_Gather(local_samples.data(), num_samples, MPI_LMP_TAGINT,
               all_samples.data(), num_samples, MPI_LMP_TAGINT,
               0, world);

    std::vector<tagint> splitters(nprocs + 1);
    if (me == 0) {
        std::sort(all_samples.begin(), all_samples.end());
        // Pick splitters that evenly divide the sorted samples
        for(int i=1; i<nprocs; i++) {
            splitters[i] = all_samples[(i * all_samples.size()) / nprocs];
        }
        splitters[0] = 0;
        splitters[nprocs] = MAXTAGINT;
    }
    MPI_Bcast(splitters.data(), nprocs + 1, MPI_LMP_TAGINT, 0, world);

    // C. Bucketing (Assign atoms to target Procs)
    std::vector<int> send_counts(nprocs, 0);
    std::vector<std::vector<SortItem>> buckets(nprocs);

    int current_bucket = 0;
    for(const auto &item : my_items) {
        // Find which processor owns this tag range
        while (current_bucket < nprocs - 1 && item.tag >= splitters[current_bucket+1]) {
            current_bucket++;
        }
        buckets[current_bucket].push_back(item);
        send_counts[current_bucket]++;
    }

    // D. Exchange Data (Alltoallv)
    std::vector<int> send_displs(nprocs, 0);
    std::vector<SortItem> send_buffer;
    send_buffer.reserve(nlocal);
    int offset = 0;
    for(int i=0; i<nprocs; i++) {
        send_displs[i] = offset;
        offset += send_counts[i];
        send_buffer.insert(send_buffer.end(), buckets[i].begin(), buckets[i].end());
    }

    std::vector<int> recv_counts(nprocs);
    MPI_Alltoall(send_counts.data(), 1, MPI_INT, recv_counts.data(), 1, MPI_INT, world);

    std::vector<int> recv_displs(nprocs, 0);
    int total_recv = 0;
    for(int i=0; i<nprocs; i++) {
        recv_displs[i] = total_recv;
        total_recv += recv_counts[i];
    }

    // Pack into flat arrays for MPI (struct SortItem is 16-24 bytes)
    // tagint=8 bytes, int=4, int=4.
    std::vector<tagint> s_tags(send_buffer.size());
    std::vector<int> s_idx(send_buffer.size()), s_rank(send_buffer.size());
    for(size_t i=0; i<send_buffer.size(); i++) {
        s_tags[i] = send_buffer[i].tag;
        s_idx[i]  = send_buffer[i].origin_index;
        s_rank[i] = send_buffer[i].origin_rank;
    }

    std::vector<tagint> r_tags(total_recv);
    std::vector<int> r_idx(total_recv), r_rank(total_recv);

    MPI_Alltoallv(s_tags.data(), send_counts.data(), send_displs.data(), MPI_LMP_TAGINT,
                  r_tags.data(), recv_counts.data(), recv_displs.data(), MPI_LMP_TAGINT, world);
    MPI_Alltoallv(s_idx.data(), send_counts.data(), send_displs.data(), MPI_INT,
                  r_idx.data(), recv_counts.data(), recv_displs.data(), MPI_INT, world);
    MPI_Alltoallv(s_rank.data(), send_counts.data(), send_displs.data(), MPI_INT,
                  r_rank.data(), recv_counts.data(), recv_displs.data(), MPI_INT, world);

    // E. Sort Received Data & Assign Global IDs
    std::vector<SortItem> sorted_items(total_recv);
    for(int i=0; i<total_recv; i++) {
        sorted_items[i] = {r_tags[i], r_idx[i], r_rank[i]};
    }
    // Final Sort: Now Proc 0 has smallest tags, Proc N has largest
    std::sort(sorted_items.begin(), sorted_items.end(),
              [](const SortItem &a, const SortItem &b) { return a.tag < b.tag; });

    // Determine Global Start ID for this processor
    long local_size = total_recv;
    long scan_offset = 0;
    MPI_Scan(&local_size, &scan_offset, 1, MPI_LONG, MPI_SUM, world);
    tagint start_id = (tagint)(scan_offset - local_size) + 1;

    // F. Return New IDs to Owners
    // Reuse buckets to send data back
    for(int i=0; i<nprocs; i++) { send_counts[i]=0; buckets[i].clear(); }

    for(int i=0; i<total_recv; i++) {
        tagint new_id = start_id + i;
        int target = sorted_items[i].origin_rank;
        // Payload: {NewID, LocalIndexOnOrigin, 0}
        buckets[target].push_back({new_id, sorted_items[i].origin_index, 0});
        send_counts[target]++;
    }

    // Re-pack
    send_buffer.clear();
    offset = 0;
    for(int i=0; i<nprocs; i++) {
        send_displs[i] = offset;
        offset += send_counts[i];
        send_buffer.insert(send_buffer.end(), buckets[i].begin(), buckets[i].end());
    }

    // Exchange counts (Return Trip)
    MPI_Alltoall(send_counts.data(), 1, MPI_INT, recv_counts.data(), 1, MPI_INT, world);
    total_recv = 0;
    for(int i=0; i<nprocs; i++) { recv_displs[i] = total_recv; total_recv += recv_counts[i]; }

    // Pack
    s_tags.resize(send_buffer.size());
    s_idx.resize(send_buffer.size());
    for(size_t i=0; i<send_buffer.size(); i++) {
        s_tags[i] = send_buffer[i].tag;  // This is the NEW ID
        s_idx[i]  = send_buffer[i].origin_index;
    }
    r_tags.resize(total_recv);
    r_idx.resize(total_recv);

    // Exchange Data (Return Trip)
    MPI_Alltoallv(s_tags.data(), send_counts.data(), send_displs.data(), MPI_LMP_TAGINT,
                  r_tags.data(), recv_counts.data(), recv_displs.data(), MPI_LMP_TAGINT, world);
    MPI_Alltoallv(s_idx.data(), send_counts.data(), send_displs.data(), MPI_INT,
                  r_idx.data(), recv_counts.data(), recv_displs.data(), MPI_INT, world);

    std::fprintf(stderr, "*** ok 1\n");
    std::fflush(stderr);

    // ---------------------------------------------------------
    // 3. APPLY NEW IDS
    // ---------------------------------------------------------
    // Create newIDs buffer (Size nmax for safety)
    double **newIDs;
    memory->create(newIDs, nmax, 1, "delete_atoms:newIDs");
    // Init to 0
    for(int i=0; i<nmax; i++) newIDs[i][0] = ubuf((tagint)0).d;

    // Fill buffer from received data
    // Now we are back on the original processor, so r_idx is the local index
    for(int i=0; i<total_recv; i++) {
        int idx = r_idx[i];
        if (idx < nmax) newIDs[idx][0] = ubuf(r_tags[i]).d;
    }
    
    std::fprintf(stderr, "*** ok 2\n");
    std::fflush(stderr);

    // ---------------------------------------------------------
    // 4. UPDATE TOPOLOGY (Standard)
    // ---------------------------------------------------------
    comm->forward_comm_array(1, newIDs);
    
    std::fprintf(stderr, "*** ok 3\n");
    std::fflush(stderr);

    if (atom->molecular != Atom::ATOMIC) {
        // --- BONDS ---
        if (atom->avec->bonds_allow) {
            for (int i = 0; i < nlocal; i++) {
                int m = 0;
                for (int j = 0; j < atom->num_bond[i]; j++) {
                   tagint oldID = atom->bond_atom[i][j];
                   int k = atom->map(oldID);
                   if (k >= 0) {
                      atom->bond_atom[i][m] = (tagint) ubuf(newIDs[k][0]).i;
                      atom->bond_type[i][m] = atom->bond_type[i][j];
                      m++;
                   }
                }
                atom->num_bond[i] = m;
            }
        }

        std::fprintf(stderr, "*** ok 4\n");
        std::fflush(stderr);

        // --- ANGLES ---
        if (atom->avec->angles_allow) {
            for (int i = 0; i < nlocal; i++) {
                int m = 0;
                for (int j = 0; j < atom->num_angle[i]; j++) {
                    int k1 = atom->map(atom->angle_atom1[i][j]);
                    int k2 = atom->map(atom->angle_atom2[i][j]);
                    int k3 = atom->map(atom->angle_atom3[i][j]);
                    if (k1 >= 0 && k2 >= 0 && k3 >= 0) {
                        atom->angle_atom1[i][m] = (tagint) ubuf(newIDs[k1][0]).i;
                        atom->angle_atom2[i][m] = (tagint) ubuf(newIDs[k2][0]).i;
                        atom->angle_atom3[i][m] = (tagint) ubuf(newIDs[k3][0]).i;
                        atom->angle_type[i][m] = atom->angle_type[i][j];
                        m++;
                    }
                }
                atom->num_angle[i] = m;
            }
        }

        // --- DIHEDRALS ---
        if (atom->avec->dihedrals_allow) {
             for (int i = 0; i < nlocal; i++) {
                int m = 0;
                for (int j = 0; j < atom->num_dihedral[i]; j++) {
                    int k1 = atom->map(atom->dihedral_atom1[i][j]);
                    int k2 = atom->map(atom->dihedral_atom2[i][j]);
                    int k3 = atom->map(atom->dihedral_atom3[i][j]);
                    int k4 = atom->map(atom->dihedral_atom4[i][j]);
                    if (k1 >= 0 && k2 >= 0 && k3 >= 0 && k4 >= 0) {
                        atom->dihedral_atom1[i][m] = (tagint) ubuf(newIDs[k1][0]).i;
                        atom->dihedral_atom2[i][m] = (tagint) ubuf(newIDs[k2][0]).i;
                        atom->dihedral_atom3[i][m] = (tagint) ubuf(newIDs[k3][0]).i;
                        atom->dihedral_atom4[i][m] = (tagint) ubuf(newIDs[k4][0]).i;
                        atom->dihedral_type[i][m] = atom->dihedral_type[i][j];
                        m++;
                    }
                }
                atom->num_dihedral[i] = m;
            }
        }

        // --- IMPROPERS ---
        if (atom->avec->impropers_allow) {
            for (int i = 0; i < nlocal; i++) {
                int m = 0;
                for (int j = 0; j < atom->num_improper[i]; j++) {
                    int k1 = atom->map(atom->improper_atom1[i][j]);
                    int k2 = atom->map(atom->improper_atom2[i][j]);
                    int k3 = atom->map(atom->improper_atom3[i][j]);
                    int k4 = atom->map(atom->improper_atom4[i][j]);
                    if (k1 >= 0 && k2 >= 0 && k3 >= 0 && k4 >= 0) {
                        atom->improper_atom1[i][m] = (tagint) ubuf(newIDs[k1][0]).i;
                        atom->improper_atom2[i][m] = (tagint) ubuf(newIDs[k2][0]).i;
                        atom->improper_atom3[i][m] = (tagint) ubuf(newIDs[k3][0]).i;
                        atom->improper_atom4[i][m] = (tagint) ubuf(newIDs[k4][0]).i;
                        atom->improper_type[i][m] = atom->improper_type[i][j];
                        m++;
                    }
                }
                atom->num_improper[i] = m;
            }
        }

        // --- SPECIAL LISTS ---
        if (atom->special) {
            for (int i = 0; i < nlocal; i++) {
                int n_total = atom->nspecial[i][0] + atom->nspecial[i][1] + atom->nspecial[i][2];
                for (int j = 0; j < n_total; j++) {
                     int k = atom->map(atom->special[i][j]);
                     if (k >= 0) atom->special[i][j] = (tagint) ubuf(newIDs[k][0]).i;
                     else atom->special[i][j] = 0;
                }
            }
        }

        std::fprintf(stderr, "*** ok 16\n");
        std::fflush(stderr);

        // --- FIXES ---
        for (auto *ifix : modify->get_fix_list()) {
            if (ifix->stores_ids) ifix->update_ids(newIDs);
        }
    }

    std::fprintf(stderr, "*** ok 17\n");
    std::fflush(stderr);

    // Assign My New Tags
    for (int i = 0; i < nlocal; i++) atom->tag[i] = (tagint) ubuf(newIDs[i][0]).i;

    memory->destroy(newIDs);
    atom->map_tag_max = -1;
    atom->map_init();
    atom->map_set();
}

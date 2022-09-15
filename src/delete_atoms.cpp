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

#include <cstring>
#include <map>
#include <utility>

using namespace LAMMPS_NS;

enum { UNKNOWN, FRACTION, COUNT };

/* ---------------------------------------------------------------------- */

DeleteAtoms::DeleteAtoms(LAMMPS *lmp) : Command(lmp) {}

/* ---------------------------------------------------------------------- */

void DeleteAtoms::command(int narg, char **arg)
{
  if (domain->box_exist == 0)
    error->all(FLERR, "Delete_atoms command before simulation box is defined");
  if (narg < 1) utils::missing_cmd_args(FLERR, "delete_atoms", error);
  if (atom->tag_enable == 0) error->all(FLERR, "Cannot use delete_atoms unless atoms have IDs");

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
    error->all(FLERR,
               "The delete_atoms 'porosity' keyword has been removed.\n"
               "Please use: delete_atoms random fraction frac exact group-ID region-ID seed\n");
  } else if (strcmp(arg[0], "variable") == 0)
    delete_variable(narg, arg);
  else
    error->all(FLERR, "Unknown delete_atoms sub-command: {}", arg[0]);

  if (allflag) {
    int igroup = group->find("all");
    if ((igroup >= 0) && modify->check_rigid_group_overlap(group->bitmask[igroup]))
      error->warning(FLERR, "Attempting to delete atoms in rigid bodies");
  } else {
    if (modify->check_rigid_list_overlap(dlist))
      error->warning(FLERR, "Attempting to delete atoms in rigid bodies");
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

  if (compress_flag) {
    if (atom->molecular == Atom::ATOMIC) {
      tagint *tag = atom->tag;
      int nlocal = atom->nlocal;
      for (int i = 0; i < nlocal; i++) tag[i] = 0;
      atom->tag_extend();
    } else if (comm->me == 0)
      error->warning(FLERR, "Ignoring 'compress yes' for molecular system");
  }

  // reset atom->natoms and also topology counts

  bigint nblocal = atom->nlocal;
  MPI_Allreduce(&nblocal, &atom->natoms, 1, MPI_LMP_BIGINT, MPI_SUM, world);

  // reset bonus data counts

  auto avec_ellipsoid = dynamic_cast<AtomVecEllipsoid *>(atom->style_match("ellipsoid"));
  auto avec_line = dynamic_cast<AtomVecLine *>(atom->style_match("line"));
  auto avec_tri = dynamic_cast<AtomVecTri *>(atom->style_match("tri"));
  auto avec_body = dynamic_cast<AtomVecBody *>(atom->style_match("body"));
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

  int igroup = group->find(arg[1]);
  if (igroup == -1) error->all(FLERR, "Could not find delete_atoms group ID {}", arg[1]);
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
  int groupbit = group->bitmask[igroup];

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) dlist[i] = 1;
}

/* ----------------------------------------------------------------------
   delete all atoms in region
------------------------------------------------------------------------- */

void DeleteAtoms::delete_region(int narg, char **arg)
{
  if (narg < 2) utils::missing_cmd_args(FLERR, "delete_atoms region", error);

  auto iregion = domain->get_region_by_id(arg[1]);
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

  int igroup1 = group->find(arg[2]);
  if (igroup1 < 0)
    error->all(FLERR, "Could not find delete_atoms overlap first group ID {}", arg[2]);
  int igroup2 = group->find(arg[3]);
  if (igroup2 < 0)
    error->all(FLERR, "Could not find delete_atoms overlap second group ID {}", arg[3]);
  options(narg - 4, &arg[4]);

  const int group1bit = group->bitmask[igroup1];
  const int group2bit = group->bitmask[igroup2];

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
  if (cut > neighbor->cutneighmin && comm->me == 0)
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

  auto list = neighbor->find_list(this);
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
      //   test on tags of I,J insures that only I or J is deleted

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

  int igroup = group->find(arg[4]);
  if (igroup == -1) error->all(FLERR, "Could not find delete_atoms random group ID {}", arg[4]);

  auto region = domain->get_region_by_id(arg[5]);
  if (!region && (strcmp(arg[5], "NULL") != 0))
    error->all(FLERR, "Could not find delete_atoms random region ID {}", arg[5]);

  int seed = utils::inumeric(FLERR, arg[6], false, lmp);
  options(narg - 7, &arg[7]);

  auto ranmars = new RanMars(lmp, seed + comm->me);

  // allocate and initialize deletion list

  int nlocal = atom->nlocal;
  memory->create(dlist, nlocal, "delete_atoms:dlist");
  for (int i = 0; i < nlocal; i++) dlist[i] = 0;

  // setup

  double **x = atom->x;
  int *mask = atom->mask;

  int groupbit = group->bitmask[igroup];
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
  // hash = for atom IDs being deleted by one processor
  // list of these IDs is sent around ring
  // at each stage of ring pass, hash is re-populated with received IDs

  hash = new std::map<tagint, int>();

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

  delete hash;
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

  hash = new std::map<tagint, int>();

  tagint *molecule = atom->molecule;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (molecule[i] == 0) continue;
    if (dlist[i] && hash->find(molecule[i]) == hash->end()) (*hash)[molecule[i]] = 1;
  }

  // list = set of unique molecule IDs from which I deleted atoms
  // pass list to all other procs via comm->ring()

  int n = hash->size();
  tagint *list;
  memory->create(list, n, "delete_atoms:list");

  n = 0;
  std::map<tagint, int>::iterator pos;
  for (pos = hash->begin(); pos != hash->end(); ++pos) list[n++] = pos->first;

  comm->ring(n, sizeof(tagint), list, 1, molring, nullptr, (void *) this);

  delete hash;
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
  auto daptr = (DeleteAtoms *) ptr;
  auto list = (tagint *) cbuf;
  std::map<tagint, int> *hash = daptr->hash;

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

  // cbuf = list of N deleted atom IDs from other proc, put them in hash

  hash->clear();
  for (int i = 0; i < nbuf; i++) (*hash)[list[i]] = 1;

  // loop over my atoms and their bond topology lists
  // if any atom in an interaction matches atom ID in hash, delete interaction

  int m, n;
  for (int i = 0; i < nlocal; i++) {
    if (num_bond) {
      m = 0;
      n = num_bond[i];
      while (m < n) {
        if (hash->find(bond_atom[i][m]) != hash->end()) {
          bond_type[i][m] = bond_type[i][n - 1];
          bond_atom[i][m] = bond_atom[i][n - 1];
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
        if (hash->find(angle_atom1[i][m]) != hash->end() ||
            hash->find(angle_atom2[i][m]) != hash->end() ||
            hash->find(angle_atom3[i][m]) != hash->end()) {
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
        if (hash->find(dihedral_atom1[i][m]) != hash->end() ||
            hash->find(dihedral_atom2[i][m]) != hash->end() ||
            hash->find(dihedral_atom3[i][m]) != hash->end() ||
            hash->find(dihedral_atom4[i][m]) != hash->end()) {
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
        if (hash->find(improper_atom1[i][m]) != hash->end() ||
            hash->find(improper_atom2[i][m]) != hash->end() ||
            hash->find(improper_atom3[i][m]) != hash->end() ||
            hash->find(improper_atom4[i][m]) != hash->end()) {
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
  auto daptr = (DeleteAtoms *) ptr;
  auto list = (tagint *) cbuf;
  int *dlist = daptr->dlist;
  std::map<tagint, int> *hash = daptr->hash;
  int nlocal = daptr->atom->nlocal;
  tagint *molecule = daptr->atom->molecule;

  // cbuf = list of N molecule IDs from other proc, put them in hash

  hash->clear();
  for (int i = 0; i < n; i++) (*hash)[list[i]] = 1;

  // loop over my atoms, if matches molecule ID in hash, delete that atom

  for (int i = 0; i < nlocal; i++)
    if (hash->find(molecule[i]) != hash->end()) dlist[i] = 1;
}

/* ----------------------------------------------------------------------
   process command options
------------------------------------------------------------------------- */

void DeleteAtoms::options(int narg, char **arg)
{
  compress_flag = 1;
  bond_flag = mol_flag = 0;

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "compress") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "delete_atoms compress", error);
      compress_flag = utils::logical(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "bond") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "delete_atoms bond", error);
      if (atom->molecular == Atom::ATOMIC)
        error->all(FLERR, "Cannot use delete_atoms bond yes for non-molecular systems");
      if (atom->molecular == Atom::TEMPLATE)
        error->all(FLERR, "Cannot use delete_atoms bond yes with atom_style template");
      bond_flag = utils::logical(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "mol") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "delete_atoms mol", error);
      if (atom->molecule_flag == 0)
        error->all(FLERR, "Delete_atoms mol yes requires atom attribute molecule");
      mol_flag = utils::logical(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else
      error->all(FLERR, "Unknown delete_atoms option: {}", arg[iarg]);
  }
}

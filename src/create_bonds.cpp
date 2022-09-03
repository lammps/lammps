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
   Contributing authors:
     Mike Salerno (NRL) added single methods
     Thomas Farmer (ISIS) added single/improper
------------------------------------------------------------------------- */

#include "create_bonds.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "group.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "special.h"

#include <cstring>

using namespace LAMMPS_NS;

enum { MANY, SBOND, SANGLE, SDIHEDRAL, SIMPROPER };

/* ---------------------------------------------------------------------- */

CreateBonds::CreateBonds(LAMMPS *lmp) : Command(lmp) {}

/* ---------------------------------------------------------------------- */

void CreateBonds::command(int narg, char **arg)
{
  if (domain->box_exist == 0)
    error->all(FLERR, "Create_bonds command before simulation box is defined");
  if (atom->tag_enable == 0) error->all(FLERR, "Cannot use create_bonds unless atoms have IDs");
  if (atom->molecular != Atom::MOLECULAR)
    error->all(FLERR, "Cannot use create_bonds with non-molecular system");

  if (narg < 4) error->all(FLERR, "Illegal create_bonds command");

  // parse args

  int style;

  int iarg = 0;
  if (strcmp(arg[0], "many") == 0) {
    style = MANY;
    if (narg != 6) error->all(FLERR, "Illegal create_bonds command");
    igroup = group->find(arg[1]);
    if (igroup == -1) error->all(FLERR, "Cannot find create_bonds group ID");
    group1bit = group->bitmask[igroup];
    igroup = group->find(arg[2]);
    if (igroup == -1) error->all(FLERR, "Cannot find create_bonds group ID");
    group2bit = group->bitmask[igroup];
    btype = utils::inumeric(FLERR, arg[3], false, lmp);
    rmin = utils::numeric(FLERR, arg[4], false, lmp);
    rmax = utils::numeric(FLERR, arg[5], false, lmp);
    if (rmin > rmax) error->all(FLERR, "Illegal create_bonds command");
    iarg = 6;
  } else if (strcmp(arg[0], "single/bond") == 0) {
    style = SBOND;
    btype = utils::inumeric(FLERR, arg[1], false, lmp);
    batom1 = utils::tnumeric(FLERR, arg[2], false, lmp);
    batom2 = utils::tnumeric(FLERR, arg[3], false, lmp);
    if (batom1 == batom2) error->all(FLERR, "Illegal create_bonds command");
    iarg = 4;
  } else if (strcmp(arg[0], "single/angle") == 0) {
    style = SANGLE;
    if (narg < 5) error->all(FLERR, "Illegal create_bonds command");
    atype = utils::inumeric(FLERR, arg[1], false, lmp);
    aatom1 = utils::tnumeric(FLERR, arg[2], false, lmp);
    aatom2 = utils::tnumeric(FLERR, arg[3], false, lmp);
    aatom3 = utils::tnumeric(FLERR, arg[4], false, lmp);
    if ((aatom1 == aatom2) || (aatom1 == aatom3) || (aatom2 == aatom3))
      error->all(FLERR, "Illegal create_bonds command");
    iarg = 5;
  } else if (strcmp(arg[0], "single/dihedral") == 0) {
    style = SDIHEDRAL;
    if (narg < 6) error->all(FLERR, "Illegal create_bonds command");
    dtype = utils::inumeric(FLERR, arg[1], false, lmp);
    datom1 = utils::tnumeric(FLERR, arg[2], false, lmp);
    datom2 = utils::tnumeric(FLERR, arg[3], false, lmp);
    datom3 = utils::tnumeric(FLERR, arg[4], false, lmp);
    datom4 = utils::tnumeric(FLERR, arg[5], false, lmp);
    if ((datom1 == datom2) || (datom1 == datom3) || (datom1 == datom4) || (datom2 == datom3) ||
        (datom2 == datom4) || (datom3 == datom4))
      error->all(FLERR, "Illegal create_bonds command");
    iarg = 6;
  } else if (strcmp(arg[0], "single/improper") == 0) {
    style = SIMPROPER;
    if (narg < 6) error->all(FLERR, "Illegal create_bonds command");
    dtype = utils::inumeric(FLERR, arg[1], false, lmp);
    datom1 = utils::tnumeric(FLERR, arg[2], false, lmp);
    datom2 = utils::tnumeric(FLERR, arg[3], false, lmp);
    datom3 = utils::tnumeric(FLERR, arg[4], false, lmp);
    datom4 = utils::tnumeric(FLERR, arg[5], false, lmp);
    if ((datom1 == datom2) || (datom1 == datom3) || (datom1 == datom4) || (datom2 == datom3) ||
        (datom2 == datom4) || (datom3 == datom4))
      error->all(FLERR, "Illegal create_bonds command");
    iarg = 6;
  } else
    error->all(FLERR, "Illegal create_bonds command");

  // optional args

  int specialflag = 1;

  while (iarg < narg) {
    if (strcmp(arg[iarg], "special") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal create_bonds command");
      specialflag = utils::logical(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else
      error->all(FLERR, "Illegal create_bonds command");
  }

  // error checks

  if (style == MANY) {
    if (btype <= 0 || btype > atom->nbondtypes)
      error->all(FLERR, "Invalid bond type in create_bonds command");
    if (specialflag == 0) error->all(FLERR, "Cannot use special no with create_bonds many");
  } else if (style == SBOND) {
    if (btype <= 0 || btype > atom->nbondtypes)
      error->all(FLERR, "Invalid bond type in create_bonds command");
  } else if (style == SANGLE) {
    if (atype <= 0 || atype > atom->nangletypes)
      error->all(FLERR, "Invalid angle type in create_bonds command");
  } else if (style == SDIHEDRAL) {
    if (dtype <= 0 || dtype > atom->ndihedraltypes)
      error->all(FLERR, "Invalid dihedral type in create_bonds command");
  } else if (style == SIMPROPER) {
    if (dtype <= 0 || dtype > atom->nimpropertypes)
      error->all(FLERR, "Invalid improper type in create_bonds command");
  }

  // invoke creation method

  if (style == MANY)
    many();
  else if (style == SBOND)
    single_bond();
  else if (style == SANGLE)
    single_angle();
  else if (style == SDIHEDRAL)
    single_dihedral();
  else if (style == SIMPROPER)
    single_improper();

  // trigger special list build

  if (specialflag) {
    Special special(lmp);
    special.build();
  }
}

/* ---------------------------------------------------------------------- */

void CreateBonds::many()
{
  double rminsq = rmin * rmin;
  double rmaxsq = rmax * rmax;

  // store state before bond creation

  bigint nbonds_previous = atom->nbonds;

  // request a full neighbor list for use by this command

  neighbor->add_request(this, "create_bonds", NeighConst::REQ_FULL);

  // init entire system since comm->borders and neighbor->build is done
  // comm::init needs neighbor::init needs pair::init needs kspace::init, etc

  lmp->init();

  // error check on cutoff
  // if no pair style, neighbor list will be empty

  if (force->pair == nullptr) error->all(FLERR, "Create_bonds requires a pair style be defined");
  if (rmax > neighbor->cutneighmax)
    error->all(FLERR, "Create_bonds max distance > neighbor cutoff");
  if (rmax > neighbor->cutneighmin && comm->me == 0)
    error->warning(FLERR, "Create_bonds max distance > minimum neighbor cutoff");

  // require special_bonds 1-2 weights = 0.0 and KSpace = nullptr
  // so that already bonded atom pairs do not appear in neighbor list
  // otherwise with newton_bond = 1,
  //   would be hard to check if I-J bond already existed
  // note that with KSpace, pair with weight = 0 could still be in neigh list

  if (force->special_lj[1] != 0.0 || force->special_coul[1] != 0.0)
    error->all(FLERR, "Create_bonds command requires special_bonds 1-2 weights be 0.0");
  if (force->kspace) error->all(FLERR, "Create_bonds command requires no kspace_style be defined");

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

  // build neighbor list this command needs based on earlier request

  auto list = neighbor->find_list(this);
  neighbor->build_one(list, 1);

  // loop over all neighs of each atom
  // compute distance between two atoms consistently on both procs
  // add bond if group and distance criteria are met
  // check that bond list does not overflow

  tagint *tag = atom->tag;
  int *mask = atom->mask;
  double **x = atom->x;
  int *num_bond = atom->num_bond;
  int **bond_type = atom->bond_type;
  tagint **bond_atom = atom->bond_atom;
  double newton_bond = force->newton_bond;
  int nlocal = atom->nlocal;

  int i, j, ii, jj, inum, jnum, flag;
  double xtmp, ytmp, ztmp, delx, dely, delz, rsq;
  int *ilist, *jlist, *numneigh, **firstneigh;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      // only consider bond creation if I,J distance between 2 cutoffs
      // compute rsq identically on both I,J loop iterations
      // if I,J tags equal, do not bond atom to itself

      if (tag[i] < tag[j]) {
        delx = xtmp - x[j][0];
        dely = ytmp - x[j][1];
        delz = ztmp - x[j][2];
      } else if (tag[i] > tag[j]) {
        delx = x[j][0] - xtmp;
        dely = x[j][1] - ytmp;
        delz = x[j][2] - ztmp;
      } else
        continue;
      rsq = delx * delx + dely * dely + delz * delz;
      if (rsq < rminsq || rsq > rmaxsq) continue;

      // only consider bond creation if igroup and jgroup match I,J atoms

      flag = 0;
      if ((mask[i] & group1bit) && (mask[j] & group2bit)) flag = 1;
      if ((mask[i] & group2bit) && (mask[j] & group1bit)) flag = 1;
      if (!flag) continue;

      // create bond, check for overflow
      // on I,J loop iterations, store with 1 or 2 atoms based on newton_bond

      if (!newton_bond || tag[i] < tag[j]) {
        if (num_bond[i] == atom->bond_per_atom)
          error->one(FLERR, "New bond exceeded bonds per atom limit of {} in create_bonds",
                     atom->bond_per_atom);
        bond_type[i][num_bond[i]] = btype;
        bond_atom[i][num_bond[i]] = tag[j];
        num_bond[i]++;
      }
    }
  }
  neighbor->init();

  // recount bonds

  bigint nbonds = 0;
  for (i = 0; i < nlocal; i++) nbonds += num_bond[i];

  MPI_Allreduce(&nbonds, &atom->nbonds, 1, MPI_LMP_BIGINT, MPI_SUM, world);
  if (!force->newton_bond) atom->nbonds /= 2;

  // print new bond count

  bigint nadd_bonds = atom->nbonds - nbonds_previous;

  if (comm->me == 0)
    utils::logmesg(lmp, "Added {} bonds, new total = {}\n", nadd_bonds, atom->nbonds);
}

/* ---------------------------------------------------------------------- */

void CreateBonds::single_bond()
{
  int m;

  // check that 2 atoms exist

  const int nlocal = atom->nlocal;
  const int idx1 = atom->map(batom1);
  const int idx2 = atom->map(batom2);

  int count = 0;
  if ((idx1 >= 0) && (idx1 < nlocal)) count++;
  if ((idx2 >= 0) && (idx2 < nlocal)) count++;

  int allcount;
  MPI_Allreduce(&count, &allcount, 1, MPI_INT, MPI_SUM, world);
  if (allcount != 2) error->all(FLERR, "Create_bonds single/bond atoms do not exist");

  // create bond once or 2x if newton_bond set

  int *num_bond = atom->num_bond;
  int **bond_type = atom->bond_type;
  tagint **bond_atom = atom->bond_atom;

  if ((m = idx1) >= 0) {
    if (num_bond[m] == atom->bond_per_atom)
      error->one(FLERR, "New bond exceeded bonds per atom in create_bonds");
    bond_type[m][num_bond[m]] = btype;
    bond_atom[m][num_bond[m]] = batom2;
    num_bond[m]++;
  }
  atom->nbonds++;

  if (force->newton_bond) return;

  if ((m = idx2) >= 0) {
    if (num_bond[m] == atom->bond_per_atom)
      error->one(FLERR, "New bond exceeded bonds per atom in create_bonds");
    bond_type[m][num_bond[m]] = btype;
    bond_atom[m][num_bond[m]] = batom1;
    num_bond[m]++;
  }
}

/* ---------------------------------------------------------------------- */

void CreateBonds::single_angle()
{
  int m;

  // check that 3 atoms exist

  const int nlocal = atom->nlocal;
  const int idx1 = atom->map(aatom1);
  const int idx2 = atom->map(aatom2);
  const int idx3 = atom->map(aatom3);

  int count = 0;
  if ((idx1 >= 0) && (idx1 < nlocal)) count++;
  if ((idx2 >= 0) && (idx2 < nlocal)) count++;
  if ((idx3 >= 0) && (idx3 < nlocal)) count++;

  int allcount;
  MPI_Allreduce(&count, &allcount, 1, MPI_INT, MPI_SUM, world);
  if (allcount != 3) error->all(FLERR, "Create_bonds single/angle atoms do not exist");

  // create angle once or 3x if newton_bond set

  int *num_angle = atom->num_angle;
  int **angle_type = atom->angle_type;
  tagint **angle_atom1 = atom->angle_atom1;
  tagint **angle_atom2 = atom->angle_atom2;
  tagint **angle_atom3 = atom->angle_atom3;

  if ((m = idx2) >= 0) {
    if (num_angle[m] == atom->angle_per_atom)
      error->one(FLERR, "New angle exceeded angles per atom in create_bonds");
    angle_type[m][num_angle[m]] = atype;
    angle_atom1[m][num_angle[m]] = aatom1;
    angle_atom2[m][num_angle[m]] = aatom2;
    angle_atom3[m][num_angle[m]] = aatom3;
    num_angle[m]++;
  }
  atom->nangles++;

  if (force->newton_bond) return;

  if ((m = idx1) >= 0) {
    if (num_angle[m] == atom->angle_per_atom)
      error->one(FLERR, "New angle exceeded angles per atom in create_bonds");
    angle_type[m][num_angle[m]] = atype;
    angle_atom1[m][num_angle[m]] = aatom1;
    angle_atom2[m][num_angle[m]] = aatom2;
    angle_atom3[m][num_angle[m]] = aatom3;
    num_angle[m]++;
  }

  if ((m = idx3) >= 0) {
    if (num_angle[m] == atom->angle_per_atom)
      error->one(FLERR, "New angle exceeded angles per atom in create_bonds");
    angle_type[m][num_angle[m]] = atype;
    angle_atom1[m][num_angle[m]] = aatom1;
    angle_atom2[m][num_angle[m]] = aatom2;
    angle_atom3[m][num_angle[m]] = aatom3;
    num_angle[m]++;
  }
}

/* ---------------------------------------------------------------------- */

void CreateBonds::single_dihedral()
{
  int m;

  // check that 4 atoms exist

  const int nlocal = atom->nlocal;
  const int idx1 = atom->map(datom1);
  const int idx2 = atom->map(datom2);
  const int idx3 = atom->map(datom3);
  const int idx4 = atom->map(datom4);

  int count = 0;
  if ((idx1 >= 0) && (idx1 < nlocal)) count++;
  if ((idx2 >= 0) && (idx2 < nlocal)) count++;
  if ((idx3 >= 0) && (idx3 < nlocal)) count++;
  if ((idx4 >= 0) && (idx4 < nlocal)) count++;

  int allcount;
  MPI_Allreduce(&count, &allcount, 1, MPI_INT, MPI_SUM, world);
  if (allcount != 4) error->all(FLERR, "Create_bonds single/dihedral atoms do not exist");

  // create bond once or 4x if newton_bond set

  int *num_dihedral = atom->num_dihedral;
  int **dihedral_type = atom->dihedral_type;
  tagint **dihedral_atom1 = atom->dihedral_atom1;
  tagint **dihedral_atom2 = atom->dihedral_atom2;
  tagint **dihedral_atom3 = atom->dihedral_atom3;
  tagint **dihedral_atom4 = atom->dihedral_atom4;

  if ((m = idx2) >= 0) {
    if (num_dihedral[m] == atom->dihedral_per_atom)
      error->one(FLERR, "New dihedral exceeded dihedrals per atom in create_bonds");
    dihedral_type[m][num_dihedral[m]] = dtype;
    dihedral_atom1[m][num_dihedral[m]] = datom1;
    dihedral_atom2[m][num_dihedral[m]] = datom2;
    dihedral_atom3[m][num_dihedral[m]] = datom3;
    dihedral_atom4[m][num_dihedral[m]] = datom4;
    num_dihedral[m]++;
  }
  atom->ndihedrals++;

  if (force->newton_bond) return;

  if ((m = idx1) >= 0) {
    if (num_dihedral[m] == atom->dihedral_per_atom)
      error->one(FLERR, "New dihedral exceeded dihedrals per atom in create_bonds");
    dihedral_type[m][num_dihedral[m]] = dtype;
    dihedral_atom1[m][num_dihedral[m]] = datom1;
    dihedral_atom2[m][num_dihedral[m]] = datom2;
    dihedral_atom3[m][num_dihedral[m]] = datom3;
    dihedral_atom4[m][num_dihedral[m]] = datom4;
    num_dihedral[m]++;
  }

  if ((m = idx3) >= 0) {
    if (num_dihedral[m] == atom->dihedral_per_atom)
      error->one(FLERR, "New dihedral exceeded dihedrals per atom in create_bonds");
    dihedral_type[m][num_dihedral[m]] = dtype;
    dihedral_atom1[m][num_dihedral[m]] = datom1;
    dihedral_atom2[m][num_dihedral[m]] = datom2;
    dihedral_atom3[m][num_dihedral[m]] = datom3;
    dihedral_atom4[m][num_dihedral[m]] = datom4;
    num_dihedral[m]++;
  }

  if ((m = idx4) >= 0) {
    if (num_dihedral[m] == atom->dihedral_per_atom)
      error->one(FLERR, "New dihedral exceeded dihedrals per atom in create_bonds");
    dihedral_type[m][num_dihedral[m]] = dtype;
    dihedral_atom1[m][num_dihedral[m]] = datom1;
    dihedral_atom2[m][num_dihedral[m]] = datom2;
    dihedral_atom3[m][num_dihedral[m]] = datom3;
    dihedral_atom4[m][num_dihedral[m]] = datom4;
    num_dihedral[m]++;
  }
}

/* ---------------------------------------------------------------------- */

void CreateBonds::single_improper()
{
  int m;

  // check that 4 atoms exist

  const int nlocal = atom->nlocal;
  const int idx1 = atom->map(datom1);
  const int idx2 = atom->map(datom2);
  const int idx3 = atom->map(datom3);
  const int idx4 = atom->map(datom4);

  int count = 0;
  if ((idx1 >= 0) && (idx1 < nlocal)) count++;
  if ((idx2 >= 0) && (idx2 < nlocal)) count++;
  if ((idx3 >= 0) && (idx3 < nlocal)) count++;
  if ((idx4 >= 0) && (idx4 < nlocal)) count++;

  int allcount;
  MPI_Allreduce(&count, &allcount, 1, MPI_INT, MPI_SUM, world);
  if (allcount != 4) error->all(FLERR, "Create_bonds single/improper atoms do not exist");

  // create bond once or 4x if newton_bond set

  int *num_improper = atom->num_improper;
  int **improper_type = atom->improper_type;
  tagint **improper_atom1 = atom->improper_atom1;
  tagint **improper_atom2 = atom->improper_atom2;
  tagint **improper_atom3 = atom->improper_atom3;
  tagint **improper_atom4 = atom->improper_atom4;

  if ((m = idx2) >= 0) {
    if (num_improper[m] == atom->improper_per_atom)
      error->one(FLERR, "New improper exceeded impropers per atom in create_bonds");
    improper_type[m][num_improper[m]] = dtype;
    improper_atom1[m][num_improper[m]] = datom1;
    improper_atom2[m][num_improper[m]] = datom2;
    improper_atom3[m][num_improper[m]] = datom3;
    improper_atom4[m][num_improper[m]] = datom4;
    num_improper[m]++;
  }
  atom->nimpropers++;

  if (force->newton_bond) return;

  if ((m = idx1) >= 0) {
    if (num_improper[m] == atom->improper_per_atom)
      error->one(FLERR, "New improper exceeded impropers per atom in create_bonds");
    improper_type[m][num_improper[m]] = dtype;
    improper_atom1[m][num_improper[m]] = datom1;
    improper_atom2[m][num_improper[m]] = datom2;
    improper_atom3[m][num_improper[m]] = datom3;
    improper_atom4[m][num_improper[m]] = datom4;
    num_improper[m]++;
  }

  if ((m = idx3) >= 0) {
    if (num_improper[m] == atom->improper_per_atom)
      error->one(FLERR, "New improper exceeded impropers per atom in create_bonds");
    improper_type[m][num_improper[m]] = dtype;
    improper_atom1[m][num_improper[m]] = datom1;
    improper_atom2[m][num_improper[m]] = datom2;
    improper_atom3[m][num_improper[m]] = datom3;
    improper_atom4[m][num_improper[m]] = datom4;
    num_improper[m]++;
  }

  if ((m = idx4) >= 0) {
    if (num_improper[m] == atom->improper_per_atom)
      error->one(FLERR, "New improper exceeded impropers per atom in create_bonds");
    improper_type[m][num_improper[m]] = dtype;
    improper_atom1[m][num_improper[m]] = datom1;
    improper_atom2[m][num_improper[m]] = datom2;
    improper_atom3[m][num_improper[m]] = datom3;
    improper_atom4[m][num_improper[m]] = datom4;
    num_improper[m]++;
  }
}
